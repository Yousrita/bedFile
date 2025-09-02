import pandas as pd
import streamlit as st
import gzip
import io
import tempfile
import os
import pybedtools
from pybedtools import BedTool

def read_bed_file(uploaded_file):
    """
    Lit un fichier BED (normal ou compressé) et retourne un DataFrame
    """
    try:
        # Vérifier si c'est un fichier compressé
        if uploaded_file.name.endswith('.gz'):
            with gzip.open(uploaded_file, 'rt') as f:
                content = f.read()
            file_obj = io.StringIO(content)
        else:
            content = uploaded_file.getvalue().decode('utf-8')
            file_obj = io.StringIO(content)
        
        # Détecter le nombre de colonnes
        first_line = content.split('\n')[0]
        num_columns = len(first_line.split('\t'))
        
        bed_cols = ['chrom', 'start', 'end', 'name', 'score', 'strand', 
                   'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 
                   'blockSizes', 'blockStarts']
        
        df = pd.read_csv(file_obj, sep='\t', header=None, 
                        names=bed_cols[:num_columns], comment='#')
        
        return df
        
    except Exception as e:
        st.error(f"Erreur lecture fichier: {str(e)}")
        return None

def load_bed_int(file1, file2):
    """
    Charge deux fichiers BED pour l'intersection
    """
    df1 = read_bed_file(file1)
    df2 = read_bed_file(file2)
    
    # Afficher la structure des fichiers pour débogage
    if df1 is not None and df2 is not None:
        st.write(f"**Structure Fichier A:** {df1.shape[1]} colonnes - {list(df1.columns)}")
        st.write(f"**Structure Fichier B:** {df2.shape[1]} colonnes - {list(df2.columns)}")
        
        if df1.shape[1] != df2.shape[1]:
            st.warning("⚠️ Les fichiers ont des structures différentes!")
            st.info("L'intersection utilisera seulement les 3 premières colonnes (chrom, start, end)")
    
    return df1, df2

def intersect_bedtools(df1, df2):
    """
    Effectue l'intersection entre deux DataFrames BED de structures différentes
    """
    try:
        # Créer des objets BedTool à partir des DataFrames avec les 3 premières colonnes
        bed1 = BedTool.from_dataframe(df1[['chrom', 'start', 'end']].copy())
        bed2 = BedTool.from_dataframe(df2[['chrom', 'start', 'end']].copy())
        
        # Effectuer l'intersection avec pybedtools
        result_bed = bed1.intersect(bed2, wa=True, wb=True)
        
        # Vérifier si il y a des résultats
        if len(result_bed) == 0:
            st.info("ℹ️ Aucune intersection trouvée entre les fichiers")
            return pd.DataFrame()
        
        # Convertir le résultat en DataFrame
        result_df = result_bed.to_dataframe()
        
        # Nommer les colonnes (pybedtools retourne 6 colonnes pour wa=True, wb=True)
        if len(result_df.columns) >= 6:
            result_df.columns = ['A_chrom', 'A_start', 'A_end', 'B_chrom', 'B_start', 'B_end']
        
        return result_df
            
    except Exception as e:
        st.error(f"Erreur lors de l'intersection: {str(e)}")
        return None

def intersect_bedtools_advanced(df1, df2, options=None):
    """
    Version avancée avec plus d'options pour l'intersection
    """
    if options is None:
        options = {}
    
    try:
        # Créer des objets BedTool à partir des DataFrames avec les 3 premières colonnes
        bed1 = BedTool.from_dataframe(df1[['chrom', 'start', 'end']].copy())
        bed2 = BedTool.from_dataframe(df2[['chrom', 'start', 'end']].copy())
        
        # Effectuer l'intersection avec les options
        result_bed = bed1.intersect(bed2, 
                                   wa=options.get('wa', False),
                                   wb=options.get('wb', False),
                                   wo=options.get('wo', False),
                                   v=options.get('v', False),
                                   f=options.get('f', 1e-9))  # Valeur par défaut très petite pour -f
        
        # Vérifier si il y a des résultats
        if len(result_bed) == 0:
            return pd.DataFrame()
        
        # Convertir le résultat en DataFrame
        result_df = result_bed.to_dataframe()
        
        # Adapter les colonnes en fonction des options
        if options.get('wo', False) and len(result_df.columns) >= 7:
            # -wo ajoute une colonne avec le nombre de paires de bases d'overlap
            result_df.columns = ['A_chrom', 'A_start', 'A_end', 'B_chrom', 'B_start', 'B_end', 'overlap']
        elif options.get('wa', False) and options.get('wb', False) and len(result_df.columns) >= 6:
            # -wa -wb : 6 colonnes
            result_df.columns = ['A_chrom', 'A_start', 'A_end', 'B_chrom', 'B_start', 'B_end']
        elif options.get('wa', False) and len(result_df.columns) >= 3:
            # -wa seulement : 3 colonnes (fichier A)
            result_df.columns = ['A_chrom', 'A_start', 'A_end']
        elif options.get('wb', False) and len(result_df.columns) >= 3:
            # -wb seulement : 3 colonnes (fichier B)
            result_df.columns = ['B_chrom', 'B_start', 'B_end']
        elif len(result_df.columns) >= 3:
            # Par défaut : 3 colonnes (fichier A)
            result_df.columns = ['A_chrom', 'A_start', 'A_end']
        
        return result_df
            
    except Exception as e:
        st.error(f"Erreur lors de l'intersection: {str(e)}")
        return None

# Interface Streamlit (identique à l'original)
def main():
    st.title("🔬 Intersection de fichiers BED")
    st.write("Cet outil permet de trouver les intersections entre deux fichiers BED")
    
    # Upload des fichiers
    file1 = st.file_uploader("Choisir le premier fichier BED (.bed ou .bed.gz)", type=['bed', 'bed.gz'])
    file2 = st.file_uploader("Choisir le deuxième fichier BED (.bed ou .bed.gz)", type=['bed', 'bed.gz'])
    
    if file1 and file2:
        df1, df2 = load_bed_int(file1, file2)
        
        if df1 is not None and df2 is not None:
            st.success("✅ Fichiers chargés avec succès!")
            
            # Options d'intersection
            st.subheader("Options d'intersection")
            col1, col2, col3 = st.columns(3)
            
            with col1:
                wa = st.checkbox("wa", value=True, 
                                help="Écrire la fonctionnalité originale A pour chaque intersection")
            with col2:
                wb = st.checkbox("wb", value=True,
                                help="Écrire la fonctionnalité originale B pour chaque intersection")
            with col3:
                wo = st.checkbox("wo", value=False,
                                help="Écrire la fonctionnalité originale A, B plus le nombre de paires de bases d'intersection")
            
            v = st.checkbox("v", value=False,
                           help="N'écrire qu'entrées dans A qui n'ont AUCUNE intersection avec B")
            
            # Option -f (fraction overlap)
            f_value = st.slider("Fraction d'overlap minimum (-f)", 
                               min_value=0.0, max_value=1.0, value=0.0, step=0.01,
                               help="Fraction minimale de recouvrement requis")
            
            # Bouton pour lancer l'intersection
            if st.button("Lancer l'intersection"):
                with st.spinner("Calcul de l'intersection en cours..."):
                    options = {
                        'wa': wa,
                        'wb': wb,
                        'wo': wo,
                        'v': v,
                        'f': f_value if f_value > 0 else 1e-9  # Éviter 0 exact pour pybedtools
                    }
                    
                    result_df = intersect_bedtools_advanced(df1, df2, options)
                    
                    if result_df is not None:
                        if not result_df.empty:
                            st.success(f"✅ Intersection terminée! {len(result_df)} régions trouvées")
                            
                            # Afficher le résultat
                            st.subheader("Résultat de l'intersection")
                            st.dataframe(result_df)
                            
                            # Télécharger le résultat
                            csv = result_df.to_csv(index=False, sep='\t')
                            st.download_button(
                                label="Télécharger le résultat",
                                data=csv,
                                file_name="intersection_result.bed",
                                mime="text/tab-separated-values"
                            )
                        else:
                            st.info("ℹ️ Aucune intersection trouvée entre les fichiers")

if __name__ == "__main__":
    main()