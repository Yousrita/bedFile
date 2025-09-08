import pandas as pd
import streamlit as st
import gzip
import io
import tempfile
import os
import pybedtools

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
            st.info("Pybedtools utilisera seulement les 3 premières colonnes (chrom, start, end) pour l'intersection")
    
    return df1, df2

def intersect_bedtools(df1, df2):
    """
    Effectue l'intersection entre deux DataFrames BED avec pybedtools
    Version basique avec wa et wb
    """
    try:
        # Créer des fichiers temporaires avec seulement les 3 premières colonnes
        with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f1:
            df1.iloc[:, :3].to_csv(f1.name, sep='\t', header=False, index=False)
            temp_file1 = f1.name
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f2:
            df2.iloc[:, :3].to_csv(f2.name, sep='\t', header=False, index=False)
            temp_file2 = f2.name
        
        # Créer les objets BedTool et effectuer l'intersection
        a = pybedtools.BedTool(temp_file1)
        b = pybedtools.BedTool(temp_file2)
        
        result = a.intersect(b, wa=True, wb=True)
        
        # Nettoyer les fichiers temporaires
        os.unlink(temp_file1)
        os.unlink(temp_file2)
        
        if len(result) == 0:
            st.info("ℹ️ Aucune intersection trouvée entre les fichiers")
            return pd.DataFrame()
        
        # Convertir en DataFrame
        result_df = result.to_dataframe(disable_auto_names=True, header=None)
        
        # Bedtools avec -wa -wb retourne TOUJOURS 6 colonnes :
        # Colonnes 1-3 : A_chrom, A_start, A_end (du fichier A)
        # Colonnes 4-6 : B_chrom, B_start, B_end (du fichier B)
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
        # Créer des fichiers temporaires avec seulement les 3 premières colonnes
        with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f1:
            df1.iloc[:, :3].to_csv(f1.name, sep='\t', header=False, index=False)
            temp_file1 = f1.name
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f2:
            df2.iloc[:, :3].to_csv(f2.name, sep='\t', header=False, index=False)
            temp_file2 = f2.name
        
        # Créer les objets BedTool
        a = pybedtools.BedTool(temp_file1)
        b = pybedtools.BedTool(temp_file2)
        
        # Configurer les options d'intersection
        intersect_kwargs = {}
        if options.get('wa', False):
            intersect_kwargs['wa'] = True
        if options.get('wb', False):
            intersect_kwargs['wb'] = True
        if options.get('wo', False):
            intersect_kwargs['wo'] = True
        if options.get('v', False):
            intersect_kwargs['v'] = True
        if options.get('f', 0) > 0:
            intersect_kwargs['f'] = options['f']
        
        # Effectuer l'intersection
        result = a.intersect(b, **intersect_kwargs)
        
        # Nettoyer les fichiers temporaires
        os.unlink(temp_file1)
        os.unlink(temp_file2)
        
        # Convertir le résultat en DataFrame
        if len(result) == 0:
            return pd.DataFrame()
        
        result_df = result.to_dataframe(disable_auto_names=True, header=None)
        
        # Adapter les colonnes en fonction des options
        if options.get('wo', False):
            # -wo ajoute une colonne avec le pourcentage d'overlap
            result_df.columns = ['A_chrom', 'A_start', 'A_end', 'B_chrom', 'B_start', 'B_end', 'overlap']
        elif options.get('wa', False) and options.get('wb', False):
            # -wa -wb : 6 colonnes
            result_df.columns = ['A_chrom', 'A_start', 'A_end', 'B_chrom', 'B_start', 'B_end']
        elif options.get('wa', False):
            # -wa seulement : 3 colonnes (fichier A)
            result_df.columns = ['A_chrom', 'A_start', 'A_end']
        elif options.get('wb', False):
            # -wb seulement : 3 colonnes (fichier B)
            result_df.columns = ['B_chrom', 'B_start', 'B_end']
        else:
            # Par défaut : 3 colonnes (fichier A)
            result_df.columns = ['A_chrom', 'A_start', 'A_end']
        
        return result_df
            
    except Exception as e:
        st.error(f"Erreur lors de l'intersection: {str(e)}")
        return None

# Interface Streamlit
st.title("🔍 Intersection de fichiers BED avec PyBedTools")
st.write("Cet outil permet de trouver les intersections entre deux fichiers BED")

# Upload des fichiers
file1 = st.file_uploader("Télécharger le premier fichier BED (A)", type=['bed', 'bed.gz'])
file2 = st.file_uploader("Télécharger le second fichier BED (B)", type=['bed', 'bed.gz'])

if file1 and file2:
    # Charger les fichiers
    df1, df2 = load_bed_int(file1, file2)
    
    if df1 is not None and df2 is not None:
        st.success("✅ Files loaded successfully")
        st.write(f"**File A:** {len(df1)} regions")
        st.write(f"**File B:** {len(df2)} regions")
        
        # Options d'intersection
        st.subheader("Options d'intersection")
        col1, col2, col3 = st.columns(3)
        with col1:
            wa = st.checkbox("Inclure régions A (-wa)", value=True)
            wb = st.checkbox("Inclure régions B (-wb)", value=True)
        with col2:
            wo = st.checkbox("Ajouter longueur overlap (-wo)", value=False)
            v = st.checkbox("Trouver non-intersections (-v)", value=False)
        with col3:
            min_overlap = st.slider("Recouvrement minimum (%)", 0, 100, 0)
        
        # Bouton pour lancer l'intersection
        if st.button("Lancer l'intersection"):
            with st.spinner("Calcul de l'intersection en cours..."):
                options = {
                    'wa': wa,
                    'wb': wb,
                    'wo': wo,
                    'v': v,
                    'f': min_overlap / 100.0
                }
                
                result_df = intersect_bedtools_advanced(df1, df2, options)
                
                if result_df is not None:
                    if len(result_df) > 0:
                        st.success(f"✅ Intersection terminée: {len(result_df)} régions trouvées")
                        st.dataframe(result_df.head(100))
                        
                        # Téléchargement des résultats
                        csv = result_df.to_csv(index=False, sep='\t')
                        st.download_button(
                            label="Télécharger les résultats",
                            data=csv,
                            file_name="intersection_results.bed",
                            mime="text/tab-separated-values"
                        )
                    else:
                        st.info("Aucune intersection trouvée avec les paramètres actuels")