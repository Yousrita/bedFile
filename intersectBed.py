import pandas as pd
import streamlit as st
import gzip
import io
import numpy as np
from collections import defaultdict

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
    Effectue l'intersection entre deux DataFrames BED
    Version basique avec wa et wb
    """
    return intersect_bedtools_advanced(df1, df2, {'wa': True, 'wb': True})

def intersect_bedtools_advanced(df1, df2, options=None):
    """
    Version avancée avec plus d'options pour l'intersection
    Implémentation Python pure optimisée
    """
    if options is None:
        options = {}
    
    try:
        # Préparer les données avec seulement les 3 premières colonnes
        df_a = df1.iloc[:, :3].copy()
        df_b = df2.iloc[:, :3].copy()
        
        df_a.columns = ['chrom', 'start', 'end']
        df_b.columns = ['chrom', 'start', 'end']
        
        # Convertir en types numériques
        df_a['start'] = pd.to_numeric(df_a['start'], errors='coerce')
        df_a['end'] = pd.to_numeric(df_a['end'], errors='coerce')
        df_b['start'] = pd.to_numeric(df_b['start'], errors='coerce')
        df_b['end'] = pd.to_numeric(df_b['end'], errors='coerce')
        
        # Supprimer les lignes avec des valeurs manquantes
        df_a = df_a.dropna()
        df_b = df_b.dropna()
        
        results = []
        
        # Optimisation: grouper par chromosome
        chrom_groups_a = df_a.groupby('chrom')
        chrom_groups_b = df_b.groupby('chrom')
        
        common_chroms = set(chrom_groups_a.groups.keys()) & set(chrom_groups_b.groups.keys())
        
        if not common_chroms:
            st.info("Aucun chromosome commun entre les fichiers")
            return pd.DataFrame()
        
        progress_bar = st.progress(0)
        total_chroms = len(common_chroms)
        
        for i, chrom in enumerate(common_chroms):
            # Mettre à jour la barre de progression
            progress_bar.progress((i + 1) / total_chroms)
            
            # Filtrer par chromosome
            regions_a = chrom_groups_a.get_group(chrom)
            regions_b = chrom_groups_b.get_group(chrom)
            
            # Convertir en arrays numpy pour performance
            a_starts = regions_a['start'].values
            a_ends = regions_a['end'].values
            b_starts = regions_b['start'].values
            b_ends = regions_b['end'].values
            
            # Vérifier les intersections de manière vectorisée
            for idx_a, (start_a, end_a) in enumerate(zip(a_starts, a_ends)):
                # Trouver les régions B qui chevauchent avec cette région A
                overlaps = (b_starts < end_a) & (b_ends > start_a)
                
                if np.any(overlaps):
                    overlapping_indices = np.where(overlaps)[0]
                    
                    for idx_b in overlapping_indices:
                        start_b, end_b = b_starts[idx_b], b_ends[idx_b]
                        
                        # Calculer l'intersection
                        overlap_start = max(start_a, start_b)
                        overlap_end = min(end_a, end_b)
                        overlap_length = overlap_end - overlap_start
                        
                        if overlap_length > 0:
                            result_row = {
                                'A_chrom': chrom,
                                'A_start': start_a,
                                'A_end': end_a,
                                'B_chrom': chrom,
                                'B_start': start_b,
                                'B_end': end_b
                            }
                            
                            if options.get('wo', False):
                                result_row['overlap'] = overlap_length
                            
                            results.append(result_row)
        
        progress_bar.empty()
        
        if not results:
            st.info("ℹ️ Aucune intersection trouvée entre les fichiers")
            return pd.DataFrame()
        
        result_df = pd.DataFrame(results)
        
        # Adapter les colonnes en fonction des options
        if options.get('wo', False):
            # -wo ajoute une colonne avec le pourcentage d'overlap
            result_df = result_df[['A_chrom', 'A_start', 'A_end', 'B_chrom', 'B_start', 'B_end', 'overlap']]
        elif options.get('wa', True) and options.get('wb', True):
            # -wa -wb : 6 colonnes
            result_df = result_df[['A_chrom', 'A_start', 'A_end', 'B_chrom', 'B_start', 'B_end']]
        elif options.get('wa', True):
            # -wa seulement : 3 colonnes (fichier A)
            result_df = result_df[['A_chrom', 'A_start', 'A_end']].drop_duplicates()
        elif options.get('wb', True):
            # -wb seulement : 3 colonnes (fichier B)
            result_df = result_df[['B_chrom', 'B_start', 'B_end']].drop_duplicates()
        else:
            # Par défaut : 3 colonnes (fichier A)
            result_df = result_df[['A_chrom', 'A_start', 'A_end']].drop_duplicates()
        
        # Option -v: trouver les non-intersections
        if options.get('v', False):
            # Pour -v, on retourne les régions de A qui n'intersectent pas B
            all_a_indices = set(range(len(df_a)))
            intersecting_a_indices = set()
            
            # Identifier les indices de A qui intersectent
            for chrom in common_chroms:
                regions_a = chrom_groups_a.get_group(chrom)
                regions_b = chrom_groups_b.get_group(chrom)
                
                a_starts = regions_a['start'].values
                a_ends = regions_a['end'].values
                b_starts = regions_b['start'].values
                b_ends = regions_b['end'].values
                
                for idx_a, (start_a, end_a) in enumerate(zip(a_starts, a_ends)):
                    overlaps = (b_starts < end_a) & (b_ends > start_a)
                    if np.any(overlaps):
                        original_idx = regions_a.index[idx_a]
                        intersecting_a_indices.add(original_idx)
            
            # Prendre les régions de A qui n'intersectent pas
            non_intersecting = df_a[~df_a.index.isin(intersecting_a_indices)]
            result_df = non_intersecting[['chrom', 'start', 'end']].copy()
            result_df.columns = ['A_chrom', 'A_start', 'A_end']
        
        return result_df
        
    except Exception as e:
        st.error(f"Erreur lors de l'intersection: {str(e)}")
        return None

# Interface Streamlit
st.title("🔍 Intersection de fichiers BED")
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