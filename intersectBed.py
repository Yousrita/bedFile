import pandas as pd
import streamlit as st
import subprocess
import tempfile
import os
from io import StringIO

# --- Vérification de l'installation de bedtools ---
def check_bedtools_installed():
    """Vérifie si bedtools est installé et accessible"""
    try:
        result = subprocess.run(['bedtools', '--version'], 
                              capture_output=True, text=True, timeout=10)
        if result.returncode == 0 and "bedtools" in result.stdout.lower():
            return True
        return False
    except (subprocess.TimeoutExpired, FileNotFoundError, subprocess.SubprocessError):
        return False

# --- Fallback manuel si bedtools échoue ---
def intersect_manual(df1, df2):
    """
    Implémentation manuelle de l'intersection comme fallback
    """
    try:
        results = []
        
        # Assurer que les colonnes sont numériques
        df1 = df1.copy()
        df2 = df2.copy()
        df1['start'] = pd.to_numeric(df1['start'], errors='coerce')
        df1['end'] = pd.to_numeric(df1['end'], errors='coerce')
        df2['start'] = pd.to_numeric(df2['start'], errors='coerce')
        df2['end'] = pd.to_numeric(df2['end'], errors='coerce')
        
        # Supprimer les lignes avec des valeurs NaN
        df1 = df1.dropna(subset=['start', 'end'])
        df2 = df2.dropna(subset=['start', 'end'])
        
        for _, row1 in df1.iterrows():
            chrom1, start1, end1 = row1['chrom'], int(row1['start']), int(row1['end'])
            
            if pd.isna(start1) or pd.isna(end1) or start1 >= end1:
                continue
            
            # Trouver les chevauchements dans df2 pour le même chromosome
            overlaps = df2[(df2['chrom'] == chrom1)].copy()
            
            overlaps['start'] = pd.to_numeric(overlaps['start'], errors='coerce')
            overlaps['end'] = pd.to_numeric(overlaps['end'], errors='coerce')
            overlaps = overlaps.dropna(subset=['start', 'end'])
            
            for _, row2 in overlaps.iterrows():
                start2, end2 = int(row2['start']), int(row2['end'])
                
                if pd.isna(start2) or pd.isna(end2) or start2 >= end2:
                    continue
                
                # Calculer le chevauchement
                overlap_start = max(start1, start2)
                overlap_end = min(end1, end2)
                
                if overlap_start < overlap_end:  # Il y a chevauchement
                    results.append({
                        'chrom_a': chrom1, 'start_a': start1, 'end_a': end1,
                        'chrom_b': row2['chrom'], 'start_b': start2, 'end_b': end2,
                        'overlap_length': overlap_end - overlap_start
                    })
        
        return pd.DataFrame(results)
        
    except Exception as e:
        st.error(f"❌ Erreur lors de l'intersection manuelle : {e}")
        return None

# --- Fonction d'intersection SIMPLIFIÉE ---
def intersect_bedtools(df1, df2):
    """
    Effectue l'intersection entre deux DataFrames BED avec l'option classique -wa -wb
    """
    # Vérifier et nettoyer les données
    df1_clean = clean_bed_data(df1)
    df2_clean = clean_bed_data(df2)
    
    if df1_clean is None or df2_clean is None:
        st.warning("Données BED invalides, utilisation du mode manuel")
        return intersect_manual(df1, df2)
    
    try:
        # Créer des fichiers temporaires avec un format BED propre
        with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f1:
            df1_clean.to_csv(f1, sep='\t', header=False, index=False)
            temp_file1 = f1.name
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f2:
            df2_clean.to_csv(f2, sep='\t', header=False, index=False)
            temp_file2 = f2.name
        
        try:
            # Commande bedtools classique : -wa -wb (écrit les deux fichiers)
            cmd = [
                'bedtools', 'intersect',
                '-a', temp_file1,
                '-b', temp_file2,
                '-wa', '-wb',
                '-sorted'
            ]
            
            # Exécuter la commande
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
            
            if result.returncode == 0:
                if result.stdout.strip():
                    # Lire le résultat en DataFrame - CORRECTION StringIO
                    result_df = pd.read_csv(StringIO(result.stdout), sep='\t', header=None)
                    
                    # Nommer les colonnes
                    n_cols_a = len(df1_clean.columns)
                    n_cols_b = len(df2_clean.columns)
                    
                    col_names = []
                    for i in range(n_cols_a):
                        col_names.append(f'A_col_{i+1}')
                    for i in range(n_cols_b):
                        col_names.append(f'B_col_{i+1}')
                    
                    result_df.columns = col_names[:len(result_df.columns)]
                    return result_df
                else:
                    return pd.DataFrame()  # Aucun résultat
            else:
                st.error(f"❌ Erreur bedtools: {result.stderr}")
                # Fallback vers le mode manuel
                return intersect_manual(df1, df2)
                
        finally:
            # Nettoyer les fichiers temporaires
            for f in [temp_file1, temp_file2]:
                if os.path.exists(f):
                    os.unlink(f)
            
    except Exception as e:
        st.error(f"❌ Erreur lors de l'intersection avec bedtools : {str(e)}")
        st.info("🔄 Utilisation du mode manuel...")
        return intersect_manual(df1, df2)

def clean_bed_data(df):
    """Nettoie et valide les données BED"""
    try:
        df_clean = df.copy()
        
        # Vérifier les colonnes requises
        required_cols = ['chrom', 'start', 'end']
        if not all(col in df_clean.columns for col in required_cols):
            return None
        
        # Convertir en numérique
        df_clean['start'] = pd.to_numeric(df_clean['start'], errors='coerce')
        df_clean['end'] = pd.to_numeric(df_clean['end'], errors='coerce')
        
        # Supprimer les lignes invalides
        df_clean = df_clean.dropna(subset=['start', 'end'])
        df_clean = df_clean[df_clean['start'] < df_clean['end']]
        
        # Trier par chromosome et position
        df_clean = df_clean.sort_values(['chrom', 'start', 'end'])
        
        return df_clean[['chrom', 'start', 'end']]
        
    except Exception as e:
        st.error(f"Erreur lors du nettoyage des données BED: {e}")
        return None

# --- Chargement de deux fichiers BED en DataFrames ---
def load_bed_int(file1, file2):
    """Charge deux fichiers BED en DataFrame"""
    try:
        # Lecture des fichiers
        df1 = pd.read_csv(file1, sep='\t', comment='#', header=None)
        df2 = pd.read_csv(file2, sep='\t', comment='#', header=None)

        # Renommage des colonnes
        df1 = df1.rename(columns={0: 'chrom', 1: 'start', 2: 'end'})
        df2 = df2.rename(columns={0: 'chrom', 1: 'start', 2: 'end'})

        return df1[['chrom', 'start', 'end']], df2[['chrom', 'start', 'end']]

    except Exception as e:
        st.error(f"❌ Erreur lors du chargement des fichiers : {e}")
        return None, None