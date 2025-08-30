import pandas as pd
import streamlit as st
import gzip
import io
import subprocess
import tempfile
import os

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
            st.info("Bedtools utilisera seulement les 3 premières colonnes (chrom, start, end) pour l'intersection")
    
    return df1, df2

def intersect_bedtools(df1, df2):
    """
    Effectue l'intersection entre deux DataFrames BED de structures différentes
    """
    try:
        # Créer des fichiers temporaires avec seulement les 3 premières colonnes
        with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f1:
            # ⚠️ CORRECTION : Utiliser seulement les 3 premières colonnes pour bedtools
            df1.iloc[:, :3].to_csv(f1.name, sep='\t', header=False, index=False)
            temp_file1 = f1.name
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f2:
            # ⚠️ CORRECTION : Utiliser seulement les 3 premières colonnes pour bedtools
            df2.iloc[:, :3].to_csv(f2.name, sep='\t', header=False, index=False)
            temp_file2 = f2.name
        
        # Exécuter bedtools intersect
        result = subprocess.run([
            'bedtools', 'intersect', 
            '-a', temp_file1, 
            '-b', temp_file2,
            '-wa', '-wb'
        ], capture_output=True, text=True)
        
        # Nettoyer les fichiers temporaires
        os.unlink(temp_file1)
        os.unlink(temp_file2)
        
        if result.returncode == 0:
            if not result.stdout.strip():
                st.info("ℹ️ Aucune intersection trouvée entre les fichiers")
                return pd.DataFrame()
            
            # Lire le résultat
            result_df = pd.read_csv(io.StringIO(result.stdout), sep='\t', header=None)
            
            # ⚠️ CORRECTION IMPORTANTE : 
            # Bedtools avec -wa -wb retourne TOUJOURS 6 colonnes :
            # Colonnes 1-3 : A_chrom, A_start, A_end (du fichier A)
            # Colonnes 4-6 : B_chrom, B_start, B_end (du fichier B)
            
            # Nommer les colonnes de base
            result_df.columns = ['A_chrom', 'A_start', 'A_end', 'B_chrom', 'B_start', 'B_end']
            
            # ⭐ OPTIONNEL : Si vous voulez garder les colonnes supplémentaires originales
            # Vous pouvez faire un merge avec les données originales plus tard
            
            return result_df
        else:
            st.error(f"Erreur bedtools: {result.stderr}")
            return None
            
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
        
        # Construire la commande bedtools
        cmd = ['bedtools', 'intersect', '-a', temp_file1, '-b', temp_file2]
        
        # Ajouter les options
        if options.get('wa', False):
            cmd.append('-wa')
        if options.get('wb', False):
            cmd.append('-wb')
        if options.get('wo', False):
            cmd.append('-wo')
        if options.get('v', False):
            cmd.append('-v')
        if options.get('f', 0) > 0:
            cmd.extend(['-f', str(options['f'])])
        
        # Exécuter bedtools intersect
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        # Nettoyer les fichiers temporaires
        os.unlink(temp_file1)
        os.unlink(temp_file2)
        
        if result.returncode == 0:
            if not result.stdout.strip():
                return pd.DataFrame()
            
            # Lire le résultat
            result_df = pd.read_csv(io.StringIO(result.stdout), sep='\t', header=None)
            
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
        else:
            st.error(f"Erreur bedtools: {result.stderr}")
            return None
            
    except Exception as e:
        st.error(f"Erreur lors de l'intersection: {str(e)}")
        return None