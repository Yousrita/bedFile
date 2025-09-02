import pandas as pd
import numpy as np
import io

def read_bed_file(uploaded_file):
    """
    Lit un fichier BED (normal ou compressé) et retourne un DataFrame
    """
    try:
        content = uploaded_file.getvalue().decode('utf-8')
        
        # Détecter le nombre de colonnes
        first_line = content.split('\n')[0]
        num_columns = len(first_line.split('\t'))
        
        bed_cols = ['chrom', 'start', 'end', 'name', 'score', 'strand', 
                   'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 
                   'blockSizes', 'blockStarts']
        
        # CORRECTION : Utiliser io.StringIO au lieu de pandas.compat.StringIO
        df = pd.read_csv(io.StringIO(content), sep='\t', header=None, 
                        names=bed_cols[:num_columns], comment='#')
        
        return df
        
    except Exception as e:
        raise Exception(f"Erreur lecture fichier: {str(e)}")

def load_bed_int(file1, file2):
    """
    Charge deux fichiers BED pour l'intersection
    """
    df1 = read_bed_file(file1)
    df2 = read_bed_file(file2)
    
    return df1, df2

def python_pure_intersect(df1, df2, options=None):
    """
    Intersection BED en pur Python - Sans pybedtools!
    """
    if options is None:
        options = {'wa': True, 'wb': True}
    
    results = []
    
    # Vérifier que les DataFrames ont les colonnes nécessaires
    required_cols = ['chrom', 'start', 'end']
    for col in required_cols:
        if col not in df1.columns:
            raise ValueError(f"DataFrame 1 doit avoir la colonne '{col}'")
        if col not in df2.columns:
            raise ValueError(f"DataFrame 2 doit avoir la colonne '{col}'")
    
    # Optimisation: travailler par chromosome
    common_chroms = set(df1['chrom']).intersection(set(df2['chrom']))
    
    if not common_chroms:
        print("⚠️ Aucun chromosome commun entre les fichiers")
        return pd.DataFrame()
    
    for chrom in common_chroms:
        df1_chrom = df1[df1['chrom'] == chrom]
        df2_chrom = df2[df2['chrom'] == chrom]
        
        # Tri pour optimisation
        df1_chrom = df1_chrom.sort_values('start').reset_index(drop=True)
        df2_chrom = df2_chrom.sort_values('start').reset_index(drop=True)
        
        for _, row1 in df1_chrom.iterrows():
            # Filtrer les régions de df2 qui pourraient chevaucher
            potential_overlaps = df2_chrom[
                (df2_chrom['start'] <= row1['end']) & 
                (df2_chrom['end'] >= row1['start'])
            ]
            
            for _, row2 in potential_overlaps.iterrows():
                overlap_start = max(row1['start'], row2['start'])
                overlap_end = min(row1['end'], row2['end'])
                
                if overlap_start < overlap_end:
                    result_row = {}
                    
                    if options.get('wa', True):
                        result_row.update({
                            'A_chrom': row1['chrom'],
                            'A_start': row1['start'],
                            'A_end': row1['end']
                        })
                    
                    if options.get('wb', True):
                        result_row.update({
                            'B_chrom': row2['chrom'],
                            'B_start': row2['start'],
                            'B_end': row2['end']
                        })
                    
                    if options.get('wo', False):
                        result_row['overlap'] = overlap_end - overlap_start
                    
                    results.append(result_row)
    
    if not results:
        return pd.DataFrame()
    
    return pd.DataFrame(results)

# Aliases pour la compatibilité avec l'ancien code
intersect_bedtools = python_pure_intersect
intersect_bedtools_advanced = python_pure_intersect

# Test local si exécuté directement
if __name__ == "__main__":
    # Créer des données de test
    test_data1 = {
        'chrom': ['chr1', 'chr1', 'chr2'],
        'start': [100, 500, 200],
        'end': [300, 800, 400]
    }
    
    test_data2 = {
        'chrom': ['chr1', 'chr1', 'chr3'],
        'start': [200, 600, 100],
        'end': [400, 900, 300]
    }
    
    df1 = pd.DataFrame(test_data1)
    df2 = pd.DataFrame(test_data2)
    
    print("Test d'intersection Python pure:")
    print("Fichier 1:")
    print(df1)
    print("\nFichier 2:")
    print(df2)
    
    # Test d'intersection
    result = python_pure_intersect(df1, df2)
    print(f"\nRésultat de l'intersection ({len(result)} overlaps):")
    print(result)