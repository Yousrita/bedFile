import pandas as pd
import streamlit as st

def sort_bed(df):
    """
    Trie un DataFrame BED par les colonnes chrom, start, end.
    """
    return df.sort_values(by=["chrom", "start", "end"])

def handle_sort_operation(df: pd.DataFrame):
    """Gère l'opération de tri, affiche le résultat et propose le téléchargement"""
    df_sorted = sort_bed(df)
    st.success("Sort Completed!")
    
    # Affichage du DataFrame trié une seule fois
    st.dataframe(df_sorted)
    
    # Proposition de téléchargement
    st.download_button(
        "💾 Download Sorted Result",
        df_sorted.to_csv(index=False, sep='\t').encode('utf-8'),
        "resultat_trie.bed",
        "text/plain"
    )

# Exemple d'utilisation (supposant que df est déjà défini)
# handle_sort_operation(df)
