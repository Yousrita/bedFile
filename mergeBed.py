import pandas as pd
import streamlit as st
from preview import load_bed, show_metrics_banner, show_data_preview  # ⬅️ CORRIGÉ ICI

def merge_bed_files(df: pd.DataFrame) -> pd.DataFrame:
    """Merges a BED dataframe with itself"""
    # Columns normalization
    cols = ['chrom', 'start', 'end']
    #df = df[cols].copy()
    df
    
    # Numeric conversion and cleaning
    df['start'] = pd.to_numeric(df['start'], errors='coerce')
    df['end'] = pd.to_numeric(df['end'], errors='coerce')
    df.dropna(inplace=True)
    
    # Interval sorting
    df = df.sort_values(cols)
    
    # Interval merging
    merged = []
    current = None
    
    for _, row in df.iterrows():
        if current is None:
            current = dict(row)
        elif (row['chrom'] == current['chrom'] and 
              row['start'] <= current['end']):
            current['end'] = max(current['end'], row['end'])
        else:
            merged.append(current)
            current = dict(row)
    
    if current:
        merged.append(current)
    
    return pd.DataFrame(merged)

def handle_merge_operation(df: pd.DataFrame):
    """Handles the merge interface"""
    st.write("---")
    st.subheader("🔀 Merge file")
    
    # Display before merging
    with st.expander("📁 Before Merge", expanded=True):
        show_metrics_banner(df)
    
    # Automatic merging
    with st.spinner("Merging in progress..."):
        result = merge_bed_files(df)
        
        # Result display
        st.success("Merge completed!")
        with st.expander("📊 After Merge", expanded=True):
            show_metrics_banner(result)
            st.dataframe(result.head(100))
        
        # Download
        st.download_button(
            "💾 Download Merged Result",
            result.to_csv(index=False, sep='\t').encode('utf-8'),
            "merged_result.bed",
            "text/plain"
        )