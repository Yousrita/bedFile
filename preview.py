import pandas as pd
import streamlit as st
from typing import Optional

def load_bed(uploaded_file) -> Optional[pd.DataFrame]:
    """
    Load a BED file with ALL its columns
    """
    try:
        # Read the file without column limit
        df = pd.read_csv(
            uploaded_file,
            sep="\t",
            header=None,
            comment='#',
            dtype=str  # ← Read everything as string first
        )
        
        # Define standard BED column names
        bed_cols = ['chrom', 'start', 'end', 'name', 'score', 'strand', 
                   'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 
                   'blockSizes', 'blockStarts']
        
        # Name available columns
        df.columns = bed_cols[:len(df.columns)]
        
        # Convert start and end to numeric
        if 'start' in df.columns:
            df['start'] = pd.to_numeric(df['start'], errors='coerce')
        if 'end' in df.columns:
            df['end'] = pd.to_numeric(df['end'], errors='coerce')
        
        # Validate
        if len(df.columns) < 3:
            st.error("❌ Invalid BED format: minimum 3 columns required")
            return None
            
        if (df["start"] > df["end"]).any():
            st.error("❌ Error: start positions > end positions detected")
            return None
            
        # REMOVED: st.success message about file loading
        return df
        
    except Exception as e:
        st.error(f"❌ Loading error: {str(e)}")
        return None

def load_bed_basic(uploaded_file) -> Optional[pd.DataFrame]:
    """
    Load only the first 3 columns (for bedtools operations)
    """
    df = load_bed(uploaded_file)
    if df is not None:
        return df[['chrom', 'start', 'end']].copy()
    return None

def show_metrics_banner(df, genome_size=3_299_210_039):
    """Display a comprehensive single-line file summary"""
    if df is None or df.empty:
        return
        
    # Calculate all metrics
    length = df['end'] - df['start']
    duplicates = df.duplicated().sum()
    null_values = df.isnull().sum().sum()
    total_bases = length.sum()
    coverage_percent = (total_bases / genome_size) * 100
    
    # Check if strand column exists
    strand_stats = df['strand'].value_counts(normalize=True).to_dict() if 'strand' in df else None
    
    # Improved stats
    stats = {
        "🧬 Genes": df['chrom'].nunique(),
        "📊 Intervals": f"{len(df):,}",
        "📏 Avg length": f"{length.mean():.0f} bp",
        "🎯 Min length": f"{length.min():,} bp", 
        "🎯 Max length": f"{length.max():,} bp",
        "🌐 Coverage": f"{coverage_percent:.1f}%",
        "🧩 Total bases": f"{total_bases/1e6:.1f} Mb",
        "🔄 Duplicates": duplicates,
        "❌ Null values": null_values
    }
    
    # Add strand metrics if available
    if strand_stats:
        stats.update({
            "➕ + strand": f"{strand_stats.get('+', 0):.1%}",
            "➖ - strand": f"{strand_stats.get('-', 0):.1%}"
        })
    
    st.write("""
    <div style='background: #f8f9fa;
                border-radius: 8px;
                padding: 12px;
                margin: 10px 0;
                border-left: 4px solid #2e86ab;
                display: flex;
                justify-content: space-between;
                align-items: center;
                flex-wrap: wrap;
                gap: 10px;'>
    """, unsafe_allow_html=True)
    
    # Create columns dynamically based on stats count
    cols = st.columns(len(stats))
    for i, (name, value) in enumerate(stats.items()):
        with cols[i]:
            st.write(f"""
            <div style='text-align: center;
                        min-width: 80px;
                        padding: 5px;
                        border-radius: 6px;
                        background: rgba(255,255,255,0.7);
                        border: 1px solid #e0e0e0;'>
                <div style='font-size: 11px; color: #555; margin-bottom: 4px;'>{name}</div>
                <div style='font-size: 13px; font-weight: 600; color: #2c3e50;'>{value}</div>
            </div>
            """, unsafe_allow_html=True)
    
    st.write("</div>", unsafe_allow_html=True)

def show_data_preview(df: pd.DataFrame) -> None:
    """Display data preview with ALL columns"""
    tab1, tab2 = st.tabs(["📋 Data Preview", "🧪 Data Quality"])
    
    with tab1:
        # Display ALL columns
        st.dataframe(df.head(100), height=300)
        st.caption(f"Displaying 100 rows out of {len(df):,} total - {len(df.columns)} columns")
        
        # Display column list
        st.write("**📝 Available columns:**")
        for i, col in enumerate(df.columns):
            st.write(f"- `{col}` ({df[col].dtype})")
    
    with tab2:
        show_data_quality(df)

def show_data_quality(df: pd.DataFrame) -> None:
    """Display data quality statistics"""
    st.subheader("🧬 Distinct Values")
    col1, col2, col3, col4 = st.columns(4)
    
    with col1:
        st.metric("Chromosomes", df["chrom"].nunique())
    with col2:
        st.metric("Total size", f"{(df['end'] - df['start']).sum():,} bp")
    with col3:
        st.metric("Null values", df.isnull().sum().sum())
    with col4:
        st.metric("Columns", len(df.columns))
    
    # Missing values analysis by column
    st.subheader("🔍 Null Values Detail")
    null_details = df.isnull().sum().to_frame("Null Values")
    null_details["Percentage"] = (null_details["Null Values"] / len(df) * 100).round(2)
    st.table(null_details)
    
    # Data types by column
    st.subheader("📊 Data Types")
    dtype_details = pd.DataFrame({
        'Column': df.columns,
        'Type': [str(dtype) for dtype in df.dtypes],
        'Unique Values': [df[col].nunique() for col in df.columns]
    })
    st.table(dtype_details)
    
    # Strand analysis if column exists
    if 'strand' in df:
        st.subheader("⚖ Strand Distribution")
        strand_dist = df['strand'].value_counts(normalize=True)
        st.bar_chart(strand_dist)
    
    # Score analysis if column exists
    if 'score' in df:
        st.subheader("📈 Score Distribution")
        st.write(f"Average score: {df['score'].mean():.2f}")
        st.write(f"Min score: {df['score'].min()}")
        st.write(f"Max score: {df['score'].max()}")