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
    """Display a comprehensive single-line file summary with colored backgrounds"""
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
    
    # Improved stats with professional color scheme
    stats = {
        "🧬 Genes": {"value": df['chrom'].nunique(), "color": "#e8f5e9"},
        "📊 Intervals": {"value": f"{len(df):,}", "color": "#e3f2fd"},
        "📏 Avg length": {"value": f"{length.mean():.0f} bp", "color": "#fff3e0"},
        "🎯 Min length": {"value": f"{length.min():,} bp", "color": "#ffebee"},
        "🎯 Max length": {"value": f"{length.max():,} bp", "color": "#f3e5f5"},
        "🌐 Coverage": {"value": f"{coverage_percent:.1f}%", "color": "#e8eaf6"},
        "🧩 Total bases": {"value": f"{total_bases/1e6:.1f} Mb", "color": "#e0f2f1"},
        "🔄 Duplicates": {"value": duplicates, "color": "#fff8e1"},
        "❌ Null values": {"value": null_values, "color": "#ffebee"}
    }
    
    # Add strand metrics if available
    if strand_stats:
        stats.update({
            "➕ + strand": {"value": f"{strand_stats.get('+', 0):.1%}", "color": "#e8f5e9"},
            "➖ - strand": {"value": f"{strand_stats.get('-', 0):.1%}", "color": "#ffebee"}
        })
    
    # Create the metrics banner
    st.write("""
    <div style='background: #f8f9fa;
                border-radius: 8px;
                padding: 15px;
                margin: 15px 0;
                border-left: 4px solid #2e86ab;
                display: flex;
                justify-content: space-between;
                align-items: center;
                flex-wrap: wrap;
                gap: 12px;'>
    """, unsafe_allow_html=True)
    
    # Create columns dynamically based on stats count
    cols = st.columns(len(stats))
    for i, (name, data) in enumerate(stats.items()):
        with cols[i]:
            st.write(f"""
            <div style='text-align: center;
                        min-width: 90px;
                        padding: 8px;
                        border-radius: 8px;
                        background: {data["color"]};
                        border: 1px solid #e0e0e0;
                        box-shadow: 0 2px 4px rgba(0,0,0,0.05);'>
                <div style='font-size: 12px; color: #455a64; margin-bottom: 5px; font-weight: 500;'>{name}</div>
                <div style='font-size: 14px; font-weight: 600; color: #263238;'>{data["value"]}</div>
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
    """Display data quality statistics with colored backgrounds"""
    st.subheader("🧬 Distinct Values")
    col1, col2, col3, col4 = st.columns(4)
    
    # Define colors for metrics
    colors = ["#e8f5e9", "#e3f2fd", "#fff3e0", "#e0f2f1"]
    
    with col1:
        st.markdown(f"""
        <div style='background-color: {colors[0]}; 
                    padding: 15px; 
                    border-radius: 10px;
                    text-align: center;
                    border: 1px solid #e0e0e0;'>
            <h3 style='margin: 0; color: #455a64;'>Chromosomes</h3>
            <p style='font-size: 24px; font-weight: bold; margin: 5px 0; color: #263238;'>{df["chrom"].nunique()}</p>
        </div>
        """, unsafe_allow_html=True)
    
    with col2:
        total_size = (df['end'] - df['start']).sum()
        st.markdown(f"""
        <div style='background-color: {colors[1]}; 
                    padding: 15px; 
                    border-radius: 10px;
                    text-align: center;
                    border: 1px solid #e0e0e0;'>
            <h3 style='margin: 0; color: #455a64;'>Total size</h3>
            <p style='font-size: 24px; font-weight: bold; margin: 5px 0; color: #263238;'>{f"{total_size:,} bp"}</p>
        </div>
        """, unsafe_allow_html=True)
    
    with col3:
        null_count = df.isnull().sum().sum()
        st.markdown(f"""
        <div style='background-color: {colors[2]}; 
                    padding: 15px; 
                    border-radius: 10px;
                    text-align: center;
                    border: 1px solid #e0e0e0;'>
            <h3 style='margin: 0; color: #455a64;'>Null values</h3>
            <p style='font-size: 24px; font-weight: bold; margin: 5px 0; color: #263238;'>{null_count}</p>
        </div>
        """, unsafe_allow_html=True)
    
    with col4:
        st.markdown(f"""
        <div style='background-color: {colors[3]}; 
                    padding: 15px; 
                    border-radius: 10px;
                    text-align: center;
                    border: 1px solid #e0e0e0;'>
            <h3 style='margin: 0; color: #455a64;'>Columns</h3>
            <p style='font-size: 24px; font-weight: bold; margin: 5px 0; color: #263238;'>{len(df.columns)}</p>
        </div>
        """, unsafe_allow_html=True)
    
    # Missing values analysis by column
    st.subheader("🔍 Null Values Detail")
    null_details = df.isnull().sum().to_frame("Null Values")
    null_details["Percentage"] = (null_details["Null Values"] / len(df) * 100).round(2)
    st.table(null_details.style.background_gradient(cmap="Reds", subset=["Null Values"]))
    
    # Data types by column
    st.subheader("📊 Data Types")
    dtype_details = pd.DataFrame({
        'Column': df.columns,
        'Type': [str(dtype) for dtype in df.dtypes],
        'Unique Values': [df[col].nunique() for col in df.columns]
    })
    st.table(dtype_details.style.background_gradient(cmap="Blues", subset=["Unique Values"]))
    
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

# Example usage (if running as main)
if __name__ == "__main__":
    st.title("BED File Analyzer")
    
    uploaded_file = st.file_uploader("Upload a BED file", type=["bed", "txt"])
    
    if uploaded_file is not None:
        df = load_bed(uploaded_file)
        if df is not None:
            show_metrics_banner(df)
            show_data_preview(df)