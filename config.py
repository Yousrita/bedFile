import streamlit as st
from pathlib import Path
import pandas as pd
from typing import Optional

def setup_page_config():
    """Configure sophisticated page settings for bioinformatics analysis"""
    st.set_page_config(
        page_title="Genomic Interval Toolkit",
        page_icon="🧬",
        layout="wide",
        initial_sidebar_state="expanded",
        menu_items={
            'Get Help': 'https://example.com',
            'Report a bug': "https://example.com",
            'About': "# BED File Processing Suite"
        }
    )

def setup_styles(file_uploaded=False):
    """Boutons avec couleur de fond une fois le fichier uploadé"""
    colors = {
        "primary": "#97999B",
        "secondary": "#CF5C78", 
        "background": "#F9FAFB",
        "accent": "#D1D5DB",
        "dark": "#679436",
        "light": "#F3F4F6",
        "preview": "#4CAF50",
        "sort": "#2196F3",
        "merge": "#FF9800",
        "intersect": "#9C27B0"
    }

    button_style = f"""
        <style>
        .stButton > button {{
            font-family: -apple-system, BlinkMacSystemFont, sans-serif;
            font-weight: 400;
            font-size: 13px;
            letter-spacing: 0.3px;
            border-radius: 6px;
            padding: 6px 12px;
            transition: all 0.2s ease;
            width: auto;
            min-width: 100px;
            margin: 2px 0;
            border: 1px solid {colors['accent']};
            box-shadow: none;
        }}
        
        .stButton > button:disabled {{
            background-color: {colors['background']} !important;
            color: {colors['secondary']} !important;
            border-color: {colors['light']} !important;
        }}
    """

    if file_uploaded:
        button_style += f"""
        .stButton > button:nth-child(1) {{
            background-color: {colors['preview']} !important;
            color: white !important;
            border-color: {colors['preview']} !important;
        }}
        
        .stButton > button:nth-child(2) {{
            background-color: {colors['sort']} !important;
            color: white !important;
            border-color: {colors['sort']} !important;
        }}
        
        .stButton > button:nth-child(3) {{
            background-color: {colors['merge']} !important;
            color: white !important;
            border-color: {colors['merge']} !important;
        }}
        
        .stButton > button:nth-child(4) {{
            background-color: {colors['intersect']} !important;
            color: white !important;
            border-color: {colors['intersect']} !important;
        }}
        
        .stButton > button:hover {{
            opacity: 0.9;
            transform: translateY(-1px);
        }}
        """
    else:
        button_style += f"""
        .stButton > button {{
            background-color: white;
            color: {colors['dark']};
        }}
        """

    st.write(button_style, unsafe_allow_html=True)

def display_header():
    """Display a premium header section with professional landing"""
    col1, col2 = st.columns([1, 3])

    with col1:
        st.write("""
        <div style='text-align: center; margin-bottom: 20px;'>
            <div style='font-size: 72px;'>🧬</div>
            <div style='font-size: 14px; color: #666;'>Bed Toolkit v1.0</div>
        </div>
        """, unsafe_allow_html=True) 

    with col2:
        st.title("BED File Explorer")
        st.write("""
        <div style='font-size: 16px; color: #555; margin-bottom: 10px;'>
            Exploring and editing genomic data
        </div>
        <div style='color: #777; font-size: 14px;'>
            <span style='background: #F0F0F0; padding: 3px 8px; border-radius: 4px; margin-right: 8px;'>BED</span>
        </div>
        """, unsafe_allow_html=True)

    st.write("---")

def display_action_buttons(disabled=True):
    """Create premium action buttons with enhanced UI"""
    st.write("""
    <div style='margin: 30px 0 20px 0;'>
        <h3 style='color: #444; font-weight: 500;'>Core Operations</h3>
    </div>
    """, unsafe_allow_html=True)

    cols = st.columns(4)
    with cols[0]:
        preview_btn = st.button("**🔍 Preview**",
                                disabled=disabled,
                                help="Comprehensive analysis of BED file contents")
    with cols[1]:
        sort_btn = st.button("**↕️ Sort**",
                             disabled=disabled,
                             help="Sort intervals by chromosome and position")
    with cols[2]:
        merge_btn = st.button("**🧬 Merge**",
                              disabled=disabled,
                              help="Advanced merging of genomic regions with configurable parameters")
    with cols[3]:
        intersect_btn = st.button("**🔀 Intersect**",
                                  disabled=disabled,
                                  help="Find intersections between genomic intervals")

    st.write("---")
    return preview_btn, sort_btn, merge_btn, intersect_btn

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
            dtype=str
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
            
        return df
        
    except Exception as e:
        st.error(f"❌ Loading error: {str(e)}")
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
    st.subheader("📋 Data Preview")
    st.dataframe(df.head(100), height=300)
    st.caption(f"Displaying 100 rows out of {len(df):,} total - {len(df.columns)} columns")
    
    # Display column list
    st.write("**📝 Available columns:**")
    for i, col in enumerate(df.columns):
        st.write(f"- `{col}` ({df[col].dtype})")

def show_data_quality(df: pd.DataFrame) -> None:
    """Display data quality statistics"""
    st.subheader("🧪 Data Quality")
    
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

def main():
    """Main application function"""
    # Configuration MUST be first
    setup_page_config()
    display_header()
    
    # File upload section
    uploaded_file = st.file_uploader("Upload a BED file", type=["bed", "txt"])
    
    # Initialize session state for file upload status
    if 'file_uploaded' not in st.session_state:
        st.session_state.file_uploaded = False
    if 'show_success_message' not in st.session_state:
        st.session_state.show_success_message = False
    
    if uploaded_file is not None:
        # Show loading state
        with st.spinner("🔄 Loading and analyzing your BED file..."):
            df = load_bed(uploaded_file)
        
        if df is not None:
            # ✅ AFFICHER LE MESSAGE DE SUCCÈS IMMÉDIATEMENT APRÈS LE CHARGEMENT
            # Utiliser un état de session pour éviter la ré-exécution
            if not st.session_state.show_success_message or st.session_state.current_file != uploaded_file.name:
                st.success("✅ **File loaded successfully!**")
                st.toast(f"📊 {uploaded_file.name} loaded with {len(df):,} regions", icon="✅")
                st.session_state.show_success_message = True
                st.session_state.current_file = uploaded_file.name
            
            st.session_state.file_uploaded = True
            
            # Apply styled buttons
            setup_styles(file_uploaded=True)
            
            # Display metrics banner
            show_metrics_banner(df)
            
            # Display action buttons
            preview_btn, sort_btn, merge_btn, intersect_btn = display_action_buttons(disabled=False)
            
            # Handle button actions
            if preview_btn:
                show_data_preview(df)
                
            if sort_btn:
                st.subheader("↕️ Sorted Data")
                sorted_df = df.sort_values(by=['chrom', 'start'])
                st.dataframe(sorted_df.head(100), height=300)
                st.success("Data sorted by chromosome and start position")
                
            if merge_btn:
                st.subheader("🧬 Merge Intervals")
                st.info("Merge functionality would be implemented here")
                
            if intersect_btn:
                st.subheader("🔀 Intersect Intervals")
                st.info("Intersect functionality would be implemented here")
                
            # Afficher l'onglet de qualité de données par défaut
            if not any([preview_btn, sort_btn, merge_btn, intersect_btn]):
                show_data_quality(df)
        else:
            st.error("❌ Failed to load the BED file. Please check the file format.")
            st.session_state.show_success_message = False
    else:
        # Reset states when no file is uploaded
        st.session_state.file_uploaded = False
        st.session_state.show_success_message = False
        st.session_state.current_file = None
        
        setup_styles(file_uploaded=False)
        display_action_buttons(disabled=True)
        st.info("👆 Please upload a BED file to begin analysis")

if __name__ == "__main__":
    main()