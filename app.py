import streamlit as st
from config import setup_page_config, display_header
from preview import load_bed, show_metrics_banner
from mergeBed import handle_merge_operation
from intersectBed import load_bed_int, intersect_bedtools
from sort import sort_bed, handle_sort_operation
import pandas as pd
import gzip
import io

# 👉 Move setup_page_config to the very beginning of the script
setup_page_config()

def setup_styles():
    """Custom CSS for button colors"""
    st.write("""
    <style>
        /* Inactive buttons */
        .stButton>button:disabled {
            background-color: #f0f2f6 !important;
            color: #9da5b1 !important;
            border: 1px solid #e2e5e9 !important;
        }
        
        /* Active green buttons */
        .stButton>button {
            background-color: #2e8b57 !important;
            color: white !important;
            border: none !important;
            transition: all 0.3s !important;
        }
        
        .stButton>button:hover {
            background-color: #3cb371 !important;
            transform: translateY(-1px);
        }
    </style>
    """, unsafe_allow_html=True)

def read_compressed_file(uploaded_file):
    """Read compressed .gz files"""
    try:
        if uploaded_file.name.endswith('.gz'):
            # Read the compressed file
            with gzip.open(uploaded_file, 'rt') as f:
                content = f.read()
            # Create a StringIO object for pandas
            file_obj = io.StringIO(content)
            return file_obj
        else:
            # Read the normal file
            content = uploaded_file.getvalue().decode('utf-8')
            file_obj = io.StringIO(content)
            return file_obj
    except Exception as e:
        st.error(f"Error reading compressed file: {str(e)}")
        return None

def main():
    # Interface setup
    setup_styles()
    display_header()
    
    # Initialize session states to keep context
    if 'intersect_mode' not in st.session_state:
        st.session_state.intersect_mode = False
    if 'file_a_uploaded' not in st.session_state:
        st.session_state.file_a_uploaded = None
    
    # File upload (main file) - ADDED .gz
    uploaded_file = st.file_uploader(
        "Upload BED file",
        type=["bed", "bed.gz", "txt", "txt.gz"],
        help="Supported formats: .bed, .bed.gz, .txt, .txt.gz"
    )
          
    # Initialize df
    df = None
    if uploaded_file is not None:
        df = load_bed(uploaded_file)
    
    # Action buttons
    cols = st.columns(4)
    with cols[0]:
        preview_btn = st.button("🔍 Preview", 
                              disabled=df is None,
                              help="Preview file contents")

    with cols[1]:
        sort_btn = st.button("🔢 Sort", 
                           disabled=df is None,
                           help="Sort genomic regions")    
    with cols[2]:
        merge_btn = st.button("🔄 Merge", 
                            disabled=df is None,
                            help="Merge overlapping intervals")
    with cols[3]:
        intersect_btn = st.button("⚡ Intersect", 
                                disabled=df is None,
                                help="Find intersecting regions")
    
    # Handle actions
    if df is not None:
        if preview_btn:
            show_metrics_banner(df)
            st.dataframe(df.head(100))
            
        if merge_btn:
            handle_merge_operation(df)

        if intersect_btn:
            # Activate intersection mode
            st.session_state.intersect_mode = True
            
        if sort_btn:
            show_metrics_banner(df)
            handle_sort_operation(df)
    
    # Section for the second file (only in intersection mode)
    if st.session_state.intersect_mode:
        st.write("---")
        st.subheader("⚡ BED Intersection")
        st.info("Search for overlapping regions between the two files")
        
        # Upload the second file in a dedicated section - ADDED .gz
        uploaded_file_a = st.file_uploader(
            "Select the second BED file", 
            type=["bed", "bed.gz", "txt", "txt.gz"], 
            key="fileA_intersect",
            help="Supported formats: .bed"
        )
        
        # Button to execute the intersection
        run_intersect_btn = st.button("🚀 Run Intersection", 
                                    disabled=uploaded_file_a is None,
                                    type="primary",
                                    use_container_width=True)
        
        if uploaded_file_a is not None:
            # Store the file in session state
            st.session_state.file_a_uploaded = uploaded_file_a
            
            # Load both files
            if uploaded_file is not None and uploaded_file_a is not None:
                df1, df2 = load_bed_int(uploaded_file, uploaded_file_a)
                
                if df1 is not None and df2 is not None:
                    st.success("✅ Files loaded successfully")
                    
                    # Display basic info
                    col_info1, col_info2 = st.columns(2)
                    with col_info1:
                        st.write(f"**File A:** {len(df1)} regions")
                    with col_info2:
                        st.write(f"**File B:** {len(df2)} regions")
                    
                    # Execute intersection when button is clicked
                    if run_intersect_btn:
                        with st.spinner("Searching for overlaps..."):
                            # Execute SIMPLIFIED intersection (without options)
                            result_df = intersect_bedtools(df1, df2)
                            
                            if result_df is not None:
                                if len(result_df) > 0:
                                    st.success(f"✅ {len(result_df)} overlaps found!")
                                    
                                    # Results display
                                    st.subheader("📊 Intersection Results")
                                    st.dataframe(result_df.head(20))
                                    
                                    # Metrics
                                    col_met1, col_met2, col_met3 = st.columns(3)
                                    with col_met1:
                                        st.metric("Regions in A", len(df1))
                                    with col_met2:
                                        st.metric("Regions in B", len(df2))
                                    with col_met3:
                                        st.metric("Overlaps", len(result_df))
                                    
                                    # Download results
                                    csv = result_df.to_csv(index=False, sep='\t')
                                    st.download_button(
                                        label="📥 Download Results",
                                        data=csv,
                                        file_name="intersection_results.bed",
                                        mime="text/plain",
                                        use_container_width=True
                                    )
                                else:
                                    st.warning("⚠️ No overlaps found between files")
        else:
            if st.session_state.file_a_uploaded is None:
                st.warning("⚠️ Please upload the second file to continue")
    
    # Button to exit intersection mode
    if st.session_state.intersect_mode:
        if st.button("❌ Exit Intersection Mode"):
            st.session_state.intersect_mode = False
            st.session_state.file_a_uploaded = None
            st.rerun()

if __name__ == "__main__":
    main()