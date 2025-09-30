import streamlit as st
from config import setup_page_config, display_header
from preview import load_bed, show_metrics_banner
from mergeBed import handle_merge_operation
from intersectBed import load_bed_int, intersect_bedtools
from sort import sort_bed, handle_sort_operation
import pandas as pd
import time

setup_page_config()

def detect_bed_problems(df):
    """Detect common BED file problems"""
    problems = {
        'chrUn': 0,
        'random': 0,
        'alt': 0,
        'scaffold': 0,
        'inverted_coords': 0,
        'other_bizarre': 0
    }
    
    # Check chromosome names
    chrom_col = df.iloc[:, 0].astype(str)
    
    problems['chrUn'] = len(chrom_col[chrom_col.str.contains('chrUn', na=False)])
    problems['random'] = len(chrom_col[chrom_col.str.contains('random', na=False)])
    problems['alt'] = len(chrom_col[chrom_col.str.contains('_alt', na=False)])
    problems['scaffold'] = len(chrom_col[chrom_col.str.contains('scaffold', na=False)])
    
    # Other bizarre patterns
    bizarre_patterns = ['_hap', '_CTG', '_v', 'KI', 'GL', 'JH', 'NT', 'NW_']
    for pattern in bizarre_patterns:
        problems['other_bizarre'] += len(chrom_col[chrom_col.str.contains(pattern, na=False)])
    
    # Check coordinates
    try:
        problems['inverted_coords'] = len(df[df.iloc[:, 1] >= df.iloc[:, 2]])
    except:
        problems['inverted_coords'] = 0
    
    problems['total'] = sum(problems.values()) - problems['total'] if 'total' in problems else sum(problems.values())
    
    return problems

def clean_bed_file(df):
    """Clean BED file by removing problematic regions"""
    # Filter out problematic chromosomes
    chrom_col = df.iloc[:, 0].astype(str)
    
    # Patterns to remove
    patterns_to_remove = [
        'chrUn', 'random', '_alt', 'scaffold', 
        '_hap', '_CTG', '_v', 'KI', 'GL', 'JH', 'NT', 'NW_'
    ]
    
    mask = ~chrom_col.str.contains('|'.join(patterns_to_remove), na=False)
    df_clean = df[mask].copy()
    
    # Fix inverted coordinates
    try:
        coord_mask = df_clean.iloc[:, 1] < df_clean.iloc[:, 2]
        df_clean = df_clean[coord_mask]
    except:
        pass
    
    return df_clean

def validate_bed_format(df):
    """Validate only BED format - check basic structure"""
    errors = []
    
    # Basic format check - minimum 3 columns
    if df.shape[1] < 3:
        errors.append("BED format requires at least 3 columns (chrom, start, end)")
        return {"is_valid": False, "errors": errors}
    
    # Check if coordinates are numeric and start < end
    try:
        # Test if first few rows have valid coordinates
        test_df = df.head(100)  # Test only first 100 rows for performance
        invalid_coords = test_df[
            (test_df.iloc[:, 1].astype(float) >= test_df.iloc[:, 2].astype(float))
        ]
        if len(invalid_coords) > 0:
            errors.append(f"Found {len(invalid_coords)} regions with start ≥ end coordinates")
            return {"is_valid": False, "errors": errors}
    except (ValueError, TypeError):
        errors.append("Coordinates must be numeric values")
        return {"is_valid": False, "errors": errors}
    
    return {"is_valid": True, "errors": errors}

def setup_styles():
    st.write("""
    <style>
        .stButton>button { 
            padding: 0.25rem 0.5rem !important; 
            font-size: 0.9em !important;
            height: auto !important;
        }
        .validation-frame {
            border: 1px solid #e0e0e0;
            border-radius: 8px;
            padding: 15px;
            margin: 10px 0;
            background-color: #fafafa;
        }
        .problem-item {
            margin: 5px 0;
            padding: 5px;
            border-radius: 3px;
        }
        .problems-section {
            background-color: #fff3cd;
            border: 1px solid #ffeaa7;
            border-radius: 5px;
            padding: 10px;
            margin: 10px 0;
        }
        .compact-warning {
            background-color: #fff3cd;
            border: 1px solid #ffeaa7;
            border-radius: 5px;
            padding: 8px 12px;
            margin: 10px 0;
            font-size: 0.9em;
        }
    </style>
    """, unsafe_allow_html=True)

def main():
    setup_styles()
    display_header()
    
    # Initialize session states
    if 'current_mode' not in st.session_state:
        st.session_state.current_mode = None
    if 'file_loaded' not in st.session_state:
        st.session_state.file_loaded = False
    if 'current_file' not in st.session_state:
        st.session_state.current_file = None
    if 'file_a_uploaded' not in st.session_state:
        st.session_state.file_a_uploaded = None
    if 'cleaned_df' not in st.session_state:
        st.session_state.cleaned_df = None
    if 'detected_problems' not in st.session_state:
        st.session_state.detected_problems = None
    
    # File upload
    uploaded_file = st.file_uploader("Upload BED file", type=["bed"])
    
    df = None
    if uploaded_file is not None:
        if st.session_state.current_file != uploaded_file.name:
            st.session_state.file_loaded = False
            st.session_state.cleaned_df = None
            st.session_state.detected_problems = None
            
        with st.spinner("Validating BED format..."):
            df = load_bed(uploaded_file)
        
        if df is not None:
            # Validate only the format
            validation_result = validate_bed_format(df)
            
            if not st.session_state.file_loaded:
                # Display simple validation result
                with st.container():
                    st.markdown('<div class="validation-frame">', unsafe_allow_html=True)
                    
                    if validation_result["is_valid"]:
                        st.success("✅ **VALID BED FILE** - Format is correct")
                        
                        # Stocker les problèmes pour la preview
                        problems = detect_bed_problems(df)
                        st.session_state.detected_problems = problems
                            
                    else:
                        st.error("❌ **INVALID BED FILE** - Format errors detected")
                        for error in validation_result["errors"]:
                            st.error(f"• {error}")
                    
                    st.session_state.cleaned_df = df
                    st.session_state.file_loaded = True
                    st.session_state.current_file = uploaded_file.name
                    st.markdown('</div>', unsafe_allow_html=True)
                
        else:
            st.error("❌ Failed to load file")
            st.session_state.file_loaded = False
    
    # Use cleaned dataframe if available
    current_df = st.session_state.cleaned_df if st.session_state.cleaned_df is not None else df
    
    # THREE SMALL BUTTONS AT THE BOTTOM
    st.write("---")
    st.markdown("### 🔧 Tools")
    
    col1, col2, col3 = st.columns(3)
    
    with col1:
        preview_btn = st.button("🔍 Preview", 
                              disabled=current_df is None,
                              help="Preview file",
                              use_container_width=True)

    with col2:
        dataviz_btn = st.button("📊 Viz", 
                              disabled=current_df is None,
                              help="Visualize data",
                              use_container_width=True)
    
    with col3:
        bedtools_btn = st.button("⚙️ Tools", 
                               disabled=current_df is None,
                               help="BED operations",
                               use_container_width=True)
    
    # Handle main mode selection
    if preview_btn:
        st.session_state.current_mode = 'preview'
        st.rerun()
        
    if dataviz_btn:
        st.session_state.current_mode = 'dataviz'
        st.rerun()
        
    if bedtools_btn:
        st.session_state.current_mode = 'bedtools'
        st.rerun()
    
    # PREVIEW MODE
    if st.session_state.current_mode == 'preview' and current_df is not None:
        st.write("---")
        st.subheader("🔍 Preview")
        
        # Afficher les métriques KPI
        show_metrics_banner(current_df)
        
        # Afficher les problèmes et bouton de nettoyage SOUS les KPI
        if st.session_state.detected_problems and any(st.session_state.detected_problems.values()):
            problems = st.session_state.detected_problems
            
            # Créer la liste des problèmes sur une ligne
            problem_list = []
            if problems['chrUn'] > 0:
                problem_list.append(f"chrUn:{problems['chrUn']}")
            if problems['random'] > 0:
                problem_list.append(f"random:{problems['random']}")
            if problems['alt'] > 0:
                problem_list.append(f"alt:{problems['alt']}")
            if problems['scaffold'] > 0:
                problem_list.append(f"scaffold:{problems['scaffold']}")
            if problems['other_bizarre'] > 0:
                problem_list.append(f"bizarre:{problems['other_bizarre']}")
            if problems['inverted_coords'] > 0:
                problem_list.append(f"inverted:{problems['inverted_coords']}")
            
            if problem_list:
                problems_text = " ⚠️ Quality issues detected: " + " | ".join(problem_list)
                
                st.markdown(f'<div class="compact-warning">{problems_text}</div>', unsafe_allow_html=True)
                
                # Bouton de nettoyage
                col1, col2 = st.columns([3, 1])
                with col1:
                    st.info("Clean the file to remove problematic regions and improve analysis quality")
                with col2:
                    if st.button("🧹 Clean File", type="primary", use_container_width=True):
                        st.session_state.cleaned_df = clean_bed_file(current_df)
                        st.session_state.detected_problems = None  # Reset après nettoyage
                        st.rerun()
        
        # Affichage du dataframe
        st.dataframe(current_df.head(50), height=300)
        
        if st.button("← Back", key="back_preview"):
            st.session_state.current_mode = None
            st.rerun()
    
    # DATAVIZ MODE
    elif st.session_state.current_mode == 'dataviz' and current_df is not None:
        st.write("---")
        st.subheader("📊 Visualization")
        
        # Compact metrics
        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric("Regions", len(current_df))
        with col2:
            st.metric("Chromosomes", current_df.iloc[:, 0].nunique())
        with col3:
            try:
                avg_len = (current_df.iloc[:, 2] - current_df.iloc[:, 1]).mean()
                st.metric("Avg Length", f"{avg_len:.0f}bp")
            except:
                st.metric("Avg Length", "N/A")
        
        st.info("Visualization features coming soon...")
        
        if st.button("← Back", key="back_dataviz"):
            st.session_state.current_mode = None
            st.rerun()
    
    # BEDTOOLS MODE
    elif st.session_state.current_mode == 'bedtools' and current_df is not None:
        st.write("---")
        st.subheader("⚙️ BED Operations")
        
        # Small sub-buttons
        col1, col2, col3 = st.columns(3)
        
        with col1:
            merge_btn = st.button("🔄 Merge", 
                                help="Merge intervals",
                                use_container_width=True)
        
        with col2:
            sort_btn = st.button("🔢 Sort", 
                               help="Sort regions",
                               use_container_width=True)
        
        with col3:
            intersect_btn = st.button("⚡ Intersect", 
                                    help="Find overlaps",
                                    use_container_width=True)
        
        # Handle operations
        if merge_btn:
            st.write("---")
            st.subheader("🔄 Merge")
            handle_merge_operation(current_df)
            
        if sort_btn:
            st.write("---")
            st.subheader("🔢 Sort")
            show_metrics_banner(current_df)
            handle_sort_operation(current_df)
            
        if intersect_btn:
            st.session_state.intersect_mode = True
            st.rerun()
        
        # Intersection sub-section
        if st.session_state.get('intersect_mode', False):
            st.write("---")
            st.subheader("⚡ Intersect")
            
            uploaded_file_a = st.file_uploader(
                "Second BED file", 
                type=["bed"],
                key="fileA_intersect"
            )
            
            run_intersect_btn = st.button("Run", 
                                        disabled=uploaded_file_a is None,
                                        use_container_width=True)
            
            if uploaded_file_a is not None:
                st.session_state.file_a_uploaded = uploaded_file_a
                
                if uploaded_file is not None and uploaded_file_a is not None:
                    df1, df2 = load_bed_int(uploaded_file, uploaded_file_a)
                    
                    if df1 is not None and df2 is not None:
                        col1, col2 = st.columns(2)
                        with col1:
                            st.metric("File A", f"{len(df1):,}")
                        with col2:
                            st.metric("File B", f"{len(df2):,}")
                        
                        if run_intersect_btn:
                            with st.spinner("Processing..."):
                                result_df = intersect_bedtools(df1, df2)
                                
                                if result_df is not None:
                                    if len(result_df) > 0:
                                        st.success(f"✅ {len(result_df):,} overlaps")
                                        
                                        with st.expander("Results"):
                                            st.dataframe(result_df.head(20))
                                        
                                        csv = result_df.to_csv(index=False, sep='\t')
                                        st.download_button(
                                            label="📥 Download",
                                            data=csv,
                                            file_name="intersection.bed",
                                            mime="text/plain",
                                            use_container_width=True
                                        )
                                    else:
                                        st.warning("⚠️ No overlaps")
            
            if st.button("← Back", key="back_intersect"):
                st.session_state.intersect_mode = False
                st.session_state.file_a_uploaded = None
                st.rerun()
        
        # Back button
        if st.button("← Back", key="back_bedtools"):
            st.session_state.current_mode = None
            st.session_state.intersect_mode = False
            st.session_state.file_a_uploaded = None
            st.rerun()

if __name__ == "__main__":
    main()