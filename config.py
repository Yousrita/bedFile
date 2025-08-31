import streamlit as st
from pathlib import Path

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
def setup_styles():
    """Boutons discrets et élégants en tons neutres"""
    colors = {
        "primary": "#97999B",  # Gris foncé subtil
        "secondary": "#CF5C78",  # Gris moyen
        "background": "#F9FAFB",  # Fond très clair
        "accent": "#D1D5DB",  # Gris clair
        "dark": "#679436",
        "light": "#F3F4F6"
    }

    st.write(f"""
    <style>
        /* Style minimaliste */
        .stButton>button {{
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
            background-color: white;
            color: {colors['dark']};
            box-shadow: none;
        }}

        /* Bouton principal */
        .stButton>button:first-child {{
            background-color: {colors['dark']};
            color: white;
            border-color: {colors['dark']};
        }}

        /* Bouton secondaire */
        .stButton>button:nth-child(2) {{
            background-color: {colors['background']};
            color: {colors['primary']};
            border-color: {colors['secondary']};
        }}

        /* Bouton tertiaire */
        .stButton>button:nth-child(3) {{
            background-color: white;
            color: {colors['primary']};
            border-color: {colors['accent']};
        }}

        /* Effets au survol */
        .stButton>button:hover {{
            background-color: {colors['light']};
            border-color: {colors['primary']};
        }}

        .stButton>button:first-child:hover {{
            background-color: {colors['primary']};
            border-color: {colors['primary']};
        }}

        .stButton>button:nth-child(2):hover {{
            background-color: {colors['accent']};
        }}

        /* État désactivé */
        .stButton>button:disabled {{
            background-color: {colors['background']} !important;
            color: {colors['secondary']} !important;
            border-color: {colors['light']} !important;
        }}
    </style>
    """, unsafe_allow_html=True)

def display_header():
    """Display a premium header section with professionalanding"""
    col1, col2 = st.columns([1, 3])

    with col1:
        st.write("""
        <div style='text-align: center; margin-bottom: 20px;'>
            <div style='font-size: 72px;'>🧬</div>
            <div style='font-size: 14px; color: #666;'>Bed Toolkit v1.0</div>
        </div>
        """, unsafe_allow_html=True) 

    with col2:
        st.title("BED File Explorer ")
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

    cols = st.columns(3)
    with cols[0]:
        preview_btn = st.button("**🔍 Data Preview**",
                                disabled=disabled,
                                help="Comprehensive analysis of BED file contents")
    with cols[1]:
        merge_btn = st.button("**🧬 Merge Intervals**",
                              disabled=disabled,
                              help="Advanced merging of genomic regions with configurable parameters")
    with cols[2]:
        export_btn = st.button("**📊 Export Analysis**",
                               disabled=disabled,
                               help="Generate professional reports and export in multiple formats")

    st.write("---")
    return preview_btn, merge_btn, export_btn