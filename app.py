import streamlit as st
from config import setup_page_config, display_header
from preview import load_bed, show_metrics_banner
from mergeBed import handle_merge_operation
from intersectBed import load_bed_int, intersect_bedtools
from sort import sort_bed, handle_sort_operation
import pandas as pd
import gzip
import io

# 👉 Déplace setup_page_config en tout début de script
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
    """Lire les fichiers compressés .gz"""
    try:
        if uploaded_file.name.endswith('.gz'):
            # Lire le fichier compressé
            with gzip.open(uploaded_file, 'rt') as f:
                content = f.read()
            # Créer un objet StringIO pour pandas
            file_obj = io.StringIO(content)
            return file_obj
        else:
            # Lire le fichier normal
            content = uploaded_file.getvalue().decode('utf-8')
            file_obj = io.StringIO(content)
            return file_obj
    except Exception as e:
        st.error(f"Erreur lecture fichier compressé: {str(e)}")
        return None

def main():
    # Interface setup
    setup_styles()
    display_header()
    
    # Initialiser les états de session pour garder le contexte
    if 'intersect_mode' not in st.session_state:
        st.session_state.intersect_mode = False
    if 'file_a_uploaded' not in st.session_state:
        st.session_state.file_a_uploaded = None
    
    # File upload (fichier principal) - AJOUT .gz
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
            # Activer le mode intersection
            st.session_state.intersect_mode = True
            
        if sort_btn:
            show_metrics_banner(df)
            handle_sort_operation(df)
    
    # Section pour le deuxième fichier (uniquement en mode intersection)
    if st.session_state.intersect_mode:
        st.write("---")
        st.subheader("⚡ Intersection BED")
        st.info("Recherche des régions qui se chevauchent entre les deux fichiers")
        
        # Uploader le deuxième fichier dans une section dédiée - AJOUT .gz
        uploaded_file_a = st.file_uploader(
            "Sélectionnez le deuxième fichier BED", 
            type=["bed", "bed.gz", "txt", "txt.gz"], 
            key="fileA_intersect",
            help="Supported formats: .bed, .bed.gz, .txt, .txt.gz"
        )
        
        # Bouton pour exécuter l'intersection
        run_intersect_btn = st.button("🚀 Lancer l'intersection", 
                                    disabled=uploaded_file_a is None,
                                    type="primary",
                                    use_container_width=True)
        
        if uploaded_file_a is not None:
            # Stocker le fichier dans l'état de session
            st.session_state.file_a_uploaded = uploaded_file_a
            
            # Charger les deux fichiers
            if uploaded_file is not None and uploaded_file_a is not None:
                df1, df2 = load_bed_int(uploaded_file, uploaded_file_a)
                
                if df1 is not None and df2 is not None:
                    st.success("✅ Fichiers chargés avec succès")
                    
                    # Afficher les infos de base
                    col_info1, col_info2 = st.columns(2)
                    with col_info1:
                        st.write(f"**Fichier A:** {len(df1)} régions")
                    with col_info2:
                        st.write(f"**Fichier B:** {len(df2)} régions")
                    
                    # Exécuter l'intersection quand le bouton est cliqué
                    if run_intersect_btn:
                        with st.spinner("Recherche des chevauchements..."):
                            # Exécuter l'intersection SIMPLIFIÉE (sans options)
                            result_df = intersect_bedtools(df1, df2)
                            
                            if result_df is not None:
                                if len(result_df) > 0:
                                    st.success(f"✅ {len(result_df)} chevauchements trouvés !")
                                    
                                    # Affichage des résultats
                                    st.subheader("📊 Résultats de l'intersection")
                                    st.dataframe(result_df.head(20))
                                    
                                    # Métriques
                                    col_met1, col_met2, col_met3 = st.columns(3)
                                    with col_met1:
                                        st.metric("Régions dans A", len(df1))
                                    with col_met2:
                                        st.metric("Régions dans B", len(df2))
                                    with col_met3:
                                        st.metric("Chevauchements", len(result_df))
                                    
                                    # Téléchargement des résultats
                                    csv = result_df.to_csv(index=False, sep='\t')
                                    st.download_button(
                                        label="📥 Télécharger les résultats",
                                        data=csv,
                                        file_name="intersection_results.bed",
                                        mime="text/plain",
                                        use_container_width=True
                                    )
                                else:
                                    st.warning("⚠️ Aucun chevauchement trouvé entre les fichiers")
        else:
            if st.session_state.file_a_uploaded is None:
                st.warning("⚠️ Merci d'uploader le deuxième fichier pour continuer.")
    
    # Bouton pour quitter le mode intersection
    if st.session_state.intersect_mode:
        if st.button("❌ Quitter le mode Intersection"):
            st.session_state.intersect_mode = False
            st.session_state.file_a_uploaded = None
            st.rerun()

if __name__ == "__main__":
    main()