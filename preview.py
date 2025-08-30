import pandas as pd
import streamlit as st
from typing import Optional

def load_bed(uploaded_file) -> Optional[pd.DataFrame]:
    """
    Charge et valide un fichier BED uploadé
    Retourne un DataFrame avec seulement les colonnes chrom/start/end
    """
    try:
        df = pd.read_csv(
            uploaded_file,
            sep="\t",
            header=None,
            comment='#',
            dtype={0: str, 1: int, 2: int}
        )
        
        if len(df.columns) < 3:
            st.error("❌ Format BED invalide : minimum 3 colonnes requises")
            return None
            
        # Ne garder que les 3 premières colonnes
        df = df.iloc[:, :3].copy()
        df.columns = ["chrom", "start", "end"]
        
        if (df["start"] > df["end"]).any():
            st.error("❌ Erreur : des positions start > end détectées")
            return None
            
        return df
    except Exception as e:
        st.error(f"❌ Erreur de chargement : {str(e)}")
        return None

def show_metrics_banner(df):
    """Display a comprehensive single-line file summary"""
    if df is None or df.empty:
        return
        
    # Calculate all metrics
    length = df['end'] - df['start']
    duplicates = df.duplicated().sum()
    null_values = df.isnull().sum().sum()
    
    # Check if strand column exists
    strand_stats = df['strand'].value_counts(normalize=True).to_dict() if 'strand' in df else None
    
    stats = {
        "🧬 Chromosomes": df['chrom'].nunique(),
        "📊 Intervals": f"{len(df):,}",
        "📏 Avg length": f"{length.mean():.0f} bp",
        "🧩 Total span": f"{length.sum()/1e6:.2f} Mb",
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
                        background: rgba(255,255,255,0.7);'>
                <div style='font-size: 11px; color: #555;'>{name}</div>
                <div style='font-size: 14px; font-weight: 600;'>{value}</div>
            </div>
            """, unsafe_allow_html=True)
    
    st.write("</div>", unsafe_allow_html=True)

def show_data_preview(df: pd.DataFrame) -> None:
    """Affiche l'aperçu des données avec onglets"""
    tab1, tab2 = st.tabs(["📋 Aperçu des données", "🧪 Qualité des données"])
    
    with tab1:
        st.dataframe(df.head(100), height=300)
    
    with tab2:
        show_data_quality(df)

def show_data_quality(df: pd.DataFrame) -> None:
    """Affiche les statistiques de qualité des données"""
    st.subheader("🧬 Valeurs Distinctes")
    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.metric("Chromosomes uniques", df["chrom"].nunique())
    
    with col2:
        st.metric("Taille totale (bp)", f"{(df['end'] - df['start']).sum():,}")
    
    with col3:
        st.metric("Valeurs nulles", df.isnull().sum().sum())
    
    # Analyse des valeurs manquantes par colonne
    st.subheader("🔍 Détail des valeurs nulles")
    null_details = df.isnull().sum().to_frame("Count")
    st.table(null_details)
    
    # Analyse des brins si colonne présente
    if 'strand' in df:
        st.subheader("⚖ Distribution des brins")
        strand_dist = df['strand'].value_counts(normalize=True)
        st.bar_chart(strand_dist)