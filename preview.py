import pandas as pd
import streamlit as st
from typing import Optional

def load_bed(uploaded_file) -> Optional[pd.DataFrame]:
    """
    Charge un fichier BED avec TOUTES ses colonnes
    """
    try:
        # Lire le fichier sans limite de colonnes
        df = pd.read_csv(
            uploaded_file,
            sep="\t",
            header=None,
            comment='#',
            dtype=str  # ← Lire tout en string d'abord
        )
        
        # Définir les noms de colonnes BED standard
        bed_cols = ['chrom', 'start', 'end', 'name', 'score', 'strand', 
                   'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 
                   'blockSizes', 'blockStarts']
        
        # Nommer les colonnes disponibles
        df.columns = bed_cols[:len(df.columns)]
        
        # Convertir start et end en numérique
        if 'start' in df.columns:
            df['start'] = pd.to_numeric(df['start'], errors='coerce')
        if 'end' in df.columns:
            df['end'] = pd.to_numeric(df['end'], errors='coerce')
        
        # Vérifier la validité
        if len(df.columns) < 3:
            st.error("❌ Format BED invalide : minimum 3 colonnes requises")
            return None
            
        if (df["start"] > df["end"]).any():
            st.error("❌ Erreur : des positions start > end détectées")
            return None
            
        st.success(f"✅ Fichier chargé : {len(df)} régions, {len(df.columns)} colonnes")
        return df
        
    except Exception as e:
        st.error(f"❌ Erreur de chargement : {str(e)}")
        return None

def load_bed_basic(uploaded_file) -> Optional[pd.DataFrame]:
    """
    Charge seulement les 3 premières colonnes (pour les opérations bedtools)
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
    
    # Stats améliorées
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
        #"📋 Columns": len(df.columns)
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
    """Affiche l'aperçu des données avec TOUTES les colonnes"""
    tab1, tab2 = st.tabs(["📋 Aperçu des données", "🧪 Qualité des données"])
    
    with tab1:
        # Afficher TOUTES les colonnes
        st.dataframe(df.head(100), height=300)
        st.caption(f"Affichage de 100 lignes sur {len(df):,} totales - {len(df.columns)} colonnes")
        
        # Afficher la liste des colonnes
        st.write("**📝 Colonnes disponibles :**")
        for i, col in enumerate(df.columns):
            st.write(f"- `{col}` ({df[col].dtype})")
    
    with tab2:
        show_data_quality(df)

def show_data_quality(df: pd.DataFrame) -> None:
    """Affiche les statistiques de qualité des données"""
    st.subheader("🧬 Valeurs Distinctes")
    col1, col2, col3, col4 = st.columns(4)
    
    with col1:
        st.metric("Chromosomes", df["chrom"].nunique())
    with col2:
        st.metric("Taille totale", f"{(df['end'] - df['start']).sum():,} bp")
    with col3:
        st.metric("Valeurs nulles", df.isnull().sum().sum())
    with col4:
        st.metric("Colonnes", len(df.columns))
    
    # Analyse des valeurs manquantes par colonne
    st.subheader("🔍 Détail des valeurs nulles")
    null_details = df.isnull().sum().to_frame("Valeurs nulles")
    null_details["Pourcentage"] = (null_details["Valeurs nulles"] / len(df) * 100).round(2)
    st.table(null_details)
    
    # Types de données par colonne
    st.subheader("📊 Types de données")
    dtype_details = pd.DataFrame({
        'Colonne': df.columns,
        'Type': [str(dtype) for dtype in df.dtypes],
        'Valeurs uniques': [df[col].nunique() for col in df.columns]
    })
    st.table(dtype_details)
    
    # Analyse des brins si colonne présente
    if 'strand' in df:
        st.subheader("⚖ Distribution des brins")
        strand_dist = df['strand'].value_counts(normalize=True)
        st.bar_chart(strand_dist)
    
    # Analyse du score si colonne présente
    if 'score' in df:
        st.subheader("📈 Distribution des scores")
        st.write(f"Score moyen : {df['score'].mean():.2f}")
        st.write(f"Score min : {df['score'].min()}")
        st.write(f"Score max : {df['score'].max()}")