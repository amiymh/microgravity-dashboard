"""PCA analysis module for RNA-seq sample-level expression data."""

import pandas as pd
import numpy as np
import plotly.graph_objects as go
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

from modules.data_loader import EARTH_COLS, SPACE_COLS

SAMPLE_COLS = EARTH_COLS + SPACE_COLS


def run_pca(df: pd.DataFrame) -> dict | None:
    """Run PCA on the 16 sample columns (8 Earth + 8 Space).

    Transposes so samples become rows and genes become columns,
    removes zero-variance genes, applies StandardScaler, then fits
    PCA with 3 components.

    Returns:
        dict with keys: pca_df, variance_explained, n_genes_used,
        n_genes_removed, pca_model, feature_names.
        None if input is empty or missing count columns.
    """
    if df is None or df.empty:
        return None

    # Check that at least some sample columns exist
    available = [c for c in SAMPLE_COLS if c in df.columns]
    if len(available) == 0:
        return None

    # Extract count matrix: genes x samples
    count_matrix = df[available].copy()
    n_total_genes = count_matrix.shape[0]

    # Transpose: samples as rows, genes as columns
    transposed = count_matrix.T  # shape: (n_samples, n_genes)

    # Remove genes (columns) with zero variance across samples
    variances = transposed.var(axis=0)
    non_zero_mask = variances > 0
    filtered = transposed.loc[:, non_zero_mask]
    n_genes_used = filtered.shape[1]
    n_genes_removed = n_total_genes - n_genes_used

    if n_genes_used == 0:
        return None

    # StandardScaler normalization
    scaler = StandardScaler()
    scaled = scaler.fit_transform(filtered)

    # PCA with 3 components (or fewer if not enough samples/genes)
    n_components = min(3, scaled.shape[0], scaled.shape[1])
    pca_model = PCA(n_components=n_components)
    components = pca_model.fit_transform(scaled)

    # Build result DataFrame
    conditions = []
    for sample in available:
        if sample in EARTH_COLS:
            conditions.append("Earth")
        else:
            conditions.append("Space")

    pca_df = pd.DataFrame(
        {
            "PC1": components[:, 0],
            "PC2": components[:, 1] if n_components >= 2 else 0.0,
            "PC3": components[:, 2] if n_components >= 3 else 0.0,
            "Sample": available,
            "Condition": conditions,
        }
    )

    variance_explained = pca_model.explained_variance_ratio_.tolist()
    # Pad to length 3 if fewer components were computed
    while len(variance_explained) < 3:
        variance_explained.append(0.0)

    # Gene names that survived filtering
    feature_names = df["Gene"].values[non_zero_mask.values] if "Gene" in df.columns else None

    return {
        "pca_df": pca_df,
        "variance_explained": variance_explained,
        "n_genes_used": n_genes_used,
        "n_genes_removed": n_genes_removed,
        "pca_model": pca_model,
        "feature_names": feature_names,
    }


def get_loadings(
    df: pd.DataFrame,
    pca_model: PCA,
    components_data: np.ndarray | None = None,
    n: int = 20,
) -> pd.DataFrame:
    """Return top N genes contributing to PC1 and PC2.

    Args:
        df: Original DataFrame (used for gene names).
        pca_model: Fitted PCA model.
        components_data: Unused, kept for API compatibility.
        n: Number of top genes to return.

    Returns:
        DataFrame with columns: Gene, PC1_loading, PC2_loading.
    """
    loadings = pca_model.components_  # shape: (n_components, n_features)

    # Identify genes that passed zero-variance filter
    available = [c for c in SAMPLE_COLS if c in df.columns]
    count_matrix = df[available].copy()
    variances = count_matrix.T.var(axis=0)
    non_zero_mask = variances > 0

    genes = df["Gene"].values[non_zero_mask.values] if "Gene" in df.columns else [
        f"Feature_{i}" for i in range(loadings.shape[1])
    ]

    pc1_loadings = loadings[0]
    pc2_loadings = loadings[1] if loadings.shape[0] >= 2 else np.zeros(loadings.shape[1])

    loading_df = pd.DataFrame(
        {
            "Gene": genes[:loadings.shape[1]],
            "PC1_loading": pc1_loadings,
            "PC2_loading": pc2_loadings,
        }
    )

    # Rank by combined absolute loading on PC1 + PC2
    loading_df["_combined"] = loading_df["PC1_loading"].abs() + loading_df["PC2_loading"].abs()
    loading_df = loading_df.nlargest(min(n, len(loading_df)), "_combined")
    loading_df = loading_df.drop(columns=["_combined"]).reset_index(drop=True)

    return loading_df


def create_pca_plot(pca_result: dict) -> go.Figure:
    """Create an interactive PC1 vs PC2 scatter plot.

    Args:
        pca_result: Output dict from run_pca().

    Returns:
        Plotly Figure object.
    """
    pca_df = pca_result["pca_df"]
    var_exp = pca_result["variance_explained"]

    color_map = {"Earth": "#2ecc71", "Space": "#e67e22"}

    fig = go.Figure()

    for condition in ["Earth", "Space"]:
        subset = pca_df[pca_df["Condition"] == condition]
        if subset.empty:
            continue
        fig.add_trace(
            go.Scatter(
                x=subset["PC1"],
                y=subset["PC2"],
                mode="markers+text",
                name=condition,
                marker=dict(color=color_map[condition], size=12, opacity=0.85),
                text=subset["Sample"],
                textposition="top center",
                textfont=dict(size=9),
                hovertemplate=(
                    "<b>%{text}</b><br>"
                    "PC1: %{x:.2f}<br>"
                    "PC2: %{y:.2f}<extra></extra>"
                ),
            )
        )

    fig.update_layout(
        title="PCA — Microgravity vs Ground Samples",
        xaxis_title=f"PC1 ({var_exp[0]*100:.1f}% variance)",
        yaxis_title=f"PC2 ({var_exp[1]*100:.1f}% variance)",
        template="plotly_white",
        legend=dict(orientation="h", yanchor="bottom", y=1.02),
        height=600,
    )

    return fig


def get_analysis_log(pca_result: dict) -> list[str]:
    """Return a human-readable analysis log documenting the PCA run.

    Args:
        pca_result: Output dict from run_pca().

    Returns:
        List of log strings.
    """
    if pca_result is None:
        return ["PCA analysis could not be performed (no valid data)."]

    pca_df = pca_result["pca_df"]
    var_exp = pca_result["variance_explained"]
    n_used = pca_result["n_genes_used"]
    n_removed = pca_result["n_genes_removed"]

    n_samples = len(pca_df)
    earth_count = len(pca_df[pca_df["Condition"] == "Earth"])
    space_count = len(pca_df[pca_df["Condition"] == "Space"])

    log = [
        f"PCA performed on {n_samples} samples ({earth_count} Earth, {space_count} Space).",
        f"Genes used: {n_used} (after removing {n_removed} zero-variance genes).",
        f"Variance explained — PC1: {var_exp[0]*100:.1f}%, PC2: {var_exp[1]*100:.1f}%, PC3: {var_exp[2]*100:.1f}%.",
    ]

    # Top 3 loading genes per component if feature names available
    pca_model = pca_result.get("pca_model")
    feature_names = pca_result.get("feature_names")
    if pca_model is not None and feature_names is not None:
        loadings = pca_model.components_
        for i, pc_label in enumerate(["PC1", "PC2", "PC3"]):
            if i >= loadings.shape[0]:
                break
            abs_loadings = np.abs(loadings[i])
            top_idx = np.argsort(abs_loadings)[::-1][:3]
            top_genes = [feature_names[j] for j in top_idx if j < len(feature_names)]
            if top_genes:
                log.append(f"Top loading genes for {pc_label}: {', '.join(top_genes)}.")

    return log
