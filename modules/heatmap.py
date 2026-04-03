"""Heatmap visualization for top differentially expressed genes."""

import pandas as pd
import numpy as np
import plotly.graph_objects as go
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist

EARTH_COLS = [
    "Earth_3D2", "Earth_4D2", "Earth_3E2", "Earth_4E2",
    "Earth_1D2", "Earth_2D2", "Earth_1E2", "Earth_2E2",
]

SPACE_COLS = [
    "Space_3D2", "Space_4D2", "Space_3E2", "Space_4E2",
    "Space_1D2", "Space_2D2", "Space_1E2", "Space_2E2",
]

ALL_SAMPLE_COLS = EARTH_COLS + SPACE_COLS


def compute_heatmap_data(
    df: pd.DataFrame,
    top_n: int = 50,
    cluster: bool = True,
) -> dict | None:
    """Select top DEGs, z-score normalize, and optionally cluster.

    Args:
        df: DataFrame with gene expression counts and metadata.
        top_n: Number of top genes to select by significance score.
        cluster: Whether to reorder genes via hierarchical clustering.

    Returns:
        Dict with z_matrix, gene_info, gene_order, sample_names, n_genes,
        or None if the input is empty.
    """
    if df.empty or "significance score" not in df.columns:
        return None

    # Ensure all sample columns are present
    available_cols = [c for c in ALL_SAMPLE_COLS if c in df.columns]
    if not available_cols:
        return None

    # Select top N genes by significance score
    n = min(top_n, len(df))
    top_genes = df.nlargest(n, "significance score").copy()

    if top_genes.empty:
        return None

    # Extract count matrix (genes x samples)
    count_matrix = top_genes[available_cols].copy()
    count_matrix.index = top_genes["Gene"].values

    # Z-score normalize per gene (row-wise): z = (x - mean) / std
    row_means = count_matrix.mean(axis=1)
    row_stds = count_matrix.std(axis=1)
    # Avoid division by zero for genes with constant expression
    row_stds = row_stds.replace(0, 1)
    z_matrix = count_matrix.sub(row_means, axis=0).div(row_stds, axis=0)

    # Determine gene order
    gene_order = list(z_matrix.index)

    if cluster and len(z_matrix) > 1:
        # Hierarchical clustering with ward linkage
        dist = pdist(z_matrix.values, metric="euclidean")
        Z = linkage(dist, method="ward")
        order_idx = leaves_list(Z)
        gene_order = [gene_order[i] for i in order_idx]
        z_matrix = z_matrix.loc[gene_order]

    # Build gene_info DataFrame
    info_cols = ["Gene", "log2FoldChange", "padj", "Direction"]
    available_info = [c for c in info_cols if c in top_genes.columns]
    gene_info = top_genes[available_info].reset_index(drop=True)

    return {
        "z_matrix": z_matrix,
        "gene_info": gene_info,
        "gene_order": gene_order,
        "sample_names": available_cols,
        "n_genes": len(gene_order),
    }


def create_heatmap_figure(heatmap_data: dict | None) -> go.Figure:
    """Create a Plotly heatmap figure from compute_heatmap_data output.

    Args:
        heatmap_data: Output dict from compute_heatmap_data, or None.

    Returns:
        Plotly Figure object.
    """
    if heatmap_data is None:
        fig = go.Figure()
        fig.update_layout(title="No data to display")
        return fig

    z_matrix = heatmap_data["z_matrix"]
    gene_order = heatmap_data["gene_order"]
    sample_names = heatmap_data["sample_names"]

    # Build hover text
    hover_text = []
    for gene in gene_order:
        row_text = []
        for sample in sample_names:
            z_val = z_matrix.loc[gene, sample]
            row_text.append(f"Gene: {gene}<br>Sample: {sample}<br>Z-score: {z_val:.2f}")
        hover_text.append(row_text)

    # Adaptive height: 15px per gene, clamped to 800
    height = min(800, max(400, len(gene_order) * 15 + 150))

    fig = go.Figure(
        data=go.Heatmap(
            z=z_matrix.values,
            x=sample_names,
            y=gene_order,
            colorscale="RdBu_r",
            zmid=0,
            text=hover_text,
            hoverinfo="text",
            colorbar=dict(title="Z-score"),
        )
    )

    fig.update_layout(
        title="Heatmap — Top Differentially Expressed Genes",
        xaxis_title="Samples",
        yaxis_title="Genes",
        template="plotly_white",
        height=height,
        yaxis=dict(autorange="reversed"),
    )

    return fig


def get_analysis_log(heatmap_data: dict | None, cluster: bool = True) -> list[str]:
    """Return a human-readable analysis log for transparency.

    Args:
        heatmap_data: Output dict from compute_heatmap_data, or None.
        cluster: Whether clustering was applied.

    Returns:
        List of strings describing the analysis steps.
    """
    if heatmap_data is None:
        return ["No heatmap data available."]

    n = heatmap_data["n_genes"]
    samples = heatmap_data["sample_names"]

    log = [
        f"Selected top {n} genes by significance score.",
        f"Selection criteria: highest significance score (|log2FC| * -log10(padj)).",
        f"Z-score normalization per gene (row): z = (x - mean) / std",
        f"Columns used: {', '.join(samples)}",
    ]

    if cluster:
        log.append("Clustering method: scipy hierarchical clustering, ward linkage")

    return log
