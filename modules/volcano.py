"""Volcano plot visualization for differential gene expression."""

import pandas as pd
import numpy as np
import plotly.graph_objects as go


def classify_genes(
    df: pd.DataFrame,
    padj_cutoff: float = 0.05,
    log2fc_cutoff: float = 1.0,
) -> pd.DataFrame:
    """Add classification column: Upregulated, Downregulated, or Not Significant."""
    result = df.copy()
    neg_log10_padj = -np.log10(result["padj"].clip(lower=1e-300))
    result["neg_log10_padj"] = neg_log10_padj

    conditions = [
        (result["padj"] <= padj_cutoff) & (result["log2FoldChange"] >= log2fc_cutoff),
        (result["padj"] <= padj_cutoff) & (result["log2FoldChange"] <= -log2fc_cutoff),
    ]
    choices = ["Upregulated", "Downregulated"]
    result["volcano_class"] = np.select(conditions, choices, default="Not Significant")
    return result


def create_volcano_plot(
    df: pd.DataFrame,
    padj_cutoff: float = 0.05,
    log2fc_cutoff: float = 1.0,
    label_top_n: int = 10,
    color_scheme: str = "Default",
) -> go.Figure:
    """Create an interactive Plotly volcano plot.

    Args:
        df: DataFrame with Gene, log2FoldChange, padj, significance score columns.
        padj_cutoff: Adjusted p-value threshold.
        log2fc_cutoff: |log2FC| threshold.
        label_top_n: Number of top genes to label.
        color_scheme: Color scheme name.

    Returns:
        Plotly Figure object.
    """
    if df.empty:
        fig = go.Figure()
        fig.update_layout(title="No data to display")
        return fig

    classified = classify_genes(df, padj_cutoff, log2fc_cutoff)

    color_map = {
        "Default": {"Upregulated": "#e74c3c", "Downregulated": "#3498db", "Not Significant": "#bdc3c7"},
        "Warm": {"Upregulated": "#ff6b6b", "Downregulated": "#4ecdc4", "Not Significant": "#95a5a6"},
        "Cool": {"Upregulated": "#e63946", "Downregulated": "#457b9d", "Not Significant": "#adb5bd"},
    }
    colors = color_map.get(color_scheme, color_map["Default"])

    fig = go.Figure()

    for cls in ["Not Significant", "Downregulated", "Upregulated"]:
        subset = classified[classified["volcano_class"] == cls]
        if subset.empty:
            continue

        hover_text = [
            f"<b>{g}</b><br>log2FC: {fc:.2f}<br>padj: {p:.2e}<br>{d}<br>Sig score: {s:.1f}"
            for g, fc, p, d, s in zip(
                subset["Gene"],
                subset["log2FoldChange"],
                subset["padj"],
                subset.get("Direction", [""] * len(subset)),
                subset.get("significance score", [0] * len(subset)),
            )
        ]

        fig.add_trace(
            go.Scattergl(
                x=subset["log2FoldChange"],
                y=subset["neg_log10_padj"],
                mode="markers",
                name=cls,
                marker=dict(color=colors[cls], size=5, opacity=0.7),
                hovertext=hover_text,
                hoverinfo="text",
            )
        )

    # Threshold lines
    fig.add_hline(y=-np.log10(padj_cutoff), line_dash="dash", line_color="grey", opacity=0.5)
    fig.add_vline(x=log2fc_cutoff, line_dash="dash", line_color="grey", opacity=0.5)
    fig.add_vline(x=-log2fc_cutoff, line_dash="dash", line_color="grey", opacity=0.5)

    # Label top N genes
    if label_top_n > 0 and "significance score" in classified.columns:
        sig = classified[classified["volcano_class"] != "Not Significant"]
        if not sig.empty:
            top = sig.nlargest(min(label_top_n, len(sig)), "significance score")
            fig.add_trace(
                go.Scattergl(
                    x=top["log2FoldChange"],
                    y=top["neg_log10_padj"],
                    mode="text",
                    text=top["Gene"],
                    textposition="top center",
                    textfont=dict(size=9),
                    showlegend=False,
                    hoverinfo="skip",
                )
            )

    fig.update_layout(
        title="Volcano Plot — Microgravity vs Ground",
        xaxis_title="log2 Fold Change (Space / Ground)",
        yaxis_title="-log10(adjusted p-value)",
        template="plotly_white",
        legend=dict(orientation="h", yanchor="bottom", y=1.02),
        height=600,
    )

    return fig


def get_top_significant(df: pd.DataFrame, n: int = 20) -> pd.DataFrame:
    """Return top N most significant genes sorted by significance score."""
    if df.empty or "significance score" not in df.columns:
        return pd.DataFrame()

    cols = ["Gene", "log2FoldChange", "padj", "foldChange", "significance score", "Direction", "Gene Type"]
    available = [c for c in cols if c in df.columns]
    return df.nlargest(min(n, len(df)), "significance score")[available].reset_index(drop=True)
