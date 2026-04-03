"""Volcano plot visualization for differential gene expression."""

import pandas as pd
import numpy as np
import plotly.graph_objects as go

try:
    from adjustText import adjust_text as _adjust_text
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    _HAS_ADJUSTTEXT = True
except ImportError:
    _HAS_ADJUSTTEXT = False


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


def _compute_adjusted_positions(x_vals, y_vals, labels):
    """Use adjustText + matplotlib to compute non-overlapping label positions."""
    if not _HAS_ADJUSTTEXT or len(x_vals) == 0:
        return None

    fig_mpl, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(x_vals, y_vals, s=1, alpha=0)
    texts = [ax.text(x, y, lab, fontsize=8) for x, y, lab in zip(x_vals, y_vals, labels)]

    _adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", color="grey", lw=0.5))

    positions = []
    for t in texts:
        pos = t.get_position()
        positions.append((pos[0], pos[1]))

    plt.close(fig_mpl)
    return positions


def create_volcano_plot(
    df: pd.DataFrame,
    padj_cutoff: float = 0.05,
    log2fc_cutoff: float = 1.0,
    label_top_n: int = 10,
    color_scheme: str = "Default",
) -> go.Figure:
    """Create an interactive Plotly volcano plot with non-overlapping gene labels."""
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

    # Label top N genes with non-overlapping positions
    if label_top_n > 0 and "significance score" in classified.columns:
        sig = classified[classified["volcano_class"] != "Not Significant"]
        if not sig.empty:
            top = sig.nlargest(min(label_top_n, len(sig)), "significance score")
            x_vals = top["log2FoldChange"].tolist()
            y_vals = top["neg_log10_padj"].tolist()
            labels = top["Gene"].tolist()

            # Try adjustText for non-overlapping positions
            adjusted = _compute_adjusted_positions(x_vals, y_vals, labels)

            if adjusted is not None:
                for (ax, ay), ox, oy, label in zip(adjusted, x_vals, y_vals, labels):
                    fig.add_annotation(
                        x=ox, y=oy, ax=ax, ay=ay,
                        text=label, font=dict(size=9),
                        showarrow=True,
                        arrowhead=0, arrowwidth=0.5, arrowcolor="grey",
                        xref="x", yref="y", axref="x", ayref="y",
                    )
            else:
                # Fallback: simple text labels
                fig.add_trace(
                    go.Scattergl(
                        x=x_vals, y=y_vals,
                        mode="text", text=labels,
                        textposition="top center",
                        textfont=dict(size=9),
                        showlegend=False, hoverinfo="skip",
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


def get_analysis_log(
    df: pd.DataFrame,
    padj_cutoff: float = 0.05,
    log2fc_cutoff: float = 1.0,
    label_top_n: int = 10,
) -> list[str]:
    """Return analysis log for the volcano plot."""
    log = []
    total = len(df) if not df.empty else 0
    log.append(f"Total genes in dataset: {total}")

    if df.empty:
        log.append("No data available for analysis.")
        return log

    classified = classify_genes(df, padj_cutoff, log2fc_cutoff)
    n_up = (classified["volcano_class"] == "Upregulated").sum()
    n_down = (classified["volcano_class"] == "Downregulated").sum()
    n_ns = (classified["volcano_class"] == "Not Significant").sum()

    log.append(f"Adjusted p-value threshold: {padj_cutoff}")
    log.append(f"|log2FoldChange| threshold: {log2fc_cutoff}")
    log.append(f"Genes passing padj threshold: {(df['padj'] <= padj_cutoff).sum()}")
    log.append(f"Genes passing |log2FC| threshold: {(df['log2FoldChange'].abs() >= log2fc_cutoff).sum()}")
    log.append(f"Upregulated (both thresholds): {n_up}")
    log.append(f"Downregulated (both thresholds): {n_down}")
    log.append(f"Not significant: {n_ns}")
    log.append(f"Top {label_top_n} genes labeled by significance score")
    log.append(f"Label method: {'adjustText (non-overlapping)' if _HAS_ADJUSTTEXT else 'Plotly text (may overlap)'}")

    return log
