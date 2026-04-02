"""Pathway enrichment analysis using Enrichr API."""

import pandas as pd
import numpy as np
import requests
import plotly.graph_objects as go

ENRICHR_URL = "https://maayanlab.cloud/Enrichr"

GENE_SET_LIBRARIES = [
    "KEGG_2021_Human",
    "GO_Biological_Process_2023",
    "GO_Molecular_Function_2023",
    "Reactome_2022",
    "WikiPathway_2023_Human",
]

# Demo results to use when API is unavailable
DEMO_RESULTS = pd.DataFrame({
    "term": [
        "TNF signaling pathway", "Apoptosis", "NF-kappa B signaling",
        "p53 signaling pathway", "MAPK signaling pathway",
        "Cytokine-cytokine receptor interaction", "IL-17 signaling pathway",
        "Cell cycle", "Toll-like receptor signaling", "Chemokine signaling",
    ],
    "pvalue": [1e-12, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2],
    "adjusted_pvalue": [1e-10, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 5e-2, 1e-1],
    "combined_score": [500, 400, 350, 300, 250, 200, 150, 100, 50, 25],
    "overlap": ["15/200", "12/161", "10/100", "8/72", "20/295",
                "18/294", "7/93", "10/124", "6/104", "5/190"],
    "genes": [
        "TNF,NFKB1,CASP3", "BCL2,BAX,CASP3", "NFKB1,TNF,RELA",
        "TP53,CDKN1A,BAX", "MAPK1,MAPK3,RAF1", "IL6,TNF,CCL2",
        "IL17A,NFKB1,TNF", "CDK1,CDK2,CCNB1", "TLR4,MYD88,NFKB1",
        "CCL2,CXCL8,CCR5",
    ],
})


def submit_gene_list(genes: list[str]) -> str | None:
    """Submit a gene list to Enrichr and return the list ID."""
    payload = {
        "list": (None, "\n".join(genes)),
        "description": (None, "microgravity DEGs"),
    }
    try:
        resp = requests.post(f"{ENRICHR_URL}/addList", files=payload, timeout=10)
        resp.raise_for_status()
        data = resp.json()
        return str(data.get("userListId"))
    except Exception:
        return None


def get_enrichment(list_id: str, gene_set: str) -> pd.DataFrame | None:
    """Fetch enrichment results for a submitted gene list."""
    try:
        resp = requests.get(
            f"{ENRICHR_URL}/enrich",
            params={"userListId": list_id, "backgroundType": gene_set},
            timeout=10,
        )
        resp.raise_for_status()
        data = resp.json()
    except Exception:
        return None

    results = data.get(gene_set, [])
    if not results:
        return pd.DataFrame()

    rows = []
    for r in results:
        rows.append({
            "rank": r[0],
            "term": r[1],
            "pvalue": r[2],
            "zscore": r[3],
            "combined_score": r[4],
            "overlap_genes": r[5],
            "adjusted_pvalue": r[6],
            "overlap": f"{len(r[5])}/{r[7]}" if len(r) > 7 else f"{len(r[5])}/unknown",
            "genes": ",".join(r[5]),
        })

    return pd.DataFrame(rows)


def run_enrichment(
    genes: list[str],
    gene_set: str = "KEGG_2021_Human",
    top_n: int = 20,
) -> tuple[pd.DataFrame, bool]:
    """Run full enrichment pipeline. Returns (results_df, is_live_data).

    If API fails, returns demo results with is_live_data=False.
    """
    if not genes:
        return pd.DataFrame(), True

    list_id = submit_gene_list(genes)
    if list_id is None:
        return DEMO_RESULTS.head(top_n).copy(), False

    results = get_enrichment(list_id, gene_set)
    if results is None or results.empty:
        return DEMO_RESULTS.head(top_n).copy(), False

    return results.head(top_n), True


def create_pathway_chart(df: pd.DataFrame, top_n: int = 20) -> go.Figure:
    """Create horizontal bar chart of pathway enrichment results."""
    if df.empty:
        fig = go.Figure()
        fig.update_layout(title="No pathway results to display")
        return fig

    plot_df = df.head(top_n).copy()
    plot_df["neg_log10_padj"] = -np.log10(plot_df["adjusted_pvalue"].clip(lower=1e-300))
    plot_df = plot_df.sort_values("neg_log10_padj", ascending=True)

    fig = go.Figure(
        go.Bar(
            y=plot_df["term"],
            x=plot_df["neg_log10_padj"],
            orientation="h",
            marker=dict(
                color=plot_df["combined_score"],
                colorscale="Viridis",
                colorbar=dict(title="Combined Score"),
            ),
            hovertext=[
                f"<b>{t}</b><br>p-adj: {p:.2e}<br>Score: {s:.1f}<br>Overlap: {o}"
                for t, p, s, o in zip(
                    plot_df["term"],
                    plot_df["adjusted_pvalue"],
                    plot_df["combined_score"],
                    plot_df["overlap"],
                )
            ],
            hoverinfo="text",
        )
    )

    fig.update_layout(
        title="Pathway Enrichment",
        xaxis_title="-log10(adjusted p-value)",
        yaxis_title="",
        template="plotly_white",
        height=max(400, len(plot_df) * 30),
        margin=dict(l=300),
    )

    return fig
