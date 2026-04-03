"""Disease cross-reference using OpenTargets Platform API."""

import pandas as pd
import numpy as np
import requests
import plotly.graph_objects as go
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OPENTARGETS_URL = "https://api.platform.opentargets.org/api/v4/graphql"


def search_disease(query: str) -> list[dict]:
    """Search OpenTargets for diseases matching a query string."""
    gql = """
    query($q: String!) {
      search(queryString: $q, entityNames: ["disease"], page: {size: 10, index: 0}) {
        hits {
          id
          name
          entity
          description
        }
      }
    }
    """
    try:
        resp = requests.post(
            OPENTARGETS_URL,
            json={"query": gql, "variables": {"q": query}},
            timeout=10,
        )
        resp.raise_for_status()
        data = resp.json()
    except Exception:
        return []

    hits = data.get("data", {}).get("search", {}).get("hits", [])
    return [{"id": h["id"], "name": h["name"], "description": h.get("description", "")} for h in hits]


def get_disease_genes(disease_id: str, min_score: float = 0.3) -> pd.DataFrame:
    """Get genes associated with a disease from OpenTargets."""
    gql = """
    query($diseaseId: String!) {
      disease(efoId: $diseaseId) {
        id
        name
        associatedTargets(page: {size: 500, index: 0}) {
          rows {
            target { id approvedSymbol }
            score
          }
        }
      }
    }
    """
    try:
        resp = requests.post(
            OPENTARGETS_URL,
            json={"query": gql, "variables": {"diseaseId": disease_id}},
            timeout=10,
        )
        resp.raise_for_status()
        data = resp.json()
    except Exception:
        return pd.DataFrame()

    disease = data.get("data", {}).get("disease")
    if not disease:
        return pd.DataFrame()

    rows = []
    for r in disease.get("associatedTargets", {}).get("rows", []):
        score = r.get("score", 0)
        if score >= min_score:
            rows.append({
                "Gene": r["target"]["approvedSymbol"],
                "Ensembl ID": r["target"]["id"],
                "Disease Association Score": score,
            })

    return pd.DataFrame(rows) if rows else pd.DataFrame(
        columns=["Gene", "Ensembl ID", "Disease Association Score"]
    )


def cross_reference(
    deg_df: pd.DataFrame,
    disease_id: str,
    min_score: float = 0.3,
    direction: str = "All",
) -> tuple[pd.DataFrame, dict]:
    """Cross-reference DEGs with disease genes.

    Returns (overlap_df, summary_dict).
    """
    if deg_df.empty or "Gene" not in deg_df.columns:
        return pd.DataFrame(), {"overlap": 0, "disease_genes": 0, "deg_total": 0}

    disease_genes = get_disease_genes(disease_id, min_score)

    if disease_genes.empty:
        return pd.DataFrame(), {
            "overlap": 0,
            "disease_genes": 0,
            "deg_total": len(deg_df),
        }

    # Filter DEGs by direction
    filtered_deg = deg_df.copy()
    if direction != "All" and "Direction" in filtered_deg.columns:
        filtered_deg = filtered_deg[filtered_deg["Direction"] == direction]

    # Find overlap
    deg_gene_set = set(filtered_deg["Gene"].dropna().str.upper())
    disease_gene_set = set(disease_genes["Gene"].str.upper())
    overlap_genes = deg_gene_set & disease_gene_set

    if not overlap_genes:
        return pd.DataFrame(), {
            "overlap": 0,
            "disease_genes": len(disease_gene_set),
            "deg_total": len(deg_gene_set),
        }

    # Build overlap table
    # Get DEG info
    deg_info = filtered_deg[filtered_deg["Gene"].str.upper().isin(overlap_genes)].copy()
    # Get disease scores
    disease_scores = disease_genes[disease_genes["Gene"].str.upper().isin(overlap_genes)].copy()
    disease_scores["Gene_upper"] = disease_scores["Gene"].str.upper()
    deg_info["Gene_upper"] = deg_info["Gene"].str.upper()

    merged = deg_info.merge(
        disease_scores[["Gene_upper", "Disease Association Score"]],
        on="Gene_upper",
        how="left",
    )
    merged = merged.drop(columns=["Gene_upper"])

    # Select display columns
    display_cols = ["Gene", "Direction", "log2FoldChange", "padj", "Disease Association Score"]
    available = [c for c in display_cols if c in merged.columns]
    result = merged[available].drop_duplicates(subset=["Gene"]).sort_values(
        "Disease Association Score", ascending=False
    )

    summary = {
        "overlap": len(overlap_genes),
        "disease_genes": len(disease_gene_set),
        "deg_total": len(deg_gene_set),
    }

    return result.reset_index(drop=True), summary


def create_venn_diagram(summary: dict, disease_name: str = "Disease") -> plt.Figure:
    """Create an approximate Venn diagram using matplotlib circles."""
    fig, ax = plt.subplots(1, 1, figsize=(6, 4))

    deg_total = summary.get("deg_total", 0)
    disease_total = summary.get("disease_genes", 0)
    overlap = summary.get("overlap", 0)

    # Draw circles
    from matplotlib.patches import Circle

    circle1 = Circle((0.35, 0.5), 0.3, alpha=0.3, color="#3498db", label="Microgravity DEGs")
    circle2 = Circle((0.65, 0.5), 0.3, alpha=0.3, color="#e74c3c", label=disease_name)
    ax.add_patch(circle1)
    ax.add_patch(circle2)

    # Labels
    ax.text(0.2, 0.5, f"{deg_total - overlap}", ha="center", va="center", fontsize=14, fontweight="bold")
    ax.text(0.5, 0.5, f"{overlap}", ha="center", va="center", fontsize=14, fontweight="bold")
    ax.text(0.8, 0.5, f"{disease_total - overlap}", ha="center", va="center", fontsize=14, fontweight="bold")

    ax.text(0.2, 0.15, "DEGs only", ha="center", fontsize=9)
    ax.text(0.5, 0.15, "Overlap", ha="center", fontsize=9)
    ax.text(0.8, 0.15, f"{disease_name}\nonly", ha="center", fontsize=9)

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(f"Microgravity DEGs vs {disease_name} Genes", fontsize=12)

    plt.tight_layout()
    return fig


def create_network_graph(overlap_df: pd.DataFrame, disease_name: str = "Disease") -> go.Figure:
    """Create a simple network visualization of overlapping genes."""
    if overlap_df.empty or "Gene" not in overlap_df.columns:
        fig = go.Figure()
        fig.update_layout(title="No overlap to display")
        return fig

    genes = overlap_df.head(30)  # Limit for readability
    n = len(genes)

    # Arrange genes in a circle around the disease node
    angles = np.linspace(0, 2 * np.pi, n, endpoint=False)
    gene_x = np.cos(angles) * 2
    gene_y = np.sin(angles) * 2

    # Edge traces (gene → center disease node)
    edge_x, edge_y = [], []
    for gx, gy in zip(gene_x, gene_y):
        edge_x.extend([0, gx, None])
        edge_y.extend([0, gy, None])

    fig = go.Figure()

    fig.add_trace(go.Scatter(
        x=edge_x, y=edge_y, mode="lines",
        line=dict(width=0.5, color="#ccc"),
        hoverinfo="none", showlegend=False,
    ))

    # Color by direction
    colors = []
    for d in genes.get("Direction", ["grey"] * n):
        if d == "Upregulated":
            colors.append("#e74c3c")
        elif d == "Downregulated":
            colors.append("#3498db")
        else:
            colors.append("#95a5a6")

    fig.add_trace(go.Scatter(
        x=gene_x.tolist(), y=gene_y.tolist(), mode="markers+text",
        text=genes["Gene"].tolist(),
        textposition="top center",
        textfont=dict(size=8),
        marker=dict(size=10, color=colors, line=dict(width=1, color="white")),
        hovertext=[
            f"<b>{g}</b><br>log2FC: {fc:.2f}" if pd.notna(fc) else f"<b>{g}</b>"
            for g, fc in zip(genes["Gene"], genes.get("log2FoldChange", [np.nan] * n))
        ],
        hoverinfo="text",
        showlegend=False,
    ))

    # Center disease node
    fig.add_trace(go.Scatter(
        x=[0], y=[0], mode="markers+text",
        text=[disease_name], textposition="bottom center",
        marker=dict(size=20, color="#f39c12", symbol="diamond",
                    line=dict(width=2, color="white")),
        showlegend=False,
    ))

    fig.update_layout(
        title=f"Gene-Disease Network: {disease_name}",
        template="plotly_white",
        showlegend=False,
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        height=500,
    )

    return fig


def get_analysis_log(
    disease_name: str = "",
    disease_id: str = "",
    min_score: float = 0.3,
    direction: str = "All",
    summary: dict | None = None,
) -> list[str]:
    """Return analysis log for disease cross-reference."""
    from datetime import datetime
    log = [
        f"Disease searched: '{disease_name}'",
        f"Disease EFO ID: {disease_id}" if disease_id else "Disease EFO ID: not resolved",
        f"Association score threshold: {min_score}",
        f"Direction filter: {direction}",
        f"API endpoint: {OPENTARGETS_URL}",
        f"Timestamp: {datetime.now().isoformat()}",
    ]
    if summary:
        log.append(f"Disease genes returned by OpenTargets: {summary.get('disease_genes', 0)}")
        log.append(f"DEGs queried: {summary.get('deg_total', 0)}")
        log.append(f"Overlapping genes found: {summary.get('overlap', 0)}")
        if summary.get('disease_genes', 0) > 0:
            pct = 100 * summary.get('overlap', 0) / summary['disease_genes']
            log.append(f"Overlap percentage: {pct:.1f}% of disease genes are DEGs")
    else:
        log.append("Cross-reference not yet run.")
    return log
