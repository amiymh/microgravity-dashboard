"""Biomarker discovery using OpenTargets API."""

import pandas as pd
import numpy as np
import requests
import plotly.graph_objects as go

OPENTARGETS_URL = "https://api.platform.opentargets.org/api/v4/graphql"


def query_opentargets_associations(gene_symbol: str, ensembl_id: str | None = None) -> list[dict]:
    """Query OpenTargets for disease associations of a gene."""
    # Try by Ensembl ID first (more reliable)
    target_id = ensembl_id if ensembl_id else gene_symbol

    query = """
    query($id: String!) {
      target(ensemblId: $id) {
        id
        approvedSymbol
        associatedDiseases(page: {size: 5, index: 0}) {
          rows {
            disease { id name }
            score
          }
        }
      }
    }
    """

    try:
        resp = requests.post(
            OPENTARGETS_URL,
            json={"query": query, "variables": {"id": target_id}},
            timeout=10,
        )
        resp.raise_for_status()
        data = resp.json()
    except Exception:
        return []

    target = data.get("data", {}).get("target")
    if not target:
        return []

    associations = []
    for row in target.get("associatedDiseases", {}).get("rows", []):
        associations.append({
            "disease": row["disease"]["name"],
            "disease_id": row["disease"]["id"],
            "score": row["score"],
        })

    return associations


def find_biomarkers(
    df: pd.DataFrame,
    min_log2fc: float = 2.0,
    min_sig_score: float = 0.0,
    max_genes: int = 50,
) -> pd.DataFrame:
    """Identify biomarker candidates among DEGs.

    Filters by thresholds, queries OpenTargets for disease associations,
    and scores candidates.
    """
    if df.empty or "Gene" not in df.columns:
        return pd.DataFrame()

    # Filter by thresholds
    candidates = df.copy()
    if "log2FoldChange" in candidates.columns:
        candidates = candidates[candidates["log2FoldChange"].abs() >= min_log2fc]
    if "significance score" in candidates.columns and min_sig_score > 0:
        candidates = candidates[candidates["significance score"] >= min_sig_score]

    if candidates.empty:
        return pd.DataFrame()

    # Limit to top genes
    if "significance score" in candidates.columns:
        candidates = candidates.nlargest(min(max_genes, len(candidates)), "significance score")

    # Query OpenTargets for disease associations
    disease_data = []
    api_success = 0

    for _, row in candidates.head(max_genes).iterrows():
        gene = row["Gene"]
        ens_id = row.get("ID", None)
        associations = query_opentargets_associations(gene, ens_id)

        if associations:
            api_success += 1
            disease_names = [a["disease"] for a in associations]
            max_score = max(a["score"] for a in associations)
        else:
            disease_names = []
            max_score = 0.0

        disease_data.append({
            "Gene": gene,
            "Disease Associations": "; ".join(disease_names) if disease_names else "None found",
            "Disease Score": max_score,
            "Known Biomarker": "Y" if max_score > 0.5 else "N",
            "Source": "OpenTargets" if associations else "N/A",
        })

    disease_df = pd.DataFrame(disease_data)

    # Merge with DEG info
    deg_cols = ["Gene", "Direction", "log2FoldChange", "padj", "significance score", "Gene Type"]
    available = [c for c in deg_cols if c in df.columns]
    deg_info = df[available].drop_duplicates(subset=["Gene"])
    result = disease_df.merge(deg_info, on="Gene", how="left")

    # Compute biomarker score
    if "log2FoldChange" in result.columns and "padj" in result.columns:
        result["Biomarker Score"] = (
            result["log2FoldChange"].abs()
            * (-np.log10(result["padj"].clip(lower=1e-300)))
            * (result["Disease Score"] + 0.1)
        )
        result = result.sort_values("Biomarker Score", ascending=False)

    return result.reset_index(drop=True)


def create_biomarker_scatter(df: pd.DataFrame) -> go.Figure:
    """Create scatter plot: log2FC vs -log10(padj), sized by disease score."""
    if df.empty or "log2FoldChange" not in df.columns:
        fig = go.Figure()
        fig.update_layout(title="No biomarker candidates to display")
        return fig

    plot_df = df.copy()
    plot_df["neg_log10_padj"] = -np.log10(plot_df["padj"].clip(lower=1e-300))
    size = plot_df.get("Disease Score", pd.Series([5] * len(plot_df)))
    plot_df["marker_size"] = (size * 30).clip(lower=5, upper=50)

    colors = plot_df.get("Direction", pd.Series(["Unknown"] * len(plot_df)))
    color_map = {"Upregulated": "#e74c3c", "Downregulated": "#3498db"}

    fig = go.Figure()

    for direction, color in color_map.items():
        mask = colors == direction
        subset = plot_df[mask]
        if subset.empty:
            continue

        fig.add_trace(
            go.Scatter(
                x=subset["log2FoldChange"],
                y=subset["neg_log10_padj"],
                mode="markers",
                name=direction,
                marker=dict(
                    color=color,
                    size=subset["marker_size"],
                    opacity=0.7,
                    line=dict(width=1, color="white"),
                ),
                text=subset["Gene"],
                hovertemplate="<b>%{text}</b><br>log2FC: %{x:.2f}<br>-log10(padj): %{y:.1f}<extra></extra>",
            )
        )

    fig.update_layout(
        title="Biomarker Candidates",
        xaxis_title="log2 Fold Change",
        yaxis_title="-log10(adjusted p-value)",
        template="plotly_white",
        height=500,
    )

    return fig
