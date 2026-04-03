"""Biomarker discovery using OpenTargets API with batched queries."""

import time
import pandas as pd
import numpy as np
import requests
import plotly.graph_objects as go

OPENTARGETS_URL = "https://api.platform.opentargets.org/api/v4/graphql"

# Rate limiting: pause between batches to respect API limits
BATCH_DELAY_SECONDS = 0.2
BATCH_SIZE = 20


def query_opentargets_batch(
    genes: list[dict],
) -> dict[str, list[dict]]:
    """Query OpenTargets for disease associations of multiple genes.

    Args:
        genes: List of dicts with 'symbol' and optionally 'ensembl_id' keys.

    Returns:
        Dict mapping gene symbol to list of disease association dicts.
    """
    results = {}

    for gene_info in genes:
        symbol = gene_info["symbol"]
        ens_id = gene_info.get("ensembl_id")
        target_id = ens_id if ens_id else symbol

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
            results[symbol] = []
            continue

        target = data.get("data", {}).get("target")
        if not target:
            results[symbol] = []
            continue

        associations = []
        for row in target.get("associatedDiseases", {}).get("rows", []):
            associations.append({
                "disease": row["disease"]["name"],
                "disease_id": row["disease"]["id"],
                "score": row["score"],
            })
        results[symbol] = associations

    return results


def find_biomarkers(
    df: pd.DataFrame,
    min_log2fc: float = 2.0,
    min_sig_score: float = 0.0,
    progress_callback=None,
) -> pd.DataFrame:
    """Identify biomarker candidates among ALL filtered DEGs.

    Processes all genes that pass the thresholds using batched API calls
    with rate limiting. No silent cutoffs.

    Args:
        df: DEG dataframe with Gene, log2FoldChange, padj, significance score, etc.
        min_log2fc: Minimum |log2FC| filter.
        min_sig_score: Minimum significance score filter.
        progress_callback: Optional callable(current, total) for progress updates.

    Returns:
        DataFrame of biomarker candidates ranked by score.
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

    # Sort by significance score for best results first
    if "significance score" in candidates.columns:
        candidates = candidates.sort_values("significance score", ascending=False)

    # Build gene list for batch querying
    gene_list = []
    seen = set()
    for _, row in candidates.iterrows():
        symbol = row["Gene"]
        if symbol in seen or pd.isna(symbol):
            continue
        seen.add(symbol)
        gene_list.append({
            "symbol": symbol,
            "ensembl_id": row.get("ID", None),
        })

    total_genes = len(gene_list)
    if total_genes == 0:
        return pd.DataFrame()

    # Process in batches with rate limiting
    all_associations = {}
    for batch_start in range(0, total_genes, BATCH_SIZE):
        batch = gene_list[batch_start : batch_start + BATCH_SIZE]
        batch_results = query_opentargets_batch(batch)
        all_associations.update(batch_results)

        if progress_callback:
            progress_callback(min(batch_start + BATCH_SIZE, total_genes), total_genes)

        # Rate limit between batches
        if batch_start + BATCH_SIZE < total_genes:
            time.sleep(BATCH_DELAY_SECONDS)

    # Build result table
    disease_data = []
    for gene_info in gene_list:
        symbol = gene_info["symbol"]
        associations = all_associations.get(symbol, [])

        if associations:
            disease_names = [a["disease"] for a in associations]
            max_score = max(a["score"] for a in associations)
        else:
            disease_names = []
            max_score = 0.0

        disease_data.append({
            "Gene": symbol,
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
