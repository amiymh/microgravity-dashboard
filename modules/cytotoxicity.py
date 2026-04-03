"""Cytotoxicity & apoptosis analysis using MSigDB hallmark gene sets."""

import pandas as pd
import numpy as np
from scipy import stats
import requests
import plotly.graph_objects as go

MSIGDB_URL = "https://www.gsea-msigdb.org/gsea/msigdb/human/download_geneset.jsp"

# MSigDB identifiers mapped to display names
HALLMARK_IDS = {
    "Apoptosis": "HALLMARK_APOPTOSIS",
    "TNF Signaling via NF-kB": "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "p53 Pathway": "HALLMARK_P53_PATHWAY",
    "Reactive Oxygen Species": "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
    "Inflammatory Response": "HALLMARK_INFLAMMATORY_RESPONSE",
}

# Minimal fallback sets used only when MSigDB API is completely unreachable
_FALLBACK_SETS = {
    "Apoptosis": ["CASP3", "CASP8", "CASP9", "BAX", "BCL2", "FAS", "TNF", "TP53", "CYCS", "APAF1"],
    "TNF Signaling via NF-kB": ["TNF", "NFKB1", "RELA", "TRAF2", "IKBKB", "IL6", "CXCL8", "CCL2"],
    "p53 Pathway": ["TP53", "MDM2", "CDKN1A", "BAX", "BBC3", "GADD45A", "ATM", "CHEK2"],
    "Reactive Oxygen Species": ["SOD1", "SOD2", "CAT", "GPX1", "PRDX1", "TXN", "NQO1", "HMOX1"],
    "Inflammatory Response": ["IL1B", "IL6", "TNF", "CXCL8", "CCL2", "TLR4", "NLRP3", "PTGS2"],
}

# Module-level cache for fetched gene sets
_gene_set_cache: dict[str, list[str]] = {}


def fetch_geneset(name: str, msigdb_id: str) -> list[str]:
    """Fetch a hallmark gene set from MSigDB API.

    Returns gene symbols list. Falls back to minimal hardcoded list on failure.
    """
    if name in _gene_set_cache:
        return _gene_set_cache[name]

    try:
        resp = requests.get(
            MSIGDB_URL,
            params={"geneSetName": msigdb_id, "fileType": "json"},
            timeout=10,
        )
        resp.raise_for_status()
        data = resp.json()
        genes = data.get(msigdb_id, {}).get("geneSymbols", [])
        if genes:
            _gene_set_cache[name] = genes
            return genes
    except Exception:
        pass

    # Fallback
    fallback = _FALLBACK_SETS.get(name, [])
    _gene_set_cache[name] = fallback
    return fallback


def get_hallmark_sets(pathways: list[str] | None = None) -> dict[str, list[str]]:
    """Fetch all requested hallmark gene sets, returns {name: [genes]}.

    Also returns a second value indicating which sets came from the live API.
    """
    selected = pathways if pathways else list(HALLMARK_IDS.keys())
    result = {}
    for name in selected:
        if name in HALLMARK_IDS:
            result[name] = fetch_geneset(name, HALLMARK_IDS[name])
    return result


def get_available_pathway_names() -> list[str]:
    """Return the list of available pathway display names."""
    return list(HALLMARK_IDS.keys())


def is_live_data(name: str) -> bool:
    """Check if a pathway's gene set came from the live API (not fallback)."""
    if name not in _gene_set_cache:
        return False
    return len(_gene_set_cache[name]) > len(_FALLBACK_SETS.get(name, []))


def compute_overlaps(df: pd.DataFrame, pathways: list[str] | None = None) -> pd.DataFrame:
    """Compute overlap between DEGs and hallmark gene sets.

    Returns a summary DataFrame with overlap counts and Fisher's exact test p-values.
    """
    if df.empty or "Gene" not in df.columns:
        return pd.DataFrame()

    hallmark_sets = get_hallmark_sets(pathways)
    if not hallmark_sets:
        return pd.DataFrame()

    deg_genes = set(df["Gene"].dropna().str.upper())

    rows = []
    for pathway, genes in hallmark_sets.items():
        pathway_genes = set(g.upper() for g in genes)
        overlap = deg_genes & pathway_genes
        overlap_genes = sorted(overlap)

        a = len(overlap)
        b = len(pathway_genes - deg_genes)
        c = len(deg_genes - pathway_genes)
        d = 20000 - a - b - c
        _, pvalue = stats.fisher_exact([[a, b], [c, d]], alternative="greater")

        rows.append({
            "Pathway": pathway,
            "Pathway Size": len(pathway_genes),
            "DEG Overlap": a,
            "Overlap %": round(100 * a / len(pathway_genes), 1) if pathway_genes else 0,
            "Fisher p-value": pvalue,
            "Overlap Genes": ", ".join(overlap_genes[:20]),
            "Source": "MSigDB API" if is_live_data(pathway) else "Fallback",
        })

    return pd.DataFrame(rows)


def get_gene_pathway_matrix(df: pd.DataFrame, pathways: list[str] | None = None) -> pd.DataFrame:
    """Build a genes x pathways matrix for heatmap, with log2FC values."""
    if df.empty or "Gene" not in df.columns:
        return pd.DataFrame()

    hallmark_sets = get_hallmark_sets(pathways)
    if not hallmark_sets:
        return pd.DataFrame()

    deg_genes = set(df["Gene"].dropna().str.upper())

    all_pathway_genes = set()
    for genes in hallmark_sets.values():
        all_pathway_genes.update(g.upper() for g in genes)

    relevant_genes = sorted(deg_genes & all_pathway_genes)
    if not relevant_genes:
        return pd.DataFrame()

    gene_fc = {}
    for _, row in df.iterrows():
        gene_fc[str(row["Gene"]).upper()] = row.get("log2FoldChange", 0)

    matrix_data = {}
    for pw, genes in hallmark_sets.items():
        pw_genes = set(g.upper() for g in genes)
        matrix_data[pw] = [gene_fc.get(g, np.nan) if g in pw_genes else np.nan for g in relevant_genes]

    matrix = pd.DataFrame(matrix_data, index=relevant_genes)
    matrix = matrix.dropna(how="all")
    return matrix


def create_overlap_bar_chart(overlap_df: pd.DataFrame) -> go.Figure:
    """Create bar chart showing pathway enrichment p-values."""
    if overlap_df.empty:
        fig = go.Figure()
        fig.update_layout(title="No pathway overlaps found")
        return fig

    plot_df = overlap_df.sort_values("Fisher p-value", ascending=False).copy()
    plot_df["neg_log10_p"] = -np.log10(plot_df["Fisher p-value"].clip(lower=1e-300))

    fig = go.Figure(
        go.Bar(
            y=plot_df["Pathway"],
            x=plot_df["neg_log10_p"],
            orientation="h",
            marker=dict(
                color=plot_df["DEG Overlap"],
                colorscale="Reds",
                colorbar=dict(title="Overlap Count"),
            ),
            text=plot_df["DEG Overlap"].astype(str) + " genes",
            textposition="auto",
        )
    )

    fig.update_layout(
        title="Cytotoxicity / Apoptosis Pathway Enrichment",
        xaxis_title="-log10(Fisher p-value)",
        template="plotly_white",
        height=max(300, len(plot_df) * 60),
        margin=dict(l=250),
    )

    return fig


def create_heatmap(matrix: pd.DataFrame) -> go.Figure:
    """Create heatmap of gene x pathway with log2FC coloring."""
    if matrix.empty:
        fig = go.Figure()
        fig.update_layout(title="No data for heatmap")
        return fig

    if len(matrix) > 50:
        matrix = matrix.head(50)

    fig = go.Figure(
        go.Heatmap(
            z=matrix.values,
            x=matrix.columns.tolist(),
            y=matrix.index.tolist(),
            colorscale="RdBu_r",
            zmid=0,
            colorbar=dict(title="log2FC"),
            hovertemplate="Gene: %{y}<br>Pathway: %{x}<br>log2FC: %{z:.2f}<extra></extra>",
        )
    )

    fig.update_layout(
        title="Gene x Pathway Heatmap (log2 Fold Change)",
        template="plotly_white",
        height=max(400, len(matrix) * 15),
        margin=dict(l=100),
    )

    return fig


def get_gene_table(df: pd.DataFrame, pathways: list[str] | None = None) -> pd.DataFrame:
    """Get table of DEGs with their pathway memberships."""
    if df.empty or "Gene" not in df.columns:
        return pd.DataFrame()

    hallmark_sets = get_hallmark_sets(pathways)
    if not hallmark_sets:
        return pd.DataFrame()

    # Pre-compute upper-case sets for each pathway
    pw_upper = {pw: set(g.upper() for g in genes) for pw, genes in hallmark_sets.items()}

    rows = []
    for _, row in df.iterrows():
        gene = str(row["Gene"]).upper()
        memberships = [pw for pw, gene_set in pw_upper.items() if gene in gene_set]

        if memberships:
            rows.append({
                "Gene": row["Gene"],
                "log2FoldChange": row.get("log2FoldChange", np.nan),
                "padj": row.get("padj", np.nan),
                "Direction": row.get("Direction", ""),
                "Pathways": "; ".join(memberships),
                "Pathway Count": len(memberships),
            })

    result = pd.DataFrame(rows)
    if not result.empty:
        result = result.sort_values("Pathway Count", ascending=False)
    return result.reset_index(drop=True)
