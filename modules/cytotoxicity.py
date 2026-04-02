"""Cytotoxicity & apoptosis analysis using curated hallmark gene sets."""

import pandas as pd
import numpy as np
from scipy import stats
import plotly.graph_objects as go

# Curated hallmark gene sets (from MSigDB literature)
HALLMARK_SETS = {
    "Apoptosis": [
        "CASP1", "CASP2", "CASP3", "CASP4", "CASP6", "CASP7", "CASP8", "CASP9", "CASP10",
        "BAX", "BAK1", "BCL2", "BCL2L1", "BCL2L11", "BID", "BIK", "BIRC2", "BIRC3",
        "XIAP", "MCL1", "CYCS", "APAF1", "DIABLO", "DFFA", "DFFB", "ENDOG",
        "AIFM1", "PARP1", "PARP2", "TNFRSF10A", "TNFRSF10B", "TNFRSF1A",
        "FAS", "FADD", "TRADD", "RIPK1", "CFLAR", "GADD45A", "GADD45B",
        "TP53", "CDKN1A", "BBC3", "PMAIP1", "ATM", "CHEK2", "PIDD1",
        "LMNA", "LMNB1", "SPTAN1", "ADD1", "SATB1", "HMGB2",
        "NFKB1", "RELA", "NFKBIA", "TNF", "TNFSF10",
        "CRADD", "DAPK1", "DAPK2", "DAP", "DAXX",
    ],
    "TNF Signaling via NF-kB": [
        "TNF", "TNFAIP3", "TNFRSF1A", "TNFRSF1B", "TRADD", "TRAF1", "TRAF2",
        "NFKB1", "NFKB2", "RELA", "RELB", "NFKBIA", "NFKBIB", "NFKBIE",
        "IKBKB", "IKBKG", "CHUK", "MAP3K7", "TAB1", "TAB2",
        "BIRC2", "BIRC3", "XIAP", "CFLAR", "RIPK1",
        "CXCL1", "CXCL2", "CXCL3", "CXCL8", "CCL2", "CCL5", "CCL20",
        "IL6", "IL1B", "IL1A", "CSF1", "CSF2", "LIF",
        "ICAM1", "VCAM1", "SELE", "MMP9", "PTGS2",
        "BCL2A1", "BCL3", "IRF1", "JUNB", "FOS", "FOSB", "ATF3",
        "NFKBIZ", "SOD2", "TNIP1", "TNIP2", "ZFP36",
    ],
    "p53 Pathway": [
        "TP53", "MDM2", "MDM4", "CDKN1A", "CDKN2A", "RB1",
        "BAX", "BBC3", "PMAIP1", "BID", "GADD45A", "GADD45B", "GADD45G",
        "SESN1", "SESN2", "ATM", "ATR", "CHEK1", "CHEK2",
        "DDB2", "XPC", "POLK", "RRM2B", "TIGAR", "SCO2",
        "PTEN", "TSC2", "DRAM1", "ZMAT3", "PERP", "EI24",
        "FAS", "TNFRSF10B", "PIDD1", "CASP1", "CASP8",
        "SERPINE1", "THBS1", "MASPIN", "GLS2", "FDXR",
        "STEAP3", "PML", "SIVA1", "TP53I3", "TP53INP1",
    ],
    "Reactive Oxygen Species": [
        "SOD1", "SOD2", "SOD3", "CAT", "GPX1", "GPX2", "GPX3", "GPX4",
        "PRDX1", "PRDX2", "PRDX3", "PRDX4", "PRDX5", "PRDX6",
        "TXN", "TXN2", "TXNRD1", "TXNRD2", "GLRX", "GLRX2",
        "GSR", "GCLC", "GCLM", "GSS", "GSTP1", "GSTA1",
        "NQO1", "HMOX1", "NFE2L2", "KEAP1", "SQSTM1",
        "NOX1", "NOX2", "NOX4", "DUOX1", "DUOX2",
        "MPO", "NOS2", "NOS3", "XDH",
    ],
    "Inflammatory Response": [
        "IL1A", "IL1B", "IL1R1", "IL1R2", "IL1RN", "IL1RAP",
        "IL6", "IL6R", "IL6ST", "IL10", "IL10RA", "IL10RB",
        "IL18", "IL18R1", "IL33", "IL1RL1",
        "TNF", "TNFRSF1A", "TNFRSF1B", "LTA", "LTB",
        "CXCL1", "CXCL2", "CXCL3", "CXCL5", "CXCL8", "CXCL10", "CXCL11",
        "CCL2", "CCL3", "CCL4", "CCL5", "CCL7", "CCL8", "CCL11", "CCL20",
        "CXCR1", "CXCR2", "CCR1", "CCR2", "CCR5",
        "PTGS2", "PTGES", "ALOX5", "LTA4H",
        "TLR1", "TLR2", "TLR4", "TLR5", "TLR6", "TLR7", "TLR8", "TLR9",
        "MYD88", "IRAK1", "IRAK4", "TRAF6",
        "NLRP3", "PYCARD", "CASP1", "IL18",
        "SELE", "SELP", "ICAM1", "VCAM1",
        "C3", "C5", "C5AR1", "CFB", "CFD",
    ],
}


def compute_overlaps(df: pd.DataFrame, pathways: list[str] | None = None) -> pd.DataFrame:
    """Compute overlap between DEGs and hallmark gene sets.

    Returns a summary DataFrame with overlap counts and Fisher's exact test p-values.
    """
    if df.empty or "Gene" not in df.columns:
        return pd.DataFrame()

    deg_genes = set(df["Gene"].dropna().str.upper())
    selected = pathways if pathways else list(HALLMARK_SETS.keys())

    rows = []
    for pathway in selected:
        if pathway not in HALLMARK_SETS:
            continue
        pathway_genes = set(g.upper() for g in HALLMARK_SETS[pathway])
        overlap = deg_genes & pathway_genes
        overlap_genes = sorted(overlap)

        # Fisher's exact test
        # 2x2 table: [in_pathway_and_deg, in_pathway_not_deg]
        #             [not_pathway_and_deg, not_pathway_not_deg]
        a = len(overlap)
        b = len(pathway_genes - deg_genes)
        c = len(deg_genes - pathway_genes)
        # Background ~20,000 human genes
        d = 20000 - a - b - c
        _, pvalue = stats.fisher_exact([[a, b], [c, d]], alternative="greater")

        rows.append({
            "Pathway": pathway,
            "Pathway Size": len(pathway_genes),
            "DEG Overlap": a,
            "Overlap %": round(100 * a / len(pathway_genes), 1) if pathway_genes else 0,
            "Fisher p-value": pvalue,
            "Overlap Genes": ", ".join(overlap_genes[:20]),
        })

    return pd.DataFrame(rows)


def get_gene_pathway_matrix(df: pd.DataFrame, pathways: list[str] | None = None) -> pd.DataFrame:
    """Build a genes × pathways matrix for heatmap, with log2FC values."""
    if df.empty or "Gene" not in df.columns:
        return pd.DataFrame()

    selected = pathways if pathways else list(HALLMARK_SETS.keys())
    deg_genes = set(df["Gene"].dropna().str.upper())

    # Find all genes that appear in any selected pathway
    all_pathway_genes = set()
    for pw in selected:
        if pw in HALLMARK_SETS:
            all_pathway_genes.update(g.upper() for g in HALLMARK_SETS[pw])

    # Genes in both DEG list and at least one pathway
    relevant_genes = sorted(deg_genes & all_pathway_genes)

    if not relevant_genes:
        return pd.DataFrame()

    # Build matrix
    gene_fc = {}
    for _, row in df.iterrows():
        gene_fc[str(row["Gene"]).upper()] = row.get("log2FoldChange", 0)

    matrix_data = {}
    for pw in selected:
        if pw not in HALLMARK_SETS:
            continue
        pw_genes = set(g.upper() for g in HALLMARK_SETS[pw])
        matrix_data[pw] = [gene_fc.get(g, np.nan) if g in pw_genes else np.nan for g in relevant_genes]

    matrix = pd.DataFrame(matrix_data, index=relevant_genes)
    # Keep only rows with at least one non-NaN
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
    """Create heatmap of gene × pathway with log2FC coloring."""
    if matrix.empty:
        fig = go.Figure()
        fig.update_layout(title="No data for heatmap")
        return fig

    # Limit to top 50 genes for readability
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
        title="Gene × Pathway Heatmap (log2 Fold Change)",
        template="plotly_white",
        height=max(400, len(matrix) * 15),
        margin=dict(l=100),
    )

    return fig


def get_gene_table(df: pd.DataFrame, pathways: list[str] | None = None) -> pd.DataFrame:
    """Get table of DEGs with their pathway memberships."""
    if df.empty or "Gene" not in df.columns:
        return pd.DataFrame()

    selected = pathways if pathways else list(HALLMARK_SETS.keys())
    rows = []

    for _, row in df.iterrows():
        gene = str(row["Gene"]).upper()
        memberships = []
        for pw in selected:
            if pw in HALLMARK_SETS and gene in set(g.upper() for g in HALLMARK_SETS[pw]):
                memberships.append(pw)

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
