"""Organ toxicity gene signature analysis.

Cross-references user DEGs against curated organ-specific toxicity gene sets
derived from ToxicoDB / TG-GATEs literature and MSigDB hallmark pathways.
Uses Fisher's exact test to determine enrichment.
"""

import pandas as pd
import numpy as np
from scipy import stats
import plotly.graph_objects as go
from datetime import datetime

# ── Organ-specific toxicity gene sets ────────────────────────────────────────
# Curated from TG-GATEs, DrugMatrix, ToxPanel (Te et al. 2016), and literature.
# Each set represents genes consistently associated with toxicity in that organ.

ORGAN_GENE_SETS: dict[str, list[str]] = {
    "Liver (Hepatotoxicity)": [
        "ALB", "CYP1A2", "CYP2E1", "CYP3A4", "CYP2B6", "CYP7A1",
        "HMGCR", "FASN", "SCD", "ACOX1", "CPT1A", "PPARA", "PPARG",
        "NR1I2", "NR1I3", "ABCB1", "ABCC2", "UGT1A1", "GSTA1", "GSTP1",
        "EPHX1", "NQO1", "HMOX1", "FMO3", "ADH1B", "ALDH2", "GGT1",
        "AFP", "GPC3", "SERPINA1", "HP", "TF", "TTR", "APOB", "APOA1",
        "FGA", "FGB", "FGG", "ASGR1", "HNF4A", "CEBPA", "FOS", "JUN",
        "MYC", "CCND1", "PCNA", "KI67", "TGFB1", "COL1A1", "ACTA2",
        "TIMP1", "MMP2", "MMP9", "IL6", "TNF", "IL1B", "CCL2",
    ],
    "Kidney (Nephrotoxicity)": [
        "KIM1", "HAVCR1", "LCN2", "NGAL", "CLU", "CYSTB", "CSTA",
        "SPP1", "TGFB1", "COL1A1", "FN1", "VIM", "ACTA2", "CD44",
        "AQP1", "AQP2", "SLC12A1", "SLC22A6", "SLC22A8", "ABCB1",
        "CYP24A1", "CYP27B1", "CYP3A5", "UMOD", "EGF", "VEGFA",
        "NPHS1", "NPHS2", "WT1", "PAX2", "SIX2", "HNF1B", "CTGF",
        "PDGFB", "MMP7", "TIMP1", "IL18", "CXCL1", "CCL2", "TNF",
        "IL6", "CASP3", "BAX", "BCL2", "SOD2", "CAT", "GPX1",
        "HMOX1", "NQO1", "NFE2L2", "KEAP1", "ATF3", "DDIT3",
    ],
    "Heart (Cardiotoxicity)": [
        "TNNT2", "TNNI3", "MYH7", "MYH6", "ACTC1", "MYL2", "MYL3",
        "MYBPC3", "TTN", "RYR2", "SCN5A", "KCNH2", "KCNQ1", "KCNJ2",
        "CACNA1C", "GJA1", "NPPA", "NPPB", "BNP", "ANP",
        "ACE", "AGT", "AGTR1", "ADRB1", "ADRB2", "EDN1",
        "VEGFA", "HIF1A", "NOS3", "PTGS2", "COX2",
        "CASP3", "BAX", "BCL2", "FAS", "TNF", "IL6", "IL1B",
        "TGFB1", "COL1A1", "COL3A1", "MMP2", "MMP9", "TIMP1",
        "SOD2", "CAT", "GPX1", "HMOX1", "NQO1", "CKM", "CKMB",
        "LDH", "LDHA", "LDHB", "GATA4", "NKX2-5", "TBX5",
    ],
    "Lung (Pulmonary Toxicity)": [
        "SFTPA1", "SFTPB", "SFTPC", "SFTPD", "SCGB1A1", "MUC5AC",
        "MUC5B", "FOXJ1", "NKX2-1", "SOX2", "TP63", "KRT5",
        "AGER", "CAV1", "AQP5", "PDPN", "VEGFA", "HIF1A", "EPAS1",
        "TGFB1", "COL1A1", "COL3A1", "FN1", "ACTA2", "MMP2",
        "MMP7", "MMP9", "TIMP1", "SERPINE1", "CTGF", "WNT3A",
        "IL6", "IL1B", "TNF", "CXCL8", "CCL2", "CCL5", "CXCL1",
        "TLR4", "NLRP3", "IL33", "TSLP", "IL13", "IL4", "IL5",
        "HMOX1", "NQO1", "SOD2", "GPX1", "CAT", "NFE2L2",
        "CASP3", "BAX", "BCL2", "FAS", "FASLG", "GADD45A",
    ],
    "Brain (Neurotoxicity)": [
        "GFAP", "ENO2", "NSE", "MAP2", "TUBB3", "SYN1", "SYP",
        "RBFOX3", "NEUN", "OLIG2", "MBP", "PLP1", "MOG",
        "SLC1A2", "SLC1A3", "SLC17A7", "SLC32A1", "GRIA1", "GRIN1",
        "GRIN2A", "GRIN2B", "GAD1", "GAD2", "TH", "SLC6A3", "SLC6A4",
        "BDNF", "NGF", "NTRK2", "APP", "PSEN1", "MAPT", "SNCA",
        "PARK2", "PINK1", "SOD1", "ALS2", "TARDBP", "FUS",
        "TNF", "IL1B", "IL6", "CCL2", "CXCL10", "TREM2", "CD68",
        "AIF1", "CASP3", "BAX", "BCL2", "HMOX1", "NQO1", "SOD2",
        "CAT", "GPX1", "NFE2L2", "NOS1", "NOS2", "PTGS2",
    ],
    "Bone Marrow (Hematotoxicity)": [
        "EPO", "EPOR", "KIT", "KITLG", "CSF3", "CSF3R", "CSF2",
        "TPO", "MPL", "GATA1", "GATA2", "TAL1", "SPI1", "PU1",
        "RUNX1", "CEBPA", "KLF1", "GFI1", "IKZF1", "PAX5",
        "HBA1", "HBB", "ALAS2", "GYPA", "CD34", "CD38", "CD71",
        "FLT3", "JAK2", "STAT5A", "STAT5B", "BCL2", "MCL1",
        "BAX", "CASP3", "CASP9", "TP53", "MDM2", "CDKN1A",
        "CCND1", "CDK4", "CDK6", "RB1", "E2F1", "MYC",
        "IL3", "IL6", "TNF", "IFNG", "TGFB1", "SPP1",
    ],
    "GI Tract (Gastrointestinal)": [
        "LGR5", "OLFM4", "ASCL2", "SOX9", "CDX2", "KRT20",
        "MUC2", "TFF3", "CHGA", "LYZ", "DEFA5", "DEFA6",
        "REG3A", "VIL1", "ALPI", "SLC5A1", "SLC2A2",
        "CLDN1", "CLDN2", "CLDN3", "OCLN", "TJP1", "CDH1",
        "TGFB1", "SMAD3", "SMAD4", "WNT3A", "CTNNB1", "APC",
        "PTGS2", "COX2", "IL6", "IL1B", "TNF", "CXCL8", "CCL2",
        "TLR4", "NOD2", "NLRP3", "IL10", "IL22", "IL17A",
        "CASP3", "BAX", "BCL2", "FAS", "PCNA", "KI67",
    ],
    "Muscle (Myotoxicity)": [
        "CKM", "CKMB", "MYH1", "MYH2", "MYH4", "MYH7",
        "ACTA1", "ACTC1", "DES", "DMD", "SGCA", "SGCB",
        "CAPN3", "TRIM63", "MURF1", "FBXO32", "ATROGIN1",
        "MSTN", "FST", "IGF1", "IGF1R", "MTOR", "AKT1",
        "FOXO1", "FOXO3", "MYOD1", "MYOG", "PAX7", "MEF2C",
        "PPARGC1A", "TFAM", "SOD2", "CAT", "GPX1",
        "IL6", "TNF", "IL1B", "TGFB1", "CTGF",
        "CASP3", "BAX", "BCL2", "BNIP3", "LAMP2",
        "MAP1LC3B", "SQSTM1", "BECN1", "ATG5", "ATG7",
    ],
}


def get_available_organs() -> list[str]:
    """Return list of available organ names."""
    return list(ORGAN_GENE_SETS.keys())


def compute_organ_overlaps(
    df: pd.DataFrame,
    organs: list[str] | None = None,
) -> pd.DataFrame:
    """Compute overlap between DEGs and organ toxicity gene sets.

    Returns a summary DataFrame with overlap counts, percentages,
    and Fisher's exact test p-values.
    """
    if df.empty or "Gene" not in df.columns:
        return pd.DataFrame()

    selected = organs if organs else list(ORGAN_GENE_SETS.keys())
    deg_genes = set(df["Gene"].dropna().str.upper())

    rows = []
    for organ in selected:
        if organ not in ORGAN_GENE_SETS:
            continue
        organ_genes = set(g.upper() for g in ORGAN_GENE_SETS[organ])
        overlap = deg_genes & organ_genes
        overlap_sorted = sorted(overlap)

        a = len(overlap)
        b = len(organ_genes - deg_genes)
        c = len(deg_genes - organ_genes)
        d = 20000 - a - b - c
        _, pvalue = stats.fisher_exact([[a, b], [c, d]], alternative="greater")

        rows.append({
            "Organ": organ,
            "Signature Size": len(organ_genes),
            "DEG Overlap": a,
            "Overlap %": round(100 * a / len(organ_genes), 1) if organ_genes else 0,
            "Fisher p-value": pvalue,
            "Significant": "Yes" if pvalue < 0.05 else "No",
            "Overlap Genes": ", ".join(overlap_sorted[:30]),
        })

    result = pd.DataFrame(rows)
    if not result.empty:
        result = result.sort_values("Fisher p-value")
    return result.reset_index(drop=True)


def get_organ_gene_table(
    df: pd.DataFrame,
    organs: list[str] | None = None,
) -> pd.DataFrame:
    """Get table of DEGs with their organ toxicity associations."""
    if df.empty or "Gene" not in df.columns:
        return pd.DataFrame()

    selected = organs if organs else list(ORGAN_GENE_SETS.keys())
    organ_upper = {
        organ: set(g.upper() for g in ORGAN_GENE_SETS[organ])
        for organ in selected if organ in ORGAN_GENE_SETS
    }

    rows = []
    for _, row in df.iterrows():
        gene = str(row["Gene"]).upper()
        memberships = [organ for organ, gene_set in organ_upper.items() if gene in gene_set]
        if memberships:
            rows.append({
                "Gene": row["Gene"],
                "log2FoldChange": row.get("log2FoldChange", np.nan),
                "padj": row.get("padj", np.nan),
                "Direction": row.get("Direction", ""),
                "Organs": "; ".join(memberships),
                "Organ Count": len(memberships),
            })

    result = pd.DataFrame(rows)
    if not result.empty:
        result = result.sort_values("Organ Count", ascending=False)
    return result.reset_index(drop=True)


def create_organ_bar_chart(overlap_df: pd.DataFrame) -> go.Figure:
    """Create horizontal bar chart showing organ toxicity enrichment."""
    if overlap_df.empty:
        fig = go.Figure()
        fig.update_layout(title="No organ toxicity overlaps found")
        return fig

    plot_df = overlap_df.sort_values("Fisher p-value", ascending=False).copy()
    plot_df["neg_log10_p"] = -np.log10(plot_df["Fisher p-value"].clip(lower=1e-300))

    colors = []
    for _, row in plot_df.iterrows():
        if row["Fisher p-value"] < 0.001:
            colors.append("#d62728")  # strong signal
        elif row["Fisher p-value"] < 0.05:
            colors.append("#ff7f0e")  # significant
        else:
            colors.append("#aec7e8")  # not significant

    fig = go.Figure(
        go.Bar(
            y=plot_df["Organ"],
            x=plot_df["neg_log10_p"],
            orientation="h",
            marker=dict(color=colors),
            text=plot_df["DEG Overlap"].astype(str) + " genes",
            textposition="auto",
            hovertemplate=(
                "<b>%{y}</b><br>"
                "-log10(p): %{x:.2f}<br>"
                "Overlap: %{text}<br>"
                "<extra></extra>"
            ),
        )
    )

    # Add significance line
    fig.add_vline(
        x=-np.log10(0.05), line_dash="dash", line_color="gray",
        annotation_text="p=0.05", annotation_position="top",
    )

    fig.update_layout(
        title="Organ Toxicity Signature Enrichment",
        xaxis_title="-log10(Fisher p-value)",
        template="plotly_white",
        height=max(350, len(plot_df) * 55),
        margin=dict(l=220),
    )

    return fig


def create_organ_heatmap(df: pd.DataFrame, organs: list[str] | None = None) -> go.Figure:
    """Create gene x organ heatmap with log2FC coloring."""
    if df.empty or "Gene" not in df.columns:
        fig = go.Figure()
        fig.update_layout(title="No data for heatmap")
        return fig

    selected = organs if organs else list(ORGAN_GENE_SETS.keys())
    deg_genes = set(df["Gene"].dropna().str.upper())

    # Find genes that appear in at least one organ set
    all_organ_genes = set()
    for organ in selected:
        if organ in ORGAN_GENE_SETS:
            all_organ_genes.update(g.upper() for g in ORGAN_GENE_SETS[organ])

    relevant = sorted(deg_genes & all_organ_genes)
    if not relevant:
        fig = go.Figure()
        fig.update_layout(title="No overlapping genes found")
        return fig

    # Cap at 60 genes for readability
    if len(relevant) > 60:
        relevant = relevant[:60]

    gene_fc = {}
    for _, row in df.iterrows():
        gene_fc[str(row["Gene"]).upper()] = row.get("log2FoldChange", 0)

    matrix_data = {}
    for organ in selected:
        if organ not in ORGAN_GENE_SETS:
            continue
        organ_genes_upper = set(g.upper() for g in ORGAN_GENE_SETS[organ])
        matrix_data[organ.split(" (")[0]] = [
            gene_fc.get(g, np.nan) if g in organ_genes_upper else np.nan
            for g in relevant
        ]

    matrix = pd.DataFrame(matrix_data, index=relevant)
    matrix = matrix.dropna(how="all")

    if matrix.empty:
        fig = go.Figure()
        fig.update_layout(title="No data for heatmap")
        return fig

    fig = go.Figure(
        go.Heatmap(
            z=matrix.values,
            x=matrix.columns.tolist(),
            y=matrix.index.tolist(),
            colorscale="RdBu_r",
            zmid=0,
            colorbar=dict(title="log2FC"),
            hovertemplate="Gene: %{y}<br>Organ: %{x}<br>log2FC: %{z:.2f}<extra></extra>",
        )
    )

    fig.update_layout(
        title="Gene × Organ Toxicity Heatmap (log2 Fold Change)",
        template="plotly_white",
        height=max(400, len(matrix) * 14),
        margin=dict(l=100, b=150),
        xaxis=dict(tickangle=-45),
    )

    return fig


def get_analysis_log(
    overlap_df: pd.DataFrame | None = None,
    organs_used: list[str] | None = None,
) -> list[str]:
    """Return analysis log for organ toxicity analysis."""
    log = [
        f"Timestamp: {datetime.now().isoformat()}",
        "Gene sets: Curated organ-specific toxicity signatures",
        "Sources: TG-GATEs, DrugMatrix, ToxPanel (Te et al. 2016), MSigDB, literature",
        "Statistical test: Fisher's exact test (one-sided, greater)",
        "Background gene universe: 20,000 human genes (approximate)",
        "Significance threshold: p < 0.05",
    ]
    selected = organs_used or list(ORGAN_GENE_SETS.keys())
    for organ in selected:
        if organ in ORGAN_GENE_SETS:
            log.append(f"Gene set '{organ}': {len(ORGAN_GENE_SETS[organ])} genes")
    if overlap_df is not None and not overlap_df.empty:
        sig = overlap_df[overlap_df["Significant"] == "Yes"]
        log.append(f"Organs with significant enrichment: {len(sig)}/{len(overlap_df)}")
        for _, row in overlap_df.iterrows():
            status = "**SIGNIFICANT**" if row["Significant"] == "Yes" else "not significant"
            log.append(
                f"  {row['Organ']}: overlap={row['DEG Overlap']}, "
                f"p={row['Fisher p-value']:.2e} ({status})"
            )
            if row.get("Overlap Genes"):
                genes_preview = row["Overlap Genes"]
                if len(genes_preview) > 100:
                    genes_preview = genes_preview[:100] + "..."
                log.append(f"    Genes: {genes_preview}")
    return log
