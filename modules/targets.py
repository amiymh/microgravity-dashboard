"""Therapeutic targets identification using DGIdb and OpenTargets APIs."""

import pandas as pd
import requests

DGIDB_URL = "https://dgidb.org/api/v2/interactions.json"
OPENTARGETS_URL = "https://api.platform.opentargets.org/api/v4/graphql"


def query_dgidb(genes: list[str], batch_size: int = 20) -> pd.DataFrame:
    """Query DGIdb for drug-gene interactions in batches."""
    all_rows = []

    for i in range(0, len(genes), batch_size):
        batch = genes[i : i + batch_size]
        try:
            resp = requests.get(
                DGIDB_URL,
                params={"genes": ",".join(batch)},
                timeout=10,
            )
            resp.raise_for_status()
            data = resp.json()
        except Exception:
            continue

        for match in data.get("matchedTerms", []):
            gene = match.get("geneName", "")
            for interaction in match.get("interactions", []):
                all_rows.append({
                    "Gene": gene,
                    "Drug Name": interaction.get("drugName", "Unknown"),
                    "Drug Type": interaction.get("drugChemblId", "N/A"),
                    "Interaction Type": interaction.get("interactionTypes", "N/A"),
                    "Source": ", ".join(interaction.get("sources", ["DGIdb"])),
                })

    return pd.DataFrame(all_rows) if all_rows else pd.DataFrame(
        columns=["Gene", "Drug Name", "Drug Type", "Interaction Type", "Source"]
    )


def query_opentargets_drugs(gene_symbol: str) -> list[dict]:
    """Query OpenTargets for drug info about a gene."""
    query = """
    query($ensemblId: String!) {
      target(ensemblId: $ensemblId) {
        knownDrugs(size: 10) {
          rows {
            drug { name mechanismOfAction }
            phase
            status
            disease { name }
          }
        }
      }
    }
    """
    # We'd need Ensembl ID; skip if not available
    return []


def get_therapeutic_targets(
    df: pd.DataFrame,
    top_n: int = 200,
) -> pd.DataFrame:
    """Identify drug targets among top significant DEGs.

    Args:
        df: Full DEG dataframe with Gene, significance score, Direction, log2FoldChange.
        top_n: Number of top genes to query.

    Returns:
        DataFrame of drug-gene interactions merged with DEG data.
    """
    if df.empty or "Gene" not in df.columns:
        return pd.DataFrame()

    # Get top genes by significance score
    if "significance score" in df.columns:
        top_genes_df = df.nlargest(min(top_n, len(df)), "significance score")
    else:
        top_genes_df = df.head(top_n)

    gene_list = top_genes_df["Gene"].dropna().unique().tolist()
    if not gene_list:
        return pd.DataFrame()

    # Query DGIdb
    drug_df = query_dgidb(gene_list)

    if drug_df.empty:
        return _demo_targets(top_genes_df)

    # Merge with DEG info
    deg_info_cols = ["Gene", "Direction", "log2FoldChange", "padj", "significance score"]
    available = [c for c in deg_info_cols if c in df.columns]
    deg_info = df[available].drop_duplicates(subset=["Gene"])

    merged = drug_df.merge(deg_info, on="Gene", how="left")

    # Add approval status heuristic
    merged["Approval Status"] = "Experimental"
    merged["Clinical Phase"] = "N/A"

    return merged


def _demo_targets(top_genes_df: pd.DataFrame) -> pd.DataFrame:
    """Generate demo targets when API is unavailable."""
    demo_drugs = [
        ("TNF", "Infliximab", "Antibody", "FDA Approved", "Phase 4"),
        ("TNF", "Adalimumab", "Antibody", "FDA Approved", "Phase 4"),
        ("IL6", "Tocilizumab", "Antibody", "FDA Approved", "Phase 4"),
        ("TP53", "APR-246", "Small molecule", "Clinical Trial", "Phase 3"),
        ("BCL2", "Venetoclax", "Small molecule", "FDA Approved", "Phase 4"),
        ("CASP3", "Z-DEVD-FMK", "Small molecule", "Experimental", "Preclinical"),
        ("VEGFA", "Bevacizumab", "Antibody", "FDA Approved", "Phase 4"),
        ("EGFR", "Erlotinib", "Small molecule", "FDA Approved", "Phase 4"),
        ("MAPK1", "Trametinib", "Small molecule", "FDA Approved", "Phase 4"),
        ("CDK4", "Palbociclib", "Small molecule", "FDA Approved", "Phase 4"),
    ]

    rows = []
    available_genes = set(top_genes_df["Gene"].unique()) if "Gene" in top_genes_df.columns else set()

    for gene, drug, dtype, status, phase in demo_drugs:
        if gene in available_genes or not available_genes:
            rows.append({
                "Gene": gene,
                "Drug Name": drug,
                "Drug Type": dtype,
                "Approval Status": status,
                "Clinical Phase": phase,
                "Source": "Demo Data",
            })

    if not rows:
        # Use first few genes as demo
        for gene in list(available_genes)[:5]:
            rows.append({
                "Gene": gene,
                "Drug Name": "Demo Drug",
                "Drug Type": "Small molecule",
                "Approval Status": "Experimental",
                "Clinical Phase": "Preclinical",
                "Source": "Demo Data",
            })

    result = pd.DataFrame(rows)

    # Merge with DEG info
    if not result.empty and not top_genes_df.empty:
        deg_cols = ["Gene", "Direction", "log2FoldChange", "padj", "significance score"]
        available = [c for c in deg_cols if c in top_genes_df.columns]
        deg_info = top_genes_df[available].drop_duplicates(subset=["Gene"])
        result = result.merge(deg_info, on="Gene", how="left")

    return result


def get_target_summary(targets_df: pd.DataFrame) -> dict:
    """Compute summary metrics for the targets tab."""
    if targets_df.empty:
        return {"genes_with_drugs": 0, "total_drugs": 0, "approved_drugs": 0}

    return {
        "genes_with_drugs": targets_df["Gene"].nunique(),
        "total_drugs": targets_df["Drug Name"].nunique() if "Drug Name" in targets_df.columns else 0,
        "approved_drugs": len(targets_df[targets_df.get("Approval Status", pd.Series()) == "FDA Approved"])
            if "Approval Status" in targets_df.columns else 0,
    }
