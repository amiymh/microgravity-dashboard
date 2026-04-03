"""Therapeutic targets identification using DGIdb and OpenTargets APIs."""

import time
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
                    "Interaction Type": interaction.get("interactionTypes", "N/A"),
                    "Source": ", ".join(interaction.get("sources", ["DGIdb"])),
                })

    return pd.DataFrame(all_rows) if all_rows else pd.DataFrame(
        columns=["Gene", "Drug Name", "Interaction Type", "Source"]
    )


def query_opentargets_drugs(ensembl_ids: list[str], batch_size: int = 10) -> pd.DataFrame:
    """Query OpenTargets for drug/clinical candidate info using Ensembl IDs.

    Only returns approval status and clinical phase when the API explicitly provides them.
    """
    all_rows = []

    for i in range(0, len(ensembl_ids), batch_size):
        batch = ensembl_ids[i : i + batch_size]

        for ens_id in batch:
            query = """
            query($id: String!) {
              target(ensemblId: $id) {
                approvedSymbol
                drugAndClinicalCandidates {
                  rows {
                    maxClinicalStage
                    drug { name drugType maximumClinicalStage }
                  }
                }
              }
            }
            """
            try:
                resp = requests.post(
                    OPENTARGETS_URL,
                    json={"query": query, "variables": {"id": ens_id}},
                    timeout=10,
                )
                resp.raise_for_status()
                data = resp.json()
            except Exception:
                continue

            target = data.get("data", {}).get("target")
            if not target:
                continue

            gene_symbol = target.get("approvedSymbol", "")
            candidates = target.get("drugAndClinicalCandidates", {})
            for row in candidates.get("rows", []):
                drug = row.get("drug", {})
                stage = row.get("maxClinicalStage", "")
                drug_stage = drug.get("maximumClinicalStage", "")

                all_rows.append({
                    "Gene": gene_symbol,
                    "Drug Name": drug.get("name", "Unknown"),
                    "Drug Type": drug.get("drugType", ""),
                    "Clinical Stage": stage or drug_stage or "",
                    "Source": "OpenTargets",
                })

        # Rate limit between batches
        if i + batch_size < len(ensembl_ids):
            time.sleep(0.2)

    return pd.DataFrame(all_rows) if all_rows else pd.DataFrame(
        columns=["Gene", "Drug Name", "Drug Type", "Clinical Stage", "Source"]
    )


def get_therapeutic_targets(
    df: pd.DataFrame,
    top_n: int = 200,
) -> pd.DataFrame:
    """Identify drug targets among top significant DEGs.

    Queries DGIdb for drug-gene interactions and OpenTargets for clinical stage data.
    Approval status is only shown when explicitly returned by an API.
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

    # Query DGIdb for drug interactions
    dgidb_df = query_dgidb(gene_list)

    # Query OpenTargets for clinical stage data (using Ensembl IDs)
    ot_df = pd.DataFrame()
    if "ID" in top_genes_df.columns:
        ensembl_ids = top_genes_df["ID"].dropna().unique().tolist()
        if ensembl_ids:
            ot_df = query_opentargets_drugs(ensembl_ids)

    # Merge DGIdb and OpenTargets results
    if dgidb_df.empty and ot_df.empty:
        return _demo_targets(top_genes_df)

    combined = pd.DataFrame()

    if not dgidb_df.empty:
        combined = dgidb_df.copy()

    if not ot_df.empty:
        if combined.empty:
            combined = ot_df.copy()
        else:
            # Add OpenTargets rows that aren't already in DGIdb results
            ot_keys = set(zip(ot_df["Gene"], ot_df["Drug Name"]))
            existing_keys = set(zip(combined["Gene"], combined["Drug Name"]))
            new_rows = ot_df[
                ot_df.apply(lambda r: (r["Gene"], r["Drug Name"]) not in existing_keys, axis=1)
            ]
            combined = pd.concat([combined, new_rows], ignore_index=True)

            # Back-fill Clinical Stage from OpenTargets into DGIdb rows
            ot_stage_map = {}
            for _, r in ot_df.iterrows():
                key = (r["Gene"], r["Drug Name"])
                if r.get("Clinical Stage"):
                    ot_stage_map[key] = r["Clinical Stage"]

            if "Clinical Stage" not in combined.columns:
                combined["Clinical Stage"] = ""
            for idx, r in combined.iterrows():
                key = (r["Gene"], r["Drug Name"])
                if not r.get("Clinical Stage") and key in ot_stage_map:
                    combined.at[idx, "Clinical Stage"] = ot_stage_map[key]

    # Merge with DEG info
    deg_info_cols = ["Gene", "Direction", "log2FoldChange", "padj", "significance score"]
    available = [c for c in deg_info_cols if c in df.columns]
    deg_info = df[available].drop_duplicates(subset=["Gene"])
    merged = combined.merge(deg_info, on="Gene", how="left")

    return merged


def _demo_targets(top_genes_df: pd.DataFrame) -> pd.DataFrame:
    """Generate demo targets when APIs are unavailable."""
    demo_drugs = [
        ("TNF", "Infliximab", "Antibody", "DGIdb"),
        ("TNF", "Adalimumab", "Antibody", "DGIdb"),
        ("IL6", "Tocilizumab", "Antibody", "DGIdb"),
        ("TP53", "APR-246", "Small molecule", "DGIdb"),
        ("BCL2", "Venetoclax", "Small molecule", "DGIdb"),
        ("CASP3", "Z-DEVD-FMK", "Small molecule", "DGIdb"),
        ("VEGFA", "Bevacizumab", "Antibody", "DGIdb"),
        ("EGFR", "Erlotinib", "Small molecule", "DGIdb"),
        ("MAPK1", "Trametinib", "Small molecule", "DGIdb"),
        ("CDK4", "Palbociclib", "Small molecule", "DGIdb"),
    ]

    rows = []
    available_genes = set(top_genes_df["Gene"].unique()) if "Gene" in top_genes_df.columns else set()

    for gene, drug, dtype, source in demo_drugs:
        if gene in available_genes or not available_genes:
            rows.append({
                "Gene": gene,
                "Drug Name": drug,
                "Drug Type": dtype,
                "Source": f"Demo ({source})",
            })

    if not rows:
        for gene in list(available_genes)[:5]:
            rows.append({
                "Gene": gene,
                "Drug Name": "Demo Drug",
                "Drug Type": "Small molecule",
                "Source": "Demo",
            })

    result = pd.DataFrame(rows)

    if not result.empty and not top_genes_df.empty:
        deg_cols = ["Gene", "Direction", "log2FoldChange", "padj", "significance score"]
        available = [c for c in deg_cols if c in top_genes_df.columns]
        deg_info = top_genes_df[available].drop_duplicates(subset=["Gene"])
        result = result.merge(deg_info, on="Gene", how="left")

    return result


def get_target_summary(targets_df: pd.DataFrame) -> dict:
    """Compute summary metrics for the targets tab."""
    if targets_df.empty:
        return {"genes_with_drugs": 0, "total_drugs": 0, "with_clinical_stage": 0}

    with_stage = 0
    if "Clinical Stage" in targets_df.columns:
        with_stage = int((targets_df["Clinical Stage"].fillna("") != "").sum())

    return {
        "genes_with_drugs": targets_df["Gene"].nunique(),
        "total_drugs": targets_df["Drug Name"].nunique() if "Drug Name" in targets_df.columns else 0,
        "with_clinical_stage": with_stage,
    }


def get_analysis_log(
    genes_queried: int = 0,
    targets_df: pd.DataFrame | None = None,
) -> list[str]:
    """Return analysis log for therapeutic targets."""
    from datetime import datetime
    log = [
        f"Genes queried: {genes_queried} (top by significance score)",
        f"DGIdb API endpoint: {DGIDB_URL}",
        f"DGIdb batch size: 20 genes per request",
        f"OpenTargets API endpoint: {OPENTARGETS_URL}",
        f"Timestamp: {datetime.now().isoformat()}",
    ]
    if targets_df is not None and not targets_df.empty:
        n_genes = targets_df["Gene"].nunique()
        n_drugs = targets_df["Drug Name"].nunique() if "Drug Name" in targets_df.columns else 0
        sources = targets_df["Source"].unique().tolist() if "Source" in targets_df.columns else []
        log.append(f"Genes with drug interactions found: {n_genes}")
        log.append(f"Unique drugs found: {n_drugs}")
        log.append(f"Data sources: {', '.join(str(s) for s in sources)}")
    else:
        log.append("No drug interactions found.")
    log.append("Clinical Stage shown only when explicitly returned by OpenTargets API.")
    return log
