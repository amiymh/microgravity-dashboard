"""Data loading, validation, and demo dataset generation."""

import pandas as pd
import numpy as np
import os

# Columns we expect in the SDEGs sheet (first 28 real columns)
REQUIRED_COLS = [
    "Gene", "ID", "padj", "log2FoldChange", "pvalue", "foldChange",
    "significance score", "Direction", "Gene Type",
]

EARTH_COLS = [
    "Earth_3D2", "Earth_4D2", "Earth_3E2", "Earth_4E2",
    "Earth_1D2", "Earth_2D2", "Earth_1E2", "Earth_2E2",
]

SPACE_COLS = [
    "Space_3D2", "Space_4D2", "Space_3E2", "Space_4E2",
    "Space_1D2", "Space_2D2", "Space_1E2", "Space_2E2",
]

# All meaningful columns in order
KEEP_COLS = [
    "Gene", "ID",
    *EARTH_COLS, *SPACE_COLS,
    "padj", "log2FoldChange", "pvalue", "stat", "foldChange",
    "log10padj", "significance score", "log2", "Direction", "Gene Type",
]

DEFAULT_PATH = os.path.expanduser(
    "~/Downloads/Supplementary Data 1 DESEQ2 normalized out and significance scores.xlsx"
)


def generate_demo_data(n: int = 50) -> pd.DataFrame:
    """Generate a synthetic demo dataset of n genes."""
    rng = np.random.default_rng(42)
    genes = [f"GENE{i}" for i in range(1, n + 1)]
    ids = [f"ENSG{str(i).zfill(11)}" for i in range(1, n + 1)]

    log2fc = rng.normal(0, 2, n)
    pvalues = 10 ** rng.uniform(-10, -1, n)
    padj_vals = np.clip(pvalues * n / np.arange(1, n + 1), 0, 1)

    earth_data = {col: rng.integers(10, 1000, n) for col in EARTH_COLS}
    space_data = {col: rng.integers(10, 1000, n) for col in SPACE_COLS}

    direction = np.where(log2fc > 0, "Upregulated", "Downregulated")
    gene_types = rng.choice(
        ["protein_coding", "lncRNA", "processed_pseudogene"], n, p=[0.7, 0.2, 0.1]
    )

    log10padj = -np.log10(np.clip(padj_vals, 1e-300, 1))
    sig_score = np.abs(log2fc) * log10padj

    df = pd.DataFrame(
        {
            "Gene": genes,
            "ID": ids,
            **earth_data,
            **space_data,
            "padj": padj_vals,
            "log2FoldChange": log2fc,
            "pvalue": pvalues,
            "stat": rng.normal(0, 5, n),
            "foldChange": 2.0 ** log2fc,
            "log10padj": log10padj,
            "significance score": sig_score,
            "log2": log2fc,
            "Direction": direction,
            "Gene Type": gene_types,
        }
    )
    return df


def _get_name(source) -> str:
    """Get filename from a path string or file-like object."""
    if isinstance(source, str):
        return source
    if hasattr(source, "name"):
        return source.name
    return ""


def _is_csv(source) -> bool:
    """Check if the source is a CSV file (by name or path)."""
    return _get_name(source).lower().endswith(".csv")


def _has_extension(source) -> bool:
    """Check if the source has a recognizable file extension."""
    name = _get_name(source).lower()
    return name.endswith(".csv") or name.endswith(".xlsx") or name.endswith(".xls")


def load_excel(
    file_path: str | None = None,
    sheet_name: str = "SDEGs",
    uploaded_file=None,
) -> pd.DataFrame:
    """Load DESeq2 data from Excel or CSV.

    Args:
        file_path: Path to xlsx/csv file on disk.
        sheet_name: Sheet to load (ignored for CSV files).
        uploaded_file: Streamlit UploadedFile object (takes priority over file_path).

    Returns:
        Cleaned DataFrame with only meaningful columns.
    """
    source = uploaded_file if uploaded_file is not None else file_path

    if source is None:
        source = DEFAULT_PATH

    try:
        if _is_csv(source):
            df = pd.read_csv(source)
        elif _has_extension(source):
            df = pd.read_excel(source, sheet_name=sheet_name, engine="openpyxl")
        else:
            # No detectable extension — try Excel first, then CSV
            try:
                df = pd.read_excel(source, sheet_name=sheet_name, engine="openpyxl")
            except Exception:
                if hasattr(source, "seek"):
                    source.seek(0)
                df = pd.read_csv(source)
    except FileNotFoundError:
        return generate_demo_data()
    except Exception:
        return generate_demo_data()

    # Clean: keep only meaningful columns that exist
    cols_to_keep = [c for c in KEEP_COLS if c in df.columns]
    df = df[cols_to_keep].copy()

    # Handle the secondary sheet which uses 'Type' instead of 'Gene Type'
    if "Type" in df.columns and "Gene Type" not in df.columns:
        df = df.rename(columns={"Type": "Gene Type"})

    # Drop rows where Gene is null
    df = df.dropna(subset=["Gene"])

    # Ensure numeric columns
    for col in ["padj", "log2FoldChange", "pvalue", "foldChange", "significance score"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    # Fill significance score if missing
    if "significance score" not in df.columns or df["significance score"].isna().all():
        if "log2FoldChange" in df.columns and "padj" in df.columns:
            df["significance score"] = (
                df["log2FoldChange"].abs() * (-np.log10(df["padj"].clip(lower=1e-300)))
            )

    # Fill Direction if missing
    if "Direction" not in df.columns or df["Direction"].isna().all():
        if "log2FoldChange" in df.columns:
            df["Direction"] = np.where(
                df["log2FoldChange"] > 0, "Upregulated", "Downregulated"
            )

    return df


def validate_columns(df: pd.DataFrame) -> list[str]:
    """Return list of missing required columns."""
    return [c for c in REQUIRED_COLS if c not in df.columns]


def filter_degs(
    df: pd.DataFrame,
    padj_threshold: float = 0.05,
    log2fc_threshold: float = 1.0,
    gene_types: list[str] | None = None,
    direction: str = "All",
) -> pd.DataFrame:
    """Apply user-selected filters to the DEG dataframe."""
    filtered = df.copy()

    if "padj" in filtered.columns:
        filtered = filtered[filtered["padj"] <= padj_threshold]

    if "log2FoldChange" in filtered.columns:
        filtered = filtered[filtered["log2FoldChange"].abs() >= log2fc_threshold]

    if gene_types and "Gene Type" in filtered.columns:
        filtered = filtered[filtered["Gene Type"].isin(gene_types)]

    if direction != "All" and "Direction" in filtered.columns:
        filtered = filtered[filtered["Direction"] == direction]

    return filtered
