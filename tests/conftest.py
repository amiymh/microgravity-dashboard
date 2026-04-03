"""Shared test fixtures — real data is the primary test source."""

import os
import sys
import pytest
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

SEARCH_PATHS = [
    os.path.join(os.path.dirname(__file__), "..", "data", "test_data.xlsx"),
    os.path.expanduser("~/Downloads/Supplementary Data 1 DESEQ2 normalized out and significance scores.xlsx"),
    os.path.expanduser("~/Downloads/Supplementary-Data-1-DESEQ2-normalized-out-and-significance-scores.xlsx"),
    "Supplementary-Data-1-DESEQ2-normalized-out-and-significance-scores.xlsx",
]

DATA_FILE = next((p for p in SEARCH_PATHS if os.path.exists(p)), None)


@pytest.fixture(scope="session")
def real_df():
    """Load the real DESeq2 SDEGs sheet. Skip if file not found."""
    if DATA_FILE is None:
        pytest.skip("Real data file not available in any search path")
    from modules.data_loader import load_excel
    return load_excel(DATA_FILE, sheet_name="SDEGs")


@pytest.fixture(scope="session")
def real_df_filtered(real_df):
    """Real data with default filters applied."""
    from modules.data_loader import filter_degs
    return filter_degs(real_df, padj_threshold=0.05, log2fc_threshold=1.0)


@pytest.fixture(scope="session")
def demo_df():
    """50-gene synthetic dataset for edge-case tests."""
    from modules.data_loader import generate_demo_data
    return generate_demo_data(50)
