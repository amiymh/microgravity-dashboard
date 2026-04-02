"""Tests for data_loader module."""

import os
import sys
import pytest
import pandas as pd
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from modules.data_loader import (
    generate_demo_data,
    load_excel,
    validate_columns,
    filter_degs,
    REQUIRED_COLS,
    DEFAULT_PATH,
)


class TestGenerateDemoData:
    def test_returns_dataframe(self):
        df = generate_demo_data()
        assert isinstance(df, pd.DataFrame)

    def test_default_50_genes(self):
        df = generate_demo_data()
        assert len(df) == 50

    def test_custom_size(self):
        df = generate_demo_data(n=20)
        assert len(df) == 20

    def test_has_required_columns(self):
        df = generate_demo_data()
        missing = validate_columns(df)
        assert missing == [], f"Missing columns: {missing}"

    def test_direction_values(self):
        df = generate_demo_data()
        assert set(df["Direction"].unique()).issubset({"Upregulated", "Downregulated"})

    def test_no_nan_in_key_fields(self):
        df = generate_demo_data()
        for col in ["Gene", "padj", "log2FoldChange", "Direction"]:
            assert df[col].isna().sum() == 0, f"NaN found in {col}"


class TestLoadExcel:
    def test_load_actual_file_sdegs(self):
        if not os.path.exists(DEFAULT_PATH):
            pytest.skip("Data file not available")
        df = load_excel(DEFAULT_PATH, sheet_name="SDEGs")
        assert len(df) > 4000
        assert "Gene" in df.columns
        assert "Direction" in df.columns

    def test_load_actual_file_all_genes(self):
        if not os.path.exists(DEFAULT_PATH):
            pytest.skip("Data file not available")
        df = load_excel(DEFAULT_PATH, sheet_name="SIGNFICANCE SCOREs Gene Types")
        assert len(df) > 14000

    def test_fallback_to_demo_on_missing_file(self):
        df = load_excel("/nonexistent/path.xlsx")
        assert isinstance(df, pd.DataFrame)
        assert len(df) == 50  # demo data

    def test_no_extra_empty_columns(self):
        if not os.path.exists(DEFAULT_PATH):
            pytest.skip("Data file not available")
        df = load_excel(DEFAULT_PATH, sheet_name="SDEGs")
        # Should not have hundreds of None/Column### columns
        assert len(df.columns) < 30


class TestValidateColumns:
    def test_valid_dataframe(self):
        df = generate_demo_data()
        assert validate_columns(df) == []

    def test_missing_columns(self):
        df = pd.DataFrame({"Gene": ["A"], "padj": [0.01]})
        missing = validate_columns(df)
        assert "log2FoldChange" in missing
        assert "Direction" in missing


class TestFilterDegs:
    def test_padj_filter(self):
        df = generate_demo_data()
        filtered = filter_degs(df, padj_threshold=0.01)
        assert all(filtered["padj"] <= 0.01)

    def test_log2fc_filter(self):
        df = generate_demo_data()
        filtered = filter_degs(df, log2fc_threshold=2.0, padj_threshold=1.0)
        assert all(filtered["log2FoldChange"].abs() >= 2.0)

    def test_direction_filter(self):
        df = generate_demo_data()
        up = filter_degs(df, direction="Upregulated", padj_threshold=1.0, log2fc_threshold=0)
        assert all(up["Direction"] == "Upregulated")

    def test_gene_type_filter(self):
        df = generate_demo_data()
        filtered = filter_degs(
            df, gene_types=["protein_coding"], padj_threshold=1.0, log2fc_threshold=0
        )
        assert all(filtered["Gene Type"] == "protein_coding")

    def test_empty_result(self):
        df = generate_demo_data()
        filtered = filter_degs(df, padj_threshold=1e-300)
        assert isinstance(filtered, pd.DataFrame)
        assert len(filtered) == 0

    def test_empty_input(self):
        df = pd.DataFrame(columns=["Gene", "padj", "log2FoldChange", "Direction", "Gene Type"])
        filtered = filter_degs(df)
        assert len(filtered) == 0
