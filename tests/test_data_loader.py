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


class TestLoadCSV:
    def test_load_csv_with_correct_columns(self, tmp_path):
        """A CSV with the right columns should load and clean correctly."""
        csv_file = tmp_path / "test_data.csv"
        demo = generate_demo_data(n=10)
        demo.to_csv(csv_file, index=False)

        df = load_excel(str(csv_file))
        assert isinstance(df, pd.DataFrame)
        assert len(df) == 10
        assert "Gene" in df.columns
        assert "Direction" in df.columns

    def test_csv_missing_columns_still_loads(self, tmp_path):
        """A CSV with some missing columns should load without crashing."""
        csv_file = tmp_path / "partial.csv"
        pd.DataFrame({
            "Gene": ["TP53", "TNF", "IL6"],
            "log2FoldChange": [2.5, -1.3, 0.8],
            "padj": [0.001, 0.01, 0.05],
        }).to_csv(csv_file, index=False)

        df = load_excel(str(csv_file))
        assert len(df) == 3
        assert "Gene" in df.columns
        # Direction should be auto-filled from log2FoldChange
        assert "Direction" in df.columns
        assert df.loc[df["Gene"] == "TP53", "Direction"].iloc[0] == "Upregulated"

    def test_csv_fallback_on_missing_file(self):
        df = load_excel("/nonexistent/path.csv")
        assert isinstance(df, pd.DataFrame)
        assert len(df) == 50  # demo fallback

    def test_csv_uploaded_file_object(self, tmp_path):
        """Simulate a Streamlit UploadedFile-like object with a .name attribute."""
        import io
        csv_file = tmp_path / "upload.csv"
        demo = generate_demo_data(n=5)
        demo.to_csv(csv_file, index=False)

        # Mimic UploadedFile: BytesIO with a .name attribute
        content = csv_file.read_bytes()
        fake_upload = io.BytesIO(content)
        fake_upload.name = "upload.csv"
        df = load_excel(uploaded_file=fake_upload)
        assert len(df) == 5

    def test_csv_and_excel_produce_same_output(self, tmp_path):
        """CSV and Excel should produce identical cleaned output for same data."""
        demo = generate_demo_data(n=15)

        csv_file = tmp_path / "data.csv"
        demo.to_csv(csv_file, index=False)

        xlsx_file = tmp_path / "data.xlsx"
        demo.to_excel(xlsx_file, index=False, sheet_name="SDEGs")

        df_csv = load_excel(str(csv_file))
        df_xlsx = load_excel(str(xlsx_file), sheet_name="SDEGs")

        assert list(df_csv.columns) == list(df_xlsx.columns)
        assert len(df_csv) == len(df_xlsx)

    def test_csv_sheet_name_ignored(self, tmp_path):
        """sheet_name parameter should be silently ignored for CSV files."""
        csv_file = tmp_path / "data.csv"
        generate_demo_data(n=3).to_csv(csv_file, index=False)

        df = load_excel(str(csv_file), sheet_name="SDEGs")
        assert len(df) == 3


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

    def test_empty_list_gene_types_returns_all(self):
        """Passing gene_types=[] should return all genes, same as None."""
        df = generate_demo_data()
        all_genes = filter_degs(df, padj_threshold=1.0, log2fc_threshold=0, gene_types=None)
        empty_list = filter_degs(df, padj_threshold=1.0, log2fc_threshold=0, gene_types=[])
        assert len(empty_list) == len(all_genes)

    def test_empty_input(self):
        df = pd.DataFrame(columns=["Gene", "padj", "log2FoldChange", "Direction", "Gene Type"])
        filtered = filter_degs(df)
        assert len(filtered) == 0
