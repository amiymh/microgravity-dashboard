"""Tests for data_loader module — real data is primary test source."""

import os
import sys
import io
import pytest
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from modules.data_loader import (
    generate_demo_data, load_excel, validate_columns, filter_degs,
)


class TestGenerateDemoData:
    def test_returns_dataframe(self):
        assert isinstance(generate_demo_data(), pd.DataFrame)

    def test_default_50_genes(self):
        assert len(generate_demo_data()) == 50

    def test_custom_size(self):
        assert len(generate_demo_data(n=20)) == 20

    def test_has_required_columns(self):
        assert validate_columns(generate_demo_data()) == []

    def test_direction_values(self):
        df = generate_demo_data()
        assert set(df["Direction"].unique()).issubset({"Upregulated", "Downregulated"})

    def test_no_nan_in_key_fields(self):
        df = generate_demo_data()
        for col in ["Gene", "padj", "log2FoldChange", "Direction"]:
            assert df[col].isna().sum() == 0


class TestLoadExcel:
    def test_load_real_file_sdegs(self, real_df):
        assert len(real_df) > 4000
        assert "Gene" in real_df.columns
        assert "Direction" in real_df.columns

    def test_load_real_file_column_count(self, real_df):
        assert len(real_df.columns) < 30

    def test_fallback_to_demo_on_missing_file(self):
        df = load_excel("/nonexistent/path.xlsx")
        assert len(df) == 50


class TestLoadCSV:
    def test_load_csv_with_correct_columns(self, tmp_path):
        csv_file = tmp_path / "test_data.csv"
        generate_demo_data(n=10).to_csv(csv_file, index=False)
        df = load_excel(str(csv_file))
        assert len(df) == 10 and "Gene" in df.columns

    def test_csv_missing_columns_still_loads(self, tmp_path):
        csv_file = tmp_path / "partial.csv"
        pd.DataFrame({
            "Gene": ["TP53", "TNF", "IL6"],
            "log2FoldChange": [2.5, -1.3, 0.8],
            "padj": [0.001, 0.01, 0.05],
        }).to_csv(csv_file, index=False)
        df = load_excel(str(csv_file))
        assert len(df) == 3 and "Direction" in df.columns
        assert df.loc[df["Gene"] == "TP53", "Direction"].iloc[0] == "Upregulated"

    def test_csv_fallback_on_missing_file(self):
        assert len(load_excel("/nonexistent/path.csv")) == 50

    def test_csv_uploaded_file_object(self, tmp_path):
        csv_file = tmp_path / "upload.csv"
        generate_demo_data(n=5).to_csv(csv_file, index=False)
        fake_upload = io.BytesIO(csv_file.read_bytes())
        fake_upload.name = "upload.csv"
        assert len(load_excel(uploaded_file=fake_upload)) == 5

    def test_csv_and_excel_produce_same_output(self, tmp_path):
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
        csv_file = tmp_path / "data.csv"
        generate_demo_data(n=3).to_csv(csv_file, index=False)
        assert len(load_excel(str(csv_file), sheet_name="SDEGs")) == 3


class TestValidateColumns:
    def test_real_data_valid(self, real_df):
        assert validate_columns(real_df) == []

    def test_missing_columns(self):
        df = pd.DataFrame({"Gene": ["A"], "padj": [0.01]})
        assert "log2FoldChange" in validate_columns(df)


class TestFilterDegs:
    def test_real_data_default_filters(self, real_df):
        filtered = filter_degs(real_df, padj_threshold=0.05, log2fc_threshold=1.0)
        assert 3000 <= len(filtered) <= 4522

    def test_padj_filter(self, real_df):
        filtered = filter_degs(real_df, padj_threshold=0.01, log2fc_threshold=0)
        assert all(filtered["padj"] <= 0.01)

    def test_log2fc_filter(self, real_df):
        filtered = filter_degs(real_df, log2fc_threshold=2.0, padj_threshold=1.0)
        assert all(filtered["log2FoldChange"].abs() >= 2.0)

    def test_direction_filter(self, real_df):
        up = filter_degs(real_df, direction="Upregulated", padj_threshold=1.0, log2fc_threshold=0)
        assert all(up["Direction"] == "Upregulated")

    def test_gene_type_filter(self, real_df):
        filtered = filter_degs(real_df, gene_types=["protein_coding"], padj_threshold=1.0, log2fc_threshold=0)
        assert all(filtered["Gene Type"] == "protein_coding")

    def test_empty_list_gene_types_returns_all(self, real_df):
        with_none = filter_degs(real_df, padj_threshold=1.0, log2fc_threshold=0, gene_types=None)
        with_empty = filter_degs(real_df, padj_threshold=1.0, log2fc_threshold=0, gene_types=[])
        assert len(with_empty) == len(with_none)

    def test_empty_result(self, demo_df):
        assert len(filter_degs(demo_df, padj_threshold=1e-300)) == 0

    def test_empty_input(self):
        df = pd.DataFrame(columns=["Gene", "padj", "log2FoldChange", "Direction", "Gene Type"])
        assert len(filter_degs(df)) == 0
