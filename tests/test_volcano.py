"""Tests for volcano plot module."""

import os
import sys
import pytest
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from modules.data_loader import generate_demo_data
from modules.volcano import (
    classify_genes, create_volcano_plot, get_top_significant,
    get_analysis_log, _HAS_ADJUSTTEXT,
)


@pytest.fixture
def demo_df():
    return generate_demo_data()


class TestClassifyGenes:
    def test_adds_volcano_class(self, demo_df):
        result = classify_genes(demo_df)
        assert "volcano_class" in result.columns
        assert "neg_log10_padj" in result.columns

    def test_classification_values(self, demo_df):
        result = classify_genes(demo_df)
        valid = {"Upregulated", "Downregulated", "Not Significant"}
        assert set(result["volcano_class"].unique()).issubset(valid)

    def test_empty_input(self):
        df = pd.DataFrame(columns=["Gene", "padj", "log2FoldChange"])
        result = classify_genes(df)
        assert len(result) == 0

    def test_strict_thresholds(self, demo_df):
        result = classify_genes(demo_df, padj_cutoff=1e-100, log2fc_cutoff=10)
        assert (result["volcano_class"] == "Not Significant").all()


class TestCreateVolcanoPlot:
    def test_returns_figure(self, demo_df):
        fig = create_volcano_plot(demo_df)
        assert fig is not None
        assert hasattr(fig, "data")

    def test_empty_input(self):
        df = pd.DataFrame(columns=["Gene", "padj", "log2FoldChange", "significance score"])
        fig = create_volcano_plot(df)
        assert fig is not None

    def test_color_schemes(self, demo_df):
        for scheme in ["Default", "Warm", "Cool", "Unknown"]:
            fig = create_volcano_plot(demo_df, color_scheme=scheme)
            assert fig is not None

    def test_no_labels(self, demo_df):
        fig = create_volcano_plot(demo_df, label_top_n=0)
        assert fig is not None

    def test_labels_with_adjusttext(self, demo_df):
        """Labels render without error whether adjustText is available or not."""
        fig = create_volcano_plot(demo_df, label_top_n=10)
        assert fig is not None

    def test_adjusttext_available(self):
        """adjustText library should be installed."""
        assert _HAS_ADJUSTTEXT is True


class TestGetTopSignificant:
    def test_returns_top_n(self, demo_df):
        result = get_top_significant(demo_df, n=10)
        assert len(result) <= 10

    def test_sorted_by_significance(self, demo_df):
        result = get_top_significant(demo_df, n=10)
        scores = result["significance score"].tolist()
        assert scores == sorted(scores, reverse=True)

    def test_empty_input(self):
        df = pd.DataFrame(columns=["Gene", "significance score"])
        result = get_top_significant(df)
        assert len(result) == 0

    def test_required_columns_present(self, demo_df):
        result = get_top_significant(demo_df)
        assert "Gene" in result.columns
        assert "log2FoldChange" in result.columns


class TestAnalysisLog:
    def test_returns_list(self, demo_df):
        log = get_analysis_log(demo_df)
        assert isinstance(log, list)
        assert len(log) > 0

    def test_all_strings(self, demo_df):
        log = get_analysis_log(demo_df)
        assert all(isinstance(s, str) for s in log)

    def test_contains_keywords(self, demo_df):
        log = get_analysis_log(demo_df)
        text = " ".join(log)
        assert "genes" in text.lower()
        assert "threshold" in text.lower()

    def test_empty_input(self):
        log = get_analysis_log(pd.DataFrame())
        assert isinstance(log, list)
        assert len(log) > 0
