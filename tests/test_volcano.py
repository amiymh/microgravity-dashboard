"""Tests for volcano plot module — uses real data as primary source."""

import os
import sys
import pytest
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from modules.volcano import (
    classify_genes, create_volcano_plot, get_top_significant,
    get_analysis_log, _HAS_ADJUSTTEXT,
)


class TestClassifyGenes:
    def test_adds_volcano_class(self, real_df):
        result = classify_genes(real_df)
        assert "volcano_class" in result.columns
        assert "neg_log10_padj" in result.columns

    def test_classification_values(self, real_df):
        result = classify_genes(real_df)
        assert set(result["volcano_class"].unique()).issubset(
            {"Upregulated", "Downregulated", "Not Significant"})

    def test_empty_input(self):
        result = classify_genes(pd.DataFrame(columns=["Gene", "padj", "log2FoldChange"]))
        assert len(result) == 0

    def test_strict_thresholds(self, real_df):
        result = classify_genes(real_df, padj_cutoff=1e-300, log2fc_cutoff=100)
        assert (result["volcano_class"] == "Not Significant").all()


class TestCreateVolcanoPlot:
    def test_returns_figure(self, real_df):
        fig = create_volcano_plot(real_df)
        assert fig is not None and hasattr(fig, "data")

    def test_empty_input(self):
        fig = create_volcano_plot(pd.DataFrame(columns=["Gene", "padj", "log2FoldChange", "significance score"]))
        assert fig is not None

    def test_color_schemes(self, demo_df):
        for scheme in ["Default", "Warm", "Cool", "Unknown"]:
            assert create_volcano_plot(demo_df, color_scheme=scheme) is not None

    def test_no_labels(self, demo_df):
        assert create_volcano_plot(demo_df, label_top_n=0) is not None

    def test_labels_with_adjusttext(self, real_df):
        assert create_volcano_plot(real_df, label_top_n=10) is not None

    def test_adjusttext_available(self):
        assert _HAS_ADJUSTTEXT is True


class TestGetTopSignificant:
    def test_returns_top_n(self, real_df):
        result = get_top_significant(real_df, n=10)
        assert len(result) == 10

    def test_sorted_by_significance(self, real_df):
        result = get_top_significant(real_df, n=10)
        scores = result["significance score"].tolist()
        assert scores == sorted(scores, reverse=True)

    def test_empty_input(self):
        assert len(get_top_significant(pd.DataFrame(columns=["Gene", "significance score"]))) == 0


class TestAnalysisLog:
    def test_returns_list(self, real_df):
        log = get_analysis_log(real_df)
        assert isinstance(log, list) and len(log) > 0

    def test_all_strings(self, real_df):
        assert all(isinstance(s, str) for s in get_analysis_log(real_df))

    def test_contains_keywords(self, real_df):
        text = " ".join(get_analysis_log(real_df))
        assert "genes" in text.lower() and "threshold" in text.lower()

    def test_empty_input(self):
        assert len(get_analysis_log(pd.DataFrame())) > 0
