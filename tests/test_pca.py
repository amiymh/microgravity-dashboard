"""Tests for PCA module — uses real data as primary source."""

import os
import sys
import pytest
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from modules.pca import run_pca, get_loadings, create_pca_plot, get_analysis_log


class TestRunPCA:
    def test_run_pca_real_data(self, real_df):
        result = run_pca(real_df)
        assert result is not None
        assert "pca_df" in result and "variance_explained" in result

    def test_pca_correct_components(self, real_df):
        result = run_pca(real_df)
        assert len(result["variance_explained"]) == 3
        assert all(v > 0 for v in result["variance_explained"])

    def test_pca_sample_count(self, real_df):
        assert len(run_pca(real_df)["pca_df"]) == 16

    def test_pca_empty_input(self):
        assert run_pca(pd.DataFrame()) is None

    def test_pca_demo_data(self, demo_df):
        assert run_pca(demo_df) is not None


class TestCreatePCAPlot:
    def test_create_pca_plot(self, real_df):
        fig = create_pca_plot(run_pca(real_df))
        assert fig is not None and hasattr(fig, "data")


class TestGetLoadings:
    def test_get_loadings(self, real_df):
        result = run_pca(real_df)
        loadings = get_loadings(real_df, result["pca_model"])
        assert isinstance(loadings, pd.DataFrame)
        assert "Gene" in loadings.columns and len(loadings) <= 20


class TestAnalysisLog:
    def test_analysis_log(self, real_df):
        log = get_analysis_log(run_pca(real_df))
        assert isinstance(log, list) and len(log) > 0

    def test_analysis_log_keywords(self, real_df):
        text = " ".join(get_analysis_log(run_pca(real_df))).lower()
        assert "variance" in text and "genes" in text and "samples" in text
