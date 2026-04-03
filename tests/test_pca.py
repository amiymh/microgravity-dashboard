"""Tests for PCA analysis module."""

import os
import sys
import pytest
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from modules.data_loader import generate_demo_data
from modules.pca import run_pca, create_pca_plot, get_loadings, get_analysis_log


@pytest.fixture
def demo_df():
    return generate_demo_data()


@pytest.fixture
def pca_result(demo_df):
    return run_pca(demo_df)


class TestRunPCA:
    def test_run_pca_demo_data(self, pca_result):
        """run_pca on demo data returns dict with expected keys."""
        assert pca_result is not None
        assert isinstance(pca_result, dict)
        expected_keys = {"pca_df", "variance_explained", "n_genes_used", "n_genes_removed"}
        assert expected_keys.issubset(pca_result.keys())

    def test_pca_correct_components(self, pca_result):
        """variance_explained has 3 values, all > 0."""
        var_exp = pca_result["variance_explained"]
        assert len(var_exp) == 3
        assert all(v > 0 for v in var_exp)

    def test_pca_sample_count(self, pca_result):
        """pca_df has 16 rows (or fewer if demo has fewer samples)."""
        pca_df = pca_result["pca_df"]
        assert len(pca_df) <= 16
        assert len(pca_df) > 0

    def test_pca_empty_input(self):
        """run_pca returns None on empty DataFrame."""
        empty_df = pd.DataFrame()
        result = run_pca(empty_df)
        assert result is None


class TestCreatePCAPlot:
    def test_create_pca_plot(self, pca_result):
        """create_pca_plot returns a Plotly figure."""
        fig = create_pca_plot(pca_result)
        assert fig is not None
        assert hasattr(fig, "data")
        assert hasattr(fig, "layout")


class TestGetLoadings:
    def test_get_loadings(self, demo_df, pca_result):
        """get_loadings returns DataFrame with Gene column, length <= 20."""
        pca_model = pca_result["pca_model"]
        loading_df = get_loadings(demo_df, pca_model, n=20)
        assert isinstance(loading_df, pd.DataFrame)
        assert "Gene" in loading_df.columns
        assert len(loading_df) <= 20


class TestAnalysisLog:
    def test_analysis_log(self, pca_result):
        """get_analysis_log returns non-empty list of strings."""
        log = get_analysis_log(pca_result)
        assert isinstance(log, list)
        assert len(log) > 0
        assert all(isinstance(entry, str) for entry in log)

    def test_analysis_log_keywords(self, pca_result):
        """Log contains 'variance', 'genes', and 'samples'."""
        log = get_analysis_log(pca_result)
        combined = " ".join(log).lower()
        assert "variance" in combined
        assert "genes" in combined
        assert "samples" in combined
