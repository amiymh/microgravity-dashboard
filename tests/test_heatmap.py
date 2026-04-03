"""Tests for heatmap module."""

import os
import sys
import pytest
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from modules.data_loader import generate_demo_data
from modules.heatmap import compute_heatmap_data, create_heatmap_figure, get_analysis_log


@pytest.fixture
def demo_df():
    return generate_demo_data()


class TestComputeHeatmapDemo:
    def test_compute_heatmap_demo(self, demo_df):
        result = compute_heatmap_data(demo_df)
        assert result is not None
        expected_keys = {"z_matrix", "gene_info", "gene_order", "sample_names", "n_genes"}
        assert set(result.keys()) == expected_keys

    def test_zscore_normalization(self, demo_df):
        result = compute_heatmap_data(demo_df, cluster=False)
        z = result["z_matrix"]
        row_means = z.mean(axis=1)
        row_stds = z.std(axis=1)
        # Each row should have approximately mean=0 and std=1
        assert np.allclose(row_means, 0, atol=1e-10)
        assert np.allclose(row_stds, 1, atol=0.15)

    def test_top_n_selection(self, demo_df):
        result = compute_heatmap_data(demo_df, top_n=5)
        assert result["gene_info"].shape[0] == 5
        assert result["n_genes"] == 5

    def test_top_n_exceeds_data(self):
        small_df = generate_demo_data(n=3)
        result = compute_heatmap_data(small_df, top_n=10)
        assert result is not None
        assert result["n_genes"] <= 3

    def test_clustering_reorders(self, demo_df):
        result_clustered = compute_heatmap_data(demo_df, top_n=20, cluster=True)
        result_unclustered = compute_heatmap_data(demo_df, top_n=20, cluster=False)
        # Clustering should produce a result without crashing
        assert result_clustered is not None
        assert result_unclustered is not None
        # Orders may differ (not guaranteed but likely with enough genes)
        # At minimum, both should have the same set of genes
        assert set(result_clustered["gene_order"]) == set(result_unclustered["gene_order"])

    def test_no_clustering(self, demo_df):
        result = compute_heatmap_data(demo_df, cluster=False)
        assert result is not None
        assert "z_matrix" in result

    def test_empty_input(self):
        empty_df = pd.DataFrame()
        result = compute_heatmap_data(empty_df)
        assert result is None

    def test_empty_with_columns(self):
        empty_df = pd.DataFrame(columns=["Gene", "significance score", "Earth_3D2", "Space_3D2"])
        result = compute_heatmap_data(empty_df)
        assert result is None


class TestCreateHeatmapFigure:
    def test_create_heatmap_figure(self, demo_df):
        heatmap_data = compute_heatmap_data(demo_df)
        fig = create_heatmap_figure(heatmap_data)
        assert fig is not None
        assert hasattr(fig, "data")
        assert len(fig.data) > 0

    def test_none_input(self):
        fig = create_heatmap_figure(None)
        assert fig is not None
        assert fig.layout.title.text == "No data to display"


class TestAnalysisLog:
    def test_analysis_log(self, demo_df):
        heatmap_data = compute_heatmap_data(demo_df)
        log = get_analysis_log(heatmap_data, cluster=True)
        assert isinstance(log, list)
        assert len(log) > 0
        log_text = " ".join(log)
        assert "z-score" in log_text.lower()

    def test_analysis_log_no_cluster(self, demo_df):
        heatmap_data = compute_heatmap_data(demo_df)
        log = get_analysis_log(heatmap_data, cluster=False)
        assert not any("ward" in line.lower() for line in log)

    def test_analysis_log_with_cluster(self, demo_df):
        heatmap_data = compute_heatmap_data(demo_df)
        log = get_analysis_log(heatmap_data, cluster=True)
        log_text = " ".join(log)
        assert "ward" in log_text.lower()

    def test_analysis_log_none_input(self):
        log = get_analysis_log(None)
        assert isinstance(log, list)
        assert len(log) > 0
