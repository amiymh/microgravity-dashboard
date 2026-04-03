"""Tests for heatmap module — uses real data as primary source."""

import os
import sys
import pytest
import pandas as pd
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from modules.heatmap import compute_heatmap_data, create_heatmap_figure, get_analysis_log


class TestComputeHeatmap:
    def test_real_data(self, real_df_filtered):
        result = compute_heatmap_data(real_df_filtered, top_n=50)
        assert result is not None
        for key in ["z_matrix", "gene_info", "gene_order", "sample_names", "n_genes"]:
            assert key in result

    def test_zscore_normalization(self, real_df_filtered):
        result = compute_heatmap_data(real_df_filtered, top_n=20)
        z = result["z_matrix"]
        for _, row in z.iterrows():
            if row.notna().sum() > 1:
                assert abs(row.mean()) < 0.5
                assert abs(row.std() - 1.0) < 0.5

    def test_top_n_selection(self, real_df_filtered):
        result = compute_heatmap_data(real_df_filtered, top_n=5)
        assert result["n_genes"] == 5

    def test_clustering(self, real_df_filtered):
        with_cluster = compute_heatmap_data(real_df_filtered, top_n=20, cluster=True)
        no_cluster = compute_heatmap_data(real_df_filtered, top_n=20, cluster=False)
        assert with_cluster is not None and no_cluster is not None
        assert set(with_cluster["gene_order"]) == set(no_cluster["gene_order"])

    def test_empty_input(self):
        assert compute_heatmap_data(pd.DataFrame()) is None

    def test_demo_data(self, demo_df):
        assert compute_heatmap_data(demo_df, top_n=10) is not None


class TestCreateHeatmapFigure:
    def test_returns_figure(self, real_df_filtered):
        fig = create_heatmap_figure(compute_heatmap_data(real_df_filtered, top_n=20))
        assert fig is not None and hasattr(fig, "data")

    def test_none_input(self):
        fig = create_heatmap_figure(None)
        assert fig is not None


class TestAnalysisLog:
    def test_analysis_log(self, real_df_filtered):
        hm = compute_heatmap_data(real_df_filtered, top_n=20)
        log = get_analysis_log(hm, True)
        assert isinstance(log, list) and len(log) > 0
        assert any("z-score" in s.lower() for s in log)

    def test_with_clustering(self, real_df_filtered):
        hm = compute_heatmap_data(real_df_filtered, top_n=20, cluster=True)
        assert any("ward" in s.lower() for s in get_analysis_log(hm, True))

    def test_none_input(self):
        assert len(get_analysis_log(None, False)) > 0
