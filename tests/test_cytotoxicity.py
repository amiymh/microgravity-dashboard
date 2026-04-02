"""Tests for cytotoxicity & apoptosis module."""

import os
import sys
import pytest
import pandas as pd
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from modules.data_loader import generate_demo_data, load_excel, DEFAULT_PATH
from modules.cytotoxicity import (
    compute_overlaps,
    get_gene_pathway_matrix,
    create_overlap_bar_chart,
    create_heatmap,
    get_gene_table,
    HALLMARK_SETS,
)


@pytest.fixture
def demo_df():
    return generate_demo_data()


@pytest.fixture
def real_df():
    if not os.path.exists(DEFAULT_PATH):
        pytest.skip("Data file not available")
    return load_excel(DEFAULT_PATH)


class TestComputeOverlaps:
    def test_returns_dataframe(self, demo_df):
        result = compute_overlaps(demo_df)
        assert isinstance(result, pd.DataFrame)

    def test_real_data_overlap(self, real_df):
        result = compute_overlaps(real_df)
        assert len(result) == 5  # All 5 pathways
        assert all(result["DEG Overlap"] > 0)

    def test_has_expected_columns(self, demo_df):
        result = compute_overlaps(demo_df)
        if not result.empty:
            for col in ["Pathway", "DEG Overlap", "Fisher p-value"]:
                assert col in result.columns

    def test_empty_input(self):
        result = compute_overlaps(pd.DataFrame())
        assert len(result) == 0

    def test_specific_pathways(self, demo_df):
        result = compute_overlaps(demo_df, pathways=["Apoptosis"])
        assert len(result) <= 1


class TestGenePathwayMatrix:
    def test_returns_dataframe(self, real_df):
        matrix = get_gene_pathway_matrix(real_df)
        assert isinstance(matrix, pd.DataFrame)
        if not matrix.empty:
            assert len(matrix.columns) <= len(HALLMARK_SETS)

    def test_empty_input(self):
        result = get_gene_pathway_matrix(pd.DataFrame())
        assert len(result) == 0


class TestCreateOverlapBarChart:
    def test_returns_figure(self, demo_df):
        overlaps = compute_overlaps(demo_df)
        fig = create_overlap_bar_chart(overlaps)
        assert fig is not None

    def test_empty_input(self):
        fig = create_overlap_bar_chart(pd.DataFrame())
        assert fig is not None


class TestCreateHeatmap:
    def test_returns_figure(self, real_df):
        matrix = get_gene_pathway_matrix(real_df)
        fig = create_heatmap(matrix)
        assert fig is not None

    def test_empty_input(self):
        fig = create_heatmap(pd.DataFrame())
        assert fig is not None


class TestGetGeneTable:
    def test_returns_dataframe(self, real_df):
        result = get_gene_table(real_df)
        assert isinstance(result, pd.DataFrame)
        if not result.empty:
            assert "Gene" in result.columns
            assert "Pathways" in result.columns

    def test_empty_input(self):
        result = get_gene_table(pd.DataFrame())
        assert len(result) == 0
