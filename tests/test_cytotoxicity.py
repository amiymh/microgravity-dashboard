"""Tests for cytotoxicity & apoptosis module."""

import os
import sys
import pytest
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from modules.data_loader import generate_demo_data, load_excel, DEFAULT_PATH
from modules.cytotoxicity import (
    fetch_geneset,
    get_hallmark_sets,
    get_available_pathway_names,
    is_live_data,
    compute_overlaps,
    get_gene_pathway_matrix,
    create_overlap_bar_chart,
    create_heatmap,
    get_gene_table,
    HALLMARK_IDS,
    _gene_set_cache,
)


@pytest.fixture
def demo_df():
    return generate_demo_data()


@pytest.fixture
def real_df():
    if not os.path.exists(DEFAULT_PATH):
        pytest.skip("Data file not available")
    return load_excel(DEFAULT_PATH)


class TestFetchGeneset:
    def test_fetches_apoptosis(self):
        genes = fetch_geneset("Apoptosis", "HALLMARK_APOPTOSIS")
        assert isinstance(genes, list)
        assert len(genes) > 10  # MSigDB has ~161 genes; fallback has 10

    def test_fetches_all_hallmarks(self):
        for name, msigdb_id in HALLMARK_IDS.items():
            genes = fetch_geneset(name, msigdb_id)
            assert len(genes) > 0, f"No genes returned for {name}"

    def test_live_data_larger_than_fallback(self):
        genes = fetch_geneset("Apoptosis", "HALLMARK_APOPTOSIS")
        # Live MSigDB returns ~161 genes, fallback is 10
        if len(genes) > 10:
            assert is_live_data("Apoptosis")

    def test_invalid_geneset_returns_fallback_or_empty(self):
        genes = fetch_geneset("FakePathway", "HALLMARK_FAKE_DOES_NOT_EXIST")
        assert isinstance(genes, list)

    def test_caching(self):
        _gene_set_cache.clear()
        genes1 = fetch_geneset("Apoptosis", "HALLMARK_APOPTOSIS")
        assert "Apoptosis" in _gene_set_cache
        genes2 = fetch_geneset("Apoptosis", "HALLMARK_APOPTOSIS")
        assert genes1 is genes2  # Same object from cache


class TestGetHallmarkSets:
    def test_returns_all_five(self):
        result = get_hallmark_sets()
        assert len(result) == 5

    def test_specific_pathways(self):
        result = get_hallmark_sets(["Apoptosis", "p53 Pathway"])
        assert len(result) == 2
        assert "Apoptosis" in result
        assert "p53 Pathway" in result

    def test_gene_sets_are_lists(self):
        result = get_hallmark_sets()
        for name, genes in result.items():
            assert isinstance(genes, list), f"{name} is not a list"
            assert all(isinstance(g, str) for g in genes)


class TestGetAvailablePathwayNames:
    def test_returns_five_names(self):
        names = get_available_pathway_names()
        assert len(names) == 5
        assert "Apoptosis" in names


class TestComputeOverlaps:
    def test_returns_dataframe(self, demo_df):
        result = compute_overlaps(demo_df)
        assert isinstance(result, pd.DataFrame)

    def test_real_data_overlap(self, real_df):
        result = compute_overlaps(real_df)
        assert len(result) == 5
        assert all(result["DEG Overlap"] > 0)
        assert "Source" in result.columns

    def test_has_expected_columns(self, demo_df):
        result = compute_overlaps(demo_df)
        if not result.empty:
            for col in ["Pathway", "DEG Overlap", "Fisher p-value", "Source"]:
                assert col in result.columns

    def test_empty_input(self):
        result = compute_overlaps(pd.DataFrame())
        assert len(result) == 0

    def test_specific_pathways(self, demo_df):
        result = compute_overlaps(demo_df, pathways=["Apoptosis"])
        assert len(result) <= 1

    def test_live_gene_counts_are_realistic(self, real_df):
        result = compute_overlaps(real_df)
        for _, row in result.iterrows():
            # MSigDB hallmark sets have 40-200 genes typically
            assert row["Pathway Size"] >= 8, f"{row['Pathway']} has only {row['Pathway Size']} genes"


class TestGenePathwayMatrix:
    def test_returns_dataframe(self, real_df):
        matrix = get_gene_pathway_matrix(real_df)
        assert isinstance(matrix, pd.DataFrame)
        if not matrix.empty:
            assert len(matrix.columns) <= 5

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
