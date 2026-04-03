"""Tests for cytotoxicity & apoptosis module — uses real data as primary source."""

import os
import sys
import pytest
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from modules.cytotoxicity import (
    fetch_geneset, get_hallmark_sets, get_available_pathway_names, is_live_data,
    compute_overlaps, get_gene_pathway_matrix, create_overlap_bar_chart,
    create_heatmap, get_gene_table, get_analysis_log,
    HALLMARK_IDS, _gene_set_cache,
)


class TestFetchGeneset:
    def test_fetches_apoptosis(self):
        genes = fetch_geneset("Apoptosis", "HALLMARK_APOPTOSIS")
        assert len(genes) > 10

    def test_fetches_all_hallmarks(self):
        for name, msigdb_id in HALLMARK_IDS.items():
            assert len(fetch_geneset(name, msigdb_id)) > 0

    def test_live_data_larger_than_fallback(self):
        genes = fetch_geneset("Apoptosis", "HALLMARK_APOPTOSIS")
        if len(genes) > 10:
            assert is_live_data("Apoptosis")

    def test_caching(self):
        _gene_set_cache.clear()
        genes1 = fetch_geneset("Apoptosis", "HALLMARK_APOPTOSIS")
        assert "Apoptosis" in _gene_set_cache
        assert fetch_geneset("Apoptosis", "HALLMARK_APOPTOSIS") is genes1


class TestGetHallmarkSets:
    def test_returns_all_five(self):
        assert len(get_hallmark_sets()) == 5

    def test_specific_pathways(self):
        result = get_hallmark_sets(["Apoptosis", "p53 Pathway"])
        assert len(result) == 2

    def test_gene_sets_are_lists(self):
        for name, genes in get_hallmark_sets().items():
            assert isinstance(genes, list)


class TestComputeOverlaps:
    def test_real_data_overlap(self, real_df_filtered):
        result = compute_overlaps(real_df_filtered)
        assert len(result) == 5
        assert all(result["DEG Overlap"] > 0)
        assert "Source" in result.columns

    def test_live_gene_counts_are_realistic(self, real_df_filtered):
        result = compute_overlaps(real_df_filtered)
        for _, row in result.iterrows():
            assert row["Pathway Size"] >= 8

    def test_empty_input(self):
        assert len(compute_overlaps(pd.DataFrame())) == 0

    def test_specific_pathways(self, demo_df):
        assert len(compute_overlaps(demo_df, pathways=["Apoptosis"])) <= 1


class TestGenePathwayMatrix:
    def test_real_data(self, real_df_filtered):
        matrix = get_gene_pathway_matrix(real_df_filtered)
        assert isinstance(matrix, pd.DataFrame)
        if not matrix.empty:
            assert len(matrix.columns) <= 5

    def test_empty_input(self):
        assert len(get_gene_pathway_matrix(pd.DataFrame())) == 0


class TestCharts:
    def test_overlap_bar_chart(self, real_df_filtered):
        assert create_overlap_bar_chart(compute_overlaps(real_df_filtered)) is not None

    def test_heatmap(self, real_df_filtered):
        matrix = get_gene_pathway_matrix(real_df_filtered)
        assert create_heatmap(matrix) is not None

    def test_empty_inputs(self):
        assert create_overlap_bar_chart(pd.DataFrame()) is not None
        assert create_heatmap(pd.DataFrame()) is not None


class TestGetGeneTable:
    def test_real_data(self, real_df_filtered):
        result = get_gene_table(real_df_filtered)
        assert isinstance(result, pd.DataFrame)
        if not result.empty:
            assert "Gene" in result.columns and "Pathways" in result.columns

    def test_empty_input(self):
        assert len(get_gene_table(pd.DataFrame())) == 0


class TestAnalysisLog:
    def test_real_data(self, real_df_filtered):
        overlaps = compute_overlaps(real_df_filtered)
        log = get_analysis_log(overlaps)
        assert isinstance(log, list) and len(log) > 0
