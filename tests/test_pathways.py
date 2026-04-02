"""Tests for pathway enrichment module."""

import os
import sys
import pytest
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from modules.pathways import (
    submit_gene_list,
    get_enrichment,
    run_enrichment,
    create_pathway_chart,
    DEMO_RESULTS,
)


TEST_GENES = ["TNF", "IL6", "TP53", "CASP3", "BCL2", "BAX", "NFKB1", "MAPK1"]


class TestSubmitGeneList:
    def test_returns_string_or_none(self):
        result = submit_gene_list(TEST_GENES)
        assert result is None or isinstance(result, str)

    def test_empty_list(self):
        result = submit_gene_list([])
        # May succeed or fail — just shouldn't crash
        assert result is None or isinstance(result, str)


class TestRunEnrichment:
    def test_returns_tuple(self):
        df, is_live = run_enrichment(TEST_GENES, top_n=5)
        assert isinstance(df, pd.DataFrame)
        assert isinstance(is_live, bool)

    def test_empty_genes(self):
        df, is_live = run_enrichment([])
        assert len(df) == 0

    def test_fallback_on_bad_gene_set(self):
        df, is_live = run_enrichment(TEST_GENES, gene_set="NONEXISTENT_DB")
        assert isinstance(df, pd.DataFrame)
        # Should either return results or demo fallback

    def test_result_has_expected_columns(self):
        df, _ = run_enrichment(TEST_GENES, top_n=5)
        if not df.empty:
            assert "term" in df.columns
            assert "pvalue" in df.columns or "adjusted_pvalue" in df.columns


class TestDemoResults:
    def test_demo_has_columns(self):
        assert "term" in DEMO_RESULTS.columns
        assert "combined_score" in DEMO_RESULTS.columns
        assert len(DEMO_RESULTS) == 10


class TestCreatePathwayChart:
    def test_returns_figure(self):
        fig = create_pathway_chart(DEMO_RESULTS)
        assert fig is not None
        assert hasattr(fig, "data")

    def test_empty_input(self):
        df = pd.DataFrame()
        fig = create_pathway_chart(df)
        assert fig is not None

    def test_top_n_limit(self):
        fig = create_pathway_chart(DEMO_RESULTS, top_n=3)
        assert fig is not None
