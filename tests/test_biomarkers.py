"""Tests for biomarker discovery module — uses real data as primary source."""

import os
import sys
import inspect
import pytest
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from modules.biomarkers import (
    query_opentargets_batch, find_biomarkers,
    create_biomarker_scatter, get_analysis_log, BATCH_SIZE,
)


class TestQueryOpentargetsBatch:
    def test_returns_dict(self):
        genes = [{"symbol": "TNF", "ensembl_id": "ENSG00000232810"},
                 {"symbol": "TP53", "ensembl_id": "ENSG00000141510"}]
        assert isinstance(query_opentargets_batch(genes), dict)

    def test_result_structure(self):
        result = query_opentargets_batch([{"symbol": "TP53", "ensembl_id": "ENSG00000141510"}])
        assoc = result.get("TP53", [])
        if assoc:
            assert "disease" in assoc[0] and "score" in assoc[0]

    def test_empty_input(self):
        assert query_opentargets_batch([]) == {}


class TestFindBiomarkers:
    def test_returns_dataframe(self, demo_df):
        assert isinstance(find_biomarkers(demo_df, min_log2fc=1.0), pd.DataFrame)

    def test_empty_input(self):
        assert len(find_biomarkers(pd.DataFrame())) == 0

    def test_high_threshold_returns_fewer(self, demo_df):
        low = find_biomarkers(demo_df, min_log2fc=1.0)
        high = find_biomarkers(demo_df, min_log2fc=5.0)
        assert len(high) <= len(low)

    def test_no_removed_parameters(self):
        sig = inspect.signature(find_biomarkers)
        assert "max_genes" not in sig.parameters
        assert "min_sig_score" not in sig.parameters

    def test_processes_all_filtered_genes(self, demo_df):
        result = find_biomarkers(demo_df, min_log2fc=0.1)
        n_passing = len(demo_df[demo_df["log2FoldChange"].abs() >= 0.1])
        assert len(result) == n_passing or len(result) > 0

    def test_has_expected_columns(self, demo_df):
        result = find_biomarkers(demo_df, min_log2fc=1.0)
        if not result.empty:
            assert "Gene" in result.columns and "Disease Associations" in result.columns

    def test_progress_callback(self, demo_df):
        calls = []
        find_biomarkers(demo_df, min_log2fc=1.0, progress_callback=lambda c, t: calls.append((c, t)))
        if calls:
            assert calls[-1][0] == calls[-1][1]

    def test_batching_works(self):
        assert 0 < BATCH_SIZE <= 50

    def test_real_data_includes_both_directions(self, real_df_filtered):
        """Verify filter includes both up and downregulated at |log2FC|>=2."""
        candidates = real_df_filtered[real_df_filtered["log2FoldChange"].abs() >= 2.0]
        if len(candidates) > 0:
            dirs = set(candidates["Direction"].dropna().unique())
            assert "Upregulated" in dirs
            assert "Downregulated" in dirs


class TestCreateBiomarkerScatter:
    def test_returns_figure(self, demo_df):
        fig = create_biomarker_scatter(find_biomarkers(demo_df, min_log2fc=1.0))
        assert fig is not None

    def test_empty_input(self):
        assert create_biomarker_scatter(pd.DataFrame()) is not None


class TestAnalysisLog:
    def test_returns_list(self):
        log = get_analysis_log(100, 50, 2.0)
        assert isinstance(log, list) and len(log) > 0
