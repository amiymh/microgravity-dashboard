"""Tests for therapeutic targets module — uses real data as primary source."""

import os
import sys
import pytest
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from modules.targets import (
    query_dgidb, query_opentargets_drugs,
    get_therapeutic_targets, get_target_summary, get_analysis_log,
)


class TestQueryDgidb:
    def test_returns_dataframe(self):
        assert isinstance(query_dgidb(["TNF", "IL6", "TP53"]), pd.DataFrame)

    def test_empty_input(self):
        assert isinstance(query_dgidb([]), pd.DataFrame)

    def test_has_expected_columns(self):
        result = query_dgidb(["TNF"])
        if not result.empty:
            assert {"Gene", "Drug Name"}.issubset(set(result.columns))

    def test_no_approval_status_column(self):
        assert "Approval Status" not in query_dgidb(["TNF"]).columns


class TestQueryOpentargetsDrugs:
    def test_returns_dataframe(self):
        assert isinstance(query_opentargets_drugs(["ENSG00000141510"]), pd.DataFrame)

    def test_empty_input(self):
        assert isinstance(query_opentargets_drugs([]), pd.DataFrame)

    def test_clinical_stage_from_api(self):
        result = query_opentargets_drugs(["ENSG00000141510"])
        if not result.empty and "Clinical Stage" in result.columns:
            for s in result["Clinical Stage"].dropna():
                if s:
                    assert s.startswith("PHASE_") or s == ""

    def test_no_approval_status_heuristic(self):
        assert "Approval Status" not in query_opentargets_drugs(["ENSG00000141510"]).columns


class TestGetTherapeuticTargets:
    def test_real_data_both_directions(self, real_df_filtered):
        result = get_therapeutic_targets(real_df_filtered, top_n=200)
        if not result.empty and "Direction" in result.columns:
            dirs = set(result["Direction"].dropna().unique())
            assert "Upregulated" in dirs
            assert "Downregulated" in dirs

    def test_empty_input(self):
        assert len(get_therapeutic_targets(pd.DataFrame())) == 0

    def test_has_gene_column(self, real_df_filtered):
        result = get_therapeutic_targets(real_df_filtered, top_n=10)
        if not result.empty:
            assert "Gene" in result.columns


class TestGetTargetSummary:
    def test_returns_dict(self, real_df_filtered):
        targets = get_therapeutic_targets(real_df_filtered, top_n=10)
        summary = get_target_summary(targets)
        assert "genes_with_drugs" in summary and "total_drugs" in summary

    def test_no_approved_drugs_key(self):
        assert "approved_drugs" not in get_target_summary(pd.DataFrame())

    def test_empty_input(self):
        assert get_target_summary(pd.DataFrame())["genes_with_drugs"] == 0


class TestAnalysisLog:
    def test_returns_list(self):
        log = get_analysis_log(200)
        assert isinstance(log, list) and len(log) > 0
        assert all(isinstance(s, str) for s in log)
