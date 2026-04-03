"""Tests for therapeutic targets module."""

import os
import sys
import pytest
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from modules.data_loader import generate_demo_data
from modules.targets import (
    query_dgidb,
    query_opentargets_drugs,
    get_therapeutic_targets,
    get_target_summary,
)


class TestQueryDgidb:
    def test_returns_dataframe(self):
        result = query_dgidb(["TNF", "IL6", "TP53"])
        assert isinstance(result, pd.DataFrame)

    def test_empty_input(self):
        result = query_dgidb([])
        assert isinstance(result, pd.DataFrame)

    def test_has_expected_columns(self):
        result = query_dgidb(["TNF"])
        expected = {"Gene", "Drug Name"}
        if not result.empty:
            assert expected.issubset(set(result.columns))

    def test_no_approval_status_column(self):
        """DGIdb should NOT have an Approval Status column — it doesn't provide that."""
        result = query_dgidb(["TNF"])
        assert "Approval Status" not in result.columns


class TestQueryOpentargetsDrugs:
    def test_returns_dataframe(self):
        # TP53 = ENSG00000141510
        result = query_opentargets_drugs(["ENSG00000141510"])
        assert isinstance(result, pd.DataFrame)

    def test_empty_input(self):
        result = query_opentargets_drugs([])
        assert isinstance(result, pd.DataFrame)

    def test_clinical_stage_from_api(self):
        """Clinical stage must come from API, not be guessed."""
        result = query_opentargets_drugs(["ENSG00000141510"])
        if not result.empty and "Clinical Stage" in result.columns:
            # Every non-empty stage should be a real API value
            stages = result["Clinical Stage"].dropna()
            stages = stages[stages != ""]
            for s in stages:
                assert s.startswith("PHASE_") or s == "", f"Unexpected stage format: {s}"

    def test_no_approval_status_heuristic(self):
        """Must not contain heuristic Approval Status column."""
        result = query_opentargets_drugs(["ENSG00000141510"])
        assert "Approval Status" not in result.columns


class TestGetTherapeuticTargets:
    def test_returns_dataframe(self):
        df = generate_demo_data()
        result = get_therapeutic_targets(df, top_n=10)
        assert isinstance(result, pd.DataFrame)

    def test_empty_input(self):
        df = pd.DataFrame()
        result = get_therapeutic_targets(df)
        assert isinstance(result, pd.DataFrame)
        assert len(result) == 0

    def test_has_gene_column(self):
        df = generate_demo_data()
        result = get_therapeutic_targets(df, top_n=10)
        if not result.empty:
            assert "Gene" in result.columns

    def test_no_fabricated_approval(self):
        """The result must not have heuristic 'Approval Status' = 'Experimental' blanket."""
        df = generate_demo_data()
        result = get_therapeutic_targets(df, top_n=10)
        if not result.empty and "Approval Status" in result.columns:
            # Should not have blanket "Experimental" — that was the old heuristic
            assert not (result["Approval Status"] == "Experimental").all()


class TestGetTargetSummary:
    def test_returns_dict(self):
        df = generate_demo_data()
        targets = get_therapeutic_targets(df, top_n=10)
        summary = get_target_summary(targets)
        assert isinstance(summary, dict)
        assert "genes_with_drugs" in summary
        assert "total_drugs" in summary
        assert "with_clinical_stage" in summary

    def test_no_approved_drugs_key(self):
        """The old 'approved_drugs' key (heuristic) should be gone."""
        summary = get_target_summary(pd.DataFrame())
        assert "approved_drugs" not in summary

    def test_empty_input(self):
        summary = get_target_summary(pd.DataFrame())
        assert summary["genes_with_drugs"] == 0
