"""Tests for therapeutic targets module."""

import os
import sys
import pytest
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from modules.data_loader import generate_demo_data
from modules.targets import query_dgidb, get_therapeutic_targets, get_target_summary


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


class TestGetTargetSummary:
    def test_returns_dict(self):
        df = generate_demo_data()
        targets = get_therapeutic_targets(df, top_n=10)
        summary = get_target_summary(targets)
        assert isinstance(summary, dict)
        assert "genes_with_drugs" in summary
        assert "total_drugs" in summary

    def test_empty_input(self):
        summary = get_target_summary(pd.DataFrame())
        assert summary["genes_with_drugs"] == 0
