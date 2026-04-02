"""Tests for biomarker discovery module."""

import os
import sys
import pytest
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from modules.data_loader import generate_demo_data
from modules.biomarkers import (
    query_opentargets_associations,
    find_biomarkers,
    create_biomarker_scatter,
)


class TestQueryOpenTargets:
    def test_returns_list(self):
        result = query_opentargets_associations("TNF", "ENSG00000232810")
        assert isinstance(result, list)

    def test_nonexistent_gene(self):
        result = query_opentargets_associations("FAKEGENE999")
        assert isinstance(result, list)

    def test_result_structure(self):
        result = query_opentargets_associations("TP53", "ENSG00000141510")
        if result:
            assert "disease" in result[0]
            assert "score" in result[0]


class TestFindBiomarkers:
    def test_returns_dataframe(self):
        df = generate_demo_data()
        result = find_biomarkers(df, min_log2fc=1.0, max_genes=5)
        assert isinstance(result, pd.DataFrame)

    def test_empty_input(self):
        result = find_biomarkers(pd.DataFrame())
        assert isinstance(result, pd.DataFrame)
        assert len(result) == 0

    def test_high_threshold_returns_fewer(self):
        df = generate_demo_data()
        low = find_biomarkers(df, min_log2fc=1.0, max_genes=5)
        high = find_biomarkers(df, min_log2fc=5.0, max_genes=5)
        assert len(high) <= len(low)

    def test_has_expected_columns(self):
        df = generate_demo_data()
        result = find_biomarkers(df, min_log2fc=1.0, max_genes=5)
        if not result.empty:
            assert "Gene" in result.columns
            assert "Disease Associations" in result.columns
            assert "Known Biomarker" in result.columns


class TestCreateBiomarkerScatter:
    def test_returns_figure(self):
        df = generate_demo_data()
        biomarkers = find_biomarkers(df, min_log2fc=1.0, max_genes=5)
        fig = create_biomarker_scatter(biomarkers)
        assert fig is not None

    def test_empty_input(self):
        fig = create_biomarker_scatter(pd.DataFrame())
        assert fig is not None
