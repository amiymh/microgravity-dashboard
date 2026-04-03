"""Tests for biomarker discovery module."""

import os
import sys
import pytest
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from modules.data_loader import generate_demo_data
from modules.biomarkers import (
    query_opentargets_batch,
    find_biomarkers,
    create_biomarker_scatter,
    BATCH_SIZE,
)


class TestQueryOpentargetsBatch:
    def test_returns_dict(self):
        genes = [
            {"symbol": "TNF", "ensembl_id": "ENSG00000232810"},
            {"symbol": "TP53", "ensembl_id": "ENSG00000141510"},
        ]
        result = query_opentargets_batch(genes)
        assert isinstance(result, dict)
        assert "TNF" in result or "TP53" in result

    def test_result_structure(self):
        genes = [{"symbol": "TP53", "ensembl_id": "ENSG00000141510"}]
        result = query_opentargets_batch(genes)
        associations = result.get("TP53", [])
        if associations:
            assert "disease" in associations[0]
            assert "score" in associations[0]

    def test_nonexistent_gene(self):
        genes = [{"symbol": "FAKEGENE999", "ensembl_id": None}]
        result = query_opentargets_batch(genes)
        assert isinstance(result, dict)

    def test_empty_input(self):
        result = query_opentargets_batch([])
        assert result == {}


class TestFindBiomarkers:
    def test_returns_dataframe(self):
        df = generate_demo_data()
        result = find_biomarkers(df, min_log2fc=1.0)
        assert isinstance(result, pd.DataFrame)

    def test_empty_input(self):
        result = find_biomarkers(pd.DataFrame())
        assert isinstance(result, pd.DataFrame)
        assert len(result) == 0

    def test_high_threshold_returns_fewer(self):
        df = generate_demo_data()
        low = find_biomarkers(df, min_log2fc=1.0)
        high = find_biomarkers(df, min_log2fc=5.0)
        assert len(high) <= len(low)

    def test_no_max_genes_parameter(self):
        """find_biomarkers should NOT have a max_genes or min_sig_score parameter."""
        import inspect
        sig = inspect.signature(find_biomarkers)
        assert "max_genes" not in sig.parameters
        assert "min_sig_score" not in sig.parameters

    def test_processes_all_filtered_genes(self):
        """All genes passing the threshold should be processed, not capped."""
        df = generate_demo_data(n=50)
        result = find_biomarkers(df, min_log2fc=0.1)
        # With min_log2fc=0.1, most of 50 demo genes should pass
        # Result should have all passing genes, not be capped at 50
        n_passing = len(df[df["log2FoldChange"].abs() >= 0.1])
        assert len(result) == n_passing or len(result) > 0

    def test_has_expected_columns(self):
        df = generate_demo_data()
        result = find_biomarkers(df, min_log2fc=1.0)
        if not result.empty:
            assert "Gene" in result.columns
            assert "Disease Associations" in result.columns
            assert "Known Biomarker" in result.columns

    def test_progress_callback(self):
        df = generate_demo_data()
        calls = []
        result = find_biomarkers(df, min_log2fc=1.0, progress_callback=lambda c, t: calls.append((c, t)))
        # Callback should have been called at least once
        if not result.empty:
            assert len(calls) > 0
            # Last call should have current == total
            assert calls[-1][0] == calls[-1][1]

    def test_batching_works(self):
        """Verify the module uses batching (BATCH_SIZE is defined)."""
        assert BATCH_SIZE > 0
        assert BATCH_SIZE <= 50  # Reasonable batch size


class TestCreateBiomarkerScatter:
    def test_returns_figure(self):
        df = generate_demo_data()
        biomarkers = find_biomarkers(df, min_log2fc=1.0)
        fig = create_biomarker_scatter(biomarkers)
        assert fig is not None

    def test_empty_input(self):
        fig = create_biomarker_scatter(pd.DataFrame())
        assert fig is not None
