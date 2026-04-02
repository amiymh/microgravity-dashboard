"""Tests for disease cross-reference module."""

import os
import sys
import pytest
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from modules.data_loader import generate_demo_data, load_excel, DEFAULT_PATH
from modules.disease import (
    search_disease,
    get_disease_genes,
    cross_reference,
    create_venn_diagram,
    create_network_graph,
)


class TestSearchDisease:
    def test_returns_list(self):
        result = search_disease("Alzheimer")
        assert isinstance(result, list)

    def test_result_structure(self):
        result = search_disease("breast cancer")
        if result:
            assert "id" in result[0]
            assert "name" in result[0]

    def test_empty_query(self):
        result = search_disease("")
        assert isinstance(result, list)


class TestGetDiseaseGenes:
    def test_returns_dataframe(self):
        # EFO_0000249 = Alzheimer disease
        result = get_disease_genes("EFO_0000249", min_score=0.5)
        assert isinstance(result, pd.DataFrame)

    def test_invalid_disease(self):
        result = get_disease_genes("FAKE_DISEASE_123")
        assert isinstance(result, pd.DataFrame)

    def test_has_gene_column(self):
        result = get_disease_genes("EFO_0000249", min_score=0.5)
        if not result.empty:
            assert "Gene" in result.columns
            assert "Disease Association Score" in result.columns


class TestCrossReference:
    def test_returns_tuple(self):
        df = generate_demo_data()
        result, summary = cross_reference(df, "EFO_0000249")
        assert isinstance(result, pd.DataFrame)
        assert isinstance(summary, dict)
        assert "overlap" in summary

    def test_empty_input(self):
        result, summary = cross_reference(pd.DataFrame(), "EFO_0000249")
        assert summary["overlap"] == 0

    def test_real_data(self):
        if not os.path.exists(DEFAULT_PATH):
            pytest.skip("Data file not available")
        df = load_excel(DEFAULT_PATH)
        result, summary = cross_reference(df, "EFO_0000249", min_score=0.5)
        assert isinstance(result, pd.DataFrame)
        assert summary["deg_total"] > 0


class TestCreateVennDiagram:
    def test_returns_figure(self):
        summary = {"overlap": 10, "disease_genes": 100, "deg_total": 4522}
        fig = create_venn_diagram(summary, "Alzheimer")
        assert fig is not None
        import matplotlib.pyplot as plt
        plt.close(fig)


class TestCreateNetworkGraph:
    def test_returns_figure(self):
        df = pd.DataFrame({
            "Gene": ["TNF", "IL6", "TP53"],
            "Direction": ["Upregulated", "Downregulated", "Upregulated"],
            "log2FoldChange": [2.1, -1.5, 3.0],
        })
        fig = create_network_graph(df, "Test Disease")
        assert fig is not None

    def test_empty_input(self):
        fig = create_network_graph(pd.DataFrame())
        assert fig is not None
