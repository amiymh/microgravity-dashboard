"""Tests for Jupyter notebook export module."""

import os
import sys
import json
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from modules.notebook import generate_notebook


class TestGenerateNotebook:
    def test_generates_without_error(self):
        data = generate_notebook(session_id="test-123", filename="test.xlsx", gene_count=100)
        assert isinstance(data, bytes)
        assert len(data) > 0

    def test_is_valid_json(self):
        data = generate_notebook()
        parsed = json.loads(data)
        assert "cells" in parsed
        assert "nbformat" in parsed

    def test_contains_expected_cell_types(self):
        data = generate_notebook(
            analysis_logs={"Volcano Plot": ["Test line 1"]},
        )
        parsed = json.loads(data)
        cell_types = [c["cell_type"] for c in parsed["cells"]]
        assert "markdown" in cell_types
        assert "code" in cell_types

    def test_non_empty_content(self):
        data = generate_notebook(session_id="abc", filename="data.xlsx", gene_count=50)
        parsed = json.loads(data)
        for cell in parsed["cells"]:
            source = cell.get("source", "")
            assert len(source) > 0 or isinstance(source, list)

    def test_analysis_logs_included(self):
        data = generate_notebook(
            analysis_logs={"PCA Plot": ["Used 16 samples", "Variance PC1: 45%"]},
        )
        text = data.decode("utf-8")
        assert "PCA Plot" in text
        assert "16 samples" in text
