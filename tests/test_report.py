"""Tests for Word methods report module."""

import os
import sys
import pytest
import io

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from modules.report import generate_report


class TestGenerateReport:
    def test_generates_without_error(self):
        data = generate_report(session_id="test-123", filename="test.xlsx", gene_count=100)
        assert isinstance(data, bytes)
        assert len(data) > 0

    def test_is_valid_docx(self):
        data = generate_report()
        from docx import Document
        doc = Document(io.BytesIO(data))
        # Should have paragraphs
        assert len(doc.paragraphs) > 0

    def test_contains_expected_headings(self):
        data = generate_report(
            analysis_logs={"Volcano Plot": ["Test line 1", "Test line 2"]},
        )
        from docx import Document
        doc = Document(io.BytesIO(data))
        headings = [p.text for p in doc.paragraphs if p.style.name.startswith("Heading")]
        heading_text = " ".join(headings)
        assert "Volcano Plot" in heading_text
        assert "Appendix A" in heading_text
        assert "Appendix C" in heading_text

    def test_gene_lists_in_appendix(self):
        data = generate_report(
            gene_lists={"Test Analysis": ["TP53", "TNF", "IL6"]},
        )
        from docx import Document
        doc = Document(io.BytesIO(data))
        full_text = " ".join(p.text for p in doc.paragraphs)
        assert "TP53" in full_text

    def test_versions_in_appendix(self):
        data = generate_report(
            versions={"Python": "3.14.0", "pandas": "2.2.0"},
        )
        from docx import Document
        doc = Document(io.BytesIO(data))
        # Check table content
        found = False
        for table in doc.tables:
            for row in table.rows:
                cells_text = [c.text for c in row.cells]
                if "Python" in cells_text:
                    found = True
        assert found

    def test_not_run_tabs_noted(self):
        data = generate_report()  # No analysis_logs
        from docx import Document
        doc = Document(io.BytesIO(data))
        full_text = " ".join(p.text for p in doc.paragraphs)
        assert "Not run" in full_text
