"""Tests for ZIP export module."""

import os
import sys
import io
import zipfile
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from modules.zip_export import generate_zip


class TestGenerateZip:
    def test_generates_without_error(self):
        data = generate_zip(session_id="test-123")
        assert isinstance(data, bytes)
        assert len(data) > 0

    def test_is_valid_zip(self):
        data = generate_zip(session_id="test-123")
        zf = zipfile.ZipFile(io.BytesIO(data))
        names = zf.namelist()
        assert len(names) > 0
        zf.close()

    def test_readme_present(self):
        data = generate_zip(session_id="test-123")
        zf = zipfile.ZipFile(io.BytesIO(data))
        names = zf.namelist()
        readme_files = [n for n in names if "README.txt" in n]
        assert len(readme_files) == 1
        readme_content = zf.read(readme_files[0]).decode("utf-8")
        assert len(readme_content) > 0
        assert "Microgravity" in readme_content
        zf.close()

    def test_report_included(self):
        fake_report = b"fake docx content"
        data = generate_zip(session_id="test", report_bytes=fake_report)
        zf = zipfile.ZipFile(io.BytesIO(data))
        names = zf.namelist()
        assert any("methods_report.docx" in n for n in names)
        zf.close()

    def test_notebook_included(self):
        fake_nb = b'{"cells": []}'
        data = generate_zip(session_id="test", notebook_bytes=fake_nb)
        zf = zipfile.ZipFile(io.BytesIO(data))
        names = zf.namelist()
        assert any("analysis_notebook.ipynb" in n for n in names)
        zf.close()

    def test_figures_included(self):
        data = generate_zip(
            session_id="test",
            figures={"volcano.png": b"png data", "volcano.svg": b"svg data"},
        )
        zf = zipfile.ZipFile(io.BytesIO(data))
        names = zf.namelist()
        assert any("volcano.png" in n for n in names)
        assert any("volcano.svg" in n for n in names)
        zf.close()

    def test_data_csvs_included(self):
        data = generate_zip(
            session_id="test",
            data_csvs={"filtered_genes.csv": b"Gene,padj\nTP53,0.01"},
        )
        zf = zipfile.ZipFile(io.BytesIO(data))
        names = zf.namelist()
        assert any("filtered_genes.csv" in n for n in names)
        zf.close()
