"""Tests for versions module."""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from modules.versions import get_session_versions, get_api_versions, format_versions_text


EXPECTED_SESSION_KEYS = [
    "python", "pandas", "numpy", "scipy", "plotly", "streamlit",
    "requests", "scikit-learn", "openpyxl", "matplotlib", "kaleido",
    "python-docx", "nbformat", "adjustText", "session_id", "session_start",
]

EXPECTED_API_KEYS = ["opentargets", "dgidb", "enrichr", "msigdb"]


class TestSessionVersions:
    def test_session_versions_keys(self):
        sv = get_session_versions()
        for key in EXPECTED_SESSION_KEYS:
            assert key in sv, f"Missing key: {key}"

    def test_session_versions_values(self):
        sv = get_session_versions()
        for key, val in sv.items():
            assert isinstance(val, str), f"{key} is not a string"
            assert len(val) > 0, f"{key} is empty"

    def test_session_id_consistent(self):
        sv1 = get_session_versions()
        sv2 = get_session_versions()
        assert sv1["session_id"] == sv2["session_id"]


class TestApiVersions:
    def test_api_versions(self):
        av = get_api_versions()
        for key in EXPECTED_API_KEYS:
            assert key in av, f"Missing key: {key}"


class TestFormatVersionsText:
    def test_format_versions_text(self):
        text = format_versions_text()
        assert isinstance(text, str)
        assert len(text) > 0
        assert "Python" in text
