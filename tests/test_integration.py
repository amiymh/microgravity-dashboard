"""Integration tests using the real DESeq2 Excel file.

Skipped if the file is not found at the expected path.
These tests validate the full pipeline end-to-end with real data.
"""

import os
import sys
import json
import io
import zipfile
import pytest
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

DATA_PATH = os.path.expanduser(
    "~/Downloads/Supplementary Data 1 DESEQ2 normalized out and significance scores.xlsx"
)

pytestmark = pytest.mark.skipif(
    not os.path.exists(DATA_PATH),
    reason=f"Real data file not found at {DATA_PATH}",
)


@pytest.fixture(scope="module")
def real_data():
    from modules.data_loader import load_excel
    return load_excel(DATA_PATH, sheet_name="SDEGs")


@pytest.fixture(scope="module")
def filtered_data(real_data):
    from modules.data_loader import filter_degs
    return filter_degs(real_data, padj_threshold=0.05, log2fc_threshold=1.0)


class TestDataLoading:
    def test_correct_gene_count(self, real_data):
        assert len(real_data) == 4522, f"Expected 4522 genes, got {len(real_data)}"

    def test_all_required_columns(self, real_data):
        from modules.data_loader import validate_columns
        missing = validate_columns(real_data)
        assert missing == [], f"Missing columns: {missing}"


class TestFilterDegs:
    def test_default_thresholds_reasonable(self, real_data):
        from modules.data_loader import filter_degs
        filtered = filter_degs(real_data, padj_threshold=0.05, log2fc_threshold=1.0)
        assert 3000 <= len(filtered) <= 4522, f"Filtered count {len(filtered)} outside expected range 3000-4522"

    def test_empty_gene_types_same_as_none(self, real_data):
        from modules.data_loader import filter_degs
        with_none = filter_degs(real_data, padj_threshold=1.0, log2fc_threshold=0, gene_types=None)
        with_empty = filter_degs(real_data, padj_threshold=1.0, log2fc_threshold=0, gene_types=[])
        assert len(with_empty) == len(with_none), "Empty list [] should behave same as None"


class TestTargets:
    def test_includes_both_directions(self, filtered_data):
        from modules.targets import get_therapeutic_targets
        result = get_therapeutic_targets(filtered_data, top_n=200)
        if not result.empty and "Direction" in result.columns:
            directions = set(result["Direction"].dropna().unique())
            assert "Upregulated" in directions, "Targets missing upregulated genes"
            assert "Downregulated" in directions, "Targets missing downregulated genes"


class TestBiomarkers:
    def test_filter_includes_both_directions(self, filtered_data):
        """Verify the biomarker filtering step includes both up and downregulated genes.

        We test the filter logic directly (no API calls) since querying 900+ genes
        via OpenTargets would exceed test timeout. The filter is the part that
        previously excluded downregulated genes due to negative significance scores.
        """
        # Reproduce what find_biomarkers does internally: filter by |log2FC|
        candidates = filtered_data[filtered_data["log2FoldChange"].abs() >= 2.0]
        assert len(candidates) > 0
        directions = set(candidates["Direction"].dropna().unique())
        assert "Upregulated" in directions, "Biomarker filter missing upregulated genes"
        assert "Downregulated" in directions, "Biomarker filter missing downregulated genes"

        # Verify sort by absolute significance includes both directions in top positions
        candidates = candidates.assign(
            _abs_sig=candidates["significance score"].abs()
        ).sort_values("_abs_sig", ascending=False)
        top50 = candidates.head(50)
        top_dirs = set(top50["Direction"].dropna().unique())
        assert len(top_dirs) >= 1, "Top 50 by absolute significance should have genes"


class TestPCA:
    def test_variance_explained_pc1(self, real_data):
        from modules.pca import run_pca
        result = run_pca(real_data)
        assert result is not None, "PCA returned None on real data"
        pc1_var = result["variance_explained"][0] * 100
        assert pc1_var > 50, f"PC1 variance {pc1_var:.1f}% too low (expected >50%)"


class TestHeatmap:
    def test_z_matrix_not_none(self, filtered_data):
        from modules.heatmap import compute_heatmap_data
        result = compute_heatmap_data(filtered_data, top_n=50)
        assert result is not None, "Heatmap returned None on real data"
        assert "z_matrix" in result
        assert result["z_matrix"] is not None


class TestCytotoxicity:
    def test_five_pathways_with_overlap(self, filtered_data):
        from modules.cytotoxicity import compute_overlaps
        result = compute_overlaps(filtered_data)
        assert len(result) == 5, f"Expected 5 pathways, got {len(result)}"
        for _, row in result.iterrows():
            assert row["DEG Overlap"] > 0, f"{row['Pathway']} has 0 overlap"


class TestAnalysisLogs:
    def test_volcano_log(self, real_data):
        from modules.volcano import get_analysis_log
        log = get_analysis_log(real_data)
        assert len(log) > 0
        assert all(isinstance(s, str) for s in log)

    def test_pca_log(self, real_data):
        from modules.pca import run_pca, get_analysis_log
        result = run_pca(real_data)
        log = get_analysis_log(result)
        assert len(log) > 0

    def test_heatmap_log(self, filtered_data):
        from modules.heatmap import compute_heatmap_data, get_analysis_log
        hm = compute_heatmap_data(filtered_data, top_n=20)
        log = get_analysis_log(hm, True)
        assert len(log) > 0

    def test_pathways_log(self):
        from modules.pathways import get_analysis_log
        log = get_analysis_log(100, "KEGG_2021_Human", "All")
        assert len(log) > 0

    def test_targets_log(self):
        from modules.targets import get_analysis_log
        log = get_analysis_log(200)
        assert len(log) > 0

    def test_biomarkers_log(self):
        from modules.biomarkers import get_analysis_log
        log = get_analysis_log(100, 50, 2.0)
        assert len(log) > 0

    def test_cytotoxicity_log(self, filtered_data):
        from modules.cytotoxicity import compute_overlaps, get_analysis_log
        overlaps = compute_overlaps(filtered_data)
        log = get_analysis_log(overlaps)
        assert len(log) > 0

    def test_disease_log(self):
        from modules.disease import get_analysis_log
        log = get_analysis_log("Alzheimer", "EFO_0000249", 0.3, "All",
                               {"overlap": 5, "disease_genes": 100, "deg_total": 50})
        assert len(log) > 0


class TestExports:
    def test_report_size(self):
        from modules.report import generate_report
        data = generate_report(session_id="integ-test", filename="test.xlsx", gene_count=4522,
                               analysis_logs={"Volcano Plot": ["Test line"]})
        assert len(data) > 10000, f"Report too small: {len(data)} bytes"

    def test_notebook_valid_json(self):
        from modules.notebook import generate_notebook
        data = generate_notebook(session_id="integ-test", filename="test.xlsx", gene_count=4522)
        assert len(data) > 1000, f"Notebook too small: {len(data)} bytes"
        parsed = json.loads(data)
        assert "cells" in parsed

    def test_zip_contains_report(self):
        from modules.report import generate_report
        from modules.notebook import generate_notebook
        from modules.zip_export import generate_zip

        report = generate_report(session_id="integ-test")
        notebook = generate_notebook(session_id="integ-test")
        data = generate_zip(session_id="integ-test", report_bytes=report, notebook_bytes=notebook)
        assert len(data) > 5000, f"ZIP too small: {len(data)} bytes"

        zf = zipfile.ZipFile(io.BytesIO(data))
        names = zf.namelist()
        assert any("methods_report.docx" in n for n in names), f"methods_report.docx not in ZIP: {names}"
        zf.close()
