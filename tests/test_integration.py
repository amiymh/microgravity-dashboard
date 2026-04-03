"""Integration tests using the real DESeq2 Excel file.

Uses fixtures from conftest.py which search multiple paths for the data file.
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


class TestDataLoading:
    def test_correct_gene_count(self, real_df):
        assert len(real_df) == 4522, f"Expected 4522 genes, got {len(real_df)}"

    def test_all_required_columns(self, real_df):
        from modules.data_loader import validate_columns
        missing = validate_columns(real_df)
        assert missing == [], f"Missing columns: {missing}"


class TestFilterDegs:
    def test_default_thresholds_reasonable(self, real_df):
        from modules.data_loader import filter_degs
        filtered = filter_degs(real_df, padj_threshold=0.05, log2fc_threshold=1.0)
        assert 3000 <= len(filtered) <= 4522, f"Filtered count {len(filtered)} outside 3000-4522"

    def test_empty_gene_types_same_as_none(self, real_df):
        from modules.data_loader import filter_degs
        with_none = filter_degs(real_df, padj_threshold=1.0, log2fc_threshold=0, gene_types=None)
        with_empty = filter_degs(real_df, padj_threshold=1.0, log2fc_threshold=0, gene_types=[])
        assert len(with_empty) == len(with_none), "Empty list [] must behave same as None"

    def test_never_zero_with_defaults(self, real_df):
        from modules.data_loader import filter_degs
        filtered = filter_degs(real_df, padj_threshold=0.05, log2fc_threshold=1.0)
        assert len(filtered) > 0, "Default filters should not produce 0 genes"
        assert len(filtered) != 16, "Got exactly 16 — likely returning samples instead of genes"


class TestGeneTypeAutoFill:
    def test_significance_sheet_has_gene_type(self):
        """SIGNIFICANCE SCOREs sheet has Type column, renamed to Gene Type."""
        from tests.conftest import DATA_FILE
        if DATA_FILE is None:
            pytest.skip("Data file not available")
        from modules.data_loader import load_excel
        df = load_excel(DATA_FILE, sheet_name="SIGNFICANCE SCOREs Gene Types")
        assert "Gene Type" in df.columns, "Type should be renamed to Gene Type"
        assert df["Gene Type"].notna().sum() > 14000, "All genes should have Gene Type"

    def test_normalized_degs_auto_fills_gene_type(self):
        """Normalized DEGs sheet has no Type/Gene Type — should auto-fill from SDEGs."""
        from tests.conftest import DATA_FILE
        if DATA_FILE is None:
            pytest.skip("Data file not available")
        from modules.data_loader import load_excel
        df = load_excel(DATA_FILE, sheet_name="Normalized DEGs")
        assert "Gene Type" in df.columns, "Gene Type should be auto-filled from SDEGs"
        assert df["Gene Type"].notna().any(), "Some genes should have Gene Type values"

    def test_missing_gene_type_graceful(self, tmp_path):
        """When no Gene Type column exists anywhere, filter should be disabled."""
        import pandas as pd
        from modules.data_loader import load_excel
        csv = tmp_path / "no_genetype.csv"
        pd.DataFrame({
            "Gene": ["A", "B"], "padj": [0.01, 0.02],
            "log2FoldChange": [1.5, -2.0],
        }).to_csv(csv, index=False)
        df = load_excel(str(csv))
        assert "Gene Type" not in df.columns  # no auto-fill possible for CSV

    def test_padj_default_filter_range(self, real_df):
        """Default padj=0.05 on SDEGs returns 3000-4522 genes."""
        from modules.data_loader import filter_degs
        filtered = filter_degs(real_df, padj_threshold=0.05, log2fc_threshold=1.0)
        assert 3000 <= len(filtered) <= 4522


class TestTargets:
    def test_includes_both_directions(self, real_df_filtered):
        from modules.targets import get_therapeutic_targets
        result = get_therapeutic_targets(real_df_filtered, top_n=200)
        if not result.empty and "Direction" in result.columns:
            directions = set(result["Direction"].dropna().unique())
            assert "Upregulated" in directions, "Targets missing upregulated genes"
            assert "Downregulated" in directions, "Targets missing downregulated genes"


class TestBiomarkers:
    def test_filter_includes_both_directions(self, real_df_filtered):
        candidates = real_df_filtered[real_df_filtered["log2FoldChange"].abs() >= 2.0]
        assert len(candidates) > 0
        directions = set(candidates["Direction"].dropna().unique())
        assert "Upregulated" in directions, "Biomarker filter missing upregulated genes"
        assert "Downregulated" in directions, "Biomarker filter missing downregulated genes"

        candidates = candidates.assign(
            _abs_sig=candidates["significance score"].abs()
        ).sort_values("_abs_sig", ascending=False)
        top50 = candidates.head(50)
        assert len(top50) > 0


class TestPCA:
    def test_variance_explained_pc1(self, real_df):
        from modules.pca import run_pca
        result = run_pca(real_df)
        assert result is not None, "PCA returned None on real data"
        pc1_var = result["variance_explained"][0] * 100
        assert pc1_var > 50, f"PC1 variance {pc1_var:.1f}% too low (expected >50%)"


class TestHeatmap:
    def test_z_matrix_not_none(self, real_df_filtered):
        from modules.heatmap import compute_heatmap_data
        result = compute_heatmap_data(real_df_filtered, top_n=50)
        assert result is not None, "Heatmap returned None on real data"
        assert "z_matrix" in result
        assert result["z_matrix"] is not None


class TestCytotoxicity:
    def test_five_pathways_with_overlap(self, real_df_filtered):
        from modules.cytotoxicity import compute_overlaps
        result = compute_overlaps(real_df_filtered)
        assert len(result) == 5, f"Expected 5 pathways, got {len(result)}"
        for _, row in result.iterrows():
            assert row["DEG Overlap"] > 0, f"{row['Pathway']} has 0 overlap"


class TestAnalysisLogs:
    def test_volcano_log(self, real_df):
        from modules.volcano import get_analysis_log
        log = get_analysis_log(real_df)
        assert len(log) > 0 and all(isinstance(s, str) for s in log)

    def test_pca_log(self, real_df):
        from modules.pca import run_pca, get_analysis_log
        result = run_pca(real_df)
        log = get_analysis_log(result)
        assert len(log) > 0

    def test_heatmap_log(self, real_df_filtered):
        from modules.heatmap import compute_heatmap_data, get_analysis_log
        hm = compute_heatmap_data(real_df_filtered, top_n=20)
        log = get_analysis_log(hm, True)
        assert len(log) > 0

    def test_pathways_log(self):
        from modules.pathways import get_analysis_log
        assert len(get_analysis_log(100, "KEGG_2021_Human", "All")) > 0

    def test_targets_log(self):
        from modules.targets import get_analysis_log
        assert len(get_analysis_log(200)) > 0

    def test_biomarkers_log(self):
        from modules.biomarkers import get_analysis_log
        assert len(get_analysis_log(100, 50, 2.0)) > 0

    def test_cytotoxicity_log(self, real_df_filtered):
        from modules.cytotoxicity import compute_overlaps, get_analysis_log
        overlaps = compute_overlaps(real_df_filtered)
        assert len(get_analysis_log(overlaps)) > 0

    def test_disease_log(self):
        from modules.disease import get_analysis_log
        assert len(get_analysis_log("Alzheimer", "EFO_0000249", 0.3, "All",
                                    {"overlap": 5, "disease_genes": 100, "deg_total": 50})) > 0


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
        assert "cells" in json.loads(data)

    def test_zip_contains_report(self):
        from modules.report import generate_report
        from modules.notebook import generate_notebook
        from modules.zip_export import generate_zip
        report = generate_report(session_id="integ-test")
        notebook = generate_notebook(session_id="integ-test")
        data = generate_zip(session_id="integ-test", report_bytes=report, notebook_bytes=notebook)
        assert len(data) > 5000, f"ZIP too small: {len(data)} bytes"
        zf = zipfile.ZipFile(io.BytesIO(data))
        assert any("methods_report.docx" in n for n in zf.namelist())
        zf.close()
