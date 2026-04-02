# Build Report — Microgravity RNA Dashboard

**Date:** 2026-04-03
**Builder:** Claude Code (Opus 4.6)

---

## 1. Files Created

```
microgravity-dashboard/
├── app.py                          # Main Streamlit entry point (6 tabs)
├── requirements.txt                # 12 Python dependencies
├── BUILD_REPORT.md                 # This file
├── README.md                       # Setup and usage instructions
├── .gitignore                      # Excludes data files, caches, secrets
├── data/
│   └── .gitkeep
├── modules/
│   ├── __init__.py
│   ├── data_loader.py              # Excel parsing, validation, demo data generation
│   ├── volcano.py                  # Interactive volcano plot (Plotly)
│   ├── pathways.py                 # Enrichr API pathway enrichment
│   ├── targets.py                  # DGIdb drug-gene interactions
│   ├── biomarkers.py               # OpenTargets biomarker discovery
│   ├── cytotoxicity.py             # Curated hallmark gene set overlap analysis
│   └── disease.py                  # OpenTargets disease cross-reference
└── tests/
    ├── test_data_loader.py         # 18 tests
    ├── test_volcano.py             # 12 tests
    ├── test_pathways.py            # 10 tests
    ├── test_targets.py             #  8 tests
    ├── test_biomarkers.py          #  9 tests
    ├── test_cytotoxicity.py        # 13 tests
    └── test_disease.py             # 12 tests
```

**Total: 20 files, 82 tests**

---

## 2. Test Results

```
============================= test session starts ==============================
platform darwin -- Python 3.14.0, pytest-9.0.2, pluggy-1.6.0

tests/test_biomarkers.py::TestQueryOpenTargets::test_returns_list PASSED
tests/test_biomarkers.py::TestQueryOpenTargets::test_nonexistent_gene PASSED
tests/test_biomarkers.py::TestQueryOpenTargets::test_result_structure PASSED
tests/test_biomarkers.py::TestFindBiomarkers::test_returns_dataframe PASSED
tests/test_biomarkers.py::TestFindBiomarkers::test_empty_input PASSED
tests/test_biomarkers.py::TestFindBiomarkers::test_high_threshold_returns_fewer PASSED
tests/test_biomarkers.py::TestFindBiomarkers::test_has_expected_columns PASSED
tests/test_biomarkers.py::TestCreateBiomarkerScatter::test_returns_figure PASSED
tests/test_biomarkers.py::TestCreateBiomarkerScatter::test_empty_input PASSED
tests/test_cytotoxicity.py::TestComputeOverlaps::test_returns_dataframe PASSED
tests/test_cytotoxicity.py::TestComputeOverlaps::test_real_data_overlap PASSED
tests/test_cytotoxicity.py::TestComputeOverlaps::test_has_expected_columns PASSED
tests/test_cytotoxicity.py::TestComputeOverlaps::test_empty_input PASSED
tests/test_cytotoxicity.py::TestComputeOverlaps::test_specific_pathways PASSED
tests/test_cytotoxicity.py::TestGenePathwayMatrix::test_returns_dataframe PASSED
tests/test_cytotoxicity.py::TestGenePathwayMatrix::test_empty_input PASSED
tests/test_cytotoxicity.py::TestCreateOverlapBarChart::test_returns_figure PASSED
tests/test_cytotoxicity.py::TestCreateOverlapBarChart::test_empty_input PASSED
tests/test_cytotoxicity.py::TestCreateHeatmap::test_returns_figure PASSED
tests/test_cytotoxicity.py::TestCreateHeatmap::test_empty_input PASSED
tests/test_cytotoxicity.py::TestGetGeneTable::test_returns_dataframe PASSED
tests/test_cytotoxicity.py::TestGetGeneTable::test_empty_input PASSED
tests/test_data_loader.py::TestGenerateDemoData::test_returns_dataframe PASSED
tests/test_data_loader.py::TestGenerateDemoData::test_default_50_genes PASSED
tests/test_data_loader.py::TestGenerateDemoData::test_custom_size PASSED
tests/test_data_loader.py::TestGenerateDemoData::test_has_required_columns PASSED
tests/test_data_loader.py::TestGenerateDemoData::test_direction_values PASSED
tests/test_data_loader.py::TestGenerateDemoData::test_no_nan_in_key_fields PASSED
tests/test_data_loader.py::TestLoadExcel::test_load_actual_file_sdegs PASSED
tests/test_data_loader.py::TestLoadExcel::test_load_actual_file_all_genes PASSED
tests/test_data_loader.py::TestLoadExcel::test_fallback_to_demo_on_missing_file PASSED
tests/test_data_loader.py::TestLoadExcel::test_no_extra_empty_columns PASSED
tests/test_data_loader.py::TestValidateColumns::test_valid_dataframe PASSED
tests/test_data_loader.py::TestValidateColumns::test_missing_columns PASSED
tests/test_data_loader.py::TestFilterDegs::test_padj_filter PASSED
tests/test_data_loader.py::TestFilterDegs::test_log2fc_filter PASSED
tests/test_data_loader.py::TestFilterDegs::test_direction_filter PASSED
tests/test_data_loader.py::TestFilterDegs::test_gene_type_filter PASSED
tests/test_data_loader.py::TestFilterDegs::test_empty_result PASSED
tests/test_data_loader.py::TestFilterDegs::test_empty_input PASSED
tests/test_disease.py::TestSearchDisease::test_returns_list PASSED
tests/test_disease.py::TestSearchDisease::test_result_structure PASSED
tests/test_disease.py::TestSearchDisease::test_empty_query PASSED
tests/test_disease.py::TestGetDiseaseGenes::test_returns_dataframe PASSED
tests/test_disease.py::TestGetDiseaseGenes::test_invalid_disease PASSED
tests/test_disease.py::TestGetDiseaseGenes::test_has_gene_column PASSED
tests/test_disease.py::TestCrossReference::test_returns_tuple PASSED
tests/test_disease.py::TestCrossReference::test_empty_input PASSED
tests/test_disease.py::TestCrossReference::test_real_data PASSED
tests/test_disease.py::TestCreateVennDiagram::test_returns_figure PASSED
tests/test_disease.py::TestCreateNetworkGraph::test_returns_figure PASSED
tests/test_disease.py::TestCreateNetworkGraph::test_empty_input PASSED
tests/test_pathways.py::TestSubmitGeneList::test_returns_string_or_none PASSED
tests/test_pathways.py::TestSubmitGeneList::test_empty_list PASSED
tests/test_pathways.py::TestRunEnrichment::test_returns_tuple PASSED
tests/test_pathways.py::TestRunEnrichment::test_empty_genes PASSED
tests/test_pathways.py::TestRunEnrichment::test_fallback_on_bad_gene_set PASSED
tests/test_pathways.py::TestRunEnrichment::test_result_has_expected_columns PASSED
tests/test_pathways.py::TestDemoResults::test_demo_has_columns PASSED
tests/test_pathways.py::TestCreatePathwayChart::test_returns_figure PASSED
tests/test_pathways.py::TestCreatePathwayChart::test_empty_input PASSED
tests/test_pathways.py::TestCreatePathwayChart::test_top_n_limit PASSED
tests/test_targets.py::TestQueryDgidb::test_returns_dataframe PASSED
tests/test_targets.py::TestQueryDgidb::test_empty_input PASSED
tests/test_targets.py::TestQueryDgidb::test_has_expected_columns PASSED
tests/test_targets.py::TestGetTherapeuticTargets::test_returns_dataframe PASSED
tests/test_targets.py::TestGetTherapeuticTargets::test_empty_input PASSED
tests/test_targets.py::TestGetTherapeuticTargets::test_has_gene_column PASSED
tests/test_targets.py::TestGetTargetSummary::test_returns_dict PASSED
tests/test_targets.py::TestGetTargetSummary::test_empty_input PASSED
tests/test_volcano.py::TestClassifyGenes::test_adds_volcano_class PASSED
tests/test_volcano.py::TestClassifyGenes::test_classification_values PASSED
tests/test_volcano.py::TestClassifyGenes::test_empty_input PASSED
tests/test_volcano.py::TestClassifyGenes::test_strict_thresholds PASSED
tests/test_volcano.py::TestCreateVolcanoPlot::test_returns_figure PASSED
tests/test_volcano.py::TestCreateVolcanoPlot::test_empty_input PASSED
tests/test_volcano.py::TestCreateVolcanoPlot::test_color_schemes PASSED
tests/test_volcano.py::TestCreateVolcanoPlot::test_no_labels PASSED
tests/test_volcano.py::TestGetTopSignificant::test_returns_top_n PASSED
tests/test_volcano.py::TestGetTopSignificant::test_sorted_by_significance PASSED
tests/test_volcano.py::TestGetTopSignificant::test_empty_input PASSED
tests/test_volcano.py::TestGetTopSignificant::test_required_columns_present PASSED

======================= 82 passed, 8 warnings in 28.30s ========================
```

---

## 3. API Connectivity

| API | Status | Notes |
|-----|--------|-------|
| **Enrichr** (maayanlab.cloud) | Working | Gene list submission + enrichment retrieval both functional |
| **DGIdb** (dgidb.org) | Working | Drug-gene interactions returned for known targets |
| **OpenTargets GraphQL** (api.platform.opentargets.org) | Working | Disease associations, gene targets, disease search all functional |

All 3 external APIs responded successfully during testing. Demo/fallback data is available for all modules in case APIs become unavailable.

---

## 4. Known Limitations

1. **OpenTargets drug details**: The therapeutic targets tab primarily uses DGIdb. OpenTargets drug mechanism queries require Ensembl IDs and are not fully integrated (approval status is heuristic-based).
2. **SDEGs sheet excess columns**: The Excel file's SDEGs sheet has ~980 columns, but only 28 contain data. The loader correctly truncates these.
3. **DisGeNET**: Requires API key registration; the app uses OpenTargets as the primary/fallback source instead.
4. **MSigDB gene sets**: Hardcoded curated gene lists (~60 genes each) rather than fetched from MSigDB API, since MSigDB requires registration. Gene lists are from published hallmark sets.
5. **Venn diagram**: Uses matplotlib circle approximation rather than a proportional area Venn.
6. **Rate limits**: No explicit rate limiting on API calls; batching (DGIdb) and session caching (Streamlit) mitigate this.

---

## 5. Data Notes

- **SDEGs sheet**: 4,522 significantly differentially expressed genes
- **All genes sheet**: 14,725 genes (all types)
- **Top upregulated**: ORM1 (log2FC=4.25, padj=1.76e-121), FGD5 (5.09, 5.29e-51)
- **Gene types present**: protein_coding, lncRNA, processed_pseudogene, and others
- **Column quirk**: SDEGs has columns 'Gene Type' and 'log2', while the all-genes sheet uses 'Type' and lacks 'log2'. The loader handles this automatically.

---

## 6. Deployment Notes

### Run Locally
```bash
pip install -r requirements.txt
streamlit run app.py
```

### Deploy to Streamlit Community Cloud
1. Push code to GitHub (done — `amiymh/microgravity-dashboard`)
2. Go to [share.streamlit.io](https://share.streamlit.io)
3. Connect the GitHub repo `amiymh/microgravity-dashboard`
4. Set main file: `app.py`
5. Deploy

The app runs in demo mode (50 synthetic genes) if no Excel file is uploaded.

---

## 7. Issues Encountered and Fixed

- **No issues encountered.** All 82 tests passed on first full run. All APIs were responsive. No import errors detected.
