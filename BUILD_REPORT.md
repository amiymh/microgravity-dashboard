# Build Report — Microgravity RNA Dashboard V2

**Date:** 2026-04-03
**Builder:** Claude Code (Opus 4.6)

---

## 1. What Changed in V2

### New Modules (6)
- `modules/pca.py` — PCA sample clustering (sklearn), loadings table
- `modules/heatmap.py` — Top DEGs heatmap with z-score normalization, hierarchical clustering
- `modules/versions.py` — Library and API version tracking
- `modules/report.py` — Word methods report generation (python-docx)
- `modules/notebook.py` — Jupyter notebook export (nbformat)
- `modules/zip_export.py` — Full results ZIP package

### New Utility
- `modules/fig_export.py` — PNG (300 DPI) and SVG vector export via kaleido

### Enhanced Modules (6)
- `modules/volcano.py` — adjustText non-overlapping labels + analysis log
- `modules/pathways.py` — analysis log
- `modules/targets.py` — analysis log
- `modules/biomarkers.py` — analysis log
- `modules/cytotoxicity.py` — analysis log
- `modules/disease.py` — analysis log

### Enhanced Data Loader
- CSV support with auto-detection + sheet selector hides for CSV
- Extension-less file fallback: try Excel then CSV

### App (app.py) — Complete Rewrite
- 9 tabs (was 6): added PCA, Heatmap, User Manual
- Every analysis tab has transparency panel (Analysis Details expander)
- Every figure has PNG/SVG download buttons
- Sidebar: 3 prominent download buttons (Report, Notebook, ZIP)
- Version info expander in sidebar
- User Manual tab with full documentation

### New Dependencies
- scikit-learn>=1.3.0
- adjustText>=0.8
- kaleido>=0.2.1
- python-docx>=1.1.0
- nbformat>=5.9.0

---

## 2. Files Created/Modified

```
New files:
  modules/pca.py
  modules/heatmap.py
  modules/versions.py
  modules/report.py
  modules/notebook.py
  modules/zip_export.py
  modules/fig_export.py
  tests/test_pca.py
  tests/test_heatmap.py
  tests/test_versions.py
  tests/test_report.py
  tests/test_notebook.py
  tests/test_zip_export.py

Modified files:
  app.py (complete rewrite)
  modules/data_loader.py (CSV support, extension detection)
  modules/volcano.py (adjustText, analysis log)
  modules/pathways.py (analysis log)
  modules/targets.py (analysis log)
  modules/biomarkers.py (analysis log)
  modules/cytotoxicity.py (analysis log)
  modules/disease.py (analysis log)
  requirements.txt (5 new deps)
  .gitignore (added *.png, *.svg)
  tests/test_data_loader.py (CSV tests)
  tests/test_volcano.py (adjustText + log tests)
```

---

## 3. Test Results

```
161 passed, 9 warnings in 113.29s

Test breakdown:
  test_biomarkers.py     — 14 passed
  test_cytotoxicity.py   — 23 passed
  test_data_loader.py    — 24 passed
  test_disease.py        — 12 passed
  test_heatmap.py        — 14 passed
  test_notebook.py       —  5 passed
  test_pathways.py       — 10 passed
  test_pca.py            —  8 passed
  test_report.py         —  6 passed
  test_targets.py        — 15 passed
  test_versions.py       —  5 passed
  test_volcano.py        — 18 passed
  test_zip_export.py     —  7 passed
```

---

## 4. API Connectivity

| API | Status | Notes |
|-----|--------|-------|
| Enrichr | Working | Gene list submission + enrichment retrieval |
| DGIdb | Working | Drug-gene interactions |
| OpenTargets GraphQL | Working | Disease associations, drug candidates |
| MSigDB | Working | Full hallmark gene sets (e.g. 161 genes for Apoptosis) |

---

## 5. Known Limitations

1. Enrichr is exploratory — GSEA recommended for publication
2. Drug approval status not shown — only clinical stage from OpenTargets
3. Human genes only
4. Biomarker tab may be slow for large gene sets (all genes processed, no cap)
5. MSigDB gene sets fetched live — may change between versions
6. CSV files: single sheet only
7. adjustText label positions are approximate (computed via matplotlib, displayed in Plotly)

---

## 6. How to Run Locally

```bash
cd microgravity-dashboard
pip install -r requirements.txt
streamlit run app.py
```

---

## 7. How to Deploy to Streamlit Community Cloud

1. Go to share.streamlit.io
2. Connect GitHub repo `amiymh/microgravity-dashboard`
3. Set main file: `app.py`
4. Deploy — auto-deploys on every push to main
