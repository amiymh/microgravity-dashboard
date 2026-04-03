# Microgravity RNA Dashboard

A multi-tool Streamlit dashboard for analyzing differentially expressed genes (DEGs) under microgravity vs ground conditions, built on DESeq2 output data.

**Live app:** [https://microgravity-rna.streamlit.app](https://microgravity-rna.streamlit.app)

## Features

Nine integrated analysis tools:

1. **Volcano Plot** — Interactive visualization with adjustText non-overlapping labels
2. **PCA Plot** — Sample clustering (Earth vs Space) with gene loadings
3. **Top DEGs Heatmap** — Z-score normalized expression heatmap with hierarchical clustering
4. **Pathway Enrichment** — Query Enrichr for KEGG, GO, Reactome, and WikiPathway enrichment
5. **Therapeutic Targets** — Identify drug-gene interactions via DGIdb + OpenTargets
6. **Biomarker Discovery** — Find potential disease biomarkers using OpenTargets
7. **Cytotoxicity & Apoptosis** — Overlap with MSigDB hallmark gene sets (live API)
8. **Disease Cross-Reference** — Search any disease and find overlap with microgravity DEGs
9. **User Manual** — Full documentation and transparency guide

Plus: Methods Report (.docx), Jupyter Notebook (.ipynb), and full Results ZIP downloads.

## Setup

```bash
# Clone
git clone https://github.com/amiymh/microgravity-dashboard.git
cd microgravity-dashboard

# Install dependencies
pip install -r requirements.txt

# Run
streamlit run app.py
```

## Loading Data

- **Upload via sidebar**: Use the file uploader to load any DESeq2 Excel file (.xlsx)
- **Expected format**: Columns must include Gene, padj, log2FoldChange, Direction, Gene Type
- **Demo mode**: If no file is uploaded, the app runs with a synthetic 50-gene dataset

The primary data file is: `Supplementary Data 1 DESEQ2 normalized out and significance scores.xlsx` (SDEGs sheet, 4,522 genes).

## Screenshot

<!-- Add screenshot after deployment -->

## Deployment

Live at: [https://microgravity-rna.streamlit.app](https://microgravity-rna.streamlit.app)

Auto-deploys on every push to `main`.

## Running Tests

```bash
pytest tests/ -v
```

## External APIs Used

- [Enrichr](https://maayanlab.cloud/Enrichr/) — Pathway enrichment (free, no key)
- [DGIdb](https://dgidb.org/) — Drug-gene interactions (free, no key)
- [OpenTargets](https://platform.opentargets.org/) — Disease associations and drug info (free, no key)

All API calls have 10-second timeouts and graceful fallbacks to demo data.
