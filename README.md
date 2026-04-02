# Microgravity RNA Dashboard

A multi-tool Streamlit dashboard for analyzing differentially expressed genes (DEGs) under microgravity vs ground conditions, built on DESeq2 output data.

## Features

Six integrated analysis tools:

1. **Volcano Plot** — Interactive visualization of differential expression (log2FC vs significance)
2. **Pathway Enrichment** — Query Enrichr for KEGG, GO, Reactome, and WikiPathway enrichment
3. **Therapeutic Targets** — Identify drug-gene interactions via DGIdb
4. **Biomarker Discovery** — Find potential disease biomarkers using OpenTargets
5. **Cytotoxicity & Apoptosis** — Overlap analysis with curated hallmark gene sets (apoptosis, TNF, p53, ROS, inflammation)
6. **Disease Cross-Reference** — Search any disease and find overlap with microgravity DEGs via OpenTargets

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

Deployed on Streamlit Community Cloud:
<!-- Add URL after deployment -->

## Running Tests

```bash
pytest tests/ -v
```

## External APIs Used

- [Enrichr](https://maayanlab.cloud/Enrichr/) — Pathway enrichment (free, no key)
- [DGIdb](https://dgidb.org/) — Drug-gene interactions (free, no key)
- [OpenTargets](https://platform.opentargets.org/) — Disease associations and drug info (free, no key)

All API calls have 10-second timeouts and graceful fallbacks to demo data.
