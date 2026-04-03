"""Microgravity RNA Dashboard V2 — Main Streamlit Application."""

import streamlit as st
import pandas as pd

from modules.data_loader import load_excel, filter_degs, validate_columns
from modules.volcano import create_volcano_plot, get_top_significant, get_analysis_log as volcano_log
from modules.pca import run_pca, get_loadings, create_pca_plot, get_analysis_log as pca_log
from modules.heatmap import compute_heatmap_data, create_heatmap_figure, get_analysis_log as heatmap_log
from modules.pathways import run_enrichment, create_pathway_chart, GENE_SET_LIBRARIES, get_analysis_log as pathways_log
from modules.targets import get_therapeutic_targets, get_target_summary, get_analysis_log as targets_log
from modules.biomarkers import find_biomarkers, create_biomarker_scatter, get_analysis_log as biomarkers_log
from modules.cytotoxicity import (
    compute_overlaps, get_gene_pathway_matrix, create_overlap_bar_chart,
    create_heatmap, get_gene_table, get_available_pathway_names,
    get_analysis_log as cyto_log,
)
from modules.disease import (
    search_disease, cross_reference, create_venn_diagram, create_network_graph,
    get_analysis_log as disease_log,
)
from modules.fig_export import add_download_buttons, fig_to_png
from modules.versions import get_session_versions, get_api_versions, format_versions_text
from modules.report import generate_report
from modules.notebook import generate_notebook
from modules.zip_export import generate_zip

st.set_page_config(
    page_title="Microgravity RNA Dashboard",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# ── Session State Init ───────────────────────────────────────────────────────

if "analysis_logs" not in st.session_state:
    st.session_state["analysis_logs"] = {}
if "figures" not in st.session_state:
    st.session_state["figures"] = {}
if "data_csvs" not in st.session_state:
    st.session_state["data_csvs"] = {}
if "gene_lists" not in st.session_state:
    st.session_state["gene_lists"] = {}

versions = get_session_versions()
session_id = versions["session_id"]

# ── Sidebar ──────────────────────────────────────────────────────────────────

st.sidebar.title("Microgravity RNA Dashboard")
st.sidebar.caption("Multi-tool bioinformatics research platform")

# File upload
uploaded_file = st.sidebar.file_uploader("Upload DESeq2 file", type=["xlsx", "csv"])

# Sheet selector — only for Excel
_is_csv = uploaded_file is not None and uploaded_file.name.lower().endswith(".csv")
if _is_csv:
    selected_sheet = "SDEGs"
else:
    selected_sheet = st.sidebar.selectbox("Sheet", ["SDEGs", "SIGNFICANCE SCOREs Gene Types"], index=0)

# Load data
@st.cache_data
def cached_load(sheet_name, _file=None):
    if _file is not None:
        return load_excel(uploaded_file=_file, sheet_name=sheet_name)
    return load_excel(sheet_name=sheet_name)

if uploaded_file:
    df_raw = cached_load(selected_sheet, uploaded_file)
    data_source = f"Uploaded {'CSV' if _is_csv else 'Excel'}"
    data_filename = uploaded_file.name
else:
    df_raw = cached_load(selected_sheet)
    data_filename = "demo_data" if len(df_raw) == 50 else "local_file"
    data_source = "Demo data" if len(df_raw) == 50 else "Local file"

missing = validate_columns(df_raw)
if missing:
    st.sidebar.warning(f"Missing columns: {', '.join(missing)}")

# Global filters
st.sidebar.markdown("---")
st.sidebar.subheader("Filters")

# Dynamic padj presets based on data range
if "padj" in df_raw.columns:
    padj_min = float(df_raw["padj"].min())
    padj_max = float(df_raw["padj"].max())
    padj_presets = {"Relaxed (0.1)": 0.1, "Standard (0.05)": 0.05, "Strict (0.01)": 0.01, "Very strict (0.001)": 0.001}
    padj_choice = st.sidebar.selectbox("padj threshold", list(padj_presets.keys()), index=1)
    padj_thresh = padj_presets[padj_choice]
    st.sidebar.caption(f"Data range: {padj_min:.2e} – {padj_max:.2e}")
else:
    padj_thresh = 0.05

log2fc_thresh = st.sidebar.slider("|log2FC| threshold", 0.0, 5.0, 1.0, 0.1)

# Gene Type filter — disabled gracefully when column missing
_has_gene_type = "Gene Type" in df_raw.columns and df_raw["Gene Type"].notna().any()
if _has_gene_type:
    gene_types = sorted(df_raw["Gene Type"].dropna().unique())
    selected_types = st.sidebar.multiselect("Gene types", gene_types, default=gene_types)
    effective_types = selected_types if selected_types else gene_types
else:
    gene_types = []
    effective_types = None
    st.sidebar.text_input("Gene types", value="Not available", disabled=True)
    st.sidebar.caption("Gene Type filtering unavailable — column not found in this sheet")

selected_direction = st.sidebar.selectbox("Direction", ["All", "Upregulated", "Downregulated"])

df = filter_degs(df_raw, padj_thresh, log2fc_thresh, effective_types if effective_types else None, selected_direction)

# Data summary
st.sidebar.markdown("---")
st.sidebar.subheader("Data Summary")
st.sidebar.metric("Total genes (filtered)", len(df))
if "Direction" in df.columns:
    c1, c2 = st.sidebar.columns(2)
    c1.metric("Up", int((df["Direction"] == "Upregulated").sum()))
    c2.metric("Down", int((df["Direction"] == "Downregulated").sum()))
st.sidebar.caption(f"Source: {data_source}")

# Download buttons
st.sidebar.markdown("---")
st.sidebar.subheader("Downloads")

report_bytes = generate_report(
    session_id=session_id, filename=data_filename, sheet_name=selected_sheet,
    gene_count=len(df_raw), analysis_logs=st.session_state.get("analysis_logs"),
    gene_lists=st.session_state.get("gene_lists"), versions=versions,
)
st.sidebar.download_button("📄 Download Methods Report (.docx)", report_bytes, "methods_report.docx",
                           "application/vnd.openxmlformats-officedocument.wordprocessingml.document")

nb_bytes = generate_notebook(
    session_id=session_id, filename=data_filename, sheet_name=selected_sheet,
    gene_count=len(df_raw), analysis_logs=st.session_state.get("analysis_logs"),
    versions=versions,
)
st.sidebar.download_button("📓 Download Jupyter Notebook (.ipynb)", nb_bytes, "analysis_notebook.ipynb",
                           "application/json")

zip_bytes = generate_zip(
    session_id=session_id, report_bytes=report_bytes, notebook_bytes=nb_bytes,
    figures=st.session_state.get("figures"), data_csvs=st.session_state.get("data_csvs"),
    analyses_run=list(st.session_state.get("analysis_logs", {}).keys()),
)
st.sidebar.download_button("📦 Download All Results (.zip)", zip_bytes, "results.zip", "application/zip")

# Version info
with st.sidebar.expander("Version Info"):
    st.text(format_versions_text())

# ── Main Tabs ────────────────────────────────────────────────────────────────

tab1, tab2, tab3, tab4, tab5, tab6, tab7, tab8, tab9 = st.tabs([
    "🌋 Volcano Plot", "📊 PCA Plot", "🔥 Top DEGs Heatmap",
    "🧪 Pathway Enrichment", "💊 Therapeutic Targets", "🔬 Biomarker Discovery",
    "☠️ Cytotoxicity & Apoptosis", "🏥 Disease Cross-Reference", "📖 User Manual",
])

# ── Tab 1: Volcano ──────────────────────────────────────────────────────────

with tab1:
    st.header("Volcano Plot")
    vc1, vc2, vc3 = st.columns(3)
    v_padj = vc1.slider("padj cutoff", 0.001, 0.1, 0.05, 0.001, key="v_padj", format="%.3f")
    v_fc = vc2.slider("|log2FC| cutoff", 0.0, 5.0, 1.0, 0.1, key="v_fc")
    v_n = vc3.slider("Label top N genes", 0, 30, 10, key="v_n")
    v_color = st.selectbox("Color scheme", ["Default", "Warm", "Cool"], key="v_color")

    fig = create_volcano_plot(df_raw, v_padj, v_fc, v_n, v_color)
    st.plotly_chart(fig, use_container_width=True)
    add_download_buttons(st, fig, "volcano_plot")
    _png = fig_to_png(fig)
    if _png:
        st.session_state["figures"]["volcano_plot.png"] = _png

    top20 = get_top_significant(df_raw, 20)
    if not top20.empty:
        st.subheader("Top 20 Most Significant Genes")
        st.dataframe(top20, use_container_width=True, hide_index=True)
    csv = df.to_csv(index=False)
    st.download_button("Download filtered gene list (CSV)", csv, "filtered_genes.csv", "text/csv")
    st.session_state["data_csvs"]["filtered_genes.csv"] = csv.encode()

    log = volcano_log(df_raw, v_padj, v_fc, v_n)
    st.session_state["analysis_logs"]["Volcano Plot"] = log
    with st.expander("🔍 Analysis Details — click to verify"):
        for line in log:
            st.markdown(f"- {line}")

# ── Tab 2: PCA ──────────────────────────────────────────────────────────────

with tab2:
    st.header("PCA Plot")
    pca_result = run_pca(df_raw)

    if pca_result:
        fig = create_pca_plot(pca_result)
        st.plotly_chart(fig, use_container_width=True)
        add_download_buttons(st, fig, "pca_plot")
        _png = fig_to_png(fig)
        if _png:
            st.session_state["figures"]["pca_plot.png"] = _png

        st.subheader("Top Gene Loadings (PC1 & PC2)")
        loadings_df = get_loadings(df_raw, pca_result["pca_model"])
        if loadings_df is not None and not loadings_df.empty:
            st.dataframe(loadings_df, use_container_width=True, hide_index=True)

        log = pca_log(pca_result)
        st.session_state["analysis_logs"]["PCA Plot"] = log
        with st.expander("🔍 Analysis Details — click to verify"):
            for line in log:
                st.markdown(f"- {line}")
    else:
        st.info("PCA requires Earth and Space count columns. Upload data with normalized counts.")

# ── Tab 3: Heatmap ──────────────────────────────────────────────────────────

with tab3:
    st.header("Top DEGs Heatmap")
    hc1, hc2 = st.columns(2)
    h_top_n = hc1.slider("Top N genes", 10, 100, 50, key="h_n")
    h_cluster = hc2.checkbox("Cluster genes", value=True, key="h_cluster")

    hm_data = compute_heatmap_data(df, top_n=h_top_n, cluster=h_cluster)
    if hm_data:
        fig = create_heatmap_figure(hm_data)
        st.plotly_chart(fig, use_container_width=True)
        add_download_buttons(st, fig, "heatmap")
        _png = fig_to_png(fig)
        if _png:
            st.session_state["figures"]["heatmap.png"] = _png

        if "gene_info" in hm_data and not hm_data["gene_info"].empty:
            st.subheader("Gene Table")
            st.dataframe(hm_data["gene_info"], use_container_width=True, hide_index=True)
            csv = hm_data["gene_info"].to_csv(index=False)
            st.download_button("Download gene table (CSV)", csv, "heatmap_genes.csv", "text/csv")

        log = heatmap_log(hm_data, h_cluster)
        st.session_state["analysis_logs"]["Top DEGs Heatmap"] = log
        with st.expander("🔍 Analysis Details — click to verify"):
            for line in log:
                st.markdown(f"- {line}")
    else:
        st.info("Heatmap requires expression count columns. Upload data with Earth/Space columns.")

# ── Tab 4: Pathway Enrichment ───────────────────────────────────────────────

with tab4:
    st.header("Pathway Enrichment")
    pc1, pc2, pc3 = st.columns(3)
    p_geneset = pc1.selectbox("Gene set library", GENE_SET_LIBRARIES)
    p_dir = pc2.selectbox("Direction", ["All", "Upregulated", "Downregulated"], key="p_dir")
    p_top = pc3.slider("Top N results", 10, 50, 20, key="p_top")

    p_df = df.copy()
    if p_dir != "All" and "Direction" in p_df.columns:
        p_df = p_df[p_df["Direction"] == p_dir]
    gene_list = p_df["Gene"].dropna().unique().tolist() if "Gene" in p_df.columns else []

    if st.button("Run Enrichment Analysis", key="run_enrichr"):
        if not gene_list:
            st.warning("No genes available.")
        else:
            with st.spinner(f"Querying Enrichr ({p_geneset})..."):
                results, is_live = run_enrichment(gene_list, p_geneset, p_top)
            if not is_live:
                st.warning("Enrichr API unavailable — showing demo results.")
            if not results.empty:
                fig = create_pathway_chart(results, p_top)
                st.plotly_chart(fig, use_container_width=True)
                add_download_buttons(st, fig, "pathways")
                _png = fig_to_png(fig)
                if _png:
                    st.session_state["figures"]["pathways.png"] = _png
                st.dataframe(results, use_container_width=True, hide_index=True)
                csv = results.to_csv(index=False)
                st.download_button("Download pathway results (CSV)", csv, "pathways.csv", "text/csv")
                st.session_state["data_csvs"]["pathway_results.csv"] = csv.encode()

                log = pathways_log(len(gene_list), p_geneset, p_dir, results, is_live)
                st.session_state["analysis_logs"]["Pathway Enrichment"] = log
                with st.expander("🔍 Analysis Details — click to verify"):
                    for line in log:
                        st.markdown(f"- {line}")
            else:
                st.info("No enrichment results found.")
    else:
        st.info("Click 'Run Enrichment Analysis' to query Enrichr.")

# ── Tab 5: Therapeutic Targets ──────────────────────────────────────────────

with tab5:
    st.header("Therapeutic Targets")
    if st.button("Search Drug-Gene Interactions", key="run_targets"):
        with st.spinner("Querying DGIdb + OpenTargets..."):
            targets = get_therapeutic_targets(df, top_n=200)
            st.session_state["targets"] = targets

    targets = st.session_state.get("targets", pd.DataFrame())
    if not targets.empty:
        summary = get_target_summary(targets)
        mc1, mc2, mc3 = st.columns(3)
        mc1.metric("Genes with drugs", summary["genes_with_drugs"])
        mc2.metric("Unique drugs", summary["total_drugs"])
        mc3.metric("With clinical stage", summary["with_clinical_stage"])
        st.dataframe(targets, use_container_width=True, hide_index=True)
        csv = targets.to_csv(index=False)
        st.download_button("Download targets (CSV)", csv, "therapeutic_targets.csv", "text/csv")
        st.session_state["data_csvs"]["therapeutic_targets.csv"] = csv.encode()

        log = targets_log(200, targets)
        st.session_state["analysis_logs"]["Therapeutic Targets"] = log
        with st.expander("🔍 Analysis Details — click to verify"):
            for line in log:
                st.markdown(f"- {line}")
    else:
        st.info("Click 'Search Drug-Gene Interactions' to query DGIdb.")

# ── Tab 6: Biomarker Discovery ──────────────────────────────────────────────

with tab6:
    st.header("Biomarker Discovery")
    b_fc = st.slider("Min |log2FC|", 0.5, 5.0, 2.0, 0.1, key="b_fc")
    st.caption("Both upregulated and downregulated genes are included.")

    if st.button("Find Biomarker Candidates", key="run_bio"):
        bar = st.progress(0, text="Querying OpenTargets...")

        def _update(cur, tot):
            bar.progress(cur / tot, text=f"Querying... {cur}/{tot} genes")

        biomarkers = find_biomarkers(df, min_log2fc=b_fc, progress_callback=_update)
        bar.empty()
        st.session_state["biomarkers"] = biomarkers

    biomarkers = st.session_state.get("biomarkers", pd.DataFrame())
    if not biomarkers.empty:
        st.dataframe(biomarkers, use_container_width=True, hide_index=True)
        fig = create_biomarker_scatter(biomarkers)
        st.plotly_chart(fig, use_container_width=True)
        add_download_buttons(st, fig, "biomarkers")
        _png = fig_to_png(fig)
        if _png:
            st.session_state["figures"]["biomarkers.png"] = _png
        csv = biomarkers.to_csv(index=False)
        st.download_button("Download biomarkers (CSV)", csv, "biomarkers.csv", "text/csv")
        st.session_state["data_csvs"]["biomarkers.csv"] = csv.encode()

        log = biomarkers_log(len(df), len(biomarkers), b_fc, biomarkers)
        st.session_state["analysis_logs"]["Biomarker Discovery"] = log
        with st.expander("🔍 Analysis Details — click to verify"):
            for line in log:
                st.markdown(f"- {line}")
    else:
        st.info("Click 'Find Biomarker Candidates' to search.")

# ── Tab 7: Cytotoxicity ─────────────────────────────────────────────────────

with tab7:
    st.header("Cytotoxicity & Apoptosis Analysis")
    sel_pw = st.multiselect("Select pathways", get_available_pathway_names(),
                            default=get_available_pathway_names(), key="cyto_pw")
    if sel_pw:
        overlaps = compute_overlaps(df, sel_pw)
        if not overlaps.empty:
            st.dataframe(overlaps, use_container_width=True, hide_index=True)
            fig = create_overlap_bar_chart(overlaps)
            st.plotly_chart(fig, use_container_width=True)
            add_download_buttons(st, fig, "cytotoxicity_bar")
            _png = fig_to_png(fig)
            if _png:
                st.session_state["figures"]["cytotoxicity_bar.png"] = _png

            matrix = get_gene_pathway_matrix(df, sel_pw)
            if not matrix.empty:
                st.subheader("Gene × Pathway Heatmap")
                hfig = create_heatmap(matrix)
                st.plotly_chart(hfig, use_container_width=True)
                add_download_buttons(st, hfig, "cytotoxicity_heatmap")
                _png = fig_to_png(hfig)
                if _png:
                    st.session_state["figures"]["cytotoxicity_heatmap.png"] = _png

            gene_tbl = get_gene_table(df, sel_pw)
            if not gene_tbl.empty:
                st.subheader("Genes in Pathways")
                st.dataframe(gene_tbl, use_container_width=True, hide_index=True)
                csv = gene_tbl.to_csv(index=False)
                st.download_button("Download gene table (CSV)", csv, "cytotoxicity_genes.csv", "text/csv")
                st.session_state["data_csvs"]["cytotoxicity_genes.csv"] = csv.encode()

            log = cyto_log(overlaps, sel_pw)
            st.session_state["analysis_logs"]["Cytotoxicity & Apoptosis"] = log
            with st.expander("🔍 Analysis Details — click to verify"):
                for line in log:
                    st.markdown(f"- {line}")
        else:
            st.info("No overlapping genes found.")
    else:
        st.info("Select at least one pathway.")

# ── Tab 8: Disease Cross-Reference ──────────────────────────────────────────

with tab8:
    st.header("Disease Cross-Reference")
    disease_query = st.text_input("Search for a disease", placeholder="e.g. Alzheimer, breast cancer")
    d_score = st.slider("Association score threshold", 0.0, 1.0, 0.3, 0.05, key="d_score")
    d_dir = st.selectbox("Direction filter", ["All", "Upregulated", "Downregulated"], key="d_dir")

    if disease_query:
        with st.spinner(f"Searching '{disease_query}'..."):
            diseases = search_disease(disease_query)
        if diseases:
            disease_names = {d["name"]: d["id"] for d in diseases}
            sel_disease = st.selectbox("Select disease", list(disease_names.keys()))
            disease_id = disease_names[sel_disease]

            if st.button("Run Cross-Reference", key="run_disease"):
                with st.spinner("Querying OpenTargets..."):
                    overlap_df, summary = cross_reference(df, disease_id, d_score, d_dir)
                    st.session_state["disease_result"] = (overlap_df, summary, sel_disease, disease_id)

            if "disease_result" in st.session_state:
                overlap_df, summary, disease_name, did = st.session_state["disease_result"]
                sc1, sc2, sc3 = st.columns(3)
                sc1.metric("Overlapping genes", summary["overlap"])
                sc2.metric("Disease genes", summary["disease_genes"])
                sc3.metric("DEGs queried", summary["deg_total"])

                if summary["overlap"] > 0:
                    venn_fig = create_venn_diagram(summary, disease_name)
                    st.pyplot(venn_fig)
                    import matplotlib.pyplot as plt
                    plt.close(venn_fig)

                    st.dataframe(overlap_df, use_container_width=True, hide_index=True)
                    net = create_network_graph(overlap_df, disease_name)
                    st.plotly_chart(net, use_container_width=True)
                    add_download_buttons(st, net, "disease_network")
                    _png = fig_to_png(net)
                    if _png:
                        st.session_state["figures"]["disease_network.png"] = _png
                    csv = overlap_df.to_csv(index=False)
                    st.download_button("Download overlap (CSV)", csv, "disease_overlap.csv", "text/csv")
                    st.session_state["data_csvs"]["disease_overlap.csv"] = csv.encode()
                else:
                    st.info(f"No overlap found with {disease_name}.")

                log = disease_log(disease_name, did, d_score, d_dir, summary)
                st.session_state["analysis_logs"]["Disease Cross-Reference"] = log
                with st.expander("🔍 Analysis Details — click to verify"):
                    for line in log:
                        st.markdown(f"- {line}")
        else:
            st.warning("No diseases found. Try a different term.")

# ── Tab 9: User Manual ──────────────────────────────────────────────────────

with tab9:
    st.header("📖 User Manual")

    st.subheader("What This Tool Does")
    st.markdown("""
This dashboard analyzes RNA sequencing data to understand how gene expression changes under
microgravity (space) compared to normal gravity (Earth). You upload your DESeq2 results,
and the tool provides interactive visualizations, pathway analysis, drug target identification,
and disease comparisons — all without writing any code.
""")

    st.subheader("Your Data")
    st.markdown("""
- **Accepted formats:** Excel (`.xlsx`) with multiple sheets, or CSV (`.csv`) for single datasets
- **Required columns:** `Gene`, `padj`, `log2FoldChange`, `Direction` — other columns are optional but enhance the analysis
- **Missing columns:** The tool auto-fills `Direction` from `log2FoldChange` and computes `significance score` if absent. Missing count columns disable PCA and Heatmap tabs.
- **Demo mode:** If no file is uploaded, the tool runs with 50 synthetic genes so you can explore every feature
""")

    st.subheader("Understanding the Sidebar")

    st.markdown("**Uploading Your File**")
    st.markdown("""
Click "Browse files" in the sidebar. Select your Excel (.xlsx) or CSV (.csv) file. The app loads it
automatically and updates all tabs. If you do not upload anything, the app runs on demo data — 50
synthetic genes — so you can explore every feature before using your own data.
""")

    st.markdown("**Sheet Selector (Excel only)**")
    st.markdown("""
If you uploaded an Excel file, you will see a sheet selector. Your file may have multiple sheets:

- **SDEGs** is the recommended sheet. It contains only the statistically significant genes and is ready
  to analyze immediately.
- **SIGNFICANCE SCOREs Gene Types** contains all genes including non-significant ones — useful when you
  want to see the full picture.

If you select a sheet that is missing certain columns (like Gene Type), the app fills them in
automatically from the SDEGs sheet if possible, or disables that filter with a note. This selector
does not appear for CSV files since CSVs have only one sheet.
""")

    st.markdown("**padj Threshold**")
    st.markdown("""
padj is the *adjusted p-value*. It measures how statistically confident we are that a gene is truly
differentially expressed and not just random noise. A lower padj means more confidence.

The app shows preset options based on your data range:
- **Relaxed (0.1)** — includes more genes, some may be borderline
- **Standard (0.05)** — the most common threshold in research papers
- **Strict (0.01)** — only genes with strong statistical evidence
- **Very strict (0.001)** — only the most confident results

Start with Standard. Tighten it if you want only the most reliable genes.
""")

    st.markdown("**|log2FC| Threshold**")
    st.markdown("""
log2 Fold Change measures how much a gene's expression changed between Earth and Space:
- A value of **1.0** means the gene doubled (or halved) its expression
- A value of **2.0** means it quadrupled (or quartered)
- A value of **0** means no change

A higher threshold shows only large expression changes. Start at 1.0 for a balanced view.
""")

    st.markdown("**Gene Types**")
    st.markdown("""
Genes in your data are classified by type:
- **protein_coding** — genes that make proteins directly (the majority of known functional genes)
- **lncRNA** — long non-coding RNA that regulates other genes without making proteins
- **processed_pseudogene** — inactive copies of genes, usually not functional

By default all types are shown. Deselect types you are not interested in.
If this filter is greyed out, it means Gene Type information was not available in your selected sheet.
""")

    st.markdown("**Direction**")
    st.markdown("""
- **Upregulated** means the gene is *more active* in Space than on Earth
- **Downregulated** means the gene is *less active* in Space
- **All** shows both

Use this filter to focus on one direction. For example, if you are looking for genes activated as a
stress response in microgravity, filter to Upregulated.
""")

    st.markdown("**Data Summary**")
    st.markdown("""
Shows how many genes passed all your current filters:
- **Total genes** is the number being analyzed right now across all tabs
- **Up** and **Down** show the split between upregulated and downregulated
- **Source** confirms whether you are using uploaded data or demo data

If the number seems too low, check your filters — especially the padj threshold.
""")

    st.markdown("**Downloads**")
    st.markdown("""
Three download buttons at the bottom of the sidebar:

- **Methods Report (.docx)** — a Word document with the full scientific methods for every analysis
  you ran in this session. Written so another scientist can reproduce your results manually, without
  using this tool. Includes exact URLs, gene lists, parameters, and expected outputs.

- **Jupyter Notebook (.ipynb)** — a notebook containing all the code and outputs from your session.
  Can be re-run by a Python user to get identical results. Serves as an auditable scientific record.

- **Download All Results (.zip)** — everything in one file: the methods report, the notebook, all
  figures as PNG (300 DPI, publication-ready) and SVG (vector), and all result tables as CSV files.
""")

    st.markdown("**Version Info**")
    st.markdown("""
Click the Version Info expander at the very bottom of the sidebar to see the exact version of every
software library and database used in this session. This is important for reproducibility — if you
publish results, include this information in your methods section so others know the exact software
environment you used.
""")

    st.subheader("How to Use Each Tab")

    st.markdown("**🌋 Volcano Plot** — Overview of all differentially expressed genes. Red dots are upregulated, blue are downregulated. Adjust thresholds to change what counts as significant.")
    st.markdown("**📊 PCA Plot** — Shows how your Earth and Space samples cluster. If Earth samples group together and Space samples group together, your experiment worked well.")
    st.markdown("**🔥 Top DEGs Heatmap** — Expression patterns of the most significant genes across all samples. Blue means low expression, red means high. Clustering groups genes with similar patterns.")
    st.markdown("**🧪 Pathway Enrichment** — Identifies biological pathways enriched in your DEGs using the Enrichr database. Click 'Run' to query the API.")
    st.markdown("**💊 Therapeutic Targets** — Finds which DEGs are known drug targets using DGIdb and OpenTargets databases.")
    st.markdown("**🔬 Biomarker Discovery** — Identifies potential disease biomarkers by querying OpenTargets for disease associations of your DEGs.")
    st.markdown("**☠️ Cytotoxicity & Apoptosis** — Checks overlap between your DEGs and known cell death/inflammation gene sets from MSigDB.")
    st.markdown("**🏥 Disease Cross-Reference** — Search any disease and find which of your DEGs are associated with it.")

    st.subheader("Understanding the Black Box — Transparency and Reproducibility")
    st.markdown("""
**Every analysis in this tool is fully transparent.** Nothing is hidden.

- **Analysis Details panel:** Every tab has a collapsible "🔍 Analysis Details" section at the bottom. Click it to see every parameter, threshold, API call, and result count.
- **Methods Report (.docx):** Download from the sidebar. Contains complete documentation of every analysis, including exact parameters, statistical methods, formulas, and step-by-step "To Reproduce Manually" instructions.
- **Jupyter Notebook (.ipynb):** Download from the sidebar. A scientific record that a Python-literate colleague can re-run to reproduce your exact results.
- **"To Reproduce Manually":** Each section in the Methods Report includes numbered instructions with exact URLs, gene lists, and parameters so you can independently verify any result outside this tool.
- **External databases used:**
  - [Enrichr](https://maayanlab.cloud/Enrichr/) — pathway enrichment
  - [DGIdb](https://dgidb.org/) — drug-gene interactions
  - [OpenTargets](https://platform.opentargets.org/) — disease associations and drug info
  - [MSigDB](https://www.gsea-msigdb.org/) — curated gene sets
- **Verification example:** To verify a Disease Cross-Reference result, go to [OpenTargets](https://platform.opentargets.org/), search for the disease, browse associated targets, and compare with the overlap genes shown in this tool.
""")

    st.subheader("Limitations")
    st.markdown("""
- **Enrichr vs GSEA:** Enrichr is a fast exploratory enrichment tool. For publication, validate results with full Gene Set Enrichment Analysis (GSEA).
- **Drug approval status:** Only clinical stage data from OpenTargets is shown. FDA approval cannot be reliably determined from free APIs alone.
- **API dependencies:** Enrichr, DGIdb, OpenTargets, and MSigDB are external services. If unavailable, fallback data is used and labeled.
- **Human genes only:** All databases are human-focused. Mouse or other organism gene symbols will return no results.
- **Biomarker processing time:** All qualifying genes are processed (no cap). Large gene sets may take several minutes. Use filters to narrow first.
- **Demo mode:** Uses 50 synthetic genes — suitable for exploring the interface, not for real analysis.
""")

    st.subheader("Getting Help")
    st.markdown("""
- **Unexpected results?** Check the Analysis Details panel for the exact parameters used
- **Verify independently:** Use the "To Reproduce Manually" section in the Methods Report
- **Technical issues:** Check that your file has the required columns and that you have internet access for API-dependent tabs
""")
