"""Microgravity RNA Dashboard — Main Streamlit Application."""

import streamlit as st
import pandas as pd

from modules.data_loader import load_excel, filter_degs, validate_columns, generate_demo_data
from modules.volcano import create_volcano_plot, get_top_significant
from modules.pathways import run_enrichment, create_pathway_chart, GENE_SET_LIBRARIES
from modules.targets import get_therapeutic_targets, get_target_summary
from modules.biomarkers import find_biomarkers, create_biomarker_scatter
from modules.cytotoxicity import (
    compute_overlaps, get_gene_pathway_matrix, create_overlap_bar_chart,
    create_heatmap, get_gene_table, get_available_pathway_names,
)
from modules.disease import (
    search_disease, cross_reference, create_venn_diagram, create_network_graph,
)

st.set_page_config(
    page_title="Microgravity RNA Dashboard",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# ── Sidebar ──────────────────────────────────────────────────────────────────

st.sidebar.title("Microgravity RNA Dashboard")
st.sidebar.markdown("Analyze differentially expressed genes under microgravity vs ground conditions.")

# File upload
uploaded_file = st.sidebar.file_uploader("Upload DESeq2 file", type=["xlsx", "csv"])

# Sheet selector
sheet_options = ["SDEGs", "SIGNFICANCE SCOREs Gene Types"]
selected_sheet = st.sidebar.selectbox("Sheet", sheet_options, index=0)

# Load data
@st.cache_data
def cached_load(sheet_name, _file=None):
    if _file is not None:
        return load_excel(uploaded_file=_file, sheet_name=sheet_name)
    return load_excel(sheet_name=sheet_name)

if uploaded_file:
    df_raw = cached_load(selected_sheet, uploaded_file)
    data_source = "Uploaded file"
else:
    df_raw = cached_load(selected_sheet)
    if len(df_raw) == 50:
        data_source = "Demo data (upload Excel for real data)"
    else:
        data_source = "Local file"

# Validate
missing = validate_columns(df_raw)
if missing:
    st.sidebar.warning(f"Missing columns: {', '.join(missing)}")

# Global filters
st.sidebar.markdown("---")
st.sidebar.subheader("Filters")

padj_thresh = st.sidebar.slider("padj threshold", 0.001, 0.1, 0.05, 0.001, format="%.3f")
log2fc_thresh = st.sidebar.slider("|log2FC| threshold", 0.0, 5.0, 1.0, 0.1)

gene_types = sorted(df_raw["Gene Type"].dropna().unique()) if "Gene Type" in df_raw.columns else []
selected_types = st.sidebar.multiselect("Gene types", gene_types, default=gene_types)

direction_opts = ["All", "Upregulated", "Downregulated"]
selected_direction = st.sidebar.selectbox("Direction", direction_opts)

# Apply filters
df = filter_degs(df_raw, padj_thresh, log2fc_thresh, selected_types or None, selected_direction)

# Summary
st.sidebar.markdown("---")
st.sidebar.subheader("Data Summary")
st.sidebar.metric("Total genes (filtered)", len(df))
if "Direction" in df.columns:
    up_count = (df["Direction"] == "Upregulated").sum()
    down_count = (df["Direction"] == "Downregulated").sum()
    col1, col2 = st.sidebar.columns(2)
    col1.metric("Upregulated", up_count)
    col2.metric("Downregulated", down_count)
st.sidebar.caption(f"Source: {data_source}")

# ── Main Area ────────────────────────────────────────────────────────────────

tab1, tab2, tab3, tab4, tab5, tab6 = st.tabs([
    "🌋 Volcano Plot",
    "🧪 Pathway Enrichment",
    "💊 Therapeutic Targets",
    "🔬 Biomarker Discovery",
    "☠️ Cytotoxicity & Apoptosis",
    "🏥 Disease Cross-Reference",
])

# ── Tab 1: Volcano Plot ─────────────────────────────────────────────────────

with tab1:
    st.header("Volcano Plot")

    vcol1, vcol2, vcol3 = st.columns(3)
    v_padj = vcol1.slider("padj cutoff", 0.001, 0.1, 0.05, 0.001, key="v_padj", format="%.3f")
    v_log2fc = vcol2.slider("|log2FC| cutoff", 0.0, 5.0, 1.0, 0.1, key="v_fc")
    v_labels = vcol3.slider("Label top N genes", 0, 30, 10, key="v_labels")
    v_color = st.selectbox("Color scheme", ["Default", "Warm", "Cool"], key="v_color")

    fig = create_volcano_plot(df_raw, v_padj, v_log2fc, v_labels, v_color)
    st.plotly_chart(fig, use_container_width=True)

    st.subheader("Top 20 Most Significant Genes")
    top20 = get_top_significant(df_raw, 20)
    if not top20.empty:
        st.dataframe(top20, use_container_width=True, hide_index=True)
    else:
        st.info("No significant genes found with current filters.")

    # Download
    csv = df.to_csv(index=False)
    st.download_button("Download filtered gene list (CSV)", csv, "filtered_genes.csv", "text/csv")

# ── Tab 2: Pathway Enrichment ───────────────────────────────────────────────

with tab2:
    st.header("Pathway Enrichment")

    pcol1, pcol2, pcol3 = st.columns(3)
    p_geneset = pcol1.selectbox("Gene set library", GENE_SET_LIBRARIES)
    p_direction = pcol2.selectbox("Direction", ["All", "Upregulated", "Downregulated"], key="p_dir")
    p_top_n = pcol3.slider("Top N results", 10, 50, 20, key="p_top")

    # Filter genes by direction
    p_df = df.copy()
    if p_direction != "All" and "Direction" in p_df.columns:
        p_df = p_df[p_df["Direction"] == p_direction]

    gene_list = p_df["Gene"].dropna().unique().tolist() if "Gene" in p_df.columns else []

    if st.button("Run Enrichment Analysis", key="run_enrichr"):
        if not gene_list:
            st.warning("No genes available for enrichment analysis.")
        else:
            with st.spinner(f"Querying Enrichr ({p_geneset})..."):
                results, is_live = run_enrichment(gene_list, p_geneset, p_top_n)

            if not is_live:
                st.warning("Enrichr API unavailable — showing demo results.")

            if not results.empty:
                chart = create_pathway_chart(results, p_top_n)
                st.plotly_chart(chart, use_container_width=True)

                st.subheader("Results Table")
                st.dataframe(results, use_container_width=True, hide_index=True)

                csv = results.to_csv(index=False)
                st.download_button("Download pathway results (CSV)", csv, "pathways.csv", "text/csv")
            else:
                st.info("No enrichment results found.")
    else:
        st.info("Click 'Run Enrichment Analysis' to query Enrichr.")

# ── Tab 3: Therapeutic Targets ──────────────────────────────────────────────

with tab3:
    st.header("Therapeutic Targets")

    if st.button("Search Drug-Gene Interactions", key="run_targets"):
        with st.spinner("Querying DGIdb for drug interactions..."):
            targets = get_therapeutic_targets(df, top_n=200)
            st.session_state["targets"] = targets

    targets = st.session_state.get("targets", pd.DataFrame())

    if not targets.empty:
        summary = get_target_summary(targets)
        mcol1, mcol2, mcol3 = st.columns(3)
        mcol1.metric("Genes with drug interactions", summary["genes_with_drugs"])
        mcol2.metric("Unique drugs found", summary["total_drugs"])
        mcol3.metric("With clinical stage data", summary["with_clinical_stage"])

        st.subheader("Drug-Gene Interactions")
        st.dataframe(targets, use_container_width=True, hide_index=True)

        csv = targets.to_csv(index=False)
        st.download_button("Download targets (CSV)", csv, "therapeutic_targets.csv", "text/csv")
    else:
        st.info("Click 'Search Drug-Gene Interactions' to query DGIdb.")

# ── Tab 4: Biomarker Discovery ──────────────────────────────────────────────

with tab4:
    st.header("Biomarker Discovery")

    bcol1, bcol2 = st.columns(2)
    b_min_fc = bcol1.slider("Min |log2FC|", 0.5, 5.0, 2.0, 0.1, key="b_fc")
    b_min_sig = bcol2.slider("Min significance score", 0.0, 100.0, 0.0, 1.0, key="b_sig")

    if st.button("Find Biomarker Candidates", key="run_biomarkers"):
        progress_bar = st.progress(0, text="Querying OpenTargets...")

        def update_progress(current, total):
            progress_bar.progress(current / total, text=f"Querying OpenTargets... {current}/{total} genes")

        biomarkers = find_biomarkers(df, min_log2fc=b_min_fc, min_sig_score=b_min_sig, progress_callback=update_progress)
        progress_bar.empty()
        st.session_state["biomarkers"] = biomarkers

    biomarkers = st.session_state.get("biomarkers", pd.DataFrame())

    if not biomarkers.empty:
        st.subheader("Ranked Biomarker Candidates")
        st.dataframe(biomarkers, use_container_width=True, hide_index=True)

        scatter = create_biomarker_scatter(biomarkers)
        st.plotly_chart(scatter, use_container_width=True)

        csv = biomarkers.to_csv(index=False)
        st.download_button("Download biomarkers (CSV)", csv, "biomarkers.csv", "text/csv")
    else:
        st.info("Click 'Find Biomarker Candidates' to search for potential biomarkers.")

# ── Tab 5: Cytotoxicity & Apoptosis ─────────────────────────────────────────

with tab5:
    st.header("Cytotoxicity & Apoptosis Analysis")

    selected_pathways = st.multiselect(
        "Select pathways",
        get_available_pathway_names(),
        default=get_available_pathway_names(),
        key="cyto_pw",
    )

    if selected_pathways:
        overlaps = compute_overlaps(df, selected_pathways)
        if not overlaps.empty:
            st.subheader("Pathway Overlap Summary")
            st.dataframe(overlaps, use_container_width=True, hide_index=True)

            bar_chart = create_overlap_bar_chart(overlaps)
            st.plotly_chart(bar_chart, use_container_width=True)

            # Heatmap
            matrix = get_gene_pathway_matrix(df, selected_pathways)
            if not matrix.empty:
                st.subheader("Gene × Pathway Heatmap")
                hmap = create_heatmap(matrix)
                st.plotly_chart(hmap, use_container_width=True)

            # Gene table
            gene_tbl = get_gene_table(df, selected_pathways)
            if not gene_tbl.empty:
                st.subheader("Genes in Cytotoxicity/Apoptosis Pathways")
                st.dataframe(gene_tbl, use_container_width=True, hide_index=True)

                csv = gene_tbl.to_csv(index=False)
                st.download_button("Download gene table (CSV)", csv, "cytotoxicity_genes.csv", "text/csv")
        else:
            st.info("No overlapping genes found.")
    else:
        st.info("Select at least one pathway to analyze.")

# ── Tab 6: Disease Cross-Reference ──────────────────────────────────────────

with tab6:
    st.header("Disease Cross-Reference")

    disease_query = st.text_input("Search for a disease", placeholder="e.g. Alzheimer, breast cancer, diabetes")
    d_score_thresh = st.slider("Association score threshold", 0.0, 1.0, 0.3, 0.05, key="d_score")
    d_direction = st.selectbox("Direction filter", ["All", "Upregulated", "Downregulated"], key="d_dir")

    if disease_query:
        with st.spinner(f"Searching for '{disease_query}'..."):
            diseases = search_disease(disease_query)

        if diseases:
            disease_names = {d["name"]: d["id"] for d in diseases}
            selected_disease = st.selectbox("Select disease", list(disease_names.keys()))
            disease_id = disease_names[selected_disease]

            if st.button("Run Cross-Reference", key="run_disease"):
                with st.spinner("Querying OpenTargets..."):
                    overlap_df, summary = cross_reference(df, disease_id, d_score_thresh, d_direction)
                    st.session_state["disease_result"] = (overlap_df, summary, selected_disease)

            if "disease_result" in st.session_state:
                overlap_df, summary, disease_name = st.session_state["disease_result"]

                st.subheader("Overlap Summary")
                scol1, scol2, scol3 = st.columns(3)
                scol1.metric("Overlapping genes", summary["overlap"])
                scol2.metric("Disease genes", summary["disease_genes"])
                scol3.metric("DEGs queried", summary["deg_total"])

                if summary["overlap"] > 0:
                    # Venn diagram
                    venn_fig = create_venn_diagram(summary, disease_name)
                    st.pyplot(venn_fig)
                    import matplotlib.pyplot as plt
                    plt.close(venn_fig)

                    st.subheader("Overlapping Genes")
                    st.dataframe(overlap_df, use_container_width=True, hide_index=True)

                    # Network graph
                    network = create_network_graph(overlap_df, disease_name)
                    st.plotly_chart(network, use_container_width=True)

                    csv = overlap_df.to_csv(index=False)
                    st.download_button("Download overlap (CSV)", csv, "disease_overlap.csv", "text/csv")
                else:
                    st.info(f"No overlapping genes found between DEGs and {disease_name}.")
        else:
            st.warning("No diseases found. Try a different search term.")
