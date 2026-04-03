"""Word methods report generator using python-docx."""

import io
from datetime import datetime
from docx import Document
from docx.shared import Pt, Inches
from docx.enum.text import WD_ALIGN_PARAGRAPH


def generate_report(
    session_id: str = "",
    filename: str = "unknown",
    sheet_name: str = "SDEGs",
    gene_count: int = 0,
    analysis_logs: dict[str, list[str]] | None = None,
    gene_lists: dict[str, list[str]] | None = None,
    versions: dict | None = None,
) -> bytes:
    """Generate a complete Word methods report.

    Args:
        session_id: UUID of the current session.
        filename: Name of the uploaded file.
        sheet_name: Which Excel sheet was used.
        gene_count: Total genes loaded.
        analysis_logs: Dict of {tab_name: [log_lines]} for each analysis run.
        gene_lists: Dict of {analysis_name: [gene_symbols]} for Appendix A.
        versions: Dict of library/API versions for Appendix C.

    Returns:
        bytes of the .docx file.
    """
    doc = Document()

    # --- Title Page ---
    title = doc.add_heading("Analysis Methods Report", level=0)
    title.alignment = WD_ALIGN_PARAGRAPH.CENTER
    doc.add_paragraph("Microgravity RNA Dashboard", style="Subtitle").alignment = WD_ALIGN_PARAGRAPH.CENTER
    doc.add_paragraph("")

    # Info table
    info_items = [
        ("Date Generated", datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
        ("File Analyzed", filename),
        ("Sheet", sheet_name),
        ("Total Genes Loaded", str(gene_count)),
        ("Session ID", session_id or "N/A"),
    ]
    table = doc.add_table(rows=len(info_items), cols=2)
    table.style = "Table Grid"
    for i, (key, val) in enumerate(info_items):
        table.cell(i, 0).text = key
        table.cell(i, 1).text = val

    doc.add_page_break()

    # --- Analysis Sections ---
    all_tabs = [
        "Volcano Plot", "PCA Plot", "Top DEGs Heatmap",
        "Pathway Enrichment", "Therapeutic Targets",
        "Biomarker Discovery", "Cytotoxicity & Apoptosis",
        "Disease Cross-Reference",
    ]

    logs = analysis_logs or {}

    for tab_name in all_tabs:
        doc.add_heading(tab_name, level=1)

        if tab_name in logs and logs[tab_name]:
            doc.add_heading("Analysis Details", level=2)
            for line in logs[tab_name]:
                doc.add_paragraph(line, style="List Bullet")

            doc.add_heading("To Reproduce Manually", level=2)
            doc.add_paragraph(
                "Refer to the Analysis Details above for exact parameters, "
                "API endpoints, and gene lists used. See Appendix A for full gene lists."
            )
        else:
            doc.add_paragraph("Not run in this session.", style="Intense Quote")

    doc.add_page_break()

    # --- Appendix A: Gene Lists ---
    doc.add_heading("Appendix A: Gene Lists", level=1)
    gene_lists = gene_lists or {}
    if gene_lists:
        for name, genes in gene_lists.items():
            doc.add_heading(name, level=2)
            doc.add_paragraph(", ".join(genes[:200]))
            if len(genes) > 200:
                doc.add_paragraph(f"... and {len(genes) - 200} more genes")
    else:
        doc.add_paragraph("No gene lists recorded in this session.")

    # --- Appendix B: API Responses ---
    doc.add_heading("Appendix B: API Response Summaries", level=1)
    doc.add_paragraph("API response details are included in each analysis section's log above.")

    # --- Appendix C: Versions ---
    doc.add_heading("Appendix C: Software and Database Versions", level=1)
    versions = versions or {}
    if versions:
        table = doc.add_table(rows=len(versions), cols=2)
        table.style = "Table Grid"
        for i, (key, val) in enumerate(versions.items()):
            table.cell(i, 0).text = str(key)
            table.cell(i, 1).text = str(val)
    else:
        doc.add_paragraph("Version information not available.")

    # --- Appendix D: References ---
    doc.add_heading("Appendix D: Method References", level=1)
    refs = [
        "DESeq2: Love MI et al. Genome Biology 2014. DOI: 10.1186/s13059-014-0550-8",
        "Enrichr: Kuleshov MV et al. Nucleic Acids Res 2016. DOI: 10.1093/nar/gkw377",
        "OpenTargets: Ochoa D et al. Nucleic Acids Res 2023. DOI: 10.1093/nar/gkac1033",
        "DGIdb: Freshour SL et al. Nucleic Acids Res 2021. DOI: 10.1093/nar/gkaa1084",
        "MSigDB: Liberzon A et al. Cell Systems 2015. DOI: 10.1016/j.cels.2015.12.004",
        "Fisher's Exact Test: Fisher RA. J Royal Stat Soc 1922.",
        "PCA: Jolliffe IT. Principal Component Analysis. Springer 2002.",
    ]
    for ref in refs:
        doc.add_paragraph(ref, style="List Number")

    # Save to bytes
    buf = io.BytesIO()
    doc.save(buf)
    buf.seek(0)
    return buf.getvalue()
