"""Full results ZIP package export."""

import io
import zipfile
from datetime import datetime


def generate_zip(
    session_id: str = "",
    report_bytes: bytes | None = None,
    notebook_bytes: bytes | None = None,
    figures: dict[str, bytes] | None = None,
    data_csvs: dict[str, bytes] | None = None,
    analyses_run: list[str] | None = None,
) -> bytes:
    """Generate a ZIP file containing all session results.

    Args:
        session_id: Session UUID for folder naming.
        report_bytes: Word report .docx bytes.
        notebook_bytes: Jupyter notebook .ipynb bytes.
        figures: Dict of {filename: PNG/SVG bytes}.
        data_csvs: Dict of {filename: CSV bytes}.
        analyses_run: List of tab names that were actually run.

    Returns:
        bytes of the ZIP file.
    """
    folder = f"results_{session_id}" if session_id else "results"
    buf = io.BytesIO()

    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        # Methods report
        if report_bytes:
            zf.writestr(f"{folder}/methods_report.docx", report_bytes)

        # Notebook
        if notebook_bytes:
            zf.writestr(f"{folder}/analysis_notebook.ipynb", notebook_bytes)

        # Figures
        figures = figures or {}
        for fname, fbytes in figures.items():
            zf.writestr(f"{folder}/figures/{fname}", fbytes)

        # Data CSVs
        data_csvs = data_csvs or {}
        for fname, fbytes in data_csvs.items():
            zf.writestr(f"{folder}/data/{fname}", fbytes)

        # README
        readme = _build_readme(session_id, analyses_run or [], figures, data_csvs, report_bytes is not None, notebook_bytes is not None)
        zf.writestr(f"{folder}/README.txt", readme)

    buf.seek(0)
    return buf.getvalue()


def _build_readme(
    session_id: str,
    analyses_run: list[str],
    figures: dict,
    data_csvs: dict,
    has_report: bool,
    has_notebook: bool,
) -> str:
    lines = [
        "Microgravity RNA Dashboard — Results Package",
        f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        f"Session ID: {session_id or 'N/A'}",
        "",
        "## Contents",
        "",
    ]

    if has_report:
        lines.append("- methods_report.docx — Full analysis methods report with reproduction instructions")
    if has_notebook:
        lines.append("- analysis_notebook.ipynb — Jupyter notebook for reproducing analyses")

    if figures:
        lines.append("")
        lines.append("## Figures")
        for fname in sorted(figures.keys()):
            lines.append(f"- figures/{fname}")

    if data_csvs:
        lines.append("")
        lines.append("## Data Files")
        for fname in sorted(data_csvs.keys()):
            lines.append(f"- data/{fname}")

    lines.append("")
    lines.append("## Analyses Run")
    if analyses_run:
        for name in analyses_run:
            lines.append(f"- {name}")
    else:
        lines.append("- No analyses were run in this session")

    not_run = [
        t for t in [
            "Volcano Plot", "PCA Plot", "Top DEGs Heatmap",
            "Pathway Enrichment", "Therapeutic Targets",
            "Biomarker Discovery", "Cytotoxicity & Apoptosis",
            "Disease Cross-Reference",
        ] if t not in (analyses_run or [])
    ]
    if not_run:
        lines.append("")
        lines.append("## Not Run (excluded from this package)")
        for name in not_run:
            lines.append(f"- {name}")

    return "\n".join(lines)
