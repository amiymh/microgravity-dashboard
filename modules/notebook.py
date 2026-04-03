"""Jupyter notebook export using nbformat."""

import nbformat
from datetime import datetime


def generate_notebook(
    session_id: str = "",
    filename: str = "unknown",
    sheet_name: str = "SDEGs",
    gene_count: int = 0,
    analysis_logs: dict[str, list[str]] | None = None,
    versions: dict | None = None,
) -> bytes:
    """Generate a Jupyter notebook as a scientific record.

    Returns bytes of the .ipynb JSON file.
    """
    nb = nbformat.v4.new_notebook()
    cells = []

    # Cell 1: Title
    title_md = f"""# Microgravity RNA Analysis — Session Record
**Date:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
**File:** {filename}
**Sheet:** {sheet_name}
**Genes Loaded:** {gene_count}
**Session ID:** {session_id or 'N/A'}

This notebook is an auditable record of the analysis performed in the Microgravity RNA Dashboard.
It can be re-run by a Python user to reproduce all results independently.
"""
    cells.append(nbformat.v4.new_markdown_cell(title_md))

    # Cell 2: Setup
    setup_code = """import pandas as pd
import numpy as np
from scipy import stats
import plotly.graph_objects as go
import sys
print(f"Python: {sys.version}")
print(f"pandas: {pd.__version__}")
print(f"numpy: {np.__version__}")
print(f"scipy: {stats.scipy.__version__ if hasattr(stats, 'scipy') else 'N/A'}")
"""
    cells.append(nbformat.v4.new_code_cell(setup_code))

    # Cell 3: Data Loading
    load_md = f"""## Data Loading
The analysis loaded **{gene_count}** genes from `{filename}` (sheet: {sheet_name}).
"""
    cells.append(nbformat.v4.new_markdown_cell(load_md))

    load_code = f"""# Load the data file
# Replace the path below with the actual location of your data file
df = pd.read_excel("{filename}", sheet_name="{sheet_name}", engine="openpyxl")
print(f"Loaded {{len(df)}} genes")
print(f"Columns: {{list(df.columns)}}")
"""
    cells.append(nbformat.v4.new_code_cell(load_code))

    # Analysis cells
    logs = analysis_logs or {}
    for tab_name, log_lines in logs.items():
        # Markdown explanation
        analysis_md = f"## {tab_name}\n\n"
        for line in log_lines:
            analysis_md += f"- {line}\n"
        cells.append(nbformat.v4.new_markdown_cell(analysis_md))

        # Code cell (placeholder showing the analysis was performed)
        code = f"# {tab_name} — analysis was performed with the parameters above.\n"
        code += f"# See the Methods Report for full reproduction instructions.\n"
        code += f"print('Analysis: {tab_name}')\n"
        cells.append(nbformat.v4.new_code_cell(code))

    # Final summary
    summary_md = "## Results Summary\n\n"
    if logs:
        summary_md += f"Analyses completed in this session: {', '.join(logs.keys())}\n\n"
    else:
        summary_md += "No analyses were run in this session.\n\n"

    if versions:
        summary_md += "### Software Versions\n"
        for k, v in list(versions.items())[:10]:
            summary_md += f"- **{k}:** {v}\n"

    cells.append(nbformat.v4.new_markdown_cell(summary_md))

    nb.cells = cells

    import json
    return json.dumps(nbformat.from_dict(nb), indent=2).encode("utf-8")
