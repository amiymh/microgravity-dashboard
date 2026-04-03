"""Publication-ready figure export utilities."""

import io


def fig_to_png(fig, width: int = 1200, height: int = 800) -> bytes | None:
    """Export a Plotly figure to PNG. Returns bytes or None."""
    if fig is None:
        return None
    try:
        return fig.to_image(format="png", width=width, height=height, scale=3)
    except Exception:
        return None


def fig_to_svg(fig, width: int = 1200, height: int = 800) -> bytes | None:
    """Export a Plotly figure to SVG. Returns bytes or None."""
    if fig is None:
        return None
    try:
        return fig.to_image(format="svg", width=width, height=height)
    except Exception:
        return None


def fig_to_html(fig) -> bytes | None:
    """Export a Plotly figure to self-contained HTML. Always works."""
    if fig is None:
        return None
    try:
        return fig.to_html(include_plotlyjs=True, full_html=True).encode("utf-8")
    except Exception:
        return None


def add_download_buttons(st_module, fig, prefix: str = "figure"):
    """Add download buttons for a Plotly figure in Streamlit."""
    col1, col2, col3 = st_module.columns(3)
    png_bytes = fig_to_png(fig)
    svg_bytes = fig_to_svg(fig)
    html_bytes = fig_to_html(fig)

    if png_bytes:
        col1.download_button(
            "Download PNG (300 DPI)",
            png_bytes,
            f"{prefix}.png",
            "image/png",
            key=f"dl_png_{prefix}",
        )
    elif html_bytes:
        col1.download_button(
            "Download HTML (interactive)",
            html_bytes,
            f"{prefix}.html",
            "text/html",
            key=f"dl_html_{prefix}",
        )

    if svg_bytes:
        col2.download_button(
            "Download SVG (vector)",
            svg_bytes,
            f"{prefix}.svg",
            "image/svg+xml",
            key=f"dl_svg_{prefix}",
        )

    if not png_bytes and not svg_bytes and html_bytes:
        col2.caption("PNG/SVG available when running locally with kaleido")
