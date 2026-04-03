"""Publication-ready figure export utilities (PNG 300 DPI, SVG)."""

import io

try:
    import kaleido  # noqa: F401 — needed by plotly for write_image
    _HAS_KALEIDO = True
except ImportError:
    _HAS_KALEIDO = False


def fig_to_png(fig, width: int = 1200, height: int = 800) -> bytes | None:
    """Export a Plotly figure to PNG at 300 DPI (scale=3).

    Returns bytes or None if kaleido is unavailable.
    """
    if not _HAS_KALEIDO or fig is None:
        return None
    try:
        return fig.to_image(format="png", width=width, height=height, scale=3)
    except Exception:
        return None


def fig_to_svg(fig, width: int = 1200, height: int = 800) -> bytes | None:
    """Export a Plotly figure to SVG vector format.

    Returns bytes or None if kaleido is unavailable.
    """
    if not _HAS_KALEIDO or fig is None:
        return None
    try:
        return fig.to_image(format="svg", width=width, height=height)
    except Exception:
        return None


def add_download_buttons(st_module, fig, prefix: str = "figure"):
    """Add PNG and SVG download buttons for a Plotly figure in Streamlit.

    Args:
        st_module: The streamlit module (pass `st` from caller).
        fig: Plotly Figure object.
        prefix: Filename prefix for downloads.
    """
    col1, col2 = st_module.columns(2)
    png_bytes = fig_to_png(fig)
    svg_bytes = fig_to_svg(fig)

    if png_bytes:
        col1.download_button(
            "Download PNG (300 DPI)",
            png_bytes,
            f"{prefix}.png",
            "image/png",
            key=f"dl_png_{prefix}",
        )
    else:
        col1.caption("PNG export unavailable (install kaleido)")

    if svg_bytes:
        col2.download_button(
            "Download SVG (vector)",
            svg_bytes,
            f"{prefix}.svg",
            "image/svg+xml",
            key=f"dl_svg_{prefix}",
        )
    else:
        col2.caption("SVG export unavailable (install kaleido)")
