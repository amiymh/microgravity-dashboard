"""Version tracking module for session and dependency auditing."""

import sys
import uuid
from datetime import datetime, timezone

# ── Session constants (generated once at module load) ─────────────────────────
_SESSION_ID = str(uuid.uuid4())
_SESSION_START = datetime.now(timezone.utc).isoformat()


def _safe_version(import_name: str, attr: str = "__version__") -> str:
    """Import *import_name* and return its version string, or 'not installed'."""
    try:
        mod = __import__(import_name)
        return str(getattr(mod, attr, "unknown"))
    except Exception:
        return "not installed"


# ── Public API ────────────────────────────────────────────────────────────────

def get_session_versions() -> dict:
    """Return a dict of runtime/library version strings plus session metadata."""
    return {
        "python": sys.version,
        "pandas": _safe_version("pandas"),
        "numpy": _safe_version("numpy"),
        "scipy": _safe_version("scipy"),
        "plotly": _safe_version("plotly"),
        "streamlit": _safe_version("streamlit"),
        "requests": _safe_version("requests"),
        "scikit-learn": _safe_version("sklearn"),
        "openpyxl": _safe_version("openpyxl"),
        "matplotlib": _safe_version("matplotlib"),
        "kaleido": _safe_version("kaleido"),
        "python-docx": _safe_version("docx"),
        "nbformat": _safe_version("nbformat"),
        "adjustText": _safe_version("adjustText"),
        "session_id": _SESSION_ID,
        "session_start": _SESSION_START,
    }


def get_api_versions() -> dict:
    """Return a dict of external-API version strings."""
    ot_version = "v4 (version check failed)"
    try:
        import requests

        resp = requests.post(
            "https://api.platform.opentargets.org/api/v4/graphql",
            json={"query": "{ meta { apiVersion { x y z } } }"},
            timeout=5,
        )
        if resp.ok:
            data = resp.json()
            v = data.get("data", {}).get("meta", {}).get("apiVersion", {})
            if v:
                ot_version = f"v{v.get('x', '?')}.{v.get('y', '?')}.{v.get('z', '?')}"
    except Exception:
        pass

    return {
        "opentargets": ot_version,
        "dgidb": "v2",
        "enrichr": "current",
        "msigdb": "current",
    }


def format_versions_text() -> str:
    """Return a human-readable multi-line summary of session and library versions."""
    sv = get_session_versions()
    av = get_api_versions()

    lines = [
        "Session Info",
        f"  Session ID : {sv['session_id']}",
        f"  Started    : {sv['session_start']}",
        "",
        "Python & Libraries",
        f"  Python       : {sv['python'].split()[0]}",
    ]

    lib_keys = [
        "pandas", "numpy", "scipy", "plotly", "streamlit",
        "requests", "scikit-learn", "openpyxl", "matplotlib",
        "kaleido", "python-docx", "nbformat", "adjustText",
    ]
    max_label = max(len(k) for k in lib_keys)
    for key in lib_keys:
        lines.append(f"  {key:<{max_label}} : {sv[key]}")

    lines.append("")
    lines.append("External APIs")
    for key, val in av.items():
        lines.append(f"  {key:<14} : {val}")

    return "\n".join(lines)
