"""SVG horizontal bar chart for per-fragment BLAST identity.

This module generates a horizontal bar chart showing the BLAST percent identity
for each fragment in a genome assembly. Bars are color-coded by identity level
to provide an immediate visual quality assessment.

Functions:
    generate_fragment_quality_svg: Build an SVG string for a fragment quality chart.
"""

from typing import Any, Dict, List, Optional

# Color thresholds for BLAST identity
_COLOR_EXCELLENT = "#2ca02c"   # green: >99%
_COLOR_GOOD = "#f0c800"        # yellow: 95-99%
_COLOR_FAIR = "#ff7f00"        # orange: 90-95%
_COLOR_POOR = "#d62728"        # red: <90%
_COLOR_UNKNOWN = "#999999"     # gray: no identity data


def _identity_color(identity: Optional[float]) -> str:
    """Return bar color based on BLAST percent identity.

    Args:
        identity: Percent identity (0-100), or None if unknown.

    Returns:
        Hex color string.
    """
    if identity is None:
        return _COLOR_UNKNOWN
    if identity > 99:
        return _COLOR_EXCELLENT
    if identity >= 95:
        return _COLOR_GOOD
    if identity >= 90:
        return _COLOR_FAIR
    return _COLOR_POOR


def _escape(text: str) -> str:
    """Escape XML special characters."""
    return (
        text.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def generate_fragment_quality_svg(
    fragments: List[Dict[str, Any]],
    genome_name: str = "",
    width: int = 600,
    height: Optional[int] = None,
) -> str:
    """Generate a horizontal bar chart SVG showing per-fragment BLAST identity.

    Args:
        fragments: List of fragment dicts with keys: number, reversed,
            blast_identity (float or None), blast_alignment_length (int or None),
            length (int).
        genome_name: Assembly name displayed as chart title.
        width: SVG width in pixels.
        height: SVG height in pixels. If None, auto-calculated from
            the number of fragments.

    Returns:
        Complete SVG document as a string.
    """
    bar_height = 28
    bar_gap = 8
    margin_top = 60
    margin_bottom = 50
    margin_left = 80
    margin_right = 40

    n_fragments = max(len(fragments), 1)
    if height is None:
        height = margin_top + n_fragments * (bar_height + bar_gap) + margin_bottom

    chart_width = width - margin_left - margin_right

    # X-axis scale: 90-100% identity
    x_min = 90.0
    x_max = 100.0
    x_range = x_max - x_min

    parts: List[str] = []

    # SVG header
    parts.append(
        f'<svg xmlns="http://www.w3.org/2000/svg" '
        f'viewBox="0 0 {width} {height}" width="{width}" height="{height}">'
    )
    parts.append("<style>")
    parts.append(
        '  text { font-family: "Helvetica Neue", Helvetica, Arial, sans-serif; }'
    )
    parts.append("</style>")

    # Background
    parts.append(f'<rect width="{width}" height="{height}" fill="white"/>')

    # Title
    title = "Fragment Quality"
    if genome_name:
        title = f"Fragment Quality - {_escape(genome_name)}"
    parts.append(
        f'<text x="{width / 2}" y="28" font-size="16" font-weight="bold" '
        f'fill="#333" text-anchor="middle">{title}</text>'
    )

    # X-axis gridlines and labels
    chart_bottom = margin_top + n_fragments * (bar_height + bar_gap)
    for pct in range(90, 101):
        x = margin_left + ((pct - x_min) / x_range) * chart_width
        # Gridline
        parts.append(
            f'<line x1="{x:.1f}" y1="{margin_top - 5}" '
            f'x2="{x:.1f}" y2="{chart_bottom}" '
            f'stroke="#e0e0e0" stroke-width="1"/>'
        )
        # Label
        parts.append(
            f'<text x="{x:.1f}" y="{chart_bottom + 18}" font-size="10" '
            f'fill="#666" text-anchor="middle">{pct}%</text>'
        )

    # X-axis title
    parts.append(
        f'<text x="{margin_left + chart_width / 2}" y="{chart_bottom + 38}" '
        f'font-size="11" fill="#666" text-anchor="middle">Percent Identity</text>'
    )

    # Bars
    for idx, frag in enumerate(fragments):
        y = margin_top + idx * (bar_height + bar_gap)
        identity = frag.get("blast_identity")
        color = _identity_color(identity)
        number = frag.get("number", "?")
        is_reversed = frag.get("reversed", False)

        # Y-axis label: fragment number with prime if reversed
        label = f"{number}'" if is_reversed else str(number)
        parts.append(
            f'<text x="{margin_left - 10}" y="{y + bar_height / 2 + 1}" '
            f'font-size="12" fill="#333" text-anchor="end" '
            f'dominant-baseline="central">{label}</text>'
        )

        if identity is not None:
            # Clamp identity to display range
            display_identity = max(identity, x_min)
            bar_width = ((display_identity - x_min) / x_range) * chart_width
            bar_width = max(bar_width, 2)  # minimum visible width

            parts.append(
                f'<rect x="{margin_left}" y="{y}" '
                f'width="{bar_width:.1f}" height="{bar_height}" '
                f'fill="{color}" rx="2"/>'
            )

            # Identity value label
            label_x = margin_left + bar_width + 5
            parts.append(
                f'<text x="{label_x:.1f}" y="{y + bar_height / 2 + 1}" '
                f'font-size="10" fill="#333" '
                f'dominant-baseline="central">{identity:.1f}%</text>'
            )
        else:
            # Unknown: gray dashed bar at half width
            dash_width = chart_width * 0.5
            parts.append(
                f'<rect x="{margin_left}" y="{y}" '
                f'width="{dash_width:.1f}" height="{bar_height}" '
                f'fill="none" stroke="{_COLOR_UNKNOWN}" '
                f'stroke-width="2" stroke-dasharray="6,4" rx="2"/>'
            )
            parts.append(
                f'<text x="{margin_left + dash_width + 5}" '
                f'y="{y + bar_height / 2 + 1}" font-size="10" fill="#999" '
                f'dominant-baseline="central">N/A</text>'
            )

    # Legend
    legend_items = [
        (_COLOR_EXCELLENT, ">99%"),
        (_COLOR_GOOD, "95-99%"),
        (_COLOR_FAIR, "90-95%"),
        (_COLOR_POOR, "<90%"),
        (_COLOR_UNKNOWN, "Unknown"),
    ]
    legend_x = margin_left
    legend_y = margin_top - 25
    for i, (color, label) in enumerate(legend_items):
        lx = legend_x + i * 90
        parts.append(
            f'<rect x="{lx}" y="{legend_y - 6}" width="12" height="12" '
            f'fill="{color}" rx="2"/>'
        )
        parts.append(
            f'<text x="{lx + 16}" y="{legend_y + 1}" font-size="9" '
            f'fill="#666" dominant-baseline="central">{label}</text>'
        )

    # Close SVG
    parts.append("</svg>")
    return "\n".join(parts)
