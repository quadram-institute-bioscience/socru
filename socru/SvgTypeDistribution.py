"""SVG bar chart for GS type frequency distribution.

This module generates a bar chart showing how frequently each genome structure
(GS) type occurs in a batch of Socru analyses. Bars can optionally be stacked
by quality level (GREEN/AMBER/RED) to show quality breakdown per type.

Functions:
    generate_type_distribution_svg: Build an SVG string for a type distribution chart.
"""

from typing import Any, Dict, List, Optional


_QUALITY_COLORS = {
    "GREEN": "#2ca02c",
    "AMBER": "#ff9900",
    "RED": "#d62728",
}

_QUALITY_ORDER = ["GREEN", "AMBER", "RED"]


def _escape(text: str) -> str:
    """Escape XML special characters."""
    return (
        text.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def generate_type_distribution_svg(
    type_counts: Dict[str, int],
    quality_counts: Optional[Dict[str, Dict[str, int]]] = None,
    width: int = 800,
    height: int = 400,
) -> str:
    """Generate a bar chart SVG of GS type frequency distribution.

    Args:
        type_counts: Dictionary mapping GS type strings to their counts,
            e.g. {"GS1.0": 15, "GS2.1": 8}.
        quality_counts: Optional dictionary mapping GS type strings to
            quality breakdowns, e.g. {"GS1.0": {"GREEN": 10, "AMBER": 3, "RED": 2}}.
            When provided, bars are stacked by quality level.
        width: SVG width in pixels.
        height: SVG height in pixels.

    Returns:
        Complete SVG document as a string.
    """
    if not type_counts:
        type_counts = {}

    margin_top = 60
    margin_bottom = 80
    margin_left = 60
    margin_right = 30

    chart_width = width - margin_left - margin_right
    chart_height = height - margin_top - margin_bottom

    # Sort types by frequency descending
    sorted_types = sorted(type_counts.items(), key=lambda x: x[1], reverse=True)
    n_bars = max(len(sorted_types), 1)
    total_count = sum(c for _, c in sorted_types)
    max_count = max((c for _, c in sorted_types), default=1)

    # Bar geometry
    bar_group_width = chart_width / n_bars
    bar_width = max(bar_group_width * 0.65, 10)
    bar_padding = (bar_group_width - bar_width) / 2

    # Y-axis: determine nice tick interval
    y_tick = _nice_tick(max_count)

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
    parts.append(
        f'<text x="{width / 2}" y="30" font-size="18" font-weight="bold" '
        f'fill="#333" text-anchor="middle">GS Type Distribution</text>'
    )

    chart_bottom = margin_top + chart_height

    # Y-axis gridlines and labels
    if max_count > 0:
        tick = y_tick
        val = 0
        while val <= max_count:
            y = chart_bottom - (val / max_count) * chart_height
            parts.append(
                f'<line x1="{margin_left}" y1="{y:.1f}" '
                f'x2="{margin_left + chart_width}" y2="{y:.1f}" '
                f'stroke="#e0e0e0" stroke-width="1"/>'
            )
            parts.append(
                f'<text x="{margin_left - 8}" y="{y:.1f}" font-size="10" '
                f'fill="#666" text-anchor="end" dominant-baseline="central">'
                f'{int(val)}</text>'
            )
            val += tick

    # Y-axis label
    parts.append(
        f'<text x="15" y="{margin_top + chart_height / 2}" font-size="12" '
        f'fill="#666" text-anchor="middle" dominant-baseline="central" '
        f'transform="rotate(-90, 15, {margin_top + chart_height / 2})">Count</text>'
    )

    # Baseline
    parts.append(
        f'<line x1="{margin_left}" y1="{chart_bottom}" '
        f'x2="{margin_left + chart_width}" y2="{chart_bottom}" '
        f'stroke="#333" stroke-width="1"/>'
    )

    # Bars
    for i, (gs_type, count) in enumerate(sorted_types):
        x = margin_left + i * bar_group_width + bar_padding

        if quality_counts and gs_type in quality_counts:
            # Stacked bar
            q_data = quality_counts[gs_type]
            y_cursor = chart_bottom
            for q_level in _QUALITY_ORDER:
                q_count = q_data.get(q_level, 0)
                if q_count <= 0:
                    continue
                seg_height = (q_count / max_count) * chart_height
                y_cursor -= seg_height
                color = _QUALITY_COLORS[q_level]
                parts.append(
                    f'<rect x="{x:.1f}" y="{y_cursor:.1f}" '
                    f'width="{bar_width:.1f}" height="{seg_height:.1f}" '
                    f'fill="{color}" rx="2"/>'
                )
        else:
            # Single color bar
            bar_h = (count / max_count) * chart_height if max_count > 0 else 0
            y = chart_bottom - bar_h
            parts.append(
                f'<rect x="{x:.1f}" y="{y:.1f}" '
                f'width="{bar_width:.1f}" height="{bar_h:.1f}" '
                f'fill="#4878cf" rx="2"/>'
            )

        # Percentage label above bar
        bar_h = (count / max_count) * chart_height if max_count > 0 else 0
        pct = (count / total_count * 100) if total_count > 0 else 0
        label_y = chart_bottom - bar_h - 6
        parts.append(
            f'<text x="{x + bar_width / 2:.1f}" y="{label_y:.1f}" '
            f'font-size="10" fill="#333" text-anchor="middle">'
            f'{pct:.0f}%</text>'
        )

        # X-axis label
        parts.append(
            f'<text x="{x + bar_width / 2:.1f}" y="{chart_bottom + 18}" '
            f'font-size="11" fill="#333" text-anchor="middle">'
            f'{_escape(gs_type)}</text>'
        )

    # Legend (only when quality_counts provided)
    if quality_counts:
        legend_x = width - margin_right - 200
        legend_y = margin_top - 15
        for i, q_level in enumerate(_QUALITY_ORDER):
            lx = legend_x + i * 70
            parts.append(
                f'<rect x="{lx}" y="{legend_y - 6}" width="12" height="12" '
                f'fill="{_QUALITY_COLORS[q_level]}" rx="2"/>'
            )
            parts.append(
                f'<text x="{lx + 16}" y="{legend_y + 1}" font-size="10" '
                f'fill="#666" dominant-baseline="central">{q_level}</text>'
            )

    # Close SVG
    parts.append("</svg>")
    return "\n".join(parts)


def _nice_tick(max_val: int) -> int:
    """Compute a nice tick interval for the Y-axis.

    Args:
        max_val: Maximum value on the axis.

    Returns:
        Integer tick interval.
    """
    if max_val <= 5:
        return 1
    if max_val <= 10:
        return 2
    if max_val <= 25:
        return 5
    if max_val <= 50:
        return 10
    if max_val <= 100:
        return 20
    # For larger values, use powers of 10
    magnitude = 10 ** (len(str(max_val)) - 1)
    return magnitude
