"""
SVG coverage pileup visualization for Socru.

This module generates a stacked panel SVG showing BLAST alignment coverage depth
across each fragment. Each fragment is rendered as a small area chart with a
filled polygon, a dashed mean-coverage line, and fragment label. Fragments with
no coverage data display a gray placeholder.

Functions:
    generate_coverage_pileup_svg: Build an SVG string for coverage pileup panels.
"""

from typing import Any, Dict, List, Optional

# ColorBrewer Paired palette (matches SvgGenomePlot)
_PALETTE = [
    "#a6cee3",  # light blue
    "#1f78b4",  # blue
    "#b2df8a",  # light green
    "#33a02c",  # green
    "#fb9a99",  # light red
    "#e31a1c",  # red
    "#fdbf6f",  # light orange
    "#ff7f00",  # orange
    "#cab2d6",  # light purple
    "#6a3d9a",  # purple
]


def _escape(text: str) -> str:
    """Escape XML special characters."""
    return (
        text.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def _downsample(depths: List[int], max_points: int) -> List[float]:
    """Downsample a coverage array to at most max_points by averaging bins.

    Args:
        depths: Raw coverage depths per base position.
        max_points: Maximum number of points in the output.

    Returns:
        List of averaged coverage values.
    """
    n = len(depths)
    if n <= max_points:
        return [float(d) for d in depths]
    bin_size = n / max_points
    result: List[float] = []
    for i in range(max_points):
        start = int(i * bin_size)
        end = int((i + 1) * bin_size)
        end = min(end, n)
        if start >= end:
            result.append(0.0)
        else:
            result.append(sum(depths[start:end]) / (end - start))
    return result


def generate_coverage_pileup_svg(
    fragments: List[Dict[str, Any]],
    genome_name: str = "",
    width: int = 800,
    height: Optional[int] = None,
) -> str:
    """Generate stacked coverage-pileup SVG panels for each fragment.

    Args:
        fragments: List of fragment dicts with keys:
            - number: Fragment identifier (int or str).
            - length: Fragment length in bp.
            - coverage_depths: List of ints (per-base coverage) or None.
        genome_name: Optional genome label for the title.
        width: SVG width in pixels.
        height: SVG height in pixels; auto-calculated if None (~120px per fragment).

    Returns:
        Complete SVG document as a string.
    """
    if not fragments:
        h = 60
        parts = [
            f'<svg xmlns="http://www.w3.org/2000/svg" '
            f'viewBox="0 0 {width} {h}" width="{width}" height="{h}">',
            '<style>text { font-family: "Helvetica Neue", Helvetica, Arial, sans-serif; }</style>',
            f'<rect width="{width}" height="{h}" fill="white"/>',
            f'<text x="{width // 2}" y="30" font-size="13" fill="#999" '
            f'text-anchor="middle" dominant-baseline="central">No fragments</text>',
            "</svg>",
        ]
        return "\n".join(parts)

    # Layout constants
    panel_height = 120
    title_height = 40
    margin_left = 80
    margin_right = 20
    margin_top = 10
    margin_bottom = 20
    plot_width = width - margin_left - margin_right
    plot_height = panel_height - margin_top - margin_bottom
    max_points = min(plot_width, 600)

    total_height = height if height is not None else title_height + panel_height * len(fragments)

    parts: List[str] = []

    # SVG header
    parts.append(
        f'<svg xmlns="http://www.w3.org/2000/svg" '
        f'viewBox="0 0 {width} {total_height}" width="{width}" height="{total_height}">'
    )
    parts.append(
        '<style>text { font-family: "Helvetica Neue", Helvetica, Arial, sans-serif; }</style>'
    )
    parts.append(f'<rect width="{width}" height="{total_height}" fill="white"/>')

    # Title
    title_text = "Coverage Pileup"
    if genome_name:
        title_text += f" - {_escape(genome_name)}"
    parts.append(
        f'<text x="{width // 2}" y="{title_height // 2 + 4}" font-size="15" '
        f'font-weight="bold" fill="#333" text-anchor="middle" '
        f'dominant-baseline="central">{title_text}</text>'
    )

    # Draw each fragment panel
    for idx, frag in enumerate(fragments):
        frag_num = frag.get("number", idx + 1)
        frag_len = frag.get("length", 0)
        depths = frag.get("coverage_depths", None)
        color = _PALETTE[idx % len(_PALETTE)]

        panel_y = title_height + idx * panel_height

        # Panel background with subtle border
        parts.append(
            f'<rect x="0" y="{panel_y}" width="{width}" height="{panel_height}" '
            f'fill="white" stroke="#eee" stroke-width="1"/>'
        )

        # Y-axis label
        parts.append(
            f'<text x="{margin_left - 12}" y="{panel_y + panel_height // 2}" '
            f'font-size="12" font-weight="bold" fill="{color}" '
            f'text-anchor="end" dominant-baseline="central">Frag {frag_num}</text>'
        )

        plot_x = margin_left
        plot_y = panel_y + margin_top

        if depths is None or len(depths) == 0:
            # No-data placeholder
            parts.append(
                f'<rect x="{plot_x}" y="{plot_y}" width="{plot_width}" '
                f'height="{plot_height}" fill="#f5f5f5" rx="4"/>'
            )
            parts.append(
                f'<text x="{plot_x + plot_width // 2}" y="{plot_y + plot_height // 2}" '
                f'font-size="12" fill="#aaa" text-anchor="middle" '
                f'dominant-baseline="central">No data</text>'
            )
            continue

        # Downsample for rendering
        sampled = _downsample(depths, max_points)
        n_points = len(sampled)
        max_depth = max(sampled) if sampled else 1.0
        if max_depth == 0:
            max_depth = 1.0
        mean_depth = sum(sampled) / n_points if n_points > 0 else 0.0

        # Build area-chart polygon points
        # Start at bottom-left, trace top, end at bottom-right
        polygon_points: List[str] = []
        polygon_points.append(f"{plot_x:.1f},{plot_y + plot_height:.1f}")
        for i, val in enumerate(sampled):
            px = plot_x + (i / max(n_points - 1, 1)) * plot_width
            py = plot_y + plot_height - (val / max_depth) * plot_height
            polygon_points.append(f"{px:.1f},{py:.1f}")
        polygon_points.append(f"{plot_x + plot_width:.1f},{plot_y + plot_height:.1f}")

        parts.append(
            f'<polygon points="{" ".join(polygon_points)}" '
            f'fill="{color}" fill-opacity="0.4" stroke="{color}" stroke-width="1"/>'
        )

        # Mean coverage dashed line
        mean_y = plot_y + plot_height - (mean_depth / max_depth) * plot_height
        parts.append(
            f'<line x1="{plot_x}" y1="{mean_y:.1f}" '
            f'x2="{plot_x + plot_width}" y2="{mean_y:.1f}" '
            f'stroke="{color}" stroke-width="1" stroke-dasharray="5,3" opacity="0.8"/>'
        )

        # Mean label
        parts.append(
            f'<text x="{plot_x + plot_width + 4}" y="{mean_y:.1f}" '
            f'font-size="9" fill="#666" dominant-baseline="central">'
            f'{mean_depth:.1f}x</text>'
        )

        # X-axis labels (0 and length)
        parts.append(
            f'<text x="{plot_x}" y="{plot_y + plot_height + 12}" '
            f'font-size="8" fill="#999" text-anchor="start">0</text>'
        )
        if frag_len > 0:
            len_label = f"{frag_len // 1000}kb" if frag_len >= 1000 else str(frag_len)
            parts.append(
                f'<text x="{plot_x + plot_width}" y="{plot_y + plot_height + 12}" '
                f'font-size="8" fill="#999" text-anchor="end">{len_label}</text>'
            )

        # Max depth label
        parts.append(
            f'<text x="{plot_x - 4}" y="{plot_y + 4}" '
            f'font-size="8" fill="#999" text-anchor="end" '
            f'dominant-baseline="hanging">{max_depth:.0f}x</text>'
        )

    # Close SVG
    parts.append("</svg>")
    return "\n".join(parts)
