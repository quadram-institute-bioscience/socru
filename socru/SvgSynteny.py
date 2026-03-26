"""
Linear synteny / fragment-order comparison SVG generator for Socru.

This module produces a multi-genome comparison diagram where each assembly is
rendered as a horizontal bar divided into colored fragment blocks.  Curved
ribbons connect matching fragments between adjacent assemblies so that
rearrangements (translocations, inversions) are immediately visible.

No external SVG libraries are required; the SVG is built as a plain Python
string, following the same approach used in SvgGenomePlot.py.

Functions:
    generate_synteny_svg: Build an SVG string for a linear synteny diagram.
"""

from typing import Any, Dict, List, Optional

# ColorBrewer Paired palette (same as SvgGenomePlot)
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

_QUALITY_COLORS: Dict[str, str] = {
    "GREEN": "#2ca02c",
    "AMBER": "#ff9900",
    "RED": "#d62728",
}


def _escape(text: str) -> str:
    """Escape XML special characters.

    Args:
        text: Raw string.

    Returns:
        XML-safe string.
    """
    return (
        text.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def _fragment_color(fragment_number: int) -> str:
    """Return palette color for a 1-based fragment number.

    Args:
        fragment_number: 1-based fragment identifier.

    Returns:
        Hex color string.
    """
    idx = (fragment_number - 1) % len(_PALETTE)
    return _PALETTE[idx]


def _ribbon_path(
    x1_left: float,
    x1_right: float,
    y1: float,
    x2_left: float,
    x2_right: float,
    y2: float,
) -> str:
    """Build a cubic-bezier ribbon path connecting two horizontal spans.

    The ribbon goes from the bottom of the upper bar to the top of the lower
    bar, with control points placed at the vertical midpoint to produce a
    smooth S-curve.

    Args:
        x1_left: Left edge x of source span (upper bar bottom).
        x1_right: Right edge x of source span.
        y1: y coordinate of upper bar bottom edge.
        x2_left: Left edge x of destination span (lower bar top).
        x2_right: Right edge x of destination span.
        y2: y coordinate of lower bar top edge.

    Returns:
        SVG path ``d`` attribute string.
    """
    cy = (y1 + y2) / 2.0
    return (
        f"M {x1_left:.2f} {y1:.2f} "
        f"C {x1_left:.2f} {cy:.2f}, {x2_left:.2f} {cy:.2f}, {x2_left:.2f} {y2:.2f} "
        f"L {x2_right:.2f} {y2:.2f} "
        f"C {x2_right:.2f} {cy:.2f}, {x1_right:.2f} {cy:.2f}, {x1_right:.2f} {y1:.2f} "
        "Z"
    )


def generate_synteny_svg(
    assemblies: List[Dict[str, Any]],
    width: int = 1000,
    height: Optional[int] = None,
    bar_height: int = 30,
    gap: int = 60,
) -> str:
    """Generate a linear multi-genome synteny comparison SVG.

    Each assembly is drawn as a horizontal bar whose coloured blocks represent
    inter-operon fragments (proportional to length).  Curved ribbons between
    adjacent bars link matching fragments, making rearrangements immediately
    visible as crossed ribbons.

    Args:
        assemblies: List of assembly dicts, each with keys:
            - name (str): Assembly label.
            - gs_type (str): Genome structure type, e.g. "GS1.0".
            - quality (str): "GREEN", "AMBER", or "RED".
            - fragments: List of fragment dicts with keys:
                - number (int): Fragment identifier (1-based).
                - reversed (bool): True if fragment is reverse-complemented.
                - length (int): Fragment length in bp.
        width: SVG width in pixels (default 1000).
        height: SVG height in pixels.  If ``None``, auto-calculated from the
            number of assemblies.
        bar_height: Height of each assembly bar in pixels (default 30).
        gap: Vertical gap between bars in pixels (default 60).

    Returns:
        Complete SVG document as a string.
    """
    if not assemblies:
        assemblies = []

    # Sort assemblies by GS type for visual clustering
    sorted_assemblies = sorted(assemblies, key=lambda a: a.get("gs_type", ""))

    n = len(sorted_assemblies)
    label_margin = 180  # left margin for labels
    right_margin = 20
    bar_width = width - label_margin - right_margin
    top_margin = 40

    if height is None:
        height = top_margin + max(n * (bar_height + gap) - gap + top_margin, 100)

    parts: List[str] = []

    # --- SVG header -----------------------------------------------------------
    parts.append(
        f'<svg xmlns="http://www.w3.org/2000/svg" '
        f'viewBox="0 0 {width} {height}" width="{width}" height="{height}">'
    )
    parts.append("<style>")
    parts.append(
        '  text { font-family: "Helvetica Neue", Helvetica, Arial, sans-serif; }'
    )
    parts.append("</style>")

    # --- Defs (hatch pattern for reversed fragments) --------------------------
    parts.append("<defs>")
    parts.append(
        '<pattern id="hatch" patternUnits="userSpaceOnUse" width="6" height="6" '
        'patternTransform="rotate(45)">'
    )
    parts.append(
        '  <line x1="0" y1="0" x2="0" y2="6" stroke="#000" '
        'stroke-width="1.5" opacity="0.35"/>'
    )
    parts.append("</pattern>")
    parts.append("</defs>")

    # --- Background -----------------------------------------------------------
    parts.append(f'<rect width="{width}" height="{height}" fill="white"/>')

    # Pre-compute block positions for every assembly so ribbons can reference
    # both the current and the next assembly.
    # block_positions[i] maps fragment_number -> (x_left, x_right)
    block_positions: List[Dict[int, tuple]] = []
    bar_tops: List[float] = []

    for i, asm in enumerate(sorted_assemblies):
        y = top_margin + i * (bar_height + gap)
        bar_tops.append(y)

        frags = asm.get("fragments", [])
        total_len = sum(f.get("length", 1) for f in frags)
        if total_len == 0:
            total_len = 1

        positions: Dict[int, tuple] = {}
        x_cursor = label_margin
        for frag in frags:
            frag_len = frag.get("length", 1)
            block_w = (frag_len / total_len) * bar_width
            fnum = frag.get("number", 0)
            positions[fnum] = (x_cursor, x_cursor + block_w)
            x_cursor += block_w
        block_positions.append(positions)

    # --- Draw ribbons between adjacent assemblies -----------------------------
    for i in range(n - 1):
        pos_upper = block_positions[i]
        pos_lower = block_positions[i + 1]
        y_upper_bottom = bar_tops[i] + bar_height
        y_lower_top = bar_tops[i + 1]

        # Find common fragment numbers
        common = set(pos_upper.keys()) & set(pos_lower.keys())
        for fnum in sorted(common):
            x1_l, x1_r = pos_upper[fnum]
            x2_l, x2_r = pos_lower[fnum]
            color = _fragment_color(fnum)
            d = _ribbon_path(x1_l, x1_r, y_upper_bottom, x2_l, x2_r, y_lower_top)
            parts.append(
                f'<path d="{d}" fill="{color}" fill-opacity="0.3" '
                f'stroke="{color}" stroke-width="0.5"/>'
            )

    # --- Draw assembly bars and labels ----------------------------------------
    for i, asm in enumerate(sorted_assemblies):
        y = bar_tops[i]
        frags = asm.get("fragments", [])
        name = asm.get("name", "")
        gs_type = asm.get("gs_type", "")
        quality = asm.get("quality", "GREEN").upper()

        # Quality badge (small colored circle)
        q_color = _QUALITY_COLORS.get(quality, "#999")
        parts.append(
            f'<circle cx="14" cy="{y + bar_height / 2:.2f}" r="5" '
            f'fill="{q_color}"/>'
        )

        # Assembly name label
        display_name = name if len(name) <= 18 else name[:15] + "..."
        parts.append(
            f'<text x="24" y="{y + bar_height / 2:.2f}" font-size="11" '
            f'fill="#333" dominant-baseline="central">'
            f"{_escape(display_name)}</text>"
        )

        # GS type label
        parts.append(
            f'<text x="{label_margin - 10}" y="{y + bar_height / 2:.2f}" '
            f'font-size="11" font-weight="bold" fill="#555" '
            f'text-anchor="end" dominant-baseline="central">'
            f"{_escape(gs_type)}</text>"
        )

        # Fragment blocks
        total_len = sum(f.get("length", 1) for f in frags)
        if total_len == 0:
            total_len = 1

        x_cursor = label_margin
        for frag in frags:
            frag_len = frag.get("length", 1)
            block_w = (frag_len / total_len) * bar_width
            fnum = frag.get("number", 0)
            color = _fragment_color(fnum)
            is_reversed = frag.get("reversed", False)

            # Filled rectangle
            parts.append(
                f'<rect x="{x_cursor:.2f}" y="{y}" '
                f'width="{block_w:.2f}" height="{bar_height}" '
                f'fill="{color}" stroke="white" stroke-width="1"/>'
            )

            # Hatch overlay for reversed fragments
            if is_reversed:
                parts.append(
                    f'<rect x="{x_cursor:.2f}" y="{y}" '
                    f'width="{block_w:.2f}" height="{bar_height}" '
                    f'fill="url(#hatch)" stroke="none"/>'
                )

            # Fragment number label (centred in block)
            mid_x = x_cursor + block_w / 2
            parts.append(
                f'<text x="{mid_x:.2f}" y="{y + bar_height / 2:.2f}" '
                f'font-size="12" font-weight="bold" fill="white" '
                f'text-anchor="middle" dominant-baseline="central" '
                f'stroke="#333" stroke-width="0.3">{fnum}</text>'
            )

            # dnaA marker (Ori) -- small red vertical line above block
            if frag.get("is_dnaA", False):
                parts.append(
                    f'<line x1="{mid_x:.2f}" y1="{y - 2:.2f}" '
                    f'x2="{mid_x:.2f}" y2="{y - 12:.2f}" '
                    f'stroke="#d62728" stroke-width="2"/>'
                )
                parts.append(
                    f'<text x="{mid_x:.2f}" y="{y - 14:.2f}" font-size="8" '
                    f'fill="#d62728" text-anchor="middle">Ori</text>'
                )

            # dif marker (Ter) -- small blue vertical line below block
            if frag.get("is_dif", False):
                parts.append(
                    f'<line x1="{mid_x:.2f}" y1="{y + bar_height + 2:.2f}" '
                    f'x2="{mid_x:.2f}" y2="{y + bar_height + 12:.2f}" '
                    f'stroke="#1f78b4" stroke-width="2"/>'
                )
                parts.append(
                    f'<text x="{mid_x:.2f}" y="{y + bar_height + 20:.2f}" '
                    f'font-size="8" fill="#1f78b4" text-anchor="middle">Ter</text>'
                )

            x_cursor += block_w

    # --- Close SVG ------------------------------------------------------------
    parts.append("</svg>")
    return "\n".join(parts)
