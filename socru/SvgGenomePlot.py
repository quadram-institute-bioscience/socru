"""
SVG circular genome diagram generator for Socru.

This module produces publication-quality SVG circular genome diagrams showing
inter-operon fragments as colored arcs, operon positions as directional markers,
origin/terminus labels, and a quality badge. No external SVG libraries are needed;
the SVG is built as a plain Python string.

Functions:
    generate_genome_svg: Build an SVG string for a circular genome diagram.
    save_genome_svg: Generate and write SVG to a file.
"""

import math
from typing import Any, Dict, List

# ColorBrewer Paired palette (colorblind-friendly, 10 colors)
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

_QUALITY_COLORS = {
    "GREEN": "#2ca02c",
    "AMBER": "#ff9900",
    "RED":   "#d62728",
}


def _polar_to_cart(cx: float, cy: float, r: float, angle_deg: float) -> tuple:
    """Convert polar coordinates to cartesian.

    Args:
        cx: Centre x.
        cy: Centre y.
        r: Radius.
        angle_deg: Angle in degrees (0 = 12 o'clock, clockwise).

    Returns:
        Tuple of (x, y).
    """
    angle_rad = math.radians(angle_deg - 90)
    return (cx + r * math.cos(angle_rad), cy + r * math.sin(angle_rad))


def _arc_path(cx: float, cy: float, r: float, start_deg: float, end_deg: float) -> str:
    """Build an SVG arc path description (stroke only, no fill).

    Args:
        cx: Centre x.
        cy: Centre y.
        r: Radius.
        start_deg: Start angle in degrees.
        end_deg: End angle in degrees.

    Returns:
        SVG path 'd' attribute string.
    """
    x1, y1 = _polar_to_cart(cx, cy, r, start_deg)
    x2, y2 = _polar_to_cart(cx, cy, r, end_deg)
    sweep = end_deg - start_deg
    large_arc = 1 if sweep > 180 else 0
    return f"M {x1:.2f} {y1:.2f} A {r:.2f} {r:.2f} 0 {large_arc} 1 {x2:.2f} {y2:.2f}"


def _annular_arc_path(
    cx: float, cy: float, r_inner: float, r_outer: float,
    start_deg: float, end_deg: float,
) -> str:
    """Build a filled annular (ring-sector) arc path.

    The path traces:
        outer arc (start -> end), line to inner arc, inner arc (end -> start), close.

    Args:
        cx: Centre x.
        cy: Centre y.
        r_inner: Inner radius.
        r_outer: Outer radius.
        start_deg: Start angle in degrees.
        end_deg: End angle in degrees.

    Returns:
        SVG path 'd' attribute string.
    """
    # Outer arc points
    ox1, oy1 = _polar_to_cart(cx, cy, r_outer, start_deg)
    ox2, oy2 = _polar_to_cart(cx, cy, r_outer, end_deg)
    # Inner arc points
    ix1, iy1 = _polar_to_cart(cx, cy, r_inner, start_deg)
    ix2, iy2 = _polar_to_cart(cx, cy, r_inner, end_deg)

    sweep = end_deg - start_deg
    large = 1 if sweep > 180 else 0

    return (
        f"M {ox1:.2f} {oy1:.2f} "
        f"A {r_outer:.2f} {r_outer:.2f} 0 {large} 1 {ox2:.2f} {oy2:.2f} "
        f"L {ix2:.2f} {iy2:.2f} "
        f"A {r_inner:.2f} {r_inner:.2f} 0 {large} 0 {ix1:.2f} {iy1:.2f} "
        "Z"
    )


def _escape(text: str) -> str:
    """Escape XML special characters."""
    return text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;").replace('"', "&quot;")


def generate_genome_svg(
    fragments: List[Dict[str, Any]],
    operons: List[Dict[str, Any]],
    genome_length: int,
    gs_type: str,
    quality: str,
    genome_name: str = "",
    width: int = 800,
    height: int = 800,
) -> str:
    """Generate a circular genome SVG diagram.

    Args:
        fragments: List of fragment dicts with keys: number, reversed, length,
            coords (list of (start, end) tuples), is_dnaA, is_dif.
        operons: List of operon dicts with keys: start, end, direction.
        genome_length: Total chromosome length in bp.
        gs_type: Genome structure type label, e.g. "GS1.0".
        quality: Quality category: "GREEN", "AMBER", or "RED".
        genome_name: Filename or label for the genome.
        width: SVG width in pixels.
        height: SVG height in pixels.

    Returns:
        Complete SVG document as a string.
    """
    cx = width / 2
    cy = height / 2

    # Radii for concentric rings
    r_outer_tick = 340  # tick marks outer
    r_outer = 320       # backbone circle
    r_arc_outer = 300   # fragment arcs outer
    r_arc_inner = 240   # fragment arcs inner
    r_operon = 220      # operon marker ring
    r_marker_outer = 355  # Ori/Ter marker lines extend to here

    parts: List[str] = []

    # --- SVG header -----------------------------------------------------------
    parts.append(
        f'<svg xmlns="http://www.w3.org/2000/svg" '
        f'viewBox="0 0 {width} {height}" width="{width}" height="{height}">'
    )
    parts.append('<style>')
    parts.append('  text { font-family: "Helvetica Neue", Helvetica, Arial, sans-serif; }')
    parts.append('</style>')

    # --- Defs (hatch pattern) -------------------------------------------------
    parts.append("<defs>")
    parts.append(
        '<pattern id="hatch" patternUnits="userSpaceOnUse" width="6" height="6" '
        'patternTransform="rotate(45)">'
    )
    parts.append('  <line x1="0" y1="0" x2="0" y2="6" stroke="#000" stroke-width="1.5" opacity="0.35"/>')
    parts.append("</pattern>")
    parts.append("</defs>")

    # --- Background -----------------------------------------------------------
    parts.append(f'<rect width="{width}" height="{height}" fill="white"/>')

    # --- Outer backbone circle ------------------------------------------------
    parts.append(
        f'<circle cx="{cx}" cy="{cy}" r="{r_outer}" '
        'fill="none" stroke="#cccccc" stroke-width="1.5"/>'
    )

    # --- Tick marks every 500 kb ----------------------------------------------
    if genome_length > 0:
        tick_interval = 500_000
        tick_count = int(genome_length / tick_interval)
        for i in range(tick_count):
            bp = i * tick_interval
            angle = (bp / genome_length) * 360.0
            x1, y1 = _polar_to_cart(cx, cy, r_outer, angle)
            x2, y2 = _polar_to_cart(cx, cy, r_outer_tick, angle)
            parts.append(
                f'<line x1="{x1:.2f}" y1="{y1:.2f}" x2="{x2:.2f}" y2="{y2:.2f}" '
                'stroke="#999" stroke-width="1"/>'
            )
            # Label
            lx, ly = _polar_to_cart(cx, cy, r_outer_tick + 14, angle)
            label = f"{bp // 1000}kb" if bp > 0 else "0"
            parts.append(
                f'<text x="{lx:.2f}" y="{ly:.2f}" font-size="9" fill="#666" '
                f'text-anchor="middle" dominant-baseline="central">{label}</text>'
            )

    # --- Fragment arcs (middle ring) ------------------------------------------
    if genome_length > 0 and fragments:
        # Determine cumulative angles for each fragment
        total_len = sum(f.get("length", 0) for f in fragments)
        if total_len == 0:
            total_len = genome_length

        gap_deg = 1.0  # gap between arcs in degrees
        total_gap = gap_deg * len(fragments)
        available = 360.0 - total_gap

        angle_cursor = 0.0
        for idx, frag in enumerate(fragments):
            frag_len = frag.get("length", 0)
            sweep = (frag_len / total_len) * available if total_len > 0 else 0
            start_angle = angle_cursor
            end_angle = angle_cursor + sweep
            color = _PALETTE[idx % len(_PALETTE)]

            # Filled annular arc
            d = _annular_arc_path(cx, cy, r_arc_inner, r_arc_outer, start_angle, end_angle)
            parts.append(f'<path d="{d}" fill="{color}" stroke="white" stroke-width="1"/>')

            # Hatch overlay for reversed fragments
            if frag.get("reversed", False):
                parts.append(f'<path d="{d}" fill="url(#hatch)" stroke="none"/>')

            # Fragment number label at midpoint
            mid_angle = (start_angle + end_angle) / 2
            lr = (r_arc_inner + r_arc_outer) / 2
            lx, ly = _polar_to_cart(cx, cy, lr, mid_angle)
            parts.append(
                f'<text x="{lx:.2f}" y="{ly:.2f}" font-size="13" font-weight="bold" '
                f'fill="white" text-anchor="middle" dominant-baseline="central" '
                f'stroke="#333" stroke-width="0.3">{frag.get("number", "")}</text>'
            )

            # Ori / Ter marker lines
            if frag.get("is_dnaA", False):
                mx1, my1 = _polar_to_cart(cx, cy, r_arc_outer + 2, mid_angle)
                mx2, my2 = _polar_to_cart(cx, cy, r_marker_outer, mid_angle)
                tx, ty = _polar_to_cart(cx, cy, r_marker_outer + 14, mid_angle)
                parts.append(
                    f'<line x1="{mx1:.2f}" y1="{my1:.2f}" x2="{mx2:.2f}" y2="{my2:.2f}" '
                    'stroke="#d62728" stroke-width="2"/>'
                )
                parts.append(
                    f'<text x="{tx:.2f}" y="{ty:.2f}" font-size="12" font-weight="bold" '
                    f'fill="#d62728" text-anchor="middle" dominant-baseline="central">Ori</text>'
                )

            if frag.get("is_dif", False):
                mx1, my1 = _polar_to_cart(cx, cy, r_arc_outer + 2, mid_angle)
                mx2, my2 = _polar_to_cart(cx, cy, r_marker_outer, mid_angle)
                tx, ty = _polar_to_cart(cx, cy, r_marker_outer + 14, mid_angle)
                parts.append(
                    f'<line x1="{mx1:.2f}" y1="{my1:.2f}" x2="{mx2:.2f}" y2="{my2:.2f}" '
                    'stroke="#1f78b4" stroke-width="2"/>'
                )
                parts.append(
                    f'<text x="{tx:.2f}" y="{ty:.2f}" font-size="12" font-weight="bold" '
                    f'fill="#1f78b4" text-anchor="middle" dominant-baseline="central">Ter</text>'
                )

            angle_cursor = end_angle + gap_deg

    # --- Operon markers (inner ring) ------------------------------------------
    if genome_length > 0 and operons:
        tri_size = 8
        for op in operons:
            mid_bp = (op.get("start", 0) + op.get("end", 0)) / 2
            angle = (mid_bp / genome_length) * 360.0
            px, py = _polar_to_cart(cx, cy, r_operon, angle)

            # Triangle pointing in transcription direction
            if op.get("direction", "forward") == "forward":
                # Clockwise-pointing triangle
                p1 = _polar_to_cart(cx, cy, r_operon + tri_size, angle)
                p2 = _polar_to_cart(cx, cy, r_operon - tri_size / 2, angle - 2.5)
                p3 = _polar_to_cart(cx, cy, r_operon - tri_size / 2, angle + 2.5)
            else:
                # Counter-clockwise-pointing triangle
                p1 = _polar_to_cart(cx, cy, r_operon - tri_size, angle)
                p2 = _polar_to_cart(cx, cy, r_operon + tri_size / 2, angle - 2.5)
                p3 = _polar_to_cart(cx, cy, r_operon + tri_size / 2, angle + 2.5)

            points = f"{p1[0]:.2f},{p1[1]:.2f} {p2[0]:.2f},{p2[1]:.2f} {p3[0]:.2f},{p3[1]:.2f}"
            parts.append(f'<polygon points="{points}" fill="#555" stroke="white" stroke-width="0.5"/>')

    # --- Centre: GS label, quality badge, genome name -------------------------
    # Quality badge
    q_color = _QUALITY_COLORS.get(quality.upper(), "#999")
    parts.append(f'<circle cx="{cx}" cy="{cy - 30}" r="10" fill="{q_color}"/>')

    # GS type label
    parts.append(
        f'<text x="{cx}" y="{cy}" font-size="22" font-weight="bold" '
        f'fill="#333" text-anchor="middle" dominant-baseline="central">'
        f'{_escape(gs_type)}</text>'
    )

    # Genome name (truncated if long)
    if genome_name:
        display_name = genome_name if len(genome_name) <= 30 else genome_name[:27] + "..."
        parts.append(
            f'<text x="{cx}" y="{cy + 28}" font-size="11" '
            f'fill="#666" text-anchor="middle" dominant-baseline="central">'
            f'{_escape(display_name)}</text>'
        )

    # Genome length label
    if genome_length > 0:
        size_label = f"{genome_length / 1_000_000:.2f} Mb"
        parts.append(
            f'<text x="{cx}" y="{cy + 48}" font-size="10" '
            f'fill="#999" text-anchor="middle" dominant-baseline="central">'
            f'{size_label}</text>'
        )

    # --- Legend ----------------------------------------------------------------
    legend_x = 20
    legend_y = height - 20 - max(len(fragments) * 20, 20) - 30
    parts.append(
        f'<text x="{legend_x}" y="{legend_y}" font-size="12" font-weight="bold" '
        f'fill="#333">Fragments</text>'
    )
    for idx, frag in enumerate(fragments):
        iy = legend_y + 18 + idx * 20
        color = _PALETTE[idx % len(_PALETTE)]
        parts.append(f'<rect x="{legend_x}" y="{iy - 8}" width="14" height="14" fill="{color}" rx="2"/>')
        rev_label = " (reversed)" if frag.get("reversed", False) else ""
        parts.append(
            f'<text x="{legend_x + 20}" y="{iy + 3}" font-size="11" fill="#333">'
            f'Fragment {frag.get("number", "?")}{rev_label}</text>'
        )

    # Hatch explanation
    hatch_y = legend_y + 18 + len(fragments) * 20 + 6
    parts.append(
        f'<rect x="{legend_x}" y="{hatch_y - 8}" width="14" height="14" '
        f'fill="#ccc" rx="2"/>'
    )
    parts.append(
        f'<rect x="{legend_x}" y="{hatch_y - 8}" width="14" height="14" '
        f'fill="url(#hatch)" rx="2"/>'
    )
    parts.append(
        f'<text x="{legend_x + 20}" y="{hatch_y + 3}" font-size="11" fill="#333">'
        f'= reversed orientation</text>'
    )

    # --- Close SVG ------------------------------------------------------------
    parts.append("</svg>")
    return "\n".join(parts)


def save_genome_svg(filename: str, **kwargs: Any) -> None:
    """Generate and save SVG to file.

    Args:
        filename: Output file path.
        **kwargs: All keyword arguments are forwarded to generate_genome_svg.
    """
    svg = generate_genome_svg(**kwargs)
    with open(filename, "w") as f:
        f.write(svg)
