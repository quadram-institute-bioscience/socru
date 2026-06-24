"""
SVG confidence heatmap for Socru.

This module generates an SVG heatmap showing per-fragment BLAST identity across
many assemblies. Rows represent assemblies (sorted by GS type then name), columns
represent fragment numbers, and cell colors encode identity percentage on a
diverging gradient from white (100%) through yellow to red (<90%).

Functions:
    generate_confidence_heatmap_svg: Build an SVG string for the heatmap.
"""

from typing import Any, Dict, List, Optional

_QUALITY_COLORS = {
    "GREEN": "#2ca02c",
    "AMBER": "#ff9900",
    "RED": "#d62728",
}


def _escape(text: str) -> str:
    """Escape XML special characters."""
    return (
        text.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def _identity_color(identity: Optional[float]) -> str:
    """Map BLAST identity percentage to a fill color.

    Color stops:
        100%    -> white (#ffffff)
        99%     -> light yellow (#ffffcc)
        97%     -> yellow (#ffff00)
        95%     -> orange (#ff9900)
        <90%    -> red (#d62728)
        None    -> gray (#dddddd)

    Args:
        identity: Percentage identity (0-100) or None for missing data.

    Returns:
        Hex color string.
    """
    if identity is None:
        return "#dddddd"

    if identity >= 100.0:
        return "#ffffff"
    elif identity >= 99.0:
        # Interpolate white -> light yellow
        t = (100.0 - identity) / 1.0  # 0 at 100%, 1 at 99%
        r = 255
        g = 255
        b = int(255 - t * (255 - 204))  # 255 -> 204
        return f"#{r:02x}{g:02x}{b:02x}"
    elif identity >= 97.0:
        # Interpolate light yellow -> yellow
        t = (99.0 - identity) / 2.0  # 0 at 99%, 1 at 97%
        r = 255
        g = 255
        b = int(204 - t * 204)  # 204 -> 0
        return f"#{r:02x}{g:02x}{b:02x}"
    elif identity >= 95.0:
        # Interpolate yellow -> orange
        t = (97.0 - identity) / 2.0  # 0 at 97%, 1 at 95%
        r = 255
        g = int(255 - t * (255 - 153))  # 255 -> 153
        b = 0
        return f"#{r:02x}{g:02x}{b:02x}"
    elif identity >= 90.0:
        # Interpolate orange -> red
        t = (95.0 - identity) / 5.0  # 0 at 95%, 1 at 90%
        r = int(255 - t * (255 - 214))  # 255 -> 214
        g = int(153 - t * (153 - 39))  # 153 -> 39
        b = int(0 + t * 40)  # 0 -> 40
        return f"#{r:02x}{g:02x}{b:02x}"
    else:
        return "#d62728"


def generate_confidence_heatmap_svg(
    assemblies: List[Dict[str, Any]],
    width: int = 800,
    height: Optional[int] = None,
    cell_size: int = 30,
) -> str:
    """Generate an SVG heatmap of per-fragment BLAST identity across assemblies.

    Args:
        assemblies: List of assembly dicts with keys:
            - name: Assembly identifier string.
            - gs_type: Genome structure type string (e.g. "GS1.0").
            - quality: Quality category ("GREEN", "AMBER", "RED").
            - fragments: List of dicts with keys:
                - number: Fragment number (int).
                - blast_identity: Percentage identity (float, 0-100) or None.
        width: SVG width in pixels.
        height: SVG height in pixels; auto-calculated if None.
        cell_size: Size of each heatmap cell in pixels.

    Returns:
        Complete SVG document as a string.
    """
    if not assemblies:
        h = 60
        parts = [
            f'<svg xmlns="http://www.w3.org/2000/svg" '
            f'viewBox="0 0 {width} {h}" width="{width}" height="{h}">',
            '<style>text { font-family: "Helvetica Neue", Helvetica, Arial, sans-serif; }</style>',
            f'<rect width="{width}" height="{h}" fill="white"/>',
            f'<text x="{width // 2}" y="30" font-size="13" fill="#999" '
            f'text-anchor="middle" dominant-baseline="central">No assemblies</text>',
            "</svg>",
        ]
        return "\n".join(parts)

    # Sort assemblies by GS type then by name
    sorted_assemblies = sorted(assemblies, key=lambda a: (a.get("gs_type", ""), a.get("name", "")))

    # Determine fragment columns from all assemblies
    all_frag_numbers: List[int] = []
    for asm in sorted_assemblies:
        for f in asm.get("fragments", []):
            num = f.get("number")
            if num is not None and num not in all_frag_numbers:
                all_frag_numbers.append(num)
    all_frag_numbers.sort()

    n_cols = len(all_frag_numbers)
    n_rows = len(sorted_assemblies)

    # Layout
    title_height = 35
    header_height = 50  # column headers
    row_label_width = 180  # left side for assembly name + quality dot
    margin_right = 20
    margin_bottom = 10

    grid_width = n_cols * cell_size
    grid_height = n_rows * cell_size

    total_width = max(width, row_label_width + grid_width + margin_right)
    total_height = height if height is not None else title_height + header_height + grid_height + margin_bottom

    grid_x = row_label_width
    grid_y = title_height + header_height

    parts: List[str] = []

    # SVG header
    parts.append(
        f'<svg xmlns="http://www.w3.org/2000/svg" '
        f'viewBox="0 0 {total_width} {total_height}" width="{total_width}" height="{total_height}">'
    )
    parts.append(
        '<style>text { font-family: "Helvetica Neue", Helvetica, Arial, sans-serif; }</style>'
    )
    parts.append(f'<rect width="{total_width}" height="{total_height}" fill="white"/>')

    # Title
    parts.append(
        f'<text x="{total_width // 2}" y="{title_height // 2 + 4}" font-size="15" '
        f'font-weight="bold" fill="#333" text-anchor="middle" '
        f'dominant-baseline="central">Fragment Identity Heatmap</text>'
    )

    # Compute column averages for header annotation
    col_averages: List[Optional[float]] = []
    for frag_num in all_frag_numbers:
        values: List[float] = []
        for asm in sorted_assemblies:
            for f in asm.get("fragments", []):
                if f.get("number") == frag_num and f.get("blast_identity") is not None:
                    values.append(f["blast_identity"])
        col_averages.append(sum(values) / len(values) if values else None)

    # Column headers
    for col_idx, frag_num in enumerate(all_frag_numbers):
        cx = grid_x + col_idx * cell_size + cell_size / 2
        # Fragment number
        parts.append(
            f'<text x="{cx:.1f}" y="{grid_y - 28}" font-size="11" font-weight="bold" '
            f'fill="#333" text-anchor="middle" dominant-baseline="central">{frag_num}</text>'
        )
        # Mean identity
        avg = col_averages[col_idx]
        avg_label = f"{avg:.1f}%" if avg is not None else "N/A"
        parts.append(
            f'<text x="{cx:.1f}" y="{grid_y - 12}" font-size="9" '
            f'fill="#666" text-anchor="middle" dominant-baseline="central">{avg_label}</text>'
        )

    # Rows
    for row_idx, asm in enumerate(sorted_assemblies):
        cy = grid_y + row_idx * cell_size + cell_size / 2

        # Quality badge (colored dot)
        quality = asm.get("quality", "").upper()
        q_color = _QUALITY_COLORS.get(quality, "#999999")
        parts.append(
            f'<circle cx="12" cy="{cy:.1f}" r="5" fill="{q_color}"/>'
        )

        # Assembly name (truncated)
        name = asm.get("name", "")
        display_name = name if len(name) <= 22 else name[:19] + "..."
        parts.append(
            f'<text x="24" y="{cy:.1f}" font-size="10" fill="#333" '
            f'dominant-baseline="central">{_escape(display_name)}</text>'
        )

        # Build lookup for this assembly's fragments
        frag_map: Dict[int, Optional[float]] = {}
        for f in asm.get("fragments", []):
            frag_map[f.get("number")] = f.get("blast_identity")

        # Cells
        for col_idx, frag_num in enumerate(all_frag_numbers):
            identity = frag_map.get(frag_num)
            color = _identity_color(identity)
            cell_x = grid_x + col_idx * cell_size
            cell_y = grid_y + row_idx * cell_size

            parts.append(
                f'<rect x="{cell_x}" y="{cell_y}" width="{cell_size}" '
                f'height="{cell_size}" fill="{color}" stroke="#eee" stroke-width="1"/>'
            )

            # Cell label
            if identity is not None:
                parts.append(
                    f'<text x="{cell_x + cell_size / 2:.1f}" '
                    f'y="{cell_y + cell_size / 2:.1f}" font-size="8" fill="#333" '
                    f'text-anchor="middle" dominant-baseline="central">'
                    f'{identity:.1f}</text>'
                )

    # Close SVG
    parts.append("</svg>")
    return "\n".join(parts)
