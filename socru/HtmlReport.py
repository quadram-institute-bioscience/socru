"""
Self-contained HTML report generator for socru analysis results.

Generates a single HTML file with inline CSS and vanilla JS that can be
shared, archived, and viewed in any modern browser without external
dependencies.

Classes:
    HtmlReport: Builds and writes the complete HTML report.
"""

from datetime import datetime
from typing import Any, Dict, List


class HtmlReport:
    """Generate a self-contained HTML report for socru analysis results.

    The report includes summary cards, a sortable results table with
    expandable per-assembly detail panels, and a QC summary section.
    All CSS and JS is inlined so the output is a single portable file.
    """

    def __init__(
        self,
        results: List[Dict[str, Any]],
        species: str,
        tool_version: str = "2.2.4",
    ) -> None:
        """Initialise the report generator.

        Args:
            results: List of dicts, each with keys:
                - genome_file, gs_type, quality, confidence_score,
                - fragment_pattern, num_operons, is_novel, qc_flags,
                - fragments (list of fragment dicts with blast stats)
                - genome_length, validation_passed
            species: Species name used in the analysis.
            tool_version: Socru version string.
        """
        self.results = results
        self.species = species
        self.tool_version = tool_version

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def generate(self) -> str:
        """Return the complete HTML report as a string."""
        return (
            self._doctype()
            + self._head()
            + "<body>\n"
            + self._header_section()
            + self._summary_cards()
            + self._results_table()
            + self._qc_summary()
            + self._footer()
            + self._inline_js()
            + "</body>\n</html>\n"
        )

    def save(self, filename: str) -> None:
        """Write the HTML report to *filename*."""
        with open(filename, "w") as fh:
            fh.write(self.generate())

    # ------------------------------------------------------------------
    # Private helpers -- HTML skeleton
    # ------------------------------------------------------------------

    @staticmethod
    def _doctype() -> str:
        return "<!DOCTYPE html>\n<html lang=\"en\">\n"

    def _head(self) -> str:
        return (
            "<head>\n"
            '<meta charset="utf-8">\n'
            '<meta name="viewport" content="width=device-width, initial-scale=1">\n'
            "<title>Socru Genome Structure Analysis Report</title>\n"
            "<style>\n" + self._css() + "</style>\n"
            "</head>\n"
        )

    # ------------------------------------------------------------------
    # CSS
    # ------------------------------------------------------------------

    @staticmethod
    def _css() -> str:
        return """\
:root {
  --green: #2da44e;
  --green-bg: #dafbe1;
  --amber: #bf8700;
  --amber-bg: #fff8c5;
  --red: #cf222e;
  --red-bg: #ffebe9;
  --border: #d0d7de;
  --bg-subtle: #f6f8fa;
  --fg: #1f2328;
  --fg-muted: #656d76;
  --font: -apple-system, BlinkMacSystemFont, 'Segoe UI', Noto Sans, Helvetica, Arial, sans-serif;
}
*, *::before, *::after { box-sizing: border-box; }
body { font-family: var(--font); color: var(--fg); margin: 0; padding: 24px; max-width: 1280px; margin: 0 auto; line-height: 1.5; }
h1, h2, h3 { margin-top: 0; }
a { color: #0969da; text-decoration: none; }
a:hover { text-decoration: underline; }

/* Header */
.report-header { border-bottom: 1px solid var(--border); padding-bottom: 16px; margin-bottom: 24px; }
.report-header h1 { font-size: 24px; font-weight: 600; }
.meta { color: var(--fg-muted); font-size: 14px; }
.meta span { margin-right: 20px; }

/* Summary cards */
.cards { display: flex; flex-wrap: wrap; gap: 16px; margin-bottom: 32px; }
.card { flex: 1 1 180px; border: 1px solid var(--border); border-radius: 6px; padding: 16px; background: var(--bg-subtle); min-width: 160px; }
.card .label { font-size: 12px; color: var(--fg-muted); text-transform: uppercase; letter-spacing: 0.5px; margin-bottom: 4px; }
.card .value { font-size: 28px; font-weight: 600; }

/* Badges */
.badge { display: inline-block; padding: 2px 8px; border-radius: 12px; font-size: 12px; font-weight: 600; line-height: 1.5; }
.badge-green { background: var(--green-bg); color: var(--green); }
.badge-amber { background: var(--amber-bg); color: var(--amber); }
.badge-red { background: var(--red-bg); color: var(--red); }
.badge-info { background: #ddf4ff; color: #0969da; }
.badge-warning { background: var(--amber-bg); color: var(--amber); }
.badge-error { background: var(--red-bg); color: var(--red); }

/* Table */
.table-wrap { overflow-x: auto; margin-bottom: 32px; }
table { width: 100%; border-collapse: collapse; font-size: 14px; }
th, td { text-align: left; padding: 8px 12px; border-bottom: 1px solid var(--border); }
th { background: var(--bg-subtle); font-weight: 600; cursor: pointer; user-select: none; white-space: nowrap; }
th:hover { background: #eaeef2; }
th .sort-arrow { margin-left: 4px; font-size: 10px; color: var(--fg-muted); }
tr.row-green { background: #f6fef9; }
tr.row-amber { background: #fffdf5; }
tr.row-red { background: #fff5f5; }
tr.clickable { cursor: pointer; }
tr.clickable:hover { background: #eaeef2; }

/* Detail panel */
.detail-panel { display: none; }
.detail-panel td { padding: 16px 12px; background: var(--bg-subtle); }
.detail-content { max-width: 900px; }
.detail-content h3 { font-size: 14px; margin: 12px 0 6px; }
.detail-content table { font-size: 13px; margin-bottom: 12px; }
.detail-content th { background: #eaeef2; }

/* Fragment visualisation */
.frag-viz { display: flex; gap: 2px; margin: 8px 0; }
.frag-block { padding: 4px 10px; border-radius: 4px; font-size: 12px; font-weight: 600; color: #fff; text-align: center; }
.frag-fwd { background: #0969da; }
.frag-rev { background: #8250df; }
.frag-unknown { background: #656d76; }

/* QC Summary */
.qc-summary { margin-bottom: 32px; }
.qc-summary table { max-width: 640px; }

/* Footer */
.report-footer { border-top: 1px solid var(--border); padding-top: 12px; font-size: 12px; color: var(--fg-muted); margin-top: 32px; }

/* Print */
@media print {
  body { padding: 0; }
  .detail-panel { display: table-row !important; }
  th { background: #eee !important; }
  tr.clickable:hover { background: inherit; }
}
"""

    # ------------------------------------------------------------------
    # Sections
    # ------------------------------------------------------------------

    def _header_section(self) -> str:
        date_str = datetime.now().strftime("%Y-%m-%d %H:%M")
        return (
            '<div class="report-header">\n'
            "  <h1>Socru Genome Structure Analysis Report</h1>\n"
            '  <div class="meta">\n'
            f'    <span>Species: <strong>{_esc(self.species)}</strong></span>\n'
            f'    <span>Date: {date_str}</span>\n'
            f'    <span>Socru v{_esc(self.tool_version)}</span>\n'
            f'    <span>Assemblies: {len(self.results)}</span>\n'
            "  </div>\n"
            "</div>\n"
        )

    def _summary_cards(self) -> str:
        total = len(self.results)
        green = sum(1 for r in self.results if r.get("quality") == "GREEN")
        amber = sum(1 for r in self.results if r.get("quality") == "AMBER")
        red = sum(1 for r in self.results if r.get("quality") == "RED")
        novel = sum(1 for r in self.results if r.get("is_novel"))
        scores = [r.get("confidence_score", 0) for r in self.results]
        avg_score = sum(scores) / len(scores) if scores else 0

        return (
            '<div class="cards">\n'
            + self._card("Total Assemblies", str(total))
            + self._card_badge("GREEN", green, "green")
            + self._card_badge("AMBER", amber, "amber")
            + self._card_badge("RED", red, "red")
            + self._card("Novel Types", str(novel))
            + self._card("Avg Confidence", f"{avg_score:.1f}")
            + "</div>\n"
        )

    @staticmethod
    def _card(label: str, value: str) -> str:
        return (
            '  <div class="card">\n'
            f'    <div class="label">{_esc(label)}</div>\n'
            f'    <div class="value">{_esc(value)}</div>\n'
            "  </div>\n"
        )

    @staticmethod
    def _card_badge(label: str, count: int, colour: str) -> str:
        return (
            '  <div class="card">\n'
            f'    <div class="label">{_esc(label)}</div>\n'
            f'    <div class="value"><span class="badge badge-{colour}">{count}</span></div>\n'
            "  </div>\n"
        )

    # ------------------------------------------------------------------
    # Results table
    # ------------------------------------------------------------------

    def _results_table(self) -> str:
        header = (
            "<h2>Results</h2>\n"
            '<div class="table-wrap">\n'
            '<table id="results-table">\n'
            "<thead><tr>\n"
            '  <th data-col="0">Assembly <span class="sort-arrow"></span></th>\n'
            '  <th data-col="1">GS Type <span class="sort-arrow"></span></th>\n'
            '  <th data-col="2">Quality <span class="sort-arrow"></span></th>\n'
            '  <th data-col="3">Confidence <span class="sort-arrow"></span></th>\n'
            '  <th data-col="4">Operons <span class="sort-arrow"></span></th>\n'
            '  <th data-col="5">Fragments <span class="sort-arrow"></span></th>\n'
            '  <th data-col="6">Novel <span class="sort-arrow"></span></th>\n'
            '  <th data-col="7">QC Flags <span class="sort-arrow"></span></th>\n'
            "</tr></thead>\n<tbody>\n"
        )

        rows = ""
        for idx, r in enumerate(self.results):
            quality = r.get("quality", "RED")
            row_class = f"row-{quality.lower()}"
            badge_class = f"badge-{quality.lower()}"
            novel_text = "Yes" if r.get("is_novel") else "No"
            flags = r.get("qc_flags", [])
            flag_count = str(len(flags)) if flags else "0"
            confidence = r.get("confidence_score", 0)

            rows += (
                f'<tr class="clickable {row_class}" data-idx="{idx}">\n'
                f"  <td>{_esc(str(r.get('genome_file', '')))}</td>\n"
                f"  <td>{_esc(str(r.get('gs_type', '')))}</td>\n"
                f'  <td><span class="badge {badge_class}">{_esc(quality)}</span></td>\n'
                f"  <td>{confidence:.1f}</td>\n"
                f"  <td>{r.get('num_operons', 0)}</td>\n"
                f"  <td>{_esc(str(r.get('fragment_pattern', '')))}</td>\n"
                f"  <td>{novel_text}</td>\n"
                f"  <td>{flag_count}</td>\n"
                "</tr>\n"
            )
            rows += self._detail_row(idx, r)

        return header + rows + "</tbody>\n</table>\n</div>\n"

    def _detail_row(self, idx: int, r: Dict[str, Any]) -> str:
        fragments = r.get("fragments", [])
        flags = r.get("qc_flags", [])
        pattern = str(r.get("fragment_pattern", ""))

        # Fragment detail table
        frag_table = ""
        if fragments:
            frag_rows = ""
            for f in fragments:
                orient = "Reversed" if f.get("reversed") else "Forward"
                identity = f.get('blast_identity')
                identity_str = f"{identity:.1f}" if identity is not None else "N/A"
                align_len = f.get('blast_alignment_length') or 0
                bit_score = f.get('blast_bit_score') or 0
                frag_rows += (
                    f"<tr>"
                    f"<td>{_esc(str(f.get('number', '')))}</td>"
                    f"<td>{orient}</td>"
                    f"<td>{identity_str}</td>"
                    f"<td>{align_len}</td>"
                    f"<td>{bit_score}</td>"
                    f"<td>{f.get('length', 0):,}</td>"
                    f"</tr>\n"
                )
            frag_table = (
                "<h3>Fragment Details</h3>\n"
                "<table><thead><tr>"
                "<th>Fragment</th><th>Orientation</th><th>BLAST Identity</th>"
                "<th>Alignment Length</th><th>Bit Score</th><th>Length</th>"
                "</tr></thead><tbody>\n"
                + frag_rows
                + "</tbody></table>\n"
            )

        # QC flags
        flags_html = ""
        if flags:
            flag_items = ""
            for fl in flags:
                sev = fl.get("severity", "info")
                badge_cls = {
                    "warning": "badge-warning",
                    "error": "badge-error",
                }.get(sev, "badge-info")
                flag_items += (
                    f'<span class="badge {badge_cls}">{_esc(fl.get("code", ""))}</span> '
                    f'{_esc(fl.get("message", ""))}<br>\n'
                )
            flags_html = "<h3>QC Flags</h3>\n" + flag_items

        # Fragment pattern visualisation
        viz = self._fragment_viz(pattern)

        return (
            f'<tr class="detail-panel" id="detail-{idx}"><td colspan="8">\n'
            '<div class="detail-content">\n'
            + viz
            + frag_table
            + flags_html
            + "</div>\n</td></tr>\n"
        )

    @staticmethod
    def _fragment_viz(pattern: str) -> str:
        parts = pattern.split()
        if not parts:
            return ""
        blocks = ""
        for p in parts:
            if p == "?":
                blocks += f'<div class="frag-block frag-unknown">?</div>'
            elif p.endswith("'"):
                blocks += f'<div class="frag-block frag-rev">{_esc(p)}</div>'
            else:
                blocks += f'<div class="frag-block frag-fwd">{_esc(p)}</div>'
        return (
            "<h3>Fragment Pattern</h3>\n"
            f'<div class="frag-viz">{blocks}</div>\n'
        )

    # ------------------------------------------------------------------
    # QC Summary
    # ------------------------------------------------------------------

    def _qc_summary(self) -> str:
        flag_counts: Dict[str, int] = {}
        for r in self.results:
            for fl in r.get("qc_flags", []):
                code = fl.get("code", "UNKNOWN")
                flag_counts[code] = flag_counts.get(code, 0) + 1

        if not flag_counts:
            return (
                '<div class="qc-summary">\n'
                "<h2>QC Summary</h2>\n"
                "<p>No QC flags were raised across any assembly.</p>\n"
                "</div>\n"
            )

        sorted_flags = sorted(flag_counts.items(), key=lambda x: -x[1])
        rows = ""
        for code, count in sorted_flags:
            rows += f"<tr><td>{_esc(code)}</td><td>{count}</td></tr>\n"

        return (
            '<div class="qc-summary">\n'
            "<h2>QC Summary</h2>\n"
            "<table><thead><tr><th>Flag</th><th>Count</th></tr></thead>\n"
            "<tbody>\n" + rows + "</tbody></table>\n"
            "</div>\n"
        )

    # ------------------------------------------------------------------
    # Footer
    # ------------------------------------------------------------------

    def _footer(self) -> str:
        return (
            '<div class="report-footer">\n'
            f"  Generated by Socru v{_esc(self.tool_version)} on "
            f"{datetime.now().strftime('%Y-%m-%d %H:%M')}\n"
            "</div>\n"
        )

    # ------------------------------------------------------------------
    # Inline JS
    # ------------------------------------------------------------------

    @staticmethod
    def _inline_js() -> str:
        return """\
<script>
(function() {
  /* Expand / collapse detail rows */
  document.querySelectorAll('tr.clickable').forEach(function(row) {
    row.addEventListener('click', function() {
      var idx = this.getAttribute('data-idx');
      var panel = document.getElementById('detail-' + idx);
      if (panel) {
        panel.style.display = panel.style.display === 'table-row' ? 'none' : 'table-row';
      }
    });
  });

  /* Column sorting */
  var table = document.getElementById('results-table');
  if (!table) return;
  var headers = table.querySelectorAll('th');
  var sortState = {};

  headers.forEach(function(th) {
    th.addEventListener('click', function() {
      var col = parseInt(this.getAttribute('data-col'), 10);
      if (isNaN(col)) return;
      var asc = sortState[col] !== 'asc';
      sortState = {};
      sortState[col] = asc ? 'asc' : 'desc';

      var tbody = table.querySelector('tbody');
      var mainRows = [];
      var detailMap = {};
      tbody.querySelectorAll('tr').forEach(function(tr) {
        if (tr.classList.contains('detail-panel')) {
          detailMap[tr.id] = tr;
        } else if (tr.classList.contains('clickable')) {
          mainRows.push(tr);
        }
      });

      mainRows.sort(function(a, b) {
        var aVal = a.children[col] ? a.children[col].textContent.trim() : '';
        var bVal = b.children[col] ? b.children[col].textContent.trim() : '';
        var aNum = parseFloat(aVal);
        var bNum = parseFloat(bVal);
        if (!isNaN(aNum) && !isNaN(bNum)) {
          return asc ? aNum - bNum : bNum - aNum;
        }
        return asc ? aVal.localeCompare(bVal) : bVal.localeCompare(aVal);
      });

      while (tbody.firstChild) tbody.removeChild(tbody.firstChild);
      mainRows.forEach(function(row) {
        tbody.appendChild(row);
        var idx = row.getAttribute('data-idx');
        var detail = detailMap['detail-' + idx];
        if (detail) tbody.appendChild(detail);
      });

      /* Update sort arrows */
      headers.forEach(function(h) {
        var arrow = h.querySelector('.sort-arrow');
        if (arrow) arrow.textContent = '';
      });
      var arrow = th.querySelector('.sort-arrow');
      if (arrow) arrow.textContent = asc ? '\\u25B2' : '\\u25BC';
    });
  });
})();
</script>
"""


def _esc(text: str) -> str:
    """Escape HTML special characters."""
    return (
        text.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )
