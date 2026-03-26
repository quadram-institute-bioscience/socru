"""Tests for the SVG confidence heatmap visualization."""

import unittest
from socru.SvgConfidenceHeatmap import generate_confidence_heatmap_svg, _identity_color


class TestSvgConfidenceHeatmap(unittest.TestCase):
    def _make_assemblies(self):
        """Helper to build 3 assemblies x 5 fragments test data."""
        return [
            {
                "name": "assembly_A",
                "gs_type": "GS1.0",
                "quality": "GREEN",
                "fragments": [
                    {"number": 1, "blast_identity": 100.0},
                    {"number": 2, "blast_identity": 99.5},
                    {"number": 3, "blast_identity": 97.0},
                    {"number": 4, "blast_identity": 95.0},
                    {"number": 5, "blast_identity": 88.0},
                ],
            },
            {
                "name": "assembly_B",
                "gs_type": "GS1.0",
                "quality": "AMBER",
                "fragments": [
                    {"number": 1, "blast_identity": 99.8},
                    {"number": 2, "blast_identity": 98.0},
                    {"number": 3, "blast_identity": 96.5},
                    {"number": 4, "blast_identity": None},
                    {"number": 5, "blast_identity": 91.0},
                ],
            },
            {
                "name": "assembly_C",
                "gs_type": "GS2.0",
                "quality": "RED",
                "fragments": [
                    {"number": 1, "blast_identity": 100.0},
                    {"number": 2, "blast_identity": 100.0},
                    {"number": 3, "blast_identity": 100.0},
                    {"number": 4, "blast_identity": 100.0},
                    {"number": 5, "blast_identity": 100.0},
                ],
            },
        ]

    def test_basic_heatmap(self):
        """Test heatmap with 3 assemblies x 5 fragments."""
        assemblies = self._make_assemblies()
        svg = generate_confidence_heatmap_svg(assemblies)
        self.assertIn("<svg", svg)
        self.assertIn("</svg>", svg)
        self.assertIn("Fragment Identity Heatmap", svg)
        self.assertIn("assembly_A", svg)
        self.assertIn("assembly_B", svg)
        self.assertIn("assembly_C", svg)
        # Should have rectangle cells
        self.assertIn("<rect", svg)

    def test_color_coding_100_percent(self):
        """100% identity maps to white."""
        color = _identity_color(100.0)
        self.assertEqual(color, "#ffffff")

    def test_color_coding_below_90(self):
        """Below 90% identity maps to red."""
        color = _identity_color(85.0)
        self.assertEqual(color, "#d62728")

    def test_color_coding_missing(self):
        """None identity maps to gray."""
        color = _identity_color(None)
        self.assertEqual(color, "#dddddd")

    def test_color_coding_gradient(self):
        """Check gradient produces distinct colors at key thresholds."""
        c100 = _identity_color(100.0)
        c99 = _identity_color(99.0)
        c97 = _identity_color(97.0)
        c95 = _identity_color(95.0)
        c90 = _identity_color(90.0)
        # All should be different
        colors = [c100, c99, c97, c95, c90]
        self.assertEqual(len(set(colors)), 5, f"Expected 5 distinct colors, got {colors}")

    def test_single_assembly(self):
        """Test with a single assembly."""
        assemblies = [
            {
                "name": "solo",
                "gs_type": "GS1.0",
                "quality": "GREEN",
                "fragments": [
                    {"number": 1, "blast_identity": 99.0},
                    {"number": 2, "blast_identity": 98.0},
                ],
            },
        ]
        svg = generate_confidence_heatmap_svg(assemblies)
        self.assertIn("<svg", svg)
        self.assertIn("solo", svg)

    def test_empty_input(self):
        """Empty assembly list should produce a valid SVG with message."""
        svg = generate_confidence_heatmap_svg([])
        self.assertIn("<svg", svg)
        self.assertIn("</svg>", svg)
        self.assertIn("No assemblies", svg)

    def test_sorted_by_gs_type(self):
        """Assemblies should be sorted by GS type then name in output."""
        assemblies = [
            {
                "name": "z_last",
                "gs_type": "GS2.0",
                "quality": "GREEN",
                "fragments": [{"number": 1, "blast_identity": 99.0}],
            },
            {
                "name": "a_first",
                "gs_type": "GS1.0",
                "quality": "GREEN",
                "fragments": [{"number": 1, "blast_identity": 99.0}],
            },
            {
                "name": "m_middle",
                "gs_type": "GS1.0",
                "quality": "AMBER",
                "fragments": [{"number": 1, "blast_identity": 95.0}],
            },
        ]
        svg = generate_confidence_heatmap_svg(assemblies)
        # a_first and m_middle (GS1.0) should appear before z_last (GS2.0)
        pos_a = svg.index("a_first")
        pos_m = svg.index("m_middle")
        pos_z = svg.index("z_last")
        self.assertLess(pos_a, pos_m)
        self.assertLess(pos_m, pos_z)

    def test_column_headers_show_fragment_numbers(self):
        """Column headers should include fragment numbers."""
        assemblies = [
            {
                "name": "test",
                "gs_type": "GS1.0",
                "quality": "GREEN",
                "fragments": [
                    {"number": 3, "blast_identity": 99.0},
                    {"number": 7, "blast_identity": 98.0},
                ],
            },
        ]
        svg = generate_confidence_heatmap_svg(assemblies)
        # Fragment numbers 3 and 7 should appear as headers
        self.assertIn(">3<", svg)
        self.assertIn(">7<", svg)

    def test_quality_badge_colors(self):
        """Quality badges should use the correct colors."""
        assemblies = [
            {
                "name": "green_asm",
                "gs_type": "GS1.0",
                "quality": "GREEN",
                "fragments": [{"number": 1, "blast_identity": 100.0}],
            },
        ]
        svg = generate_confidence_heatmap_svg(assemblies)
        self.assertIn("#2ca02c", svg)  # GREEN color


if __name__ == "__main__":
    unittest.main()
