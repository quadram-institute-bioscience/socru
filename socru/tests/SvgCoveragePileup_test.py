"""Tests for the SVG coverage pileup visualization."""

import unittest
from socru.SvgCoveragePileup import generate_coverage_pileup_svg


class TestSvgCoveragePileup(unittest.TestCase):
    def test_basic_svg_with_coverage_data(self):
        """Test SVG generation with valid coverage depth arrays."""
        fragments = [
            {
                "number": 1,
                "length": 100,
                "coverage_depths": [1, 2, 3, 4, 5, 4, 3, 2, 1, 0] * 10,
            },
            {
                "number": 2,
                "length": 80,
                "coverage_depths": [3, 3, 3, 4, 4, 4, 5, 5, 5, 2] * 8,
            },
        ]
        svg = generate_coverage_pileup_svg(fragments, genome_name="test_genome")
        self.assertIn("<svg", svg)
        self.assertIn("</svg>", svg)
        self.assertIn("Frag 1", svg)
        self.assertIn("Frag 2", svg)
        self.assertIn("test_genome", svg)
        self.assertIn("Coverage Pileup", svg)
        # Should contain polygon for area chart
        self.assertIn("<polygon", svg)
        # Should contain dashed mean line
        self.assertIn("stroke-dasharray", svg)

    def test_no_coverage_data(self):
        """Fragments with None coverage should show a placeholder."""
        fragments = [
            {"number": 1, "length": 1000, "coverage_depths": None},
            {"number": 2, "length": 800, "coverage_depths": None},
        ]
        svg = generate_coverage_pileup_svg(fragments)
        self.assertIn("<svg", svg)
        self.assertIn("No data", svg)
        # Should not contain polygon since there is no data
        self.assertNotIn("<polygon", svg)

    def test_single_fragment(self):
        """Test with exactly one fragment."""
        fragments = [
            {
                "number": 5,
                "length": 50,
                "coverage_depths": [2, 4, 6, 8, 10, 8, 6, 4, 2, 0],
            },
        ]
        svg = generate_coverage_pileup_svg(fragments)
        self.assertIn("<svg", svg)
        self.assertIn("Frag 5", svg)
        self.assertIn("<polygon", svg)

    def test_multiple_fragments_mixed(self):
        """Test with a mix of fragments: some with data, some without."""
        fragments = [
            {"number": 1, "length": 100, "coverage_depths": [1, 2, 3, 2, 1]},
            {"number": 2, "length": 80, "coverage_depths": None},
            {"number": 3, "length": 120, "coverage_depths": [5, 5, 5, 5, 5]},
        ]
        svg = generate_coverage_pileup_svg(fragments)
        self.assertIn("Frag 1", svg)
        self.assertIn("Frag 2", svg)
        self.assertIn("Frag 3", svg)
        self.assertIn("No data", svg)
        # Two fragments with data produce two polygons
        self.assertEqual(svg.count("<polygon"), 2)

    def test_empty_input(self):
        """Empty fragment list should produce a valid SVG with message."""
        svg = generate_coverage_pileup_svg([])
        self.assertIn("<svg", svg)
        self.assertIn("</svg>", svg)
        self.assertIn("No fragments", svg)

    def test_zero_coverage(self):
        """All-zero coverage should still render without errors."""
        fragments = [
            {"number": 1, "length": 50, "coverage_depths": [0, 0, 0, 0, 0]},
        ]
        svg = generate_coverage_pileup_svg(fragments)
        self.assertIn("<svg", svg)
        self.assertIn("<polygon", svg)

    def test_empty_depths_list(self):
        """Empty list for coverage_depths treated same as None."""
        fragments = [
            {"number": 1, "length": 50, "coverage_depths": []},
        ]
        svg = generate_coverage_pileup_svg(fragments)
        self.assertIn("No data", svg)

    def test_custom_dimensions(self):
        """Test with custom width and height."""
        fragments = [
            {"number": 1, "length": 100, "coverage_depths": [1, 2, 3]},
        ]
        svg = generate_coverage_pileup_svg(fragments, width=1200, height=500)
        self.assertIn('width="1200"', svg)
        self.assertIn('height="500"', svg)

    def test_large_coverage_array_downsampled(self):
        """Large coverage arrays should be downsampled for rendering."""
        depths = list(range(10000))
        fragments = [
            {"number": 1, "length": 10000, "coverage_depths": depths},
        ]
        svg = generate_coverage_pileup_svg(fragments)
        self.assertIn("<svg", svg)
        self.assertIn("<polygon", svg)

    def test_mean_depth_label_present(self):
        """Mean depth annotation should appear in the SVG."""
        fragments = [
            {"number": 1, "length": 50, "coverage_depths": [10, 10, 10, 10, 10]},
        ]
        svg = generate_coverage_pileup_svg(fragments)
        # Mean of all 10s is 10.0
        self.assertIn("10.0x", svg)


if __name__ == "__main__":
    unittest.main()
