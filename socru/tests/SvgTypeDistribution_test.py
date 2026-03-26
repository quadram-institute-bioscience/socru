"""Tests for the SVG GS type distribution bar chart generator."""

import unittest
from socru.SvgTypeDistribution import generate_type_distribution_svg


class TestSvgTypeDistributionBasic(unittest.TestCase):
    def test_basic_svg_generation(self):
        """Test that valid SVG is produced with expected structure."""
        type_counts = {"GS1.0": 15, "GS2.1": 8, "GS3.0": 3}
        svg = generate_type_distribution_svg(type_counts)
        self.assertIn("<svg", svg)
        self.assertIn("</svg>", svg)
        self.assertIn("GS Type Distribution", svg)
        self.assertIn("GS1.0", svg)
        self.assertIn("GS2.1", svg)
        self.assertIn("GS3.0", svg)

    def test_percentage_labels(self):
        """Percentage labels should appear on bars."""
        type_counts = {"GS1.0": 10, "GS2.0": 10}
        svg = generate_type_distribution_svg(type_counts)
        self.assertIn("50%", svg)

    def test_viewbox_dimensions(self):
        """SVG viewBox should match requested dimensions."""
        svg = generate_type_distribution_svg({"GS1.0": 1}, width=1000, height=500)
        self.assertIn('viewBox="0 0 1000 500"', svg)


class TestSvgTypeDistributionSingleType(unittest.TestCase):
    def test_single_type(self):
        """Chart should work with a single GS type."""
        svg = generate_type_distribution_svg({"GS1.0": 5})
        self.assertIn("<svg", svg)
        self.assertIn("GS1.0", svg)
        self.assertIn("100%", svg)


class TestSvgTypeDistributionMultipleTypes(unittest.TestCase):
    def test_many_types(self):
        """Chart should handle many different GS types."""
        type_counts = {f"GS{i}.0": i for i in range(1, 8)}
        svg = generate_type_distribution_svg(type_counts)
        self.assertIn("<svg", svg)
        for i in range(1, 8):
            self.assertIn(f"GS{i}.0", svg)

    def test_empty_counts(self):
        """Should handle empty type_counts dict."""
        svg = generate_type_distribution_svg({})
        self.assertIn("<svg", svg)
        self.assertIn("</svg>", svg)


class TestSvgTypeDistributionQualityStacked(unittest.TestCase):
    def test_stacked_bars_with_quality(self):
        """When quality_counts provided, bars should be stacked with quality colors."""
        type_counts = {"GS1.0": 10, "GS2.1": 5}
        quality_counts = {
            "GS1.0": {"GREEN": 7, "AMBER": 2, "RED": 1},
            "GS2.1": {"GREEN": 3, "AMBER": 1, "RED": 1},
        }
        svg = generate_type_distribution_svg(type_counts, quality_counts=quality_counts)
        self.assertIn("<svg", svg)
        # Quality colors should be present
        self.assertIn("#2ca02c", svg)   # green
        self.assertIn("#ff9900", svg)   # amber
        self.assertIn("#d62728", svg)   # red
        # Legend should appear
        self.assertIn("GREEN", svg)
        self.assertIn("AMBER", svg)
        self.assertIn("RED", svg)

    def test_no_legend_without_quality(self):
        """Legend should not appear when quality_counts is None."""
        svg = generate_type_distribution_svg({"GS1.0": 5})
        # GREEN/AMBER/RED legend labels should not be present
        # (the word GREEN could appear elsewhere, so check for legend pattern)
        self.assertNotIn("GREEN</text>", svg)


if __name__ == "__main__":
    unittest.main()
