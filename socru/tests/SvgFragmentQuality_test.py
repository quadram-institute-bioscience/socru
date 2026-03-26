"""Tests for the SVG fragment quality bar chart generator."""

import unittest
from socru.SvgFragmentQuality import generate_fragment_quality_svg


class TestSvgFragmentQualityBasic(unittest.TestCase):
    def test_basic_svg_generation(self):
        """Test that valid SVG is produced with expected structure."""
        fragments = [
            {"number": 1, "reversed": False, "blast_identity": 99.5, "blast_alignment_length": 1000000, "length": 1000000},
            {"number": 2, "reversed": True, "blast_identity": 97.2, "blast_alignment_length": 800000, "length": 800000},
            {"number": 3, "reversed": False, "blast_identity": 92.1, "blast_alignment_length": 600000, "length": 600000},
        ]
        svg = generate_fragment_quality_svg(fragments, genome_name="test_assembly")
        self.assertIn("<svg", svg)
        self.assertIn("</svg>", svg)
        self.assertIn("test_assembly", svg)
        self.assertIn("Fragment Quality", svg)

    def test_fragment_labels_present(self):
        """Fragment numbers should appear as Y-axis labels."""
        fragments = [
            {"number": 1, "reversed": False, "blast_identity": 99.0, "length": 1000},
            {"number": 2, "reversed": False, "blast_identity": 98.0, "length": 1000},
        ]
        svg = generate_fragment_quality_svg(fragments)
        # Labels should be present (as text content, not just any "1" in SVG)
        self.assertIn(">1<", svg)
        self.assertIn(">2<", svg)

    def test_reversed_fragment_has_prime(self):
        """Reversed fragments should have prime symbol in label."""
        fragments = [
            {"number": 3, "reversed": True, "blast_identity": 95.0, "length": 1000},
        ]
        svg = generate_fragment_quality_svg(fragments)
        self.assertIn("3'", svg)


class TestSvgFragmentQualityColors(unittest.TestCase):
    def test_excellent_identity_green(self):
        """Identity >99% should produce green bars."""
        fragments = [
            {"number": 1, "reversed": False, "blast_identity": 99.8, "length": 1000},
        ]
        svg = generate_fragment_quality_svg(fragments)
        self.assertIn("#2ca02c", svg)  # green

    def test_good_identity_yellow(self):
        """Identity 95-99% should produce yellow bars."""
        fragments = [
            {"number": 1, "reversed": False, "blast_identity": 97.0, "length": 1000},
        ]
        svg = generate_fragment_quality_svg(fragments)
        self.assertIn("#f0c800", svg)  # yellow

    def test_fair_identity_orange(self):
        """Identity 90-95% should produce orange bars."""
        fragments = [
            {"number": 1, "reversed": False, "blast_identity": 92.0, "length": 1000},
        ]
        svg = generate_fragment_quality_svg(fragments)
        self.assertIn("#ff7f00", svg)  # orange

    def test_poor_identity_red(self):
        """Identity <90% should produce red bars."""
        fragments = [
            {"number": 1, "reversed": False, "blast_identity": 85.0, "length": 1000},
        ]
        svg = generate_fragment_quality_svg(fragments)
        self.assertIn("#d62728", svg)  # red


class TestSvgFragmentQualityUnknown(unittest.TestCase):
    def test_unknown_identity_dashed(self):
        """Fragments with no identity should show dashed gray bars."""
        fragments = [
            {"number": 1, "reversed": False, "blast_identity": None, "length": 1000},
        ]
        svg = generate_fragment_quality_svg(fragments)
        self.assertIn("stroke-dasharray", svg)
        self.assertIn("#999999", svg)
        self.assertIn("N/A", svg)

    def test_mixed_known_and_unknown(self):
        """Chart should handle a mix of known and unknown identities."""
        fragments = [
            {"number": 1, "reversed": False, "blast_identity": 99.5, "length": 1000},
            {"number": 2, "reversed": False, "blast_identity": None, "length": 1000},
            {"number": 3, "reversed": True, "blast_identity": 88.0, "length": 1000},
        ]
        svg = generate_fragment_quality_svg(fragments)
        self.assertIn("<svg", svg)
        self.assertIn("N/A", svg)
        self.assertIn("99.5%", svg)
        self.assertIn("88.0%", svg)

    def test_empty_fragments(self):
        """Should handle empty fragment list without error."""
        svg = generate_fragment_quality_svg([])
        self.assertIn("<svg", svg)
        self.assertIn("</svg>", svg)

    def test_custom_dimensions(self):
        """Custom width and height should appear in viewBox."""
        svg = generate_fragment_quality_svg([], width=900, height=500)
        self.assertIn('viewBox="0 0 900 500"', svg)


if __name__ == "__main__":
    unittest.main()
