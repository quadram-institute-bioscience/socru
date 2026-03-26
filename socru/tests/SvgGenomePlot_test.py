"""Tests for the SVG circular genome diagram generator."""

import unittest
from socru.SvgGenomePlot import generate_genome_svg


class TestSvgGenomePlot(unittest.TestCase):
    def test_basic_svg_generation(self):
        """Test that SVG is generated with correct structure."""
        fragments = [
            {'number': 1, 'reversed': False, 'length': 1000000, 'coords': [(0, 1000000)], 'is_dnaA': True, 'is_dif': False},
            {'number': 2, 'reversed': True, 'length': 800000, 'coords': [(1005000, 1805000)], 'is_dnaA': False, 'is_dif': False},
            {'number': 3, 'reversed': False, 'length': 600000, 'coords': [(1810000, 2410000)], 'is_dnaA': False, 'is_dif': True},
        ]
        operons = [
            {'start': 1000000, 'end': 1005000, 'direction': 'forward'},
            {'start': 1805000, 'end': 1810000, 'direction': 'reverse'},
        ]
        svg = generate_genome_svg(
            fragments=fragments,
            operons=operons,
            genome_length=2500000,
            gs_type="GS1.0",
            quality="GREEN",
            genome_name="test_genome"
        )
        self.assertIn('<svg', svg)
        self.assertIn('</svg>', svg)
        self.assertIn('GS1.0', svg)
        self.assertIn('test_genome', svg)
        self.assertIn('Ori', svg)
        self.assertIn('Ter', svg)

    def test_reversed_fragment_has_hatch(self):
        """Reversed fragments should have a hatch pattern."""
        fragments = [
            {'number': 1, 'reversed': True, 'length': 1000000, 'coords': [(0, 1000000)], 'is_dnaA': True, 'is_dif': False},
        ]
        svg = generate_genome_svg(fragments=fragments, operons=[], genome_length=1000000, gs_type="GS1.1", quality="AMBER")
        self.assertIn('pattern', svg.lower())  # SVG hatch pattern

    def test_empty_fragments(self):
        """Should handle empty fragment list gracefully."""
        svg = generate_genome_svg(fragments=[], operons=[], genome_length=0, gs_type="GS0.0", quality="RED")
        self.assertIn('<svg', svg)

    def test_seven_fragments(self):
        """Test with 7 fragments (typical E. coli)."""
        fragments = [
            {'number': i, 'reversed': (i % 3 == 0), 'length': 500000 + i * 100000,
             'coords': [(i * 600000, i * 600000 + 500000 + i * 100000)],
             'is_dnaA': (i == 1), 'is_dif': (i == 4)}
            for i in range(1, 8)
        ]
        operons = [
            {'start': i * 600000 - 5000, 'end': i * 600000, 'direction': 'forward' if i % 2 == 0 else 'reverse'}
            for i in range(1, 7)
        ]
        svg = generate_genome_svg(fragments=fragments, operons=operons, genome_length=4800000, gs_type="GS1.3", quality="GREEN", genome_name="E_coli_K12")
        self.assertIn('E_coli_K12', svg)
        # Should have 7 different colored arcs
        for i in range(1, 8):
            self.assertIn(str(i), svg)

    def test_quality_colors(self):
        """Each quality level should produce the appropriate badge color."""
        for quality, expected_color in [("GREEN", "#2ca02c"), ("AMBER", "#ff9900"), ("RED", "#d62728")]:
            svg = generate_genome_svg(
                fragments=[{'number': 1, 'reversed': False, 'length': 100, 'coords': [(0, 100)], 'is_dnaA': False, 'is_dif': False}],
                operons=[], genome_length=100, gs_type="GS1.0", quality=quality,
            )
            self.assertIn(expected_color, svg, f"Expected color {expected_color} for quality {quality}")

    def test_tick_marks(self):
        """Tick marks should appear for a large genome."""
        svg = generate_genome_svg(
            fragments=[{'number': 1, 'reversed': False, 'length': 2000000, 'coords': [(0, 2000000)], 'is_dnaA': False, 'is_dif': False}],
            operons=[], genome_length=2000000, gs_type="GS1.0", quality="GREEN",
        )
        self.assertIn('500kb', svg)
        self.assertIn('1000kb', svg)

    def test_viewbox_dimensions(self):
        """SVG viewBox should match requested dimensions."""
        svg = generate_genome_svg(fragments=[], operons=[], genome_length=0, gs_type="GS0.0", quality="RED", width=900, height=900)
        self.assertIn('viewBox="0 0 900 900"', svg)

    def test_genome_length_label(self):
        """Genome size in Mb should appear in the centre."""
        svg = generate_genome_svg(
            fragments=[{'number': 1, 'reversed': False, 'length': 4600000, 'coords': [(0, 4600000)], 'is_dnaA': False, 'is_dif': False}],
            operons=[], genome_length=4600000, gs_type="GS1.0", quality="GREEN",
        )
        self.assertIn('4.60 Mb', svg)

    def test_long_genome_name_truncated(self):
        """Very long genome names should be truncated."""
        long_name = "a" * 50
        svg = generate_genome_svg(
            fragments=[], operons=[], genome_length=0, gs_type="GS0.0", quality="RED", genome_name=long_name,
        )
        self.assertIn('...', svg)
        # Full name should NOT appear
        self.assertNotIn(long_name, svg)


if __name__ == '__main__':
    unittest.main()
