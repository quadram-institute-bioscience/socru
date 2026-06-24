"""Tests for the linear synteny comparison SVG generator."""

import unittest

from socru.SvgSynteny import _fragment_color, _ribbon_path, generate_synteny_svg


def _make_assembly(name, gs_type, fragments, quality="GREEN"):
    """Helper to build an assembly dict."""
    return {
        "name": name,
        "gs_type": gs_type,
        "quality": quality,
        "fragments": fragments,
    }


def _make_frags(*specs):
    """Build a fragment list from (number, reversed, length) tuples."""
    return [
        {"number": n, "reversed": r, "length": length}
        for n, r, length in specs
    ]


class TestSvgSyntenyBasicStructure(unittest.TestCase):
    """SVG output should always contain fundamental structural elements."""

    def test_svg_tags_present(self):
        assemblies = [
            _make_assembly("A", "GS1.0", _make_frags((1, False, 500), (2, False, 500))),
        ]
        svg = generate_synteny_svg(assemblies)
        self.assertIn("<svg", svg)
        self.assertIn("</svg>", svg)

    def test_contains_rect_elements(self):
        assemblies = [
            _make_assembly("A", "GS1.0", _make_frags((1, False, 500), (2, False, 500))),
        ]
        svg = generate_synteny_svg(assemblies)
        self.assertIn("<rect", svg)

    def test_contains_text_elements(self):
        assemblies = [
            _make_assembly("A", "GS1.0", _make_frags((1, False, 500))),
        ]
        svg = generate_synteny_svg(assemblies)
        self.assertIn("<text", svg)

    def test_empty_assemblies(self):
        svg = generate_synteny_svg([])
        self.assertIn("<svg", svg)
        self.assertIn("</svg>", svg)


class TestSingleAssembly(unittest.TestCase):
    """A single assembly should render a bar with no ribbons."""

    def test_single_assembly_no_path(self):
        assemblies = [
            _make_assembly("single", "GS1.0", _make_frags((1, False, 1000), (2, False, 2000))),
        ]
        svg = generate_synteny_svg(assemblies)
        # No connecting ribbons expected (path elements come from ribbons)
        # The only <path> would be from ribbons; with one assembly there are none.
        self.assertNotIn("<path", svg)

    def test_single_assembly_labels(self):
        assemblies = [
            _make_assembly("my_genome", "GS2.1", _make_frags((1, False, 500))),
        ]
        svg = generate_synteny_svg(assemblies)
        self.assertIn("my_genome", svg)
        self.assertIn("GS2.1", svg)


class TestTwoAssembliesSameOrder(unittest.TestCase):
    """Two assemblies with identical fragment order produce parallel ribbons."""

    def setUp(self):
        frags = _make_frags((1, False, 500), (2, False, 500), (3, False, 500))
        self.assemblies = [
            _make_assembly("A", "GS1.0", frags),
            _make_assembly("B", "GS1.0", list(frags)),  # same order
        ]
        self.svg = generate_synteny_svg(self.assemblies)

    def test_has_ribbons(self):
        self.assertIn("<path", self.svg)

    def test_has_both_labels(self):
        self.assertIn("A", self.svg)
        self.assertIn("B", self.svg)

    def test_fragment_numbers_shown(self):
        for num in (1, 2, 3):
            self.assertIn(f">{num}<", self.svg)

    def test_rect_count(self):
        """Should have rects for fragments (3 per assembly = 6) plus background."""
        count = self.svg.count("<rect")
        # At least 6 fragment rects + 1 background
        self.assertGreaterEqual(count, 7)


class TestTwoAssembliesDifferentOrder(unittest.TestCase):
    """Two assemblies with rearranged fragments produce crossing ribbons."""

    def setUp(self):
        self.assemblies = [
            _make_assembly("ref", "GS1.0", _make_frags((1, False, 500), (2, False, 500), (3, False, 500))),
            _make_assembly("rearranged", "GS2.0", _make_frags((3, False, 500), (1, False, 500), (2, False, 500))),
        ]
        self.svg = generate_synteny_svg(self.assemblies)

    def test_has_ribbons(self):
        self.assertIn("<path", self.svg)

    def test_ribbon_count(self):
        """Should have 3 ribbons (one per shared fragment)."""
        count = self.svg.count("<path")
        self.assertEqual(count, 3)


class TestReversedFragments(unittest.TestCase):
    """Reversed fragments should show diagonal hatching."""

    def test_hatch_pattern_defined(self):
        assemblies = [
            _make_assembly("A", "GS1.0", _make_frags((1, True, 500))),
        ]
        svg = generate_synteny_svg(assemblies)
        self.assertIn('id="hatch"', svg)

    def test_reversed_fragment_uses_hatch(self):
        assemblies = [
            _make_assembly("A", "GS1.0", _make_frags((1, True, 500), (2, False, 500))),
        ]
        svg = generate_synteny_svg(assemblies)
        self.assertIn('url(#hatch)', svg)

    def test_forward_only_no_hatch_fill(self):
        assemblies = [
            _make_assembly("A", "GS1.0", _make_frags((1, False, 500), (2, False, 500))),
        ]
        svg = generate_synteny_svg(assemblies)
        # Hatch pattern is still defined in defs, but never used as fill
        fill_count = svg.count('fill="url(#hatch)"')
        self.assertEqual(fill_count, 0)


class TestFiveAssemblies(unittest.TestCase):
    """Five or more assemblies should all render correctly."""

    def setUp(self):
        self.assemblies = [
            _make_assembly(f"asm_{i}", f"GS{i}.0", _make_frags(
                (1, False, 400), (2, i % 2 == 0, 300),
                (3, False, 500), (4, False, 200),
            ))
            for i in range(1, 6)
        ]
        self.svg = generate_synteny_svg(self.assemblies)

    def test_all_labels_present(self):
        for i in range(1, 6):
            self.assertIn(f"asm_{i}", self.svg)

    def test_ribbon_groups(self):
        """Should have ribbons between 4 adjacent pairs, 4 fragments each = 16."""
        count = self.svg.count("<path")
        self.assertEqual(count, 16)

    def test_height_auto_calculated(self):
        self.assertIn("viewBox", self.svg)


class TestSortingByGsType(unittest.TestCase):
    """Assemblies should be sorted by GS type."""

    def test_gs_type_order_in_svg(self):
        assemblies = [
            _make_assembly("z_last", "GS3.0", _make_frags((1, False, 100))),
            _make_assembly("a_first", "GS1.0", _make_frags((1, False, 100))),
            _make_assembly("m_mid", "GS2.0", _make_frags((1, False, 100))),
        ]
        svg = generate_synteny_svg(assemblies)
        # GS1.0 should appear before GS2.0 which appears before GS3.0
        idx1 = svg.index("GS1.0")
        idx2 = svg.index("GS2.0")
        idx3 = svg.index("GS3.0")
        self.assertLess(idx1, idx2)
        self.assertLess(idx2, idx3)


class TestQualityBadges(unittest.TestCase):
    """Quality badges should use the correct colours."""

    def test_green_badge(self):
        svg = generate_synteny_svg([
            _make_assembly("g", "GS1.0", _make_frags((1, False, 100)), quality="GREEN"),
        ])
        self.assertIn("#2ca02c", svg)

    def test_amber_badge(self):
        svg = generate_synteny_svg([
            _make_assembly("a", "GS1.0", _make_frags((1, False, 100)), quality="AMBER"),
        ])
        self.assertIn("#ff9900", svg)

    def test_red_badge(self):
        svg = generate_synteny_svg([
            _make_assembly("r", "GS1.0", _make_frags((1, False, 100)), quality="RED"),
        ])
        self.assertIn("#d62728", svg)


class TestDnaADifMarkers(unittest.TestCase):
    """dnaA (Ori) and dif (Ter) markers should render when flagged."""

    def test_ori_marker(self):
        frags = [{"number": 1, "reversed": False, "length": 500, "is_dnaA": True, "is_dif": False}]
        svg = generate_synteny_svg([_make_assembly("x", "GS1.0", frags)])
        self.assertIn("Ori", svg)

    def test_ter_marker(self):
        frags = [{"number": 1, "reversed": False, "length": 500, "is_dnaA": False, "is_dif": True}]
        svg = generate_synteny_svg([_make_assembly("x", "GS1.0", frags)])
        self.assertIn("Ter", svg)

    def test_no_markers_by_default(self):
        frags = [{"number": 1, "reversed": False, "length": 500}]
        svg = generate_synteny_svg([_make_assembly("x", "GS1.0", frags)])
        self.assertNotIn("Ori", svg)
        self.assertNotIn("Ter", svg)


class TestCustomDimensions(unittest.TestCase):
    """Custom width/height/bar_height/gap should be honoured."""

    def test_custom_width(self):
        svg = generate_synteny_svg(
            [_make_assembly("w", "GS1.0", _make_frags((1, False, 100)))],
            width=1200,
        )
        self.assertIn('width="1200"', svg)

    def test_custom_height(self):
        svg = generate_synteny_svg(
            [_make_assembly("h", "GS1.0", _make_frags((1, False, 100)))],
            height=500,
        )
        self.assertIn('height="500"', svg)

    def test_viewbox_matches(self):
        svg = generate_synteny_svg(
            [_make_assembly("v", "GS1.0", _make_frags((1, False, 100)))],
            width=800, height=400,
        )
        self.assertIn('viewBox="0 0 800 400"', svg)


class TestHelperFunctions(unittest.TestCase):
    """Unit tests for internal helpers."""

    def test_fragment_color_wraps(self):
        """Palette should wrap around for large fragment numbers."""
        c1 = _fragment_color(1)
        c11 = _fragment_color(11)
        self.assertEqual(c1, c11)

    def test_ribbon_path_returns_string(self):
        path = _ribbon_path(10, 50, 100, 20, 60, 200)
        self.assertIsInstance(path, str)
        self.assertIn("M", path)
        self.assertIn("C", path)
        self.assertIn("Z", path)

    def test_fragment_color_first(self):
        self.assertEqual(_fragment_color(1), "#a6cee3")


if __name__ == "__main__":
    unittest.main()
