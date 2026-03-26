"""Golden-file / snapshot tests for socru visualization and report outputs.

These tests ensure output stability for SVG, JSON, and HTML formats. If the
output format changes unexpectedly, the test catches the regression. On first
run, snapshot files are created automatically; subsequent runs compare against
the saved snapshots.
"""

import json
import os
import unittest
from datetime import datetime
from unittest.mock import patch

SNAPSHOT_DIR = os.path.join(os.path.dirname(__file__), 'data', 'snapshots')


class SnapshotTestBase(unittest.TestCase):
    """Base class for snapshot/golden-file tests."""

    def assert_snapshot(self, content, snapshot_name):
        """Compare content against a saved snapshot. Create snapshot if missing."""
        snapshot_path = os.path.join(SNAPSHOT_DIR, snapshot_name)
        if not os.path.exists(snapshot_path):
            os.makedirs(os.path.dirname(snapshot_path), exist_ok=True)
            with open(snapshot_path, 'w') as f:
                f.write(content)
            self.skipTest(f"Snapshot created: {snapshot_name}. Re-run to verify.")
        with open(snapshot_path) as f:
            expected = f.read()
        self.assertEqual(content, expected, f"Output differs from snapshot {snapshot_name}")


class TestSvgGenomePlotSnapshot(SnapshotTestBase):
    """Snapshot test for the circular genome SVG diagram."""

    def test_genome_plot_snapshot(self):
        from socru.SvgGenomePlot import generate_genome_svg

        fragments = [
            {'number': 1, 'reversed': False, 'length': 1500000,
             'coords': [(0, 1500000)], 'is_dnaA': True, 'is_dif': False},
            {'number': 2, 'reversed': True, 'length': 1200000,
             'coords': [(1505000, 2705000)], 'is_dnaA': False, 'is_dif': False},
            {'number': 3, 'reversed': False, 'length': 900000,
             'coords': [(2710000, 3610000)], 'is_dnaA': False, 'is_dif': True},
        ]
        operons = [
            {'start': 1500000, 'end': 1505000, 'direction': 'forward'},
            {'start': 2705000, 'end': 2710000, 'direction': 'reverse'},
        ]
        svg = generate_genome_svg(
            fragments=fragments,
            operons=operons,
            genome_length=3700000,
            gs_type="GS1.1",
            quality="GREEN",
            genome_name="Test_genome",
        )
        self.assertIn('<svg', svg)
        self.assertIn('</svg>', svg)
        self.assertIn('GS1.1', svg)
        self.assert_snapshot(svg, 'genome_plot.svg')


class TestSvgSyntenySnapshot(SnapshotTestBase):
    """Snapshot test for the linear synteny comparison SVG."""

    def test_synteny_snapshot(self):
        from socru.SvgSynteny import generate_synteny_svg

        assemblies = [
            {
                'name': 'ref', 'gs_type': 'GS1.0', 'quality': 'GREEN',
                'fragments': [
                    {'number': 1, 'reversed': False, 'length': 1000000},
                    {'number': 2, 'reversed': False, 'length': 800000},
                    {'number': 3, 'reversed': False, 'length': 600000},
                ],
            },
            {
                'name': 'query', 'gs_type': 'GS2.1', 'quality': 'AMBER',
                'fragments': [
                    {'number': 1, 'reversed': True, 'length': 1000000},
                    {'number': 3, 'reversed': False, 'length': 600000},
                    {'number': 2, 'reversed': False, 'length': 800000},
                ],
            },
        ]
        svg = generate_synteny_svg(assemblies=assemblies)
        self.assertIn('<svg', svg)
        self.assertIn('</svg>', svg)
        self.assert_snapshot(svg, 'synteny.svg')


class TestSvgFragmentQualitySnapshot(SnapshotTestBase):
    """Snapshot test for the fragment quality bar chart SVG."""

    def test_fragment_quality_snapshot(self):
        from socru.SvgFragmentQuality import generate_fragment_quality_svg

        fragments = [
            {'number': 1, 'reversed': False, 'blast_identity': 99.8,
             'blast_alignment_length': 50000, 'length': 50000},
            {'number': 2, 'reversed': True, 'blast_identity': 96.2,
             'blast_alignment_length': 40000, 'length': 45000},
            {'number': 3, 'reversed': False, 'blast_identity': 88.5,
             'blast_alignment_length': 30000, 'length': 35000},
        ]
        svg = generate_fragment_quality_svg(fragments=fragments, genome_name="Test")
        self.assertIn('<svg', svg)
        self.assertIn('Fragment Quality', svg)
        self.assert_snapshot(svg, 'fragment_quality.svg')


class TestSvgTypeDistributionSnapshot(SnapshotTestBase):
    """Snapshot test for the GS type distribution bar chart SVG."""

    def test_type_distribution_snapshot(self):
        from socru.SvgTypeDistribution import generate_type_distribution_svg

        svg = generate_type_distribution_svg({'GS1.0': 15, 'GS2.1': 8, 'GS1.3': 3})
        self.assertIn('<svg', svg)
        self.assertIn('GS Type Distribution', svg)
        self.assert_snapshot(svg, 'type_distribution.svg')


class TestSvgConfidenceHeatmapSnapshot(SnapshotTestBase):
    """Snapshot test for the confidence heatmap SVG."""

    def test_confidence_heatmap_snapshot(self):
        from socru.SvgConfidenceHeatmap import generate_confidence_heatmap_svg

        assemblies = [
            {
                'name': f'sample_{i}',
                'gs_type': 'GS1.0',
                'quality': 'GREEN',
                'fragments': [
                    {'number': j, 'blast_identity': 99.0 - i - j * 0.5}
                    for j in range(1, 4)
                ],
            }
            for i in range(3)
        ]
        svg = generate_confidence_heatmap_svg(assemblies=assemblies)
        self.assertIn('<svg', svg)
        self.assertIn('Fragment Identity Heatmap', svg)
        self.assert_snapshot(svg, 'confidence_heatmap.svg')


class TestJsonSnapshots(SnapshotTestBase):
    """Snapshot test for JSON serialization of AnalysisResult."""

    def test_analysis_result_json_structure(self):
        from socru.AnalysisResult import AnalysisResult, FragmentResult, OperonResult

        result = AnalysisResult(
            genome_file='test.fasta',
            genome_length=4800000,
            is_circular=True,
            num_operons=7,
            gs_type='GS1.0',
            quality='GREEN',
            is_novel=False,
            fragment_pattern='1 2 3 4 5 6 7',
            orientation_binary=0,
            confidence_score=97.5,
            validation_passed=True,
            fragments=[
                FragmentResult(
                    number=i,
                    reversed=False,
                    is_dnaA=(i == 1),
                    is_dif=(i == 4),
                    length=500000 + i * 50000,
                    coords=[(i * 600000, i * 600000 + 500000)],
                )
                for i in range(1, 8)
            ],
            operons=[
                OperonResult(
                    start=i * 600000 - 5000,
                    end=i * 600000,
                    direction='forward',
                )
                for i in range(1, 7)
            ],
            qc_flags=[],
        )
        j = result.to_json(indent=2)

        # Verify structure, not just snapshot equality
        parsed = json.loads(j)
        self.assertEqual(parsed['gs_type'], 'GS1.0')
        self.assertEqual(parsed['genome_length'], 4800000)
        self.assertTrue(parsed['is_circular'])
        self.assertEqual(len(parsed['fragments']), 7)
        self.assertEqual(len(parsed['operons']), 6)
        self.assertEqual(parsed['confidence_score'], 97.5)
        self.assertTrue(parsed['validation_passed'])
        self.assertEqual(parsed['fragments'][0]['number'], 1)
        self.assertTrue(parsed['fragments'][0]['is_dnaA'])
        self.assertFalse(parsed['fragments'][0]['is_dif'])

        # Snapshot the full JSON output
        self.assert_snapshot(j, 'analysis_result.json')


class TestHtmlSnapshots(SnapshotTestBase):
    """Snapshot test for the HTML report generator."""

    @patch('socru.HtmlReport.datetime')
    def test_html_report_structure(self, mock_datetime):
        # Pin datetime so the report is deterministic across runs
        fixed_dt = datetime(2025, 1, 15, 10, 30, 0)
        mock_datetime.now.return_value = fixed_dt
        mock_datetime.side_effect = lambda *a, **kw: datetime(*a, **kw)

        from socru.HtmlReport import HtmlReport

        results = [
            {
                'genome_file': 'test.fasta',
                'gs_type': 'GS1.0',
                'quality': 'GREEN',
                'confidence_score': 97.5,
                'fragment_pattern': '1 2 3',
                'num_operons': 3,
                'is_novel': False,
                'qc_flags': [],
                'fragments': [
                    {
                        'number': 1,
                        'reversed': False,
                        'blast_identity': 99.5,
                        'blast_alignment_length': 50000,
                        'blast_bit_score': 90000,
                        'length': 50000,
                    },
                ],
                'genome_length': 3000000,
                'validation_passed': True,
            },
        ]
        report = HtmlReport(results, species='Test_species')
        html = report.generate()

        # Basic structure checks
        self.assertIn('<!DOCTYPE html>', html)
        self.assertIn('Socru Genome Structure Analysis Report', html)
        self.assertIn('Test_species', html)
        self.assertIn('GS1.0', html)
        self.assertIn('test.fasta', html)

        self.assert_snapshot(html, 'report.html')


if __name__ == '__main__':
    unittest.main()
