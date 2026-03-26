"""Integration tests for output format modules using real AnalysisResult data.

These tests verify that all output modules (JSON, HTML, SVG generators, batch
statistics) work correctly with realistic AnalysisResult data, exercising the
serialization and rendering paths end-to-end.
"""

import json
import os
import tempfile
import shutil
import unittest

from socru.AnalysisResult import AnalysisResult, FragmentResult, OperonResult, QCFlag
from socru.BatchStats import BatchStats
from socru.HtmlReport import HtmlReport
from socru.SvgConfidenceHeatmap import generate_confidence_heatmap_svg
from socru.SvgFragmentQuality import generate_fragment_quality_svg
from socru.SvgGenomePlot import generate_genome_svg
from socru.SvgSynteny import generate_synteny_svg
from socru.SvgTypeDistribution import generate_type_distribution_svg


def _make_fragment(number, reversed_frag=False, is_dnaA=False, is_dif=False,
                   length=100000, identity=99.5, bit_score=5000.0,
                   alignment_length=95000):
    """Helper to create a realistic FragmentResult."""
    return FragmentResult(
        number=number,
        reversed=reversed_frag,
        is_dnaA=is_dnaA,
        is_dif=is_dif,
        length=length,
        coords=[(0, length)],
        blast_identity=identity,
        blast_alignment_length=alignment_length,
        blast_bit_score=bit_score,
        blast_e_value=0.0,
        blast_mismatches=int(alignment_length * (100 - identity) / 100) if identity else 0,
        blast_subject=str(number),
    )


def _make_operon(start, end, direction='forward'):
    """Helper to create an OperonResult."""
    return OperonResult(start=start, end=end, direction=direction)


def _make_analysis_result(genome_file='genome1.fa', gs_type='GS1.0',
                          quality='GREEN', is_novel=False,
                          confidence_score=95.0, genome_length=4800000,
                          num_fragments=7):
    """Build a realistic AnalysisResult."""
    fragments = []
    frag_length = genome_length // (num_fragments + 1)
    for i in range(1, num_fragments + 1):
        fragments.append(_make_fragment(
            number=i,
            reversed_frag=(i % 3 == 0),
            is_dnaA=(i == 1),
            is_dif=(i == 4),
            length=frag_length,
            identity=99.5 - (i * 0.1),
            bit_score=5000.0 - (i * 50),
            alignment_length=frag_length - 1000,
        ))

    operons = []
    for i in range(num_fragments + 1):
        start = i * (frag_length + 5000)
        operons.append(_make_operon(
            start=start,
            end=start + 5000,
            direction='forward' if i % 2 == 0 else 'reverse',
        ))

    fragment_pattern = ' '.join(
        str(f.number) + ("'" if f.reversed else '') for f in fragments
    )

    qc_flags = []
    if quality != 'GREEN':
        qc_flags.append(QCFlag(
            code='LOW_IDENTITY',
            severity='warning',
            message='One or more fragments have low BLAST identity',
        ))

    result = AnalysisResult(
        genome_file=genome_file,
        genome_length=genome_length,
        is_circular=True,
        num_operons=len(operons),
        gs_type=gs_type,
        quality=quality,
        is_novel=is_novel,
        fragment_pattern=fragment_pattern,
        orientation_binary=0b1010101,
        confidence_score=confidence_score,
        fragments=fragments,
        operons=operons,
        qc_flags=qc_flags,
        validation_passed=True,
        operon_direction_string='Valid\t+ - + - + - + -',
    )

    if is_novel:
        result.novelty_assessment = {
            'nearest_known_type': 'GS1.0',
            'edit_distance': 2,
            'is_likely_real': True,
            'assessment': 'likely_real',
            'reasoning': 'Single inversion compared to GS1.0',
            'fragment_differences': ['Fragment 3: inverted'],
        }

    return result


def _make_assembly_dict(name, gs_type, quality, num_fragments=7):
    """Build assembly dict in the format expected by SVG generators."""
    fragments = []
    for i in range(1, num_fragments + 1):
        fragments.append({
            'number': i,
            'reversed': (i % 3 == 0),
            'length': 500000 + i * 10000,
            'blast_identity': 99.5 - (i * 0.2),
            'blast_alignment_length': 490000,
            'is_dnaA': (i == 1),
            'is_dif': (i == 4),
        })
    return {
        'name': name,
        'gs_type': gs_type,
        'quality': quality,
        'fragments': fragments,
    }


class TestAnalysisResultJsonRoundtrip(unittest.TestCase):
    """Test AnalysisResult serialization and deserialization."""

    def test_to_json_roundtrip(self):
        """to_json -> json.loads -> verify structure matches."""
        result = _make_analysis_result()
        json_str = result.to_json()
        parsed = json.loads(json_str)

        self.assertEqual(parsed['gs_type'], 'GS1.0')
        self.assertEqual(parsed['quality'], 'GREEN')
        self.assertEqual(parsed['genome_length'], 4800000)
        self.assertIsInstance(parsed['confidence_score'], (int, float))
        self.assertIsInstance(parsed['fragments'], list)
        self.assertEqual(len(parsed['fragments']), 7)
        self.assertIsInstance(parsed['operons'], list)
        self.assertTrue(len(parsed['operons']) > 0)

    def test_to_json_contains_all_keys(self):
        """Verify JSON contains all expected keys."""
        result = _make_analysis_result()
        parsed = json.loads(result.to_json())
        expected_keys = {
            'genome_file', 'genome_length', 'is_circular', 'num_operons',
            'gs_type', 'quality', 'is_novel', 'fragment_pattern',
            'orientation_binary', 'confidence_score', 'fragments', 'operons',
            'qc_flags', 'validation_passed', 'operon_direction_string',
            'novelty_assessment',
        }
        self.assertTrue(expected_keys.issubset(set(parsed.keys())),
                        "Missing keys: {}".format(expected_keys - set(parsed.keys())))

    def test_to_json_fragment_detail(self):
        """Verify fragment data is correctly serialized."""
        result = _make_analysis_result()
        parsed = json.loads(result.to_json())
        frag = parsed['fragments'][0]
        self.assertEqual(frag['number'], 1)
        self.assertTrue(frag['is_dnaA'])
        self.assertFalse(frag['reversed'])
        self.assertIsNotNone(frag['blast_identity'])
        self.assertIsNotNone(frag['blast_bit_score'])

    def test_to_dict_equals_json_loads(self):
        """to_dict() and json.loads(to_json()) should have the same keys and scalar values.

        Note: tuples in to_dict() (e.g. coords) become lists in JSON, so we
        compare the JSON-roundtripped versions for equality.
        """
        result = _make_analysis_result()
        from_dict = result.to_dict()
        from_json = json.loads(result.to_json())
        # Round-trip to_dict through JSON to normalize tuples to lists
        normalized_dict = json.loads(json.dumps(from_dict))
        self.assertEqual(normalized_dict, from_json)

    def test_novel_result_has_novelty_assessment(self):
        """A novel result should have novelty_assessment in JSON."""
        result = _make_analysis_result(is_novel=True, quality='AMBER',
                                       confidence_score=60.0)
        parsed = json.loads(result.to_json())
        self.assertTrue(parsed['is_novel'])
        self.assertIsNotNone(parsed['novelty_assessment'])
        self.assertIn('assessment', parsed['novelty_assessment'])


class TestHtmlReportFromAnalysisResults(unittest.TestCase):
    """Test HTML report generation from real AnalysisResult data."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def test_html_report_from_single_result(self):
        """Feed one real result to HtmlReport, verify output."""
        result = _make_analysis_result()
        results_data = [result.to_dict()]
        report = HtmlReport(results_data, species='Salmonella_enterica')
        html = report.generate()

        self.assertIn('<!DOCTYPE html>', html)
        self.assertIn('Salmonella_enterica', html)
        self.assertIn('GS1.0', html)
        self.assertIn('genome1.fa', html)

    def test_html_report_from_multiple_results(self):
        """Feed multiple real results to HtmlReport, verify output."""
        results = [
            _make_analysis_result(genome_file='genome1.fa', gs_type='GS1.0', quality='GREEN'),
            _make_analysis_result(genome_file='genome2.fa', gs_type='GS2.0', quality='AMBER',
                                  confidence_score=70.0),
            _make_analysis_result(genome_file='genome3.fa', gs_type='GS1.0', quality='RED',
                                  confidence_score=30.0, is_novel=True),
        ]
        results_data = [r.to_dict() for r in results]
        report = HtmlReport(results_data, species='Salmonella_enterica')
        html = report.generate()

        self.assertIn('genome1.fa', html)
        self.assertIn('genome2.fa', html)
        self.assertIn('genome3.fa', html)
        self.assertIn('GS1.0', html)
        self.assertIn('GS2.0', html)

    def test_html_report_save_to_file(self):
        """Verify HtmlReport.save() writes a valid file."""
        result = _make_analysis_result()
        results_data = [result.to_dict()]
        report = HtmlReport(results_data, species='Test_species')
        path = os.path.join(self.tmpdir, 'report.html')
        report.save(path)

        self.assertTrue(os.path.exists(path))
        with open(path) as fh:
            content = fh.read()
        self.assertIn('<!DOCTYPE html>', content)
        self.assertIn('</html>', content)


class TestSvgGenomePlotFromAnalysisResult(unittest.TestCase):
    """Test SVG genome plot generation from AnalysisResult fragment data."""

    def test_svg_genome_plot_basic(self):
        """Feed real fragment data to SVG generator, verify output."""
        result = _make_analysis_result()
        fragments = [f.__dict__ for f in result.fragments]
        operons = [o.__dict__ for o in result.operons]
        svg = generate_genome_svg(
            fragments=fragments,
            operons=operons,
            genome_length=result.genome_length,
            gs_type=result.gs_type,
            quality=result.quality,
            genome_name=result.genome_file,
        )
        self.assertIn('<svg', svg)
        self.assertIn('</svg>', svg)
        self.assertIn('GS1.0', svg)

    def test_svg_genome_plot_contains_fragments(self):
        """Verify the SVG contains fragment number labels."""
        result = _make_analysis_result(num_fragments=5)
        fragments = [f.__dict__ for f in result.fragments]
        operons = [o.__dict__ for o in result.operons]
        svg = generate_genome_svg(
            fragments=fragments,
            operons=operons,
            genome_length=result.genome_length,
            gs_type=result.gs_type,
            quality=result.quality,
        )
        # Each fragment number should appear somewhere in the SVG
        for i in range(1, 6):
            self.assertIn('>{}<'.format(i), svg,
                          "SVG should contain fragment label {}".format(i))

    def test_svg_genome_plot_with_quality_badge(self):
        """Verify quality badge color is present."""
        result = _make_analysis_result(quality='AMBER', confidence_score=60.0)
        fragments = [f.__dict__ for f in result.fragments]
        operons = [o.__dict__ for o in result.operons]
        svg = generate_genome_svg(
            fragments=fragments,
            operons=operons,
            genome_length=result.genome_length,
            gs_type=result.gs_type,
            quality=result.quality,
        )
        # AMBER quality color
        self.assertIn('#ff9900', svg)


class TestSyntenyFromMultipleResults(unittest.TestCase):
    """Test synteny SVG generation from batch results."""

    def test_synteny_from_two_assemblies(self):
        """Generate synteny SVG from 2 assembly results."""
        assemblies = [
            _make_assembly_dict('genome1.fa', 'GS1.0', 'GREEN'),
            _make_assembly_dict('genome2.fa', 'GS2.0', 'AMBER'),
        ]
        svg = generate_synteny_svg(assemblies)
        self.assertIn('<svg', svg)
        self.assertIn('</svg>', svg)
        self.assertIn('GS1.0', svg)
        self.assertIn('GS2.0', svg)

    def test_synteny_from_many_assemblies(self):
        """Generate synteny SVG from 5 assemblies."""
        assemblies = [
            _make_assembly_dict('g{}.fa'.format(i), 'GS{}.0'.format(i), 'GREEN')
            for i in range(1, 6)
        ]
        svg = generate_synteny_svg(assemblies)
        self.assertIn('<svg', svg)
        # All assemblies should be represented
        for i in range(1, 6):
            self.assertIn('GS{}.0'.format(i), svg)

    def test_synteny_empty_assemblies(self):
        """generate_synteny_svg should handle empty list gracefully."""
        svg = generate_synteny_svg([])
        self.assertIn('<svg', svg)
        self.assertIn('</svg>', svg)


class TestBatchStatsFromAnalysisResults(unittest.TestCase):
    """Test batch statistics computation from real results."""

    def _make_batch_results(self):
        """Build a batch of realistic analysis results."""
        results = [
            _make_analysis_result(genome_file='g1.fa', gs_type='GS1.0',
                                  quality='GREEN', confidence_score=95.0),
            _make_analysis_result(genome_file='g2.fa', gs_type='GS1.0',
                                  quality='GREEN', confidence_score=92.0),
            _make_analysis_result(genome_file='g3.fa', gs_type='GS2.0',
                                  quality='AMBER', confidence_score=70.0),
            _make_analysis_result(genome_file='g4.fa', gs_type='GS1.0',
                                  quality='GREEN', confidence_score=98.0),
            _make_analysis_result(genome_file='g5.fa', gs_type='GS3.0',
                                  quality='RED', confidence_score=25.0,
                                  is_novel=True),
        ]
        return [r.to_dict() for r in results]

    def test_type_distribution(self):
        """Verify type_distribution returns correct counts."""
        result_dicts = self._make_batch_results()
        stats = BatchStats(result_dicts)
        dist = stats.type_distribution()
        self.assertEqual(dist['GS1.0'], 3)
        self.assertEqual(dist['GS2.0'], 1)
        self.assertEqual(dist['GS3.0'], 1)

    def test_quality_summary(self):
        """Verify quality_summary returns correct counts."""
        result_dicts = self._make_batch_results()
        stats = BatchStats(result_dicts)
        summary = stats.quality_summary()
        self.assertEqual(summary['GREEN'], 3)
        self.assertEqual(summary['AMBER'], 1)
        self.assertEqual(summary['RED'], 1)

    def test_mean_confidence(self):
        """Verify mean_confidence is calculated correctly."""
        result_dicts = self._make_batch_results()
        stats = BatchStats(result_dicts)
        mean = stats.mean_confidence()
        self.assertIsNotNone(mean)
        expected = (95.0 + 92.0 + 70.0 + 98.0 + 25.0) / 5
        self.assertAlmostEqual(mean, expected, places=1)

    def test_flag_summary(self):
        """Verify flag_summary aggregates QC flags."""
        result_dicts = self._make_batch_results()
        stats = BatchStats(result_dicts)
        flags = stats.flag_summary()
        # AMBER and RED results should have LOW_IDENTITY flags
        self.assertIn('LOW_IDENTITY', flags)
        self.assertEqual(flags['LOW_IDENTITY'], 2)  # AMBER + RED

    def test_outlier_assemblies(self):
        """Verify outlier detection works."""
        result_dicts = self._make_batch_results()
        stats = BatchStats(result_dicts)
        outliers = stats.outlier_assemblies()
        # All genomes have same length, so outliers based on confidence
        self.assertIsInstance(outliers, list)


class TestFragmentQualitySvg(unittest.TestCase):
    """Test fragment quality SVG generation from real result data."""

    def test_fragment_quality_basic(self):
        """Generate fragment quality bars from real data."""
        result = _make_analysis_result()
        fragments = [
            {
                'number': f.number,
                'reversed': f.reversed,
                'blast_identity': f.blast_identity,
                'blast_alignment_length': f.blast_alignment_length,
                'length': f.length,
            }
            for f in result.fragments
        ]
        svg = generate_fragment_quality_svg(fragments, genome_name='test_genome.fa')
        self.assertIn('<svg', svg)
        self.assertIn('</svg>', svg)
        self.assertIn('test_genome.fa', svg)

    def test_fragment_quality_with_missing_identity(self):
        """Handle fragments with no BLAST identity gracefully."""
        fragments = [
            {'number': 1, 'reversed': False, 'blast_identity': 99.8, 'length': 100000},
            {'number': 2, 'reversed': False, 'blast_identity': None, 'length': 80000},
            {'number': 3, 'reversed': True, 'blast_identity': 95.0, 'length': 120000},
        ]
        svg = generate_fragment_quality_svg(fragments)
        self.assertIn('<svg', svg)
        self.assertIn('N/A', svg)  # Missing identity should show N/A


class TestConfidenceHeatmapSvg(unittest.TestCase):
    """Test confidence heatmap SVG generation from batch results."""

    def test_heatmap_from_multiple_assemblies(self):
        """Generate heatmap from multiple assembly dicts."""
        assemblies = [
            _make_assembly_dict('genome1.fa', 'GS1.0', 'GREEN'),
            _make_assembly_dict('genome2.fa', 'GS2.0', 'AMBER'),
            _make_assembly_dict('genome3.fa', 'GS1.0', 'RED'),
        ]
        svg = generate_confidence_heatmap_svg(assemblies)
        self.assertIn('<svg', svg)
        self.assertIn('</svg>', svg)
        self.assertIn('Heatmap', svg)

    def test_heatmap_empty_assemblies(self):
        """Handle empty assembly list gracefully."""
        svg = generate_confidence_heatmap_svg([])
        self.assertIn('<svg', svg)
        self.assertIn('No assemblies', svg)

    def test_heatmap_single_assembly(self):
        """Handle single assembly."""
        assemblies = [_make_assembly_dict('genome1.fa', 'GS1.0', 'GREEN')]
        svg = generate_confidence_heatmap_svg(assemblies)
        self.assertIn('<svg', svg)
        self.assertIn('</svg>', svg)


class TestTypeDistributionSvg(unittest.TestCase):
    """Test type distribution bar chart SVG."""

    def test_type_distribution_basic(self):
        """Generate type distribution from batch stats."""
        results = [
            _make_analysis_result(gs_type='GS1.0').to_dict(),
            _make_analysis_result(gs_type='GS1.0').to_dict(),
            _make_analysis_result(gs_type='GS2.0').to_dict(),
        ]
        stats = BatchStats(results)
        type_dist = stats.type_distribution()
        svg = generate_type_distribution_svg(type_dist)
        self.assertIn('<svg', svg)
        self.assertIn('</svg>', svg)
        self.assertIn('GS1.0', svg)
        self.assertIn('GS2.0', svg)

    def test_type_distribution_empty(self):
        """Handle empty type counts gracefully."""
        svg = generate_type_distribution_svg({})
        self.assertIn('<svg', svg)
        self.assertIn('</svg>', svg)


if __name__ == '__main__':
    unittest.main()
