"""Tests for cross-module serialization consistency.

Verifies that all dataclasses serialize correctly to JSON and that
SocruConfig fields match the CLI argparse options.
"""

import argparse
import json
import types
import unittest
from dataclasses import fields as dc_fields

from socru.AnalysisResult import AnalysisResult, FragmentResult, OperonResult, QCFlag
from socru.SocruConfig import SocruConfig


class TestFullAnalysisResultSerializes(unittest.TestCase):
    """Verify that a fully-populated AnalysisResult round-trips through JSON."""

    def test_full_analysis_result_serializes(self):
        fragments = [
            FragmentResult(
                number=1, reversed=False, is_dnaA=True, is_dif=False,
                length=50000, coords=[[0, 50000]],
                blast_identity=99.8, blast_alignment_length=49900,
                blast_bit_score=1500.0, blast_e_value=0.0,
                blast_mismatches=2, blast_subject='1',
            ),
            FragmentResult(
                number=2, reversed=True, is_dnaA=False, is_dif=False,
                length=40000, coords=[[50000, 90000]],
                blast_identity=98.5, blast_alignment_length=39500,
                blast_bit_score=1200.0, blast_e_value=0.0,
                blast_mismatches=8, blast_subject='2',
            ),
            FragmentResult(
                number=3, reversed=False, is_dnaA=False, is_dif=True,
                length=60000, coords=[[90000, 150000]],
                blast_identity=99.9, blast_alignment_length=59900,
                blast_bit_score=1800.0, blast_e_value=0.0,
                blast_mismatches=1, blast_subject='3',
            ),
        ]
        operons = [
            OperonResult(start=100, end=5000, direction='forward'),
            OperonResult(start=50100, end=55000, direction='reverse'),
            OperonResult(start=90100, end=95000, direction='forward'),
        ]
        qc_flags = [
            QCFlag(code='LOW_IDENTITY', severity='warning',
                   message='Fragment 2 identity below threshold',
                   details='subject=2'),
            QCFlag(code='NOVEL_ORIENTATION', severity='warning',
                   message='New orientation combination'),
        ]
        novelty = {
            'nearest_known_type': 'GS1.0',
            'edit_distance': 1,
            'assessment': 'likely_biological',
            'confidence': 'high',
            'mean_identity': 99.4,
        }

        result = AnalysisResult(
            genome_file='sample.fasta',
            genome_length=4800000,
            is_circular=True,
            num_operons=3,
            gs_type='GS2.1',
            quality='AMBER',
            is_novel=True,
            fragment_pattern="1 2' 3",
            orientation_binary=2,
            confidence_score=78.5,
            fragments=fragments,
            operons=operons,
            qc_flags=qc_flags,
            validation_passed=True,
            operon_direction_string="Valid\t--> 1 <-- 2' --> 3(Dif)",
            novelty_assessment=novelty,
        )

        # Serialize to JSON
        json_str = result.to_json()
        self.assertIsInstance(json_str, str)

        # Parse back
        parsed = json.loads(json_str)

        # Verify top-level fields
        self.assertEqual(parsed['genome_file'], 'sample.fasta')
        self.assertEqual(parsed['genome_length'], 4800000)
        self.assertTrue(parsed['is_circular'])
        self.assertEqual(parsed['num_operons'], 3)
        self.assertEqual(parsed['gs_type'], 'GS2.1')
        self.assertEqual(parsed['quality'], 'AMBER')
        self.assertTrue(parsed['is_novel'])
        self.assertAlmostEqual(parsed['confidence_score'], 78.5)
        self.assertTrue(parsed['validation_passed'])

        # Verify nested structures
        self.assertEqual(len(parsed['fragments']), 3)
        self.assertEqual(len(parsed['operons']), 3)
        self.assertEqual(len(parsed['qc_flags']), 2)

        # Verify fragment details survived round-trip
        frag0 = parsed['fragments'][0]
        self.assertEqual(frag0['number'], 1)
        self.assertTrue(frag0['is_dnaA'])
        self.assertAlmostEqual(frag0['blast_identity'], 99.8)
        self.assertEqual(frag0['blast_subject'], '1')

        # Verify operons
        self.assertEqual(parsed['operons'][1]['direction'], 'reverse')

        # Verify QC flags
        self.assertEqual(parsed['qc_flags'][0]['code'], 'LOW_IDENTITY')
        self.assertEqual(parsed['qc_flags'][0]['details'], 'subject=2')

        # Verify novelty assessment
        self.assertIsNotNone(parsed['novelty_assessment'])
        self.assertEqual(parsed['novelty_assessment']['nearest_known_type'], 'GS1.0')
        self.assertEqual(parsed['novelty_assessment']['edit_distance'], 1)


class TestAnalysisResultWithNoveltySerializes(unittest.TestCase):
    """Verify that novelty_assessment dict serializes correctly."""

    def test_analysis_result_with_novelty_serializes(self):
        novelty = {
            'nearest_known_type': 'GS1.0',
            'edit_distance': 2,
            'assessment': 'possible_artifact',
            'confidence': 'low',
            'mean_identity': 85.0,
            'details': ['low identity on fragment 2', 'unusual rearrangement'],
        }

        result = AnalysisResult(
            genome_file='novel.fasta',
            genome_length=5000000,
            is_circular=True,
            num_operons=7,
            gs_type='GS0.0',
            quality='RED',
            is_novel=True,
            fragment_pattern="1 2 ? 4 5 6 7",
            orientation_binary=0,
            confidence_score=25.0,
            novelty_assessment=novelty,
        )

        json_str = result.to_json()
        parsed = json.loads(json_str)

        self.assertEqual(parsed['novelty_assessment']['edit_distance'], 2)
        self.assertEqual(parsed['novelty_assessment']['assessment'], 'possible_artifact')
        self.assertIsInstance(parsed['novelty_assessment']['details'], list)
        self.assertEqual(len(parsed['novelty_assessment']['details']), 2)


class TestEmptyAnalysisResultSerializes(unittest.TestCase):
    """Verify that a minimal AnalysisResult with defaults serializes."""

    def test_empty_analysis_result_serializes(self):
        result = AnalysisResult(
            genome_file='empty.fasta',
            genome_length=0,
            is_circular=False,
            num_operons=0,
            gs_type='GS0.0',
            quality='RED',
            is_novel=True,
            fragment_pattern='',
            orientation_binary=0,
            confidence_score=0.0,
        )

        json_str = result.to_json()
        parsed = json.loads(json_str)

        self.assertEqual(parsed['genome_file'], 'empty.fasta')
        self.assertEqual(parsed['fragments'], [])
        self.assertEqual(parsed['operons'], [])
        self.assertEqual(parsed['qc_flags'], [])
        self.assertIsNone(parsed['novelty_assessment'])
        self.assertTrue(parsed['validation_passed'])
        self.assertEqual(parsed['operon_direction_string'], '')


class TestSocruConfigHasAllCliOptions(unittest.TestCase):
    """Verify SocruConfig fields match the argparse options in the CLI script.

    This test constructs the same argparse parser as the CLI and checks that
    every argparse dest has a corresponding SocruConfig field, and that
    from_options() maps them all.
    """

    def _build_parser(self):
        """Replicate the argparse setup from scripts/socru."""
        parser = argparse.ArgumentParser()
        parser.add_argument('species', type=str)
        parser.add_argument('input_files', nargs='+', type=str)
        parser.add_argument('--db_dir', '-d', type=str)
        parser.add_argument('--threads', '-t', type=int, default=1)
        parser.add_argument('--output_file', '-o', type=str)
        parser.add_argument('--output_plot_file', '-p', type=str, default='genome_structure.pdf')
        parser.add_argument('--novel_profiles', '-n', type=str, default='profile.txt.novel')
        parser.add_argument('--new_fragments', '-f', type=str, default='novel_fragments.fa')
        parser.add_argument('--top_blast_hits', '-b', type=str)
        parser.add_argument('--output_html', type=str)
        parser.add_argument('--output_operon_directions_file', '-r', type=str, default='operon_directions.txt')
        parser.add_argument('--output_svg', '-s', type=str)
        parser.add_argument('--output_json', '-j', type=str)
        parser.add_argument('--output_dir', type=str)
        parser.add_argument('--max_bases_from_ends', '-m', type=int)
        parser.add_argument('--not_circular', '-c', action='store_true', default=False)
        parser.add_argument('--min_bit_score', type=int, default=100)
        parser.add_argument('--min_alignment_length', type=int, default=100)
        parser.add_argument('--debug', action='store_true', default=False)
        parser.add_argument('--verbose', '-v', action='store_true', default=False)
        parser.add_argument('--version', action='version', version='test')
        return parser

    def test_socru_config_has_all_cli_options(self):
        """Every CLI option dest (except debug/version) should map to a SocruConfig field."""
        parser = self._build_parser()

        # Collect argparse dest names, excluding meta-options
        skip = {'debug', 'version'}
        cli_dests = set()
        for action in parser._actions:
            if action.dest in skip or action.dest == 'help':
                continue
            cli_dests.add(action.dest)

        # Collect SocruConfig field names
        config_fields = {f.name for f in dc_fields(SocruConfig)}

        # Every CLI dest should be a SocruConfig field
        missing_from_config = cli_dests - config_fields
        self.assertEqual(
            missing_from_config, set(),
            f"CLI options not in SocruConfig: {missing_from_config}"
        )

    def test_from_options_maps_all_cli_dests(self):
        """from_options() should map every CLI dest to the config."""
        parser = self._build_parser()

        # Parse with all options provided
        args = parser.parse_args([
            'Escherichia_coli', 'test.fa',
            '--db_dir', '/db',
            '--threads', '4',
            '--output_file', 'out.tsv',
            '--output_plot_file', 'plot.pdf',
            '--novel_profiles', 'novel.txt',
            '--new_fragments', 'novel.fa',
            '--top_blast_hits', 'blast.txt',
            '--output_html', 'report.html',
            '--output_operon_directions_file', 'dirs.txt',
            '--output_svg', 'genome.svg',
            '--output_json', 'results.json',
            '--output_dir', '/batch',
            '--max_bases_from_ends', '500',
            '--not_circular',
            '--min_bit_score', '200',
            '--min_alignment_length', '150',
            '--verbose',
        ])

        config = SocruConfig.from_options(args)

        # Verify key fields that correspond to output options
        self.assertEqual(config.output_json, 'results.json')
        self.assertEqual(config.output_svg, 'genome.svg')
        self.assertEqual(config.output_html, 'report.html')
        self.assertEqual(config.output_dir, '/batch')
        self.assertEqual(config.output_file, 'out.tsv')
        self.assertEqual(config.db_dir, '/db')
        self.assertEqual(config.threads, 4)
        self.assertEqual(config.species, 'Escherichia_coli')
        self.assertEqual(config.input_files, ['test.fa'])
        self.assertEqual(config.max_bases_from_ends, 500)
        self.assertTrue(config.not_circular)
        self.assertEqual(config.min_bit_score, 200)
        self.assertEqual(config.min_alignment_length, 150)
        self.assertTrue(config.verbose)


if __name__ == '__main__':
    unittest.main()
