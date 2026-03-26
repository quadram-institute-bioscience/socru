"""Comprehensive end-to-end integration tests for the Socru analysis pipeline.

These tests exercise the full pipeline: FASTA input -> barrnap -> BLAST -> type
assignment -> structured result output, verifying that all components work
together correctly.
"""

import json
import os
import shutil
import tempfile
import unittest

from socru.Socru import Socru
from socru.SocruConfig import SocruConfig, SocruCreateConfig
from socru.SocruCreate import SocruCreate

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data', 'socru')
create_data_dir = os.path.join(test_modules_dir, 'data', 'create')


class TestEndToEndAnalysis(unittest.TestCase):
    """Full pipeline tests: FASTA in, structured results out."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.test_fasta = os.path.join(data_dir, 'test.fa')

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def _make_config(self, **overrides):
        """Build a SocruConfig with sensible defaults for testing."""
        defaults = dict(
            input_files=[self.test_fasta],
            species='Salmonella_enterica',
            min_bit_score=1000,
            min_alignment_length=1000,
            threads=1,
            output_file=os.path.join(self.tmpdir, 'output.txt'),
            novel_profiles=os.path.join(self.tmpdir, 'novel'),
            new_fragments=os.path.join(self.tmpdir, 'newfrag.fa'),
            top_blast_hits=os.path.join(self.tmpdir, 'blast'),
            output_plot_file=os.path.join(self.tmpdir, 'plot.png'),
            output_operon_directions_file=os.path.join(self.tmpdir, 'directions'),
        )
        defaults.update(overrides)
        return SocruConfig(**defaults)

    def _run_analysis(self, **overrides):
        """Run a full Socru analysis and return the Socru instance."""
        config = self._make_config(**overrides)
        s = Socru(config)
        s.run()
        return s

    # ------------------------------------------------------------------
    # 1. Basic result production
    # ------------------------------------------------------------------

    def test_single_genome_produces_result(self):
        """Run Socru on a test genome and verify it returns a non-empty result string."""
        s = self._run_analysis()
        output_path = os.path.join(self.tmpdir, 'output.txt')
        self.assertTrue(os.path.exists(output_path))
        with open(output_path) as fh:
            content = fh.read().strip()
        self.assertTrue(len(content) > 0, "Output file should not be empty")
        # Should contain the input filename and a GS type
        self.assertIn('test.fa', content)

    # ------------------------------------------------------------------
    # 2. AnalysisResult is populated
    # ------------------------------------------------------------------

    def test_analysis_result_is_populated(self):
        """Run and check AnalysisResult has gs_type, quality, confidence_score, fragments, operons."""
        s = self._run_analysis()
        self.assertEqual(len(s.analysis_results), 1)
        result = s.analysis_results[0]
        self.assertTrue(len(result.gs_type) > 0, "gs_type should be non-empty")
        self.assertIn(result.quality, ('GREEN', 'AMBER', 'RED'))
        self.assertIsInstance(result.confidence_score, (int, float))
        self.assertIsInstance(result.fragments, list)
        self.assertTrue(len(result.fragments) > 0, "Should have at least one fragment")
        self.assertIsInstance(result.operons, list)
        self.assertTrue(len(result.operons) > 0, "Should have at least one operon")

    # ------------------------------------------------------------------
    # 3. JSON output is valid
    # ------------------------------------------------------------------

    def test_json_output_is_valid(self):
        """Use output_json, verify file is valid JSON with expected structure."""
        json_path = os.path.join(self.tmpdir, 'results.json')
        s = self._run_analysis(output_json=json_path)
        self.assertTrue(os.path.exists(json_path))
        with open(json_path) as fh:
            data = json.load(fh)
        self.assertIn('results', data)
        self.assertIsInstance(data['results'], list)
        self.assertTrue(len(data['results']) > 0)
        first = data['results'][0]
        for key in ('gs_type', 'quality', 'confidence_score', 'fragments', 'operons',
                     'genome_file', 'genome_length', 'is_novel', 'qc_flags'):
            self.assertIn(key, first, "JSON result should contain key '{}'".format(key))

    # ------------------------------------------------------------------
    # 4. SVG output is valid
    # ------------------------------------------------------------------

    def test_svg_output_is_valid(self):
        """Use output_svg, verify file contains <svg tag."""
        svg_path = os.path.join(self.tmpdir, 'genome.svg')
        s = self._run_analysis(output_svg=svg_path)
        self.assertTrue(os.path.exists(svg_path))
        with open(svg_path) as fh:
            content = fh.read()
        self.assertIn('<svg', content)
        self.assertIn('</svg>', content)

    # ------------------------------------------------------------------
    # 5. HTML output is valid
    # ------------------------------------------------------------------

    def test_html_output_is_valid(self):
        """Use output_html, verify file contains <!DOCTYPE html>."""
        html_path = os.path.join(self.tmpdir, 'report.html')
        s = self._run_analysis(output_html=html_path)
        self.assertTrue(os.path.exists(html_path))
        with open(html_path) as fh:
            content = fh.read()
        self.assertIn('<!DOCTYPE html>', content)
        self.assertIn('Socru', content)

    # ------------------------------------------------------------------
    # 6. Confidence score range
    # ------------------------------------------------------------------

    def test_confidence_score_range(self):
        """Verify confidence_score is between 0 and 100."""
        s = self._run_analysis()
        result = s.analysis_results[0]
        self.assertGreaterEqual(result.confidence_score, 0)
        self.assertLessEqual(result.confidence_score, 100)

    # ------------------------------------------------------------------
    # 7. QC flags are a list
    # ------------------------------------------------------------------

    def test_qc_flags_are_list(self):
        """Verify qc_flags is a list."""
        s = self._run_analysis()
        result = s.analysis_results[0]
        self.assertIsInstance(result.qc_flags, list)
        # Each flag should have code, severity, message
        for flag in result.qc_flags:
            self.assertTrue(hasattr(flag, 'code'))
            self.assertTrue(hasattr(flag, 'severity'))
            self.assertTrue(hasattr(flag, 'message'))

    # ------------------------------------------------------------------
    # 8. Fragment count matches species
    # ------------------------------------------------------------------

    def test_fragment_count_matches_species(self):
        """For the test genome against Salmonella_enterica, verify fragment count is reasonable."""
        s = self._run_analysis()
        result = s.analysis_results[0]
        # The test.fa genome is a partial/synthetic genome; it should produce
        # at least one fragment and the number should be consistent
        self.assertGreater(len(result.fragments), 0,
                           "Should have at least one fragment")
        # Salmonella enterica profile has 7 columns, but the test genome
        # may only have a subset of operons detected
        self.assertLessEqual(len(result.fragments), 7,
                             "Should not exceed the species fragment count")

    # ------------------------------------------------------------------
    # 9. Quality is a valid value
    # ------------------------------------------------------------------

    def test_quality_is_valid_value(self):
        """Verify quality is one of GREEN, AMBER, RED."""
        s = self._run_analysis()
        result = s.analysis_results[0]
        self.assertIn(result.quality, ('GREEN', 'AMBER', 'RED'))

    # ------------------------------------------------------------------
    # 10. Context manager cleanup
    # ------------------------------------------------------------------

    def test_context_manager_cleanup(self):
        """Verify temp dirs are cleaned up after using Socru as a context manager."""
        config = self._make_config()
        with Socru(config) as s:
            s.run()
            # Capture temp dirs created during analysis
            temp_dirs = list(s.dirs_to_cleanup)
            self.assertTrue(len(temp_dirs) > 0, "Should have temp dirs during analysis")
            # Verify they exist during analysis
            for d in temp_dirs:
                self.assertTrue(os.path.exists(d))

        # After exiting context, temp dirs should be cleaned up
        for d in temp_dirs:
            self.assertFalse(os.path.exists(d),
                             "Temp dir {} should have been cleaned up".format(d))

    # ------------------------------------------------------------------
    # 11. SocruConfig dataclass works directly
    # ------------------------------------------------------------------

    def test_config_dataclass_works(self):
        """Use SocruConfig directly (not argparse) to run analysis."""
        json_path = os.path.join(self.tmpdir, 'config_test.json')
        config = SocruConfig(
            input_files=[self.test_fasta],
            species='Salmonella_enterica',
            min_bit_score=1000,
            min_alignment_length=1000,
            threads=1,
            output_file=os.path.join(self.tmpdir, 'config_output.txt'),
            output_json=json_path,
            novel_profiles=os.path.join(self.tmpdir, 'novel2'),
            new_fragments=os.path.join(self.tmpdir, 'newfrag2.fa'),
            output_plot_file=os.path.join(self.tmpdir, 'plot2.png'),
            output_operon_directions_file=os.path.join(self.tmpdir, 'directions2'),
        )
        with Socru(config) as s:
            s.run()
        self.assertTrue(os.path.exists(json_path))
        with open(json_path) as fh:
            data = json.load(fh)
        self.assertIn('results', data)
        self.assertTrue(len(data['results']) > 0)

    # ------------------------------------------------------------------
    # 12. Novelty assessment on novel result
    # ------------------------------------------------------------------

    def test_novelty_assessment_on_novel(self):
        """If a result is novel, verify novelty_assessment is populated."""
        s = self._run_analysis()
        result = s.analysis_results[0]
        if result.is_novel:
            self.assertIsNotNone(result.novelty_assessment,
                                 "Novel results should have novelty_assessment populated")
            self.assertIsInstance(result.novelty_assessment, dict)
            # Should have key fields from NoveltyDetector
            self.assertIn('assessment', result.novelty_assessment)
            self.assertIn('nearest_known_type', result.novelty_assessment)
            self.assertIn('edit_distance', result.novelty_assessment)
        else:
            # Known types may or may not have novelty_assessment
            # Just verify the field exists
            # (it may be None for known types, which is fine)
            self.assertTrue(hasattr(result, 'novelty_assessment'))

    # ------------------------------------------------------------------
    # 13. Batch JSON includes batch_stats
    # ------------------------------------------------------------------

    def test_batch_json_includes_batch_stats(self):
        """Run on multiple files, verify batch_stats in JSON output."""
        test_fasta_1 = os.path.join(data_dir, 'test.fa')
        test_fasta_2 = os.path.join(data_dir, 'wrapped.fa')
        json_path = os.path.join(self.tmpdir, 'batch.json')
        config = self._make_config(
            input_files=[test_fasta_1, test_fasta_2],
            output_json=json_path,
        )
        s = Socru(config)
        s.run()

        self.assertTrue(os.path.exists(json_path))
        with open(json_path) as fh:
            data = json.load(fh)
        self.assertIn('results', data)

        # batch_stats should be present when there are multiple results
        # (only if both files produce results)
        if len(data['results']) > 1:
            self.assertIn('batch_stats', data)
            batch_stats = data['batch_stats']
            for key in ('type_distribution', 'quality_summary', 'mean_confidence',
                         'total_assemblies'):
                self.assertIn(key, batch_stats,
                              "batch_stats should contain '{}'".format(key))
            self.assertEqual(batch_stats['total_assemblies'], len(data['results']))

    # ------------------------------------------------------------------
    # Additional: wrapped genome still produces valid output
    # ------------------------------------------------------------------

    def test_wrapped_genome_produces_output(self):
        """The wrapped.fa test case should produce valid output or gracefully handle."""
        wrapped_fasta = os.path.join(data_dir, 'wrapped.fa')
        config = self._make_config(input_files=[wrapped_fasta])
        s = Socru(config)
        s.run()
        # Should not raise; output file should exist
        output_path = os.path.join(self.tmpdir, 'output.txt')
        self.assertTrue(os.path.exists(output_path))

    # ------------------------------------------------------------------
    # AnalysisResult.to_dict() and to_json()
    # ------------------------------------------------------------------

    def test_analysis_result_to_dict_complete(self):
        """Verify to_dict() contains all expected keys."""
        s = self._run_analysis()
        result = s.analysis_results[0]
        d = result.to_dict()
        expected_keys = [
            'genome_file', 'genome_length', 'is_circular', 'num_operons',
            'gs_type', 'quality', 'is_novel', 'fragment_pattern',
            'orientation_binary', 'confidence_score', 'fragments', 'operons',
            'qc_flags', 'validation_passed', 'operon_direction_string',
        ]
        for key in expected_keys:
            self.assertIn(key, d, "to_dict() should contain '{}'".format(key))

    def test_analysis_result_genome_length_positive(self):
        """Verify genome_length is a positive integer."""
        s = self._run_analysis()
        result = s.analysis_results[0]
        self.assertIsInstance(result.genome_length, int)
        self.assertGreater(result.genome_length, 0)

    def test_fragment_results_have_blast_data(self):
        """Verify matched fragments have BLAST alignment data."""
        s = self._run_analysis()
        result = s.analysis_results[0]
        matched = [f for f in result.fragments if str(f.number) != '?']
        self.assertTrue(len(matched) > 0, "Should have at least one matched fragment")
        for frag in matched:
            self.assertIsNotNone(frag.blast_identity,
                                 "Matched fragment should have blast_identity")
            self.assertIsNotNone(frag.blast_bit_score,
                                 "Matched fragment should have blast_bit_score")
            self.assertIsNotNone(frag.blast_alignment_length,
                                 "Matched fragment should have blast_alignment_length")

    def test_operon_results_have_coordinates(self):
        """Verify operon results have start, end, and direction."""
        s = self._run_analysis()
        result = s.analysis_results[0]
        for operon in result.operons:
            self.assertIsInstance(operon.start, int)
            self.assertIsInstance(operon.end, int)
            self.assertIn(operon.direction, ('forward', 'reverse'))
            self.assertGreater(operon.end, operon.start)


class TestEndToEndCreate(unittest.TestCase):
    """End-to-end tests for the SocruCreate database creation workflow."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.test_fasta = os.path.join(create_data_dir, 'test.fa')
        self.dnaa_fasta = os.path.join(create_data_dir, 'dnaA.fa.gz')
        self.dif_fasta = os.path.join(create_data_dir, 'dif.fa.gz')

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    # ------------------------------------------------------------------
    # 1. Create database produces fragment files
    # ------------------------------------------------------------------

    def test_create_database_produces_fragments(self):
        """Run socru_create, verify fragment FASTA files created."""
        output_dir = os.path.join(self.tmpdir, 'new_db')
        config = SocruCreateConfig(
            input_file=self.test_fasta,
            output_directory=output_dir,
            threads=1,
            dnaa_fasta=self.dnaa_fasta,
            dif_fasta=self.dif_fasta,
        )
        creator = SocruCreate(config)
        creator.run()

        # Should have at least one numbered fragment file
        fragment_files = [f for f in os.listdir(output_dir) if f.endswith('.fa') and f[0].isdigit()]
        self.assertTrue(len(fragment_files) > 0, "Should create fragment FASTA files")

        # Each fragment file should be non-empty and contain FASTA data
        for ff in fragment_files:
            path = os.path.join(output_dir, ff)
            self.assertTrue(os.path.getsize(path) > 0,
                            "Fragment file {} should not be empty".format(ff))
            with open(path) as fh:
                first_line = fh.readline()
            self.assertTrue(first_line.startswith('>'),
                            "Fragment file {} should start with '>'".format(ff))

    # ------------------------------------------------------------------
    # 2. Create database produces profile
    # ------------------------------------------------------------------

    def test_create_database_produces_profile(self):
        """Verify profile.txt created."""
        output_dir = os.path.join(self.tmpdir, 'profile_db')
        config = SocruCreateConfig(
            input_file=self.test_fasta,
            output_directory=output_dir,
            threads=1,
            dnaa_fasta=self.dnaa_fasta,
            dif_fasta=self.dif_fasta,
        )
        creator = SocruCreate(config)
        creator.run()

        profile_path = os.path.join(output_dir, 'profile.txt')
        self.assertTrue(os.path.exists(profile_path),
                        "profile.txt should be created")
        with open(profile_path) as fh:
            content = fh.read()
        self.assertTrue(len(content.strip()) > 0, "profile.txt should not be empty")
        # Should contain at least a header line and one profile
        lines = content.strip().split('\n')
        self.assertTrue(len(lines) >= 2,
                        "profile.txt should have header + at least one profile")

    # ------------------------------------------------------------------
    # 3. Created database is usable for analysis
    # ------------------------------------------------------------------

    def test_created_database_is_usable(self):
        """Create a DB, then use it to analyze the same genome."""
        output_dir = os.path.join(self.tmpdir, 'usable_db')
        config = SocruCreateConfig(
            input_file=self.test_fasta,
            output_directory=output_dir,
            threads=1,
            dnaa_fasta=self.dnaa_fasta,
            dif_fasta=self.dif_fasta,
        )
        creator = SocruCreate(config)
        creator.run()

        # Now use the created database to analyze the same genome
        results_path = os.path.join(self.tmpdir, 'analysis_output.txt')
        analysis_config = SocruConfig(
            input_files=[self.test_fasta],
            db_dir=output_dir,
            min_bit_score=100,
            min_alignment_length=100,
            threads=1,
            output_file=results_path,
            novel_profiles=os.path.join(self.tmpdir, 'novel'),
            new_fragments=os.path.join(self.tmpdir, 'newfrag.fa'),
            output_plot_file=os.path.join(self.tmpdir, 'plot.png'),
            output_operon_directions_file=os.path.join(self.tmpdir, 'directions'),
        )
        with Socru(analysis_config) as s:
            s.run()

        self.assertTrue(os.path.exists(results_path))
        with open(results_path) as fh:
            content = fh.read().strip()
        self.assertTrue(len(content) > 0, "Analysis output should not be empty")
        # The genome should produce a GS1.0 result since it was the reference
        self.assertIn('GS1.0', content,
                      "Analyzing reference genome against its own DB should yield GS1.0")

    def test_create_with_config_dataclass(self):
        """Verify SocruCreateConfig works for database creation."""
        output_dir = os.path.join(self.tmpdir, 'config_db')
        config = SocruCreateConfig(
            input_file=self.test_fasta,
            output_directory=output_dir,
            threads=1,
            dnaa_fasta=self.dnaa_fasta,
            dif_fasta=self.dif_fasta,
            verbose=False,
        )
        with SocruCreate(config) as creator:
            creator.run()

        self.assertTrue(os.path.exists(os.path.join(output_dir, 'profile.txt')))
        fragment_files = [f for f in os.listdir(output_dir) if f.endswith('.fa') and f[0].isdigit()]
        self.assertTrue(len(fragment_files) >= 2,
                        "Should have at least 2 fragment files")


if __name__ == '__main__':
    unittest.main()
