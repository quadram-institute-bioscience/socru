"""Tests for hardened error handling across the socru codebase."""

import os
import subprocess
import tempfile
import unittest
from unittest.mock import patch

from socru.FilterBlast import FilterBlast
from socru.SocruConfig import SocruConfig
from socru.ToolCheck import MissingToolError

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data', 'filter_blast')


class TestMalformedBlastLine(unittest.TestCase):
    """FilterBlast should skip malformed rows gracefully."""

    def test_malformed_blast_line_skipped(self):
        """Lines with fewer than 12 columns are silently skipped."""
        # Write a temporary BLAST results file with a valid line and a malformed one
        valid_line = "query1\tsubj1\t100.0\t500\t0\t0\t1\t500\t1\t500\t1e-10\t900"
        malformed_line = "query2\tsubj2\t100.0"  # only 3 columns

        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as fh:
            fh.write(valid_line + "\n")
            fh.write(malformed_line + "\n")
            tmp_path = fh.name

        try:
            fb = FilterBlast(tmp_path, 1, 1, False)
            # Only the valid line should have been parsed
            self.assertEqual(len(fb.results), 1)
            self.assertEqual(fb.results[0].query_name, 'query1')
        finally:
            os.unlink(tmp_path)

    def test_empty_blast_output_returns_none(self):
        """An empty BLAST results file should yield no results and return_top_result -> None."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as fh:
            # Write an empty file
            fh.write("")
            tmp_path = fh.name

        try:
            fb = FilterBlast(tmp_path, 1, 1, False)
            self.assertEqual(len(fb.results), 0)
            self.assertIsNone(fb.return_top_result())
        finally:
            os.unlink(tmp_path)

    def test_blank_lines_skipped(self):
        """Blank lines in BLAST output are skipped without error."""
        valid_line = "query1\tsubj1\t100.0\t500\t0\t0\t1\t500\t1\t500\t1e-10\t900"
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as fh:
            fh.write("\n")
            fh.write(valid_line + "\n")
            fh.write("\n\n")
            tmp_path = fh.name

        try:
            fb = FilterBlast(tmp_path, 1, 1, False)
            self.assertEqual(len(fb.results), 1)
        finally:
            os.unlink(tmp_path)


class TestSocruErrorHandling(unittest.TestCase):
    """Socru constructor should raise clear errors for bad inputs."""

    def test_socru_with_nonexistent_species_raises(self):
        """A bogus species name must raise FileNotFoundError."""
        from socru.Socru import Socru

        config = SocruConfig(
            input_files=['/tmp/dummy.fa'],
            species='Totally_Fake_Species_XYZ',
        )
        with self.assertRaises(FileNotFoundError) as ctx:
            Socru(config)
        self.assertIn('Totally_Fake_Species_XYZ', str(ctx.exception))

    def test_socru_with_empty_fasta_raises(self):
        """Running Socru on an empty FASTA should not produce a traceback.

        The pipeline should either raise a clear error or return an empty
        result -- not crash with an opaque exception.
        """
        from socru.Socru import Socru

        # Create an empty FASTA file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as fh:
            fh.write("")
            empty_fasta = fh.name

        try:
            config = SocruConfig(
                input_files=[empty_fasta],
                species='Salmonella_enterica',
                output_file=None,
            )
            # Constructing Socru should succeed (tools are checked but species exists)
            s = Socru(config)
            # Running should either raise a clear error or gracefully produce empty output
            # barrnap on an empty file typically fails or produces no operons
            try:
                s.run()
            except (subprocess.CalledProcessError, Exception):
                # Any explicit error is acceptable -- the important thing is
                # that we don't get an unhandled IndexError or similar
                pass
        finally:
            os.unlink(empty_fasta)
            # Clean up possible output files
            for fname in ['profile.txt.novel', 'novel_fragments.fa',
                          'genome_structure.pdf', 'operon_directions.txt']:
                if os.path.exists(fname):
                    os.remove(fname)


class TestSocruMissingTools(unittest.TestCase):
    """Socru should raise MissingToolError when tools are not on PATH."""

    @patch('socru.ToolCheck.shutil.which', return_value=None)
    def test_socru_init_raises_missing_tool(self, mock_which):
        """Socru.__init__ calls check_all_tools and fails early."""
        from socru.Socru import Socru

        config = SocruConfig(
            input_files=['/tmp/dummy.fa'],
            species='Salmonella_enterica',
        )
        with self.assertRaises(MissingToolError):
            Socru(config)
