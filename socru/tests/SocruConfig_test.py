"""Tests for SocruConfig and SocruCreateConfig dataclasses."""

import os
import types
import unittest
from unittest.mock import patch

from socru.SocruConfig import SocruConfig, SocruCreateConfig


class TestSocruConfigDefaults(unittest.TestCase):
    """Verify that SocruConfig has sensible defaults for every field."""

    def test_default_construction(self):
        cfg = SocruConfig()
        self.assertEqual(cfg.input_files, [])
        self.assertIsNone(cfg.db_dir)
        self.assertEqual(cfg.species, "")
        self.assertIsNone(cfg.output_file)
        self.assertIsNone(cfg.output_json)
        self.assertIsNone(cfg.output_svg)
        self.assertIsNone(cfg.output_html)
        self.assertEqual(cfg.output_plot_file, "genome_structure.pdf")
        self.assertEqual(cfg.output_operon_directions_file, "operon_directions.txt")
        self.assertEqual(cfg.novel_profiles, "profile.txt.novel")
        self.assertEqual(cfg.new_fragments, "novel_fragments.fa")
        self.assertIsNone(cfg.top_blast_hits)
        self.assertEqual(cfg.min_bit_score, 100)
        self.assertEqual(cfg.min_alignment_length, 100)
        self.assertIsNone(cfg.max_bases_from_ends)
        self.assertEqual(cfg.threads, 1)
        self.assertFalse(cfg.not_circular)
        self.assertFalse(cfg.verbose)
        self.assertIsNone(cfg.output_dir)

    def test_custom_construction(self):
        cfg = SocruConfig(
            input_files=["a.fa", "b.fa"],
            species="Escherichia_coli",
            threads=4,
            min_bit_score=200,
            not_circular=True,
        )
        self.assertEqual(cfg.input_files, ["a.fa", "b.fa"])
        self.assertEqual(cfg.species, "Escherichia_coli")
        self.assertEqual(cfg.threads, 4)
        self.assertEqual(cfg.min_bit_score, 200)
        self.assertTrue(cfg.not_circular)

    def test_default_factory_independence(self):
        """Each instance should get its own list, not a shared reference."""
        cfg1 = SocruConfig()
        cfg2 = SocruConfig()
        cfg1.input_files.append("x.fa")
        self.assertEqual(cfg2.input_files, [])


class TestSocruConfigFromOptions(unittest.TestCase):
    """Verify from_options() maps an argparse-like namespace correctly."""

    def _make_namespace(self, **overrides):
        """Build a SimpleNamespace mimicking argparse output."""
        defaults = dict(
            input_files=["test.fa"],
            db_dir="/some/db",
            species="Salmonella_enterica",
            output_file="out.tsv",
            output_json="out.json",
            output_svg="out.svg",
            output_html="out.html",
            output_plot_file="plot.pdf",
            output_operon_directions_file="dirs.txt",
            novel_profiles="novel.txt",
            new_fragments="novel.fa",
            top_blast_hits="blast.txt",
            min_bit_score=150,
            min_alignment_length=200,
            max_bases_from_ends=500,
            threads=8,
            not_circular=True,
            verbose=True,
            output_dir="/batch/output",
        )
        defaults.update(overrides)
        return types.SimpleNamespace(**defaults)

    def test_all_fields_mapped(self):
        ns = self._make_namespace()
        cfg = SocruConfig.from_options(ns)

        self.assertEqual(cfg.input_files, ["test.fa"])
        self.assertEqual(cfg.db_dir, "/some/db")
        self.assertEqual(cfg.species, "Salmonella_enterica")
        self.assertEqual(cfg.output_file, "out.tsv")
        self.assertEqual(cfg.output_json, "out.json")
        self.assertEqual(cfg.output_svg, "out.svg")
        self.assertEqual(cfg.output_html, "out.html")
        self.assertEqual(cfg.output_plot_file, "plot.pdf")
        self.assertEqual(cfg.output_operon_directions_file, "dirs.txt")
        self.assertEqual(cfg.novel_profiles, "novel.txt")
        self.assertEqual(cfg.new_fragments, "novel.fa")
        self.assertEqual(cfg.top_blast_hits, "blast.txt")
        self.assertEqual(cfg.min_bit_score, 150)
        self.assertEqual(cfg.min_alignment_length, 200)
        self.assertEqual(cfg.max_bases_from_ends, 500)
        self.assertEqual(cfg.threads, 8)
        self.assertTrue(cfg.not_circular)
        self.assertTrue(cfg.verbose)
        self.assertEqual(cfg.output_dir, "/batch/output")

    def test_missing_optional_attrs_use_defaults(self):
        """from_options should tolerate a namespace that lacks optional attrs."""
        ns = types.SimpleNamespace(
            input_files=["x.fa"],
            species="Salmonella_enterica",
            min_bit_score=100,
            min_alignment_length=100,
            threads=1,
            not_circular=False,
            verbose=False,
            novel_profiles="novel.txt",
            new_fragments="novel.fa",
            max_bases_from_ends=None,
            output_plot_file="plot.pdf",
            output_operon_directions_file="dirs.txt",
        )
        cfg = SocruConfig.from_options(ns)
        self.assertIsNone(cfg.db_dir)
        self.assertIsNone(cfg.output_file)
        self.assertIsNone(cfg.output_json)
        self.assertIsNone(cfg.output_svg)
        self.assertIsNone(cfg.output_html)
        self.assertIsNone(cfg.top_blast_hits)
        self.assertIsNone(cfg.output_dir)


class TestSocruCreateConfigDefaults(unittest.TestCase):
    """Verify that SocruCreateConfig has sensible defaults."""

    def test_default_construction(self):
        cfg = SocruCreateConfig()
        self.assertEqual(cfg.input_file, "")
        self.assertEqual(cfg.output_directory, "")
        self.assertIsNone(cfg.fragment_order)
        self.assertEqual(cfg.threads, 1)
        self.assertIsNone(cfg.dnaa_fasta)
        self.assertIsNone(cfg.dif_fasta)
        self.assertFalse(cfg.verbose)
        self.assertIsNone(cfg.max_bases_from_ends)

    def test_custom_construction(self):
        cfg = SocruCreateConfig(
            input_file="ref.fa",
            output_directory="/out",
            threads=4,
            fragment_order="1-2-3-4-5-6-7",
        )
        self.assertEqual(cfg.input_file, "ref.fa")
        self.assertEqual(cfg.output_directory, "/out")
        self.assertEqual(cfg.threads, 4)
        self.assertEqual(cfg.fragment_order, "1-2-3-4-5-6-7")


class TestSocruCreateConfigFromOptions(unittest.TestCase):
    """Verify from_options() maps an argparse-like namespace correctly."""

    def test_all_fields_mapped(self):
        ns = types.SimpleNamespace(
            input_file="ref.fa",
            output_directory="/db_out",
            fragment_order="1-3-2-4",
            threads=2,
            dnaa_fasta="/path/dnaA.fa",
            dif_fasta="/path/dif.fa",
            verbose=True,
            max_bases_from_ends=300,
        )
        cfg = SocruCreateConfig.from_options(ns)
        self.assertEqual(cfg.input_file, "ref.fa")
        self.assertEqual(cfg.output_directory, "/db_out")
        self.assertEqual(cfg.fragment_order, "1-3-2-4")
        self.assertEqual(cfg.threads, 2)
        self.assertEqual(cfg.dnaa_fasta, "/path/dnaA.fa")
        self.assertEqual(cfg.dif_fasta, "/path/dif.fa")
        self.assertTrue(cfg.verbose)
        self.assertEqual(cfg.max_bases_from_ends, 300)

    def test_missing_optional_attrs_use_defaults(self):
        ns = types.SimpleNamespace(
            input_file="ref.fa",
            output_directory="/db_out",
            threads=1,
            verbose=False,
            max_bases_from_ends=None,
        )
        cfg = SocruCreateConfig.from_options(ns)
        self.assertIsNone(cfg.fragment_order)
        self.assertIsNone(cfg.dnaa_fasta)
        self.assertIsNone(cfg.dif_fasta)


class TestSocruAcceptsSocruConfig(unittest.TestCase):
    """Verify that the Socru class accepts a SocruConfig directly."""

    @patch('socru.ToolCheck.shutil.which', return_value='/usr/bin/fake')
    def test_socru_init_with_config(self, _mock_which):
        """Socru.__init__ should work when given a SocruConfig object."""
        test_modules_dir = os.path.dirname(os.path.realpath(__file__))
        data_dir = os.path.join(test_modules_dir, "data", "socru")
        test_fa = os.path.join(data_dir, "test.fa")

        if not os.path.exists(test_fa):
            self.skipTest("Test data not available")

        from socru.Socru import Socru

        cfg = SocruConfig(
            input_files=[test_fa],
            species="Salmonella_enterica",
            min_bit_score=1000,
            min_alignment_length=1000,
            threads=1,
            output_file="cfg_output_file",
            novel_profiles="cfg_novel",
            new_fragments="cfg_newfrag.fa",
            top_blast_hits="cfg_blast",
            output_plot_file="cfg_plot.png",
            output_operon_directions_file="cfg_directions",
        )

        # Instantiation should succeed without error
        g = Socru(cfg)
        self.assertEqual(g.input_files, [test_fa])
        self.assertEqual(g.min_bit_score, 1000)
        self.assertTrue(g.is_circular)

        # Clean up any output files if they were created
        for fn in [
            "cfg_output_file",
            "cfg_novel",
            "cfg_newfrag.fa",
            "cfg_blast",
            "cfg_plot.png",
            "cfg_directions",
        ]:
            if os.path.exists(fn):
                os.remove(fn)

    @patch('socru.ToolCheck.shutil.which', return_value='/usr/bin/fake')
    def test_socru_init_with_legacy_namespace(self, _mock_which):
        """Socru.__init__ should still work with a plain namespace object."""
        test_modules_dir = os.path.dirname(os.path.realpath(__file__))
        data_dir = os.path.join(test_modules_dir, "data", "socru")
        test_fa = os.path.join(data_dir, "test.fa")

        if not os.path.exists(test_fa):
            self.skipTest("Test data not available")

        from socru.Socru import Socru

        ns = types.SimpleNamespace(
            input_files=[test_fa],
            species="Salmonella_enterica",
            db_dir=None,
            min_bit_score=1000,
            min_alignment_length=1000,
            threads=1,
            output_file="ns_output_file",
            novel_profiles="ns_novel",
            new_fragments="ns_newfrag.fa",
            max_bases_from_ends=None,
            top_blast_hits="ns_blast",
            output_plot_file="ns_plot.png",
            output_operon_directions_file="ns_directions",
            not_circular=False,
            verbose=False,
        )

        g = Socru(ns)
        self.assertEqual(g.input_files, [test_fa])
        self.assertTrue(g.is_circular)

        for fn in [
            "ns_output_file",
            "ns_novel",
            "ns_newfrag.fa",
            "ns_blast",
            "ns_plot.png",
            "ns_directions",
        ]:
            if os.path.exists(fn):
                os.remove(fn)


if __name__ == "__main__":
    unittest.main()
