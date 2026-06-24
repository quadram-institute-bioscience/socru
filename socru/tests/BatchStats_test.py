"""Tests for the BatchStats module."""

import unittest

from socru.BatchStats import BatchStats


def _make_result(
    genome_file="genome.fasta",
    gs_type="GS1.0",
    quality="GREEN",
    confidence_score=95.0,
    genome_length=4800000,
    qc_flags=None,
):
    """Helper to build a minimal AnalysisResult-style dict."""
    return {
        "genome_file": genome_file,
        "gs_type": gs_type,
        "quality": quality,
        "confidence_score": confidence_score,
        "genome_length": genome_length,
        "qc_flags": qc_flags or [],
    }


class TestTypeDistribution(unittest.TestCase):
    def test_counts_sorted_descending(self):
        results = [
            _make_result(gs_type="GS1.0"),
            _make_result(gs_type="GS1.0"),
            _make_result(gs_type="GS1.0"),
            _make_result(gs_type="GS2.1"),
            _make_result(gs_type="GS2.1"),
            _make_result(gs_type="GS3.0"),
        ]
        dist = BatchStats(results).type_distribution()
        self.assertEqual(dist, {"GS1.0": 3, "GS2.1": 2, "GS3.0": 1})
        # Verify ordering
        keys = list(dist.keys())
        self.assertEqual(keys, ["GS1.0", "GS2.1", "GS3.0"])

    def test_single_type(self):
        results = [_make_result(gs_type="GS1.0") for _ in range(5)]
        dist = BatchStats(results).type_distribution()
        self.assertEqual(dist, {"GS1.0": 5})

    def test_empty_results(self):
        dist = BatchStats([]).type_distribution()
        self.assertEqual(dist, {})


class TestQualitySummary(unittest.TestCase):
    def test_all_three_levels(self):
        results = [
            _make_result(quality="GREEN"),
            _make_result(quality="GREEN"),
            _make_result(quality="AMBER"),
            _make_result(quality="RED"),
            _make_result(quality="RED"),
            _make_result(quality="RED"),
        ]
        summary = BatchStats(results).quality_summary()
        self.assertEqual(summary, {"GREEN": 2, "AMBER": 1, "RED": 3})

    def test_all_green(self):
        results = [_make_result(quality="GREEN") for _ in range(4)]
        summary = BatchStats(results).quality_summary()
        self.assertEqual(summary, {"GREEN": 4, "AMBER": 0, "RED": 0})

    def test_empty_results(self):
        summary = BatchStats([]).quality_summary()
        self.assertEqual(summary, {"GREEN": 0, "AMBER": 0, "RED": 0})


class TestMeanConfidence(unittest.TestCase):
    def test_mean_calculation(self):
        results = [
            _make_result(confidence_score=90.0),
            _make_result(confidence_score=80.0),
            _make_result(confidence_score=100.0),
        ]
        mean = BatchStats(results).mean_confidence()
        self.assertAlmostEqual(mean, 90.0)

    def test_single_result(self):
        results = [_make_result(confidence_score=75.5)]
        mean = BatchStats(results).mean_confidence()
        self.assertAlmostEqual(mean, 75.5)

    def test_empty_results(self):
        mean = BatchStats([]).mean_confidence()
        self.assertIsNone(mean)

    def test_none_scores_excluded(self):
        results = [
            _make_result(confidence_score=80.0),
            _make_result(confidence_score=None),
            _make_result(confidence_score=100.0),
        ]
        mean = BatchStats(results).mean_confidence()
        self.assertAlmostEqual(mean, 90.0)


class TestFlagSummary(unittest.TestCase):
    def test_aggregates_flags(self):
        results = [
            _make_result(qc_flags=[
                {"code": "LOW_IDENTITY", "severity": "warning", "message": "Low"},
                {"code": "MISSING_FRAGMENT", "severity": "error", "message": "Missing"},
            ]),
            _make_result(qc_flags=[
                {"code": "LOW_IDENTITY", "severity": "warning", "message": "Low"},
            ]),
            _make_result(qc_flags=[]),
        ]
        flags = BatchStats(results).flag_summary()
        self.assertEqual(flags, {"LOW_IDENTITY": 2, "MISSING_FRAGMENT": 1})

    def test_no_flags(self):
        results = [_make_result(), _make_result()]
        flags = BatchStats(results).flag_summary()
        self.assertEqual(flags, {})

    def test_empty_results(self):
        flags = BatchStats([]).flag_summary()
        self.assertEqual(flags, {})


class TestOutlierAssemblies(unittest.TestCase):
    def test_detects_genome_length_outlier(self):
        results = [
            _make_result(genome_file="normal1.fasta", genome_length=4800000, confidence_score=95.0),
            _make_result(genome_file="normal2.fasta", genome_length=4850000, confidence_score=96.0),
            _make_result(genome_file="normal3.fasta", genome_length=4820000, confidence_score=94.0),
            _make_result(genome_file="normal4.fasta", genome_length=4810000, confidence_score=95.0),
            _make_result(genome_file="normal5.fasta", genome_length=4830000, confidence_score=95.0),
            _make_result(genome_file="outlier.fasta", genome_length=2000000, confidence_score=95.0),
        ]
        outliers = BatchStats(results).outlier_assemblies(z_threshold=2.0)
        self.assertIn("outlier.fasta", outliers)
        self.assertNotIn("normal1.fasta", outliers)

    def test_detects_confidence_outlier(self):
        results = [
            _make_result(genome_file="normal1.fasta", genome_length=4800000, confidence_score=95.0),
            _make_result(genome_file="normal2.fasta", genome_length=4800000, confidence_score=96.0),
            _make_result(genome_file="normal3.fasta", genome_length=4800000, confidence_score=94.0),
            _make_result(genome_file="normal4.fasta", genome_length=4800000, confidence_score=95.0),
            _make_result(genome_file="normal5.fasta", genome_length=4800000, confidence_score=95.0),
            _make_result(genome_file="bad.fasta", genome_length=4800000, confidence_score=10.0),
        ]
        outliers = BatchStats(results).outlier_assemblies(z_threshold=2.0)
        self.assertIn("bad.fasta", outliers)

    def test_no_outliers_in_uniform_data(self):
        results = [
            _make_result(genome_file=f"g{i}.fasta", genome_length=4800000, confidence_score=95.0)
            for i in range(10)
        ]
        outliers = BatchStats(results).outlier_assemblies()
        self.assertEqual(outliers, [])

    def test_single_result_no_outliers(self):
        results = [_make_result()]
        outliers = BatchStats(results).outlier_assemblies()
        self.assertEqual(outliers, [])

    def test_empty_results(self):
        outliers = BatchStats([]).outlier_assemblies()
        self.assertEqual(outliers, [])


if __name__ == "__main__":
    unittest.main()
