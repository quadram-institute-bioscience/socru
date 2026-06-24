"""Tests for the QCFlags module."""

import unittest

from socru.AnalysisResult import AnalysisResult, FragmentResult, OperonResult
from socru.QCFlags import generate_qc_flags


def _make_perfect_result():
    """Create an AnalysisResult with no issues (should produce no flags)."""
    fragments = [
        FragmentResult(
            number=i, reversed=False, is_dnaA=(i == 3), is_dif=(i == 1),
            length=50000, coords=[[0, 50000]],
            blast_identity=100.0, blast_alignment_length=50000,
            blast_bit_score=1000.0, blast_e_value=0.0,
            blast_mismatches=0, blast_subject=str(i),
        )
        for i in range(1, 8)
    ]
    return AnalysisResult(
        genome_file='perfect.fasta',
        genome_length=4800000,
        is_circular=True,
        num_operons=7,
        gs_type='GS1.0',
        quality='GREEN',
        is_novel=False,
        fragment_pattern='1\t2\t3\t4\t5\t6\t7',
        orientation_binary=0,
        confidence_score=100.0,
        fragments=fragments,
        operons=[OperonResult(start=i * 100, end=i * 100 + 5000, direction='forward') for i in range(7)],
        validation_passed=True,
    )


class TestQCFlags(unittest.TestCase):
    """Tests for QC flag generation."""

    def test_no_flags_for_perfect_result(self):
        """A perfect result should produce no QC flags."""
        result = _make_perfect_result()
        flags = generate_qc_flags(result, expected_fragment_count=7)
        self.assertEqual(len(flags), 0)

    def test_low_identity_flag(self):
        """LOW_IDENTITY should trigger when identity is below 95%."""
        result = _make_perfect_result()
        result.fragments[2] = FragmentResult(
            number=3, reversed=False, is_dnaA=True, is_dif=False,
            length=50000, coords=[[0, 50000]],
            blast_identity=90.0, blast_alignment_length=50000,
            blast_bit_score=800.0, blast_e_value=0.0,
            blast_mismatches=50, blast_subject='3',
        )
        flags = generate_qc_flags(result, expected_fragment_count=7)
        flag_codes = [f.code for f in flags]
        self.assertIn('LOW_IDENTITY', flag_codes)

    def test_low_identity_not_triggered_at_95(self):
        """LOW_IDENTITY should NOT trigger when identity is exactly 95%."""
        result = _make_perfect_result()
        result.fragments[0] = FragmentResult(
            number=1, reversed=False, is_dnaA=False, is_dif=True,
            length=50000, coords=[[0, 50000]],
            blast_identity=95.0, blast_alignment_length=50000,
            blast_bit_score=900.0, blast_e_value=0.0,
            blast_mismatches=10, blast_subject='1',
        )
        flags = generate_qc_flags(result, expected_fragment_count=7)
        flag_codes = [f.code for f in flags]
        self.assertNotIn('LOW_IDENTITY', flag_codes)

    def test_missing_fragment_flag(self):
        """MISSING_FRAGMENT should trigger for '?' fragments."""
        result = _make_perfect_result()
        result.fragments[4] = FragmentResult(
            number='?', reversed=False, is_dnaA=False, is_dif=False,
            length=10000, coords=[[0, 10000]],
        )
        flags = generate_qc_flags(result, expected_fragment_count=7)
        flag_codes = [f.code for f in flags]
        self.assertIn('MISSING_FRAGMENT', flag_codes)

    def test_short_alignment_flag(self):
        """SHORT_ALIGNMENT should trigger when coverage is below 80%."""
        result = _make_perfect_result()
        result.fragments[1] = FragmentResult(
            number=2, reversed=False, is_dnaA=False, is_dif=False,
            length=50000, coords=[[0, 50000]],
            blast_identity=100.0, blast_alignment_length=30000,  # 60% coverage
            blast_bit_score=600.0, blast_e_value=0.0,
            blast_mismatches=0, blast_subject='2',
        )
        flags = generate_qc_flags(result, expected_fragment_count=7)
        flag_codes = [f.code for f in flags]
        self.assertIn('SHORT_ALIGNMENT', flag_codes)

    def test_unexpected_operon_count_flag(self):
        """UNEXPECTED_OPERON_COUNT should trigger when count differs."""
        result = _make_perfect_result()
        result.num_operons = 5
        flags = generate_qc_flags(result, expected_fragment_count=7)
        flag_codes = [f.code for f in flags]
        self.assertIn('UNEXPECTED_OPERON_COUNT', flag_codes)

    def test_invalid_arrangement_flag(self):
        """INVALID_ARRANGEMENT should trigger when validation failed."""
        result = _make_perfect_result()
        result.validation_passed = False
        result.quality = 'RED'
        flags = generate_qc_flags(result, expected_fragment_count=7)
        flag_codes = [f.code for f in flags]
        self.assertIn('INVALID_ARRANGEMENT', flag_codes)

    def test_novel_orientation_flag(self):
        """NOVEL_ORIENTATION should trigger for novel AMBER results."""
        result = _make_perfect_result()
        result.is_novel = True
        result.quality = 'AMBER'
        flags = generate_qc_flags(result, expected_fragment_count=7)
        flag_codes = [f.code for f in flags]
        self.assertIn('NOVEL_ORIENTATION', flag_codes)

    def test_novel_pattern_flag(self):
        """NOVEL_PATTERN should trigger for novel RED results."""
        result = _make_perfect_result()
        result.is_novel = True
        result.quality = 'RED'
        result.validation_passed = False
        flags = generate_qc_flags(result, expected_fragment_count=7)
        flag_codes = [f.code for f in flags]
        self.assertIn('NOVEL_PATTERN', flag_codes)

    def test_low_confidence_flag(self):
        """LOW_CONFIDENCE should trigger when score is below 50."""
        result = _make_perfect_result()
        result.confidence_score = 30.0
        flags = generate_qc_flags(result, expected_fragment_count=7)
        flag_codes = [f.code for f in flags]
        self.assertIn('LOW_CONFIDENCE', flag_codes)

    def test_small_chromosome_flag(self):
        """SMALL_CHROMOSOME should trigger when length is below 500,000 bp."""
        result = _make_perfect_result()
        result.genome_length = 300000
        flags = generate_qc_flags(result, expected_fragment_count=7)
        flag_codes = [f.code for f in flags]
        self.assertIn('SMALL_CHROMOSOME', flag_codes)

    def test_small_chromosome_not_triggered_at_500k(self):
        """SMALL_CHROMOSOME should NOT trigger at exactly 500,000 bp."""
        result = _make_perfect_result()
        result.genome_length = 500000
        flags = generate_qc_flags(result, expected_fragment_count=7)
        flag_codes = [f.code for f in flags]
        self.assertNotIn('SMALL_CHROMOSOME', flag_codes)


if __name__ == '__main__':
    unittest.main()
