"""Tests for the AnalysisResult dataclass and its serialization."""

import json
import unittest

from socru.AnalysisResult import AnalysisResult, FragmentResult, OperonResult, QCFlag


class TestFragmentResult(unittest.TestCase):
    """Tests for FragmentResult creation and field access."""

    def test_create_with_blast_data(self):
        fr = FragmentResult(
            number=3,
            reversed=False,
            is_dnaA=True,
            is_dif=False,
            length=50000,
            coords=[[100, 50100]],
            blast_identity=99.5,
            blast_alignment_length=49800,
            blast_bit_score=1200.0,
            blast_e_value=0.0,
            blast_mismatches=5,
            blast_subject='3',
        )
        self.assertEqual(fr.number, 3)
        self.assertTrue(fr.is_dnaA)
        self.assertAlmostEqual(fr.blast_identity, 99.5)

    def test_create_unknown_fragment(self):
        fr = FragmentResult(
            number='?',
            reversed=False,
            is_dnaA=False,
            is_dif=False,
            length=1000,
            coords=[[0, 1000]],
        )
        self.assertEqual(str(fr.number), '?')
        self.assertIsNone(fr.blast_identity)


class TestOperonResult(unittest.TestCase):
    """Tests for OperonResult creation."""

    def test_create(self):
        op = OperonResult(start=100, end=5000, direction='forward')
        self.assertEqual(op.start, 100)
        self.assertEqual(op.direction, 'forward')


class TestQCFlag(unittest.TestCase):
    """Tests for QCFlag creation."""

    def test_create_with_details(self):
        flag = QCFlag(
            code='LOW_IDENTITY',
            severity='warning',
            message='Identity below threshold',
            details='fragment 2',
        )
        self.assertEqual(flag.code, 'LOW_IDENTITY')
        self.assertEqual(flag.severity, 'warning')
        self.assertEqual(flag.details, 'fragment 2')

    def test_create_without_details(self):
        flag = QCFlag(code='MISSING_FRAGMENT', severity='error', message='Missing')
        self.assertIsNone(flag.details)


class TestAnalysisResult(unittest.TestCase):
    """Tests for AnalysisResult creation, serialization, and round-trip."""

    def _make_result(self):
        fragments = [
            FragmentResult(
                number=1, reversed=False, is_dnaA=False, is_dif=True,
                length=30000, coords=[[0, 30000]],
                blast_identity=100.0, blast_alignment_length=30000,
                blast_bit_score=900.0, blast_e_value=0.0,
                blast_mismatches=0, blast_subject='1',
            ),
            FragmentResult(
                number=2, reversed=True, is_dnaA=False, is_dif=False,
                length=40000, coords=[[30000, 70000]],
                blast_identity=98.0, blast_alignment_length=39000,
                blast_bit_score=800.0, blast_e_value=0.0,
                blast_mismatches=10, blast_subject='2',
            ),
            FragmentResult(
                number=3, reversed=False, is_dnaA=True, is_dif=False,
                length=50000, coords=[[70000, 120000]],
                blast_identity=99.5, blast_alignment_length=49800,
                blast_bit_score=1200.0, blast_e_value=0.0,
                blast_mismatches=5, blast_subject='3',
            ),
        ]
        operons = [
            OperonResult(start=100, end=5000, direction='forward'),
            OperonResult(start=30100, end=35000, direction='reverse'),
        ]
        qc_flags = [
            QCFlag(code='NOVEL_ORIENTATION', severity='warning', message='New orientation'),
        ]
        return AnalysisResult(
            genome_file='test.fasta',
            genome_length=4800000,
            is_circular=True,
            num_operons=2,
            gs_type='GS1.3',
            quality='GREEN',
            is_novel=False,
            fragment_pattern='1\t2\'\t3',
            orientation_binary=2,
            confidence_score=95.0,
            fragments=fragments,
            operons=operons,
            qc_flags=qc_flags,
            validation_passed=True,
            operon_direction_string='Valid\t--> 1 <-- 2\' --> 3(Ori)',
        )

    def test_creation(self):
        result = self._make_result()
        self.assertEqual(result.gs_type, 'GS1.3')
        self.assertEqual(result.quality, 'GREEN')
        self.assertEqual(len(result.fragments), 3)
        self.assertEqual(len(result.operons), 2)

    def test_to_dict(self):
        result = self._make_result()
        d = result.to_dict()
        self.assertIsInstance(d, dict)
        self.assertEqual(d['genome_file'], 'test.fasta')
        self.assertEqual(d['genome_length'], 4800000)
        self.assertEqual(len(d['fragments']), 3)
        self.assertEqual(d['fragments'][0]['number'], 1)
        self.assertEqual(d['qc_flags'][0]['code'], 'NOVEL_ORIENTATION')

    def test_to_json(self):
        result = self._make_result()
        json_str = result.to_json()
        parsed = json.loads(json_str)
        self.assertEqual(parsed['gs_type'], 'GS1.3')
        self.assertAlmostEqual(parsed['confidence_score'], 95.0)

    def test_to_dict_roundtrip(self):
        """Verify that to_dict produces output that can be serialized and deserialized."""
        result = self._make_result()
        d = result.to_dict()
        json_str = json.dumps(d)
        parsed = json.loads(json_str)
        self.assertEqual(parsed['genome_file'], result.genome_file)
        self.assertEqual(parsed['gs_type'], result.gs_type)
        self.assertEqual(len(parsed['fragments']), len(result.fragments))
        self.assertEqual(len(parsed['operons']), len(result.operons))
        self.assertEqual(len(parsed['qc_flags']), len(result.qc_flags))

    def test_default_empty_lists(self):
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
        self.assertEqual(result.fragments, [])
        self.assertEqual(result.operons, [])
        self.assertEqual(result.qc_flags, [])
        self.assertTrue(result.validation_passed)


if __name__ == '__main__':
    unittest.main()
