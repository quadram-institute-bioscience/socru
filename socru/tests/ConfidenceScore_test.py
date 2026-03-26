"""Tests for the ConfidenceScore module."""

import unittest

from socru.AnalysisResult import FragmentResult
from socru.ConfidenceScore import calculate_confidence


def _make_fragment(number, length=50000, identity=100.0, alignment_length=50000):
    """Helper to create a FragmentResult with BLAST data."""
    return FragmentResult(
        number=number,
        reversed=False,
        is_dnaA=False,
        is_dif=False,
        length=length,
        coords=[[0, length]],
        blast_identity=identity,
        blast_alignment_length=alignment_length,
        blast_bit_score=1000.0,
        blast_e_value=0.0,
        blast_mismatches=0,
        blast_subject=str(number),
    )


def _make_unknown_fragment():
    """Helper to create an unknown ('?') FragmentResult."""
    return FragmentResult(
        number='?',
        reversed=False,
        is_dnaA=False,
        is_dif=False,
        length=10000,
        coords=[[0, 10000]],
    )


class TestConfidenceScore(unittest.TestCase):
    """Tests for confidence score calculation."""

    def test_perfect_score(self):
        """All identities 100%, exact match, valid => 100."""
        frags = [_make_fragment(i) for i in range(1, 8)]
        score = calculate_confidence(frags, 'GREEN', True, is_novel=False)
        self.assertAlmostEqual(score, 100.0)

    def test_worst_case(self):
        """No BLAST results, novel, invalid => 0."""
        frags = [_make_unknown_fragment() for _ in range(3)]
        score = calculate_confidence(frags, 'RED', False, is_novel=True)
        self.assertAlmostEqual(score, 0.0)

    def test_penalty_for_unknown_fragments(self):
        """Unknown fragments should subtract 10 points each."""
        good_frags = [_make_fragment(i) for i in range(1, 6)]
        score_without_unknown = calculate_confidence(
            good_frags, 'GREEN', True, is_novel=False,
        )

        # Add one unknown fragment
        frags_with_unknown = good_frags + [_make_unknown_fragment()]
        score_with_unknown = calculate_confidence(
            frags_with_unknown, 'GREEN', True, is_novel=False,
        )
        # The unknown fragment lowers the mean identity and coverage (since it
        # has no BLAST data), AND applies the 10-point penalty.
        self.assertLess(score_with_unknown, score_without_unknown)

    def test_penalty_magnitude(self):
        """Two unknowns should score lower than one unknown."""
        frags_one = [_make_fragment(1), _make_unknown_fragment()]
        frags_two = [_make_fragment(1), _make_unknown_fragment(), _make_unknown_fragment()]
        score_one = calculate_confidence(frags_one, 'AMBER', True, is_novel=True)
        score_two = calculate_confidence(frags_two, 'AMBER', True, is_novel=True)
        self.assertLess(score_two, score_one)

    def test_low_identity_reduces_score(self):
        """Fragments with lower identity should give a lower score."""
        perfect = [_make_fragment(i, identity=100.0) for i in range(1, 4)]
        low_id = [_make_fragment(i, identity=80.0) for i in range(1, 4)]
        score_perfect = calculate_confidence(perfect, 'GREEN', True)
        score_low = calculate_confidence(low_id, 'GREEN', True)
        self.assertGreater(score_perfect, score_low)

    def test_validation_failure_reduces_score(self):
        """Invalid arrangement should reduce score by the validation weight."""
        frags = [_make_fragment(i) for i in range(1, 4)]
        score_valid = calculate_confidence(frags, 'GREEN', True)
        score_invalid = calculate_confidence(frags, 'RED', False)
        self.assertGreater(score_valid, score_invalid)

    def test_score_clamped_to_zero(self):
        """Score should never go below zero even with many penalties."""
        frags = [_make_unknown_fragment() for _ in range(15)]
        score = calculate_confidence(frags, 'RED', False, is_novel=True)
        self.assertEqual(score, 0.0)

    def test_score_clamped_to_hundred(self):
        """Score should never exceed 100."""
        frags = [_make_fragment(i) for i in range(1, 4)]
        score = calculate_confidence(frags, 'GREEN', True)
        self.assertLessEqual(score, 100.0)

    def test_empty_fragments(self):
        """Empty fragment list should give 0 for identity/coverage components."""
        score = calculate_confidence([], 'GREEN', True)
        # profile_match (GREEN) = 25 + validation = 20 = 45
        self.assertAlmostEqual(score, 45.0)

    def test_partial_coverage_reduces_score(self):
        """Fragments with alignment shorter than length should score lower."""
        full = [_make_fragment(1, length=50000, alignment_length=50000)]
        partial = [_make_fragment(1, length=50000, alignment_length=25000)]
        score_full = calculate_confidence(full, 'GREEN', True)
        score_partial = calculate_confidence(partial, 'GREEN', True)
        self.assertGreater(score_full, score_partial)


if __name__ == '__main__':
    unittest.main()
