import unittest
from socru.NoveltyDetector import assess_novelty, NoveltyAssessment


class TestNoveltyDetector(unittest.TestCase):

    def test_exact_match_distance_zero(self):
        """Query matching a known profile exactly should have distance 0."""
        query = ["1", "2", "3", "4"]
        known = [["1", "2", "3", "4"], ["1", "3", "2", "4"]]
        result = assess_novelty(query, known, confidence_score=90, blast_identities=[99, 98, 97, 99])
        self.assertEqual(result.edit_distance, 0)
        self.assertIsInstance(result, NoveltyAssessment)

    def test_single_inversion_distance_one(self):
        """A single inversion should give edit distance 1."""
        query = ["1", "2'", "3", "4"]
        known = [["1", "2", "3", "4"]]
        result = assess_novelty(query, known, confidence_score=85, blast_identities=[99, 98, 97, 96])
        self.assertEqual(result.edit_distance, 1)
        self.assertTrue(len(result.fragment_differences) > 0)

    def test_translocated_fragments(self):
        """Translocated fragments should increase edit distance."""
        query = ["1", "3", "2", "4"]
        known = [["1", "2", "3", "4"]]
        result = assess_novelty(query, known, confidence_score=85, blast_identities=[99, 98, 97, 96])
        self.assertGreater(result.edit_distance, 0)

    def test_likely_real_assessment(self):
        """High confidence + high BLAST identities -> likely_real."""
        query = ["1", "2'", "3", "4"]
        known = [["1", "2", "3", "4"]]
        result = assess_novelty(
            query, known,
            confidence_score=95,
            blast_identities=[99.5, 98.0, 97.0, 96.0],
        )
        self.assertEqual(result.assessment, "likely_real")
        self.assertTrue(result.is_likely_real)

    def test_possibly_artifactual_low_confidence(self):
        """Low confidence -> possibly_artifactual."""
        query = ["1", "2'", "3", "4"]
        known = [["1", "2", "3", "4"]]
        result = assess_novelty(
            query, known,
            confidence_score=30,
            blast_identities=[99, 98, 97, 96],
        )
        self.assertEqual(result.assessment, "possibly_artifactual")
        self.assertFalse(result.is_likely_real)

    def test_possibly_artifactual_low_blast(self):
        """Low BLAST identity -> possibly_artifactual."""
        query = ["1", "2'", "3", "4"]
        known = [["1", "2", "3", "4"]]
        result = assess_novelty(
            query, known,
            confidence_score=90,
            blast_identities=[99, 85, 97, 96],
        )
        self.assertEqual(result.assessment, "possibly_artifactual")
        self.assertFalse(result.is_likely_real)

    def test_uncertain_assessment(self):
        """Middling confidence and BLAST -> uncertain."""
        query = ["1", "2'", "3", "4"]
        known = [["1", "2", "3", "4"]]
        result = assess_novelty(
            query, known,
            confidence_score=65,
            blast_identities=[93, 92, 91, 94],
        )
        self.assertEqual(result.assessment, "uncertain")
        self.assertFalse(result.is_likely_real)

    def test_missing_fragment_question_mark(self):
        """Profiles with '?' for missing fragments should still work."""
        query = ["1", "?", "3", "4"]
        known = [["1", "2", "3", "4"]]
        result = assess_novelty(
            query, known,
            confidence_score=30,
            blast_identities=[99, 0, 97, 96],
        )
        self.assertGreater(result.edit_distance, 0)
        self.assertEqual(result.assessment, "possibly_artifactual")

    def test_no_known_profiles(self):
        """Empty known profiles list should return uncertain."""
        query = ["1", "2", "3"]
        result = assess_novelty(query, [], confidence_score=90, blast_identities=[99, 98, 97])
        self.assertEqual(result.assessment, "uncertain")
        self.assertEqual(result.nearest_known_type, "none")

    def test_nearest_known_type_label(self):
        """Should identify the correct nearest profile by index."""
        query = ["1", "2", "3"]
        known = [
            ["1", "3", "2"],   # distance > 0
            ["1", "2", "3"],   # exact match, distance 0
        ]
        result = assess_novelty(query, known, confidence_score=90, blast_identities=[99, 98, 97])
        self.assertEqual(result.nearest_known_type, "profile_1")
        self.assertEqual(result.edit_distance, 0)


if __name__ == "__main__":
    unittest.main()
