import unittest
from socru.RearrangementDistance import rearrangement_distance, pairwise_distance_matrix


class TestRearrangementDistance(unittest.TestCase):

    def test_identical_profiles(self):
        """Identical profiles should have distance 0."""
        result = rearrangement_distance(["1", "2", "3"], ["1", "2", "3"])
        self.assertEqual(result["distance"], 0)
        self.assertEqual(result["inversions"], 0)
        self.assertEqual(result["translocations"], 0)
        self.assertEqual(result["details"], [])

    def test_single_inversion(self):
        """One fragment flipped should yield 1 inversion, 0 translocations."""
        result = rearrangement_distance(["1", "2", "3"], ["1", "2'", "3"])
        self.assertEqual(result["inversions"], 1)
        self.assertEqual(result["translocations"], 0)
        self.assertEqual(result["distance"], 1)

    def test_single_translocation(self):
        """Swapping two fragments (same orientation) should yield translocations."""
        result = rearrangement_distance(["1", "2", "3"], ["1", "3", "2"])
        self.assertEqual(result["translocations"], 2)
        self.assertEqual(result["inversions"], 0)
        self.assertEqual(result["distance"], 2)

    def test_combined_inversion_and_translocation(self):
        """Inversion + translocation counted together."""
        # Fragment 2 moves from pos 1->2 and flips; fragment 3 moves from pos 2->1
        result = rearrangement_distance(["1", "2", "3"], ["1", "3", "2'"])
        self.assertGreater(result["inversions"], 0)
        self.assertGreater(result["translocations"], 0)
        self.assertEqual(
            result["distance"],
            result["inversions"] + result["translocations"] + result["missing"],
        )

    def test_symmetry(self):
        """distance(A, B) should equal distance(B, A)."""
        a = ["1", "2'", "3", "4"]
        b = ["1", "3", "2", "4'"]
        r1 = rearrangement_distance(a, b)
        r2 = rearrangement_distance(b, a)
        self.assertEqual(r1["distance"], r2["distance"])
        self.assertEqual(r1["inversions"], r2["inversions"])
        self.assertEqual(r1["translocations"], r2["translocations"])

    def test_missing_fragment(self):
        """Profiles of different lengths with a missing fragment."""
        result = rearrangement_distance(["1", "2", "3", "4"], ["1", "2", "3"])
        self.assertGreater(result["missing"], 0)
        self.assertGreater(result["distance"], 0)

    def test_pairwise_distance_matrix(self):
        """Matrix with 3 profiles should contain all 9 pairs."""
        profiles = {
            "GS1.0": ["1", "2", "3"],
            "GS2.0": ["1", "3", "2"],
            "GS3.0": ["1", "2'", "3"],
        }
        matrix = pairwise_distance_matrix(profiles)

        # Should have n*n entries (including self-pairs)
        self.assertEqual(len(matrix), 9)

        # Self-distances should be 0
        for name in profiles:
            self.assertEqual(matrix[(name, name)]["distance"], 0)

        # Symmetry
        for a in profiles:
            for b in profiles:
                self.assertEqual(
                    matrix[(a, b)]["distance"], matrix[(b, a)]["distance"]
                )

    def test_all_inverted(self):
        """Every fragment inverted should count all as inversions."""
        result = rearrangement_distance(
            ["1", "2", "3"], ["1'", "2'", "3'"]
        )
        self.assertEqual(result["inversions"], 3)
        self.assertEqual(result["translocations"], 0)
        self.assertEqual(result["distance"], 3)


if __name__ == "__main__":
    unittest.main()
