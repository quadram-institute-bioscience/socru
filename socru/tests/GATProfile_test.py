import unittest

from socru.GATProfile import GATProfile


class TestGATProfile(unittest.TestCase):

    def test_gat_profile(self):
        p = GATProfile(False, fragments = ['1','5','7','3'])
        self.assertTrue(p.is_profile_in_correct_orientation())

    def test_inverted(self):
        p = GATProfile(False, fragments = ['1','5','7','3\''])
        self.assertFalse(p.is_profile_in_correct_orientation())
        self.assertEqual(["1'", '3', "7'", "5'"] , p.invert_fragments())

        p.orientate_for_dnaA()
        self.assertTrue(p.is_profile_in_correct_orientation())

    def test_reorientate_list_to_start_with_one(self):
        p = GATProfile(False, fragments = [])
        self.assertEqual(['1','3','7','2','5','6'],  p.reorientate_list_to_start_with_one(['2','5','6','1','3','7']))

    def test_reorientate(self):
        p = GATProfile(False, fragments = ['1','5','7','3\''])
        self.assertEqual(['1','5','7','3'],  p.orientationless_fragments())

    def test_orientation_binary(self):
        p = GATProfile(False, fragments = ['1\'','5','7','3'])
        self.assertEqual(1,  p.orientation_binary())

        p = GATProfile(False, fragments = ['1','5','7','2\''])
        self.assertEqual(2,  p.orientation_binary())

        p = GATProfile(False, fragments = ['1','3\'','7','2'])
        self.assertEqual(4,  p.orientation_binary())

        p = GATProfile(False, fragments = ['1','3','4\'','2'])
        self.assertEqual(8,  p.orientation_binary())

        p = GATProfile(False, fragments = ['1\'','3\'','4\'','2\''])
        self.assertEqual(15,  p.orientation_binary())

    def test_mutable_default_fragments_not_shared(self):
        """Two GATProfile instances without explicit fragments should not share the same list."""
        g1 = GATProfile(False)
        g2 = GATProfile(False)
        g1.fragments.append("1")
        self.assertEqual(g2.fragments, [])

    # ---- Edge case tests for empty, single, all-reversed, unknown fragments ----

    def test_empty_fragments_list(self):
        p = GATProfile(False, fragments=[])
        self.assertEqual(p.fragment_str(), '')
        self.assertEqual(p.orientation_binary(), 0)
        self.assertEqual(p.orientationless_fragments(), [])
        self.assertTrue(p.is_profile_in_correct_orientation())

    def test_single_fragment(self):
        p = GATProfile(False, fragments=['1'])
        self.assertEqual(p.fragment_str(), '1')
        self.assertEqual(p.orientation_binary(), 0)
        self.assertEqual(p.orientationless_fragments(), ['1'])

    def test_single_fragment_reversed(self):
        p = GATProfile(False, fragments=["1'"])
        self.assertEqual(p.orientation_binary(), 1)
        self.assertEqual(p.orientationless_fragments(), ['1'])

    def test_all_fragments_reversed(self):
        p = GATProfile(False, fragments=["1'", "2'", "3'", "4'"])
        self.assertEqual(p.orientation_binary(), 15)
        self.assertEqual(p.orientationless_fragments(), ['1', '2', '3', '4'])
        inverted = p.invert_fragments()
        # After inversion, all primes should be removed and order reversed
        for frag in inverted:
            self.assertNotIn("'", frag)

    def test_profile_with_unknown_question_mark_fragments(self):
        """Fragments marked '?' should pass through without crashing."""
        p = GATProfile(False, fragments=['1', '?', '3'])
        self.assertEqual(p.fragment_str(), '1\t?\t3')
        # '?' does not match digit-prime pattern, so orientation_binary ignores it
        self.assertEqual(p.orientation_binary(), 0)
        # orientationless also passes '?' through unchanged
        self.assertEqual(p.orientationless_fragments(), ['1', '?', '3'])

    def test_invert_fragments_twice_returns_original(self):
        """Inverting a profile twice should return the original fragment list."""
        original = ['1', '2', '3', '4', '5']
        p = GATProfile(False, fragments=list(original))
        inverted_once = p.invert_fragments()
        p2 = GATProfile(False, fragments=list(inverted_once))
        inverted_twice = p2.invert_fragments()
        self.assertEqual(inverted_twice, original)

    def test_invert_fragments_twice_with_mixed_orientations(self):
        """Double inversion restores original even with mixed orientations."""
        original = ['1', "2'", '3', "4'", '5']
        p = GATProfile(False, fragments=list(original))
        inverted_once = p.invert_fragments()
        p2 = GATProfile(False, fragments=list(inverted_once))
        inverted_twice = p2.invert_fragments()
        self.assertEqual(inverted_twice, original)

    def test_orientation_binary_deterministic(self):
        """orientation_binary should return the same value for identical fragments."""
        frags = ['1', "2'", '3', "5'", '7']
        p1 = GATProfile(False, fragments=list(frags))
        p2 = GATProfile(False, fragments=list(frags))
        self.assertEqual(p1.orientation_binary(), p2.orientation_binary())

    def test_does_the_profile_match_true(self):
        p1 = GATProfile(False, fragments=['1', '2', '3'])
        p2 = GATProfile(False, fragments=['1', '2', '3'])
        self.assertTrue(p1.does_the_profile_match(p2))

    def test_does_the_profile_match_false(self):
        p1 = GATProfile(False, fragments=['1', '2', '3'])
        p2 = GATProfile(False, fragments=['1', '3', '2'])
        self.assertFalse(p1.does_the_profile_match(p2))

    def test_order_extraction(self):
        p = GATProfile(False, gat_number='2.5')
        self.assertEqual(p.order(), 2)

    def test_order_extraction_no_dot(self):
        p = GATProfile(False, gat_number=0)
        self.assertEqual(p.order(), 0)

    # ---- Hypothesis stubs (hypothesis not available) ----
    # If hypothesis becomes available, the following property-based tests should
    # be implemented:
    #
    # @given(st.lists(st.sampled_from(['1','2','3','4','5',"1'","2'","3'","4'","5'"]), min_size=1, max_size=5))
    # def test_property_invert_twice_identity(self, frags):
    #     """Inverting twice should yield original (property-based)."""
    #     p = GATProfile(False, fragments=list(frags))
    #     inv = p.invert_fragments()
    #     p2 = GATProfile(False, fragments=list(inv))
    #     self.assertEqual(p2.invert_fragments(), frags)
    #
    # @given(st.lists(st.sampled_from(['1','2','3',"1'","2'","3'"]), min_size=0, max_size=5))
    # def test_property_orientation_binary_deterministic(self, frags):
    #     """orientation_binary is deterministic for any fragment list."""
    #     p1 = GATProfile(False, fragments=list(frags))
    #     p2 = GATProfile(False, fragments=list(frags))
    #     self.assertEqual(p1.orientation_binary(), p2.orientation_binary())
