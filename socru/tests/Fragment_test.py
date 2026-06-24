import unittest

from socru.Fragment import Fragment


class TestFragment(unittest.TestCase):

    def test_fragment_seq_name_one_coord(self):
        f = Fragment([[10, 123]])
        self.assertEqual(str(f), '0 10_123')

    def test_fragment_seq_name_two_coords(self):
        f = Fragment([[10, 123], [555, 999]])
        self.assertEqual(str(f), '0 10_123__555_999')

    def test_num_bases_empty_sequence(self):
        f = Fragment([[10, 123]])
        self.assertEqual(f.num_bases(), 0)

    def test_num_bases_with_sequence(self):
        f = Fragment([[10, 123]], sequence="ACGTACGT")
        self.assertEqual(f.num_bases(), 8)

    def test_num_bases_long_sequence(self):
        seq = "A" * 5000
        f = Fragment([[1, 5000]], sequence=seq)
        self.assertEqual(f.num_bases(), 5000)

    def test_output_filename_default(self):
        f = Fragment([[10, 123]])
        self.assertEqual(f.output_filename(), '0.fa')

    def test_output_filename_numbered(self):
        f = Fragment([[10, 123]], number=7)
        self.assertEqual(f.output_filename(), '7.fa')

    def test_operon_direction_str_forward_not_reversed(self):
        f = Fragment([[1, 100]], number=3, operon_forward_start=True, reversed_frag=False)
        self.assertEqual(f.operon_direction_str(), '--> 3')

    def test_operon_direction_str_reverse_start_not_reversed(self):
        f = Fragment([[1, 100]], number=3, operon_forward_start=False, reversed_frag=False)
        self.assertEqual(f.operon_direction_str(), '<-- 3')

    def test_operon_direction_str_forward_reversed_frag(self):
        f = Fragment([[1, 100]], number=3, operon_forward_start=True, reversed_frag=True)
        self.assertEqual(f.operon_direction_str(), "--> 3'")

    def test_operon_direction_str_reverse_start_reversed_frag(self):
        f = Fragment([[1, 100]], number=5, operon_forward_start=False, reversed_frag=True)
        self.assertEqual(f.operon_direction_str(), "<-- 5'")

    def test_operon_direction_str_with_dnaA(self):
        f = Fragment([[1, 100]], number=3, dna_A=True)
        self.assertEqual(f.operon_direction_str(), '--> 3(Ori)')

    def test_operon_direction_str_with_dif(self):
        f = Fragment([[1, 100]], number=1, dif=True)
        self.assertEqual(f.operon_direction_str(), '--> 1(Ter)')

    def test_operon_direction_str_with_both_markers(self):
        f = Fragment([[1, 100]], number=3, dna_A=True, dif=True)
        self.assertEqual(f.operon_direction_str(), '--> 3(Ori)(Ter)')

    def test_operon_direction_str_reversed_with_markers(self):
        f = Fragment([[1, 100]], number=3, reversed_frag=True, dna_A=True, dif=True)
        self.assertEqual(f.operon_direction_str(), "--> 3'(Ori)(Ter)")

    def test_multiple_coordinate_ranges(self):
        f = Fragment([[1, 100], [500, 600], [800, 900]], number=2)
        self.assertEqual(str(f), '2 1_100__500_600__800_900')

    def test_reversed_frag_default_false(self):
        f = Fragment([[1, 100]])
        self.assertFalse(f.reversed_frag)

    def test_reversed_frag_set_true(self):
        f = Fragment([[1, 100]], reversed_frag=True)
        self.assertTrue(f.reversed_frag)

    def test_reversed_frag_can_be_toggled(self):
        f = Fragment([[1, 100]], reversed_frag=False)
        self.assertFalse(f.reversed_frag)
        f.reversed_frag = True
        self.assertTrue(f.reversed_frag)

    def test_dna_A_default_false(self):
        f = Fragment([[1, 100]])
        self.assertFalse(f.dna_A)

    def test_dif_default_false(self):
        f = Fragment([[1, 100]])
        self.assertFalse(f.dif)

    def test_operon_forward_defaults(self):
        f = Fragment([[1, 100]])
        self.assertTrue(f.operon_forward_start)
        self.assertTrue(f.operon_forward_end)
