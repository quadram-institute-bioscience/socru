import tempfile
import unittest

from Bio.Seq import Seq

from socru.Fragment import Fragment
from socru.FragmentFiles import FragmentFiles


class TestFragmentFiles(unittest.TestCase):

    def test_reversed_fragment_sets_reversed_frag_attribute(self):
        """When fragment_order specifies a reversed fragment (e.g. "3'"), reversed_frag should be True."""
        # Create 3 fragments with different sizes so the largest-first logic is predictable.
        # Fragment sizes: 100, 50, 30 - so order after largest-first is [f1, f2, f3].
        f1 = Fragment([[0, 100]], sequence=Seq("A" * 100))
        f2 = Fragment([[0, 50]], sequence=Seq("C" * 50))
        f3 = Fragment([[0, 30]], sequence=Seq("G" * 30))

        with tempfile.TemporaryDirectory() as tmpdir:
            # fragment_order "1-2-3'" means third fragment should be reversed
            ff = FragmentFiles([f1, f2, f3], tmpdir, False, fragment_order="1-2-3'")
            ordered = ff.ordered_fragments

            self.assertFalse(ordered[0].reversed_frag)
            self.assertFalse(ordered[1].reversed_frag)
            self.assertTrue(ordered[2].reversed_frag)
