import os
import unittest
from unittest.mock import patch

from socru.FilterBlast import FilterBlast

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data','filter_blast')

class TestFilterBlast(unittest.TestCase):

    def test_filter_blast_valid(self):
        fb = FilterBlast(os.path.join(data_dir, 'blast_results_shrink'), 1, 1, False)

        self.assertEqual(len(fb.pileup_fragment(4)), 51494)
        pileup = fb.pileup_fragment(4)
        self.assertEqual(fb.calc_coverage_threshold(10, pileup, 50), 9)
        self.assertEqual(fb.identify_regions( 4, 10000),[[7229, 30663]])

    def test_return_top_result(self):
        """Verify return_top_result returns the first filtered result without re-filtering."""
        fb = FilterBlast(os.path.join(data_dir, 'blast_results_shrink'), 1, 1, False)
        top = fb.return_top_result()
        self.assertIsNotNone(top)
        # The top result should be the same as filter_results()[0]
        expected = fb.filter_results()[0]
        self.assertEqual(str(top), str(expected))

    def test_return_top_result_calls_filter_once(self):
        """Verify return_top_result does not call filter_results more than once."""
        fb = FilterBlast(os.path.join(data_dir, 'blast_results_shrink'), 1, 1, False)
        with patch.object(fb, 'filter_results', wraps=fb.filter_results) as mock_filter:
            fb.return_top_result()
            self.assertEqual(mock_filter.call_count, 1)
