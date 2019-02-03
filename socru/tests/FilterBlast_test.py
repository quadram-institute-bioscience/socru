import unittest
import os
import shutil
import filecmp
from socru.FilterBlast  import FilterBlast

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data','filter_blast')

class TestFilterBlast(unittest.TestCase):
   
    def test_filter_blast_valid(self):
        fb = FilterBlast(os.path.join(data_dir, 'blast_results_shrink'), 1, 1)

        self.assertEqual(len(fb.pileup_fragment(4)), 51494)
        pileup = fb.pileup_fragment(4)
        self.assertEqual(fb.calc_coverage_threshold(10, pileup, 50), 9)
        self.assertEqual(fb.identify_regions( 4, 10000),[[7229, 30663]])