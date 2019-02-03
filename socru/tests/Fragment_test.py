import unittest
import os
import shutil
from socru.Fragment  import Fragment

class TestFragment(unittest.TestCase):
   
    def test_fragment_seq_name_one_coord(self):
        f = Fragment([[10,123]])
        self.assertEqual(str(f), '0 10_123')

    def test_fragment_seq_name_two_coords(self):
        f = Fragment([[10,123],[555,999]])
        self.assertEqual(str(f), '0 10_123__555_999')
        