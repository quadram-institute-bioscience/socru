import unittest
import os
import shutil
from socru.GATProfile  import GATProfile

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
        
    def test_reorientate(self):
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