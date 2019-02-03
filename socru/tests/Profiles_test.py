import unittest
import os
import shutil
from socru.Profiles  import Profiles


test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data','profiles')

class TestProfiles(unittest.TestCase):
   
    def test_profiles(self):
        p = Profiles(os.path.join(data_dir, 'profile.txt'))
        self.assertEqual(7, len(p.gats))
        self.assertEqual("1	2	3	4	5	6	7", str(p.gats[0]))
        
        p.gats[0].fragments = p.gats[0].invert_fragments()
        self.assertEqual("1'	7'	6'	5'	4'	3'	2'", str(p.gats[0]))
        
        self.assertEqual(5, p.dnaA_fragment_number)
        