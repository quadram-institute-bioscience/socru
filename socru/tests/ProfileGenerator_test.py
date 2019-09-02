import unittest
import os
import shutil
import filecmp
from socru.ProfileGenerator  import ProfileGenerator

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data','profile_generator')

class TestProfileGenerator(unittest.TestCase):
   
    def test_profile_generator(self):
        p = ProfileGenerator( os.path.join(data_dir,'database'), 7, os.path.join(data_dir, 'dnaA.fa.gz'),os.path.join(data_dir, 'dif.fa.gz'), 1, 'reffile.fa', False)
        p.write_output_file()
        
        self.assertTrue(os.path.exists(os.path.join(data_dir,'database', 'profile.txt')))
        self.assertTrue(filecmp.cmp( os.path.join(data_dir,'database', 'profile.txt'), os.path.join(data_dir, 'expected_profile.txt')))
        
        self.assertTrue(os.path.exists(os.path.join(data_dir,'database', 'profile.txt.yml')))
        
        os.remove(os.path.join(data_dir,'database', 'profile.txt'))
        os.remove(os.path.join(data_dir,'database', 'profile.txt.yml'))
