import unittest
import os
import shutil
import filecmp
from socru.SocruUpdateProfile  import SocruUpdateProfile

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data','update_profile')

class TestOptions:
    def __init__(self, socru_output_filename, profile_filename, output_file):
        self.socru_output_filename = socru_output_filename
        self.profile_filename = profile_filename
        self.output_file = output_file
        self.verbose = False
 
class TestSocruUpdateProfile(unittest.TestCase):

    def test_socru_lookup(self):
        if os.path.exists('output_profile.txt'):
            os.remove('output_profile.txt')
            
        s = SocruUpdateProfile(TestOptions(os.path.join(data_dir, 'yersinia'),os.path.join(data_dir, 'profile.txt'), 'output_profile.txt'))
        s.run()
        self.assertTrue(os.path.exists('output_profile.txt'))
        self.assertTrue(filecmp.cmp('output_profile.txt', os.path.join(data_dir, 'expected_output_profile.txt')))
    
        os.remove('output_profile.txt')
        
    def test_same_update_twice(self):
        if os.path.exists('output_profile.txt'):
            os.remove('output_profile.txt')
            
        s = SocruUpdateProfile(TestOptions(os.path.join(data_dir, 'yersinia'),os.path.join(data_dir, 'profile_full.txt'), 'output_profile.txt'))
        s.run()
        
        self.assertTrue(os.path.exists('output_profile.txt'))
        self.assertTrue(filecmp.cmp('output_profile.txt', os.path.join(data_dir, 'expected_output_profile.txt')))
    
        os.remove('output_profile.txt')

