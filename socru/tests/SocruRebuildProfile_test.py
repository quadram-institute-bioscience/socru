import unittest
import os
import shutil
import filecmp
from socru.SocruRebuildProfile  import SocruRebuildProfile

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data','rebuild_profile')

class TestOptions:
    def __init__(self, profile_filename, output_file,prefix):
        self.profile_filename = profile_filename
        self.output_file = output_file
        self.prefix = prefix

 
class TestSocruUpdateProfile(unittest.TestCase):

    def test_socru_lookup(self):
        if os.path.exists('output_profile.txt'):
            os.remove('output_profile.txt')
            
        s = SocruRebuildProfile(TestOptions(os.path.join(data_dir, 'expanded_se_profile.txt'), 'output_profile.txt', 'GS'))
        s.run()
        
        self.assertTrue(os.path.exists('output_profile.txt'))
        self.assertTrue(filecmp.cmp('output_profile.txt', os.path.join(data_dir, 'expected_output_profile.txt')))
    
        os.remove('output_profile.txt')
        
