import unittest
import os
import shutil
import filecmp
from socru.SocruShrinkDatabase  import SocruShrinkDatabase

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data','shrink_database')

class TestOptions:
    def __init__(self, blast_results, input_database, output_database, min_fragment_size):
        self.blast_results = blast_results
        self.input_database = input_database
        self.output_database = output_database
        self.min_fragment_size = min_fragment_size
 
class TestSocruShrinkDatabase(unittest.TestCase):

    def test_shrink_database(self):     
        if os.path.exists('output_dir'):
             shutil.rmtree('output_dir')       
        s = SocruShrinkDatabase(
                TestOptions(
                    os.path.join(data_dir, 'blast_results'),
                    os.path.join(data_dir, 'input_database'), 
                    'output_dir', 
                    100))
        self.assertEqual(s.run(), "input_database	0.249	99.75")
        shutil.rmtree('output_dir')
        
