import unittest
import os
import shutil
from socru.Blast  import Blast

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data','blast')

class TestBlast(unittest.TestCase):
   
    def test_blast_forward(self):
        b = Blast(os.path.join(data_dir, 'database'), 1)
        results = b.run_blast(os.path.join(data_dir,'query_7.fa'))
        self.assertTrue(os.path.exists(results))
        
    def test_blast_reverse(self):
        b = Blast(os.path.join(data_dir, 'database'), 1)
        results = b.run_blast(os.path.join(data_dir,'query_7_reverse.fa'))
        self.assertTrue(os.path.exists(results))        
    