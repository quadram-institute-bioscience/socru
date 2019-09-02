import unittest
import os
import shutil
from socru.Dif  import Dif

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data','dif')

class TestDif(unittest.TestCase):
   
    def test_dif_compressed(self):
        d = Dif(os.path.join(data_dir,'dif.fa.gz'), os.path.join(data_dir,'database'), 1 , False)
        
        self.assertFalse(d.dif_orientation)
        self.assertTrue(int(d.fragment_with_dif) >= 2)
    
        