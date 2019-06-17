import unittest
import os
import shutil
from socru.Results  import Results


test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data','results')

class TestResults(unittest.TestCase):
   
    def test_results(self):
        r = Results(os.path.join(data_dir, 'yersinia'), False)
        
        valid_patterns  = [str(p)for p in r.filter(7)]
        self.assertEqual(valid_patterns, ["1'\t7'\t6'\t5'\t3\t4\t2'", 
                                          "1'\t2\t3\t4\t5\t6\t7", 
                                          "1'\t7'\t3\t4\t5\t6\t2'", 
                                          "1\t7'\t3\t4\t5\t6\t2'"])


