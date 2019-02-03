import unittest
import os
import shutil
from socru.SocruLookup  import SocruLookup

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data','lookup')

class TestOptions:
    def __init__(self, fragments, db_dir):
        self.fragments = fragments
        self.db_dir = db_dir


class TestSocruLookup(unittest.TestCase):

    def test_socru_lookup(self):
        g = SocruLookup(TestOptions('1-2-3-4-5-6-7', data_dir))
        self.assertEqual('GS1.0', g.calc_type())
        
        g = SocruLookup(TestOptions("1'-2-3-4'-5-6-7", data_dir))
        self.assertEqual('GS0.9', g.calc_type())
