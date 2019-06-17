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
        self.verbose = False

class TestSocruLookup(unittest.TestCase):

    def test_socru_lookup(self):
        g = SocruLookup(TestOptions('1-2-3-4-5-6-7', data_dir))
        self.assertEqual('GS1.0', g.calc_type())
        
    def test_different_orientation(self):
        g = SocruLookup(TestOptions("1'-2-3-4'-5-6-7", data_dir))
        self.assertEqual('GS1.9', g.calc_type())
    
    def test_inverted(self):
        g = SocruLookup(TestOptions("1-7'-6'-5'-4'-3'-2", data_dir))
        self.assertEqual('GS1.124', g.calc_type())

    def test_unknown_lookup(self):
        g = SocruLookup(TestOptions("1'-5-2-3'-6-7", data_dir))
        self.assertEqual('GS0.5', g.calc_type())