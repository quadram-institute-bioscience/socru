import unittest
import os
import shutil
from socru.SocruLookup  import SocruLookup

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data','lookup')
salmonella_dir = os.path.join(data_dir, 'salmonella')

class TestOptions:
    def __init__(self, fragments, db_dir):
        self.fragments = fragments
        self.db_dir = db_dir
        self.verbose = False

class TestSocruLookup(unittest.TestCase):

    def test_socru_lookup(self):
        g = SocruLookup(TestOptions('1-2-3-4-5-6-7', data_dir))
        self.assertEqual('GS1.0', g.calc_type())
        self.assertEqual('GREEN', g.calc_quality())
        
    def test_different_orientation(self):
        g = SocruLookup(TestOptions("1'-2-3-4'-5-6-7", data_dir))
        self.assertEqual('GS1.9', g.calc_type())
        self.assertEqual('RED', g.calc_quality())
    
    def test_inverted_invalid(self):
        g = SocruLookup(TestOptions("1-7'-6'-5'-4'-3'-2", data_dir))
        self.assertEqual('GS1.3', g.calc_type())
        self.assertEqual('RED', g.calc_quality())
   
    def test_unknown_lookup(self):
        g = SocruLookup(TestOptions("1'-5-2-3'-6-7", data_dir))
        self.assertEqual('GS0.114', g.calc_type())
        self.assertEqual('RED', g.calc_quality())
        
        
    def test_salmonella_gs(self):
        valid_gs = [ "1-2-3-4-5-6-7",
        "1'-2-3-4-5-6-7",
        "1'-7'-6'-5'-4'-3-2'",
        "1-7'-3-5-6-4-2'",
        "1'-7'-3-5-6-4-2'",
        "1'-7'-3-4-6-5-2'",
        "1-2-3-6-4-5-7",
        "1'-2-3-6-4-5-7",
        "1'-7'-2-3-4-5-6",
        "1-7'-3-5-4-6-2'",
        "1'-7'-3-5-4-6-2'",
        "1-2-4'-5'-6'-3-7",
        "1-7'-3-4-5-6-2'",
        "1'-7'-3-4-5-6-2'",
        "1-3-4-5-6-7-2'",
        "1'-3-4-5-6-7-2'",
        "1-7'-4'-3-5-6-2'",
        "1-2-3-5-6-4-7",
        "1'-2-3-5-6-4-7",
        "1-7'-4'-6'-5'-3-2'",
        "1'-2-3-4-6-5-7",
        "1'-5'-2-3-6-4-7",
        "1-2-3-6-5-4-7",
        "1'-2-3-6-5-4-7",
        "1'-2-5'-3-6-4-7",
        "1'-2-6'-3-5-4-7",
        "1-2-3-5-4-6-7",
        "1-7'-6'-4'-5'-3-2'",
        "1-7'-5'-3-4-6-2'",
        "1'-7'-5'-3-4-6-2'",
        "1'-7'-3-6-4-5-2'"]
        
        for gs in valid_gs:
            g = SocruLookup(TestOptions(gs, salmonella_dir))
            self.assertEqual('GREEN', g.calc_quality())
        