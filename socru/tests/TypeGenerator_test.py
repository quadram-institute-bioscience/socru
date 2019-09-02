import unittest
import os
import shutil
from socru.TypeGenerator  import TypeGenerator
from socru.Profiles  import Profiles
from socru.GATProfile import GATProfile

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data','type_generator')

class TestTypeGenerator(unittest.TestCase):
   
    def test_type_generator_valid_novel(self):
        profiles = Profiles( os.path.join(data_dir, 'initial_profiles.txt'), False)
        p = GATProfile(False, fragments = ['1', '7\'', '6\'', '4\'', '5\'', '3', '2\''])
        tg = TypeGenerator(profiles, p, False, True, prefix = 'GS')
        self.assertEqual(2, tg.calculate_orientationless_order())
        self.assertEqual('AMBER', tg.quality)
        self.assertEqual('GS2.122', tg.calculate_type())
        
    def test_type_generator_valid_seen(self):
        profiles = Profiles( os.path.join(data_dir, 'initial_profiles.txt'), False)
        p = GATProfile(False, fragments = ['1', '2', '3', '5', '4', '6', '7'])
        tg = TypeGenerator(profiles, p, False, True, prefix = 'GS')
        self.assertEqual(2, tg.calculate_orientationless_order())
        self.assertEqual('GREEN', tg.quality)
        self.assertEqual('GS2.0', tg.calculate_type())
        
    def test_type_generator_invalid(self):
        profiles = Profiles( os.path.join(data_dir, 'initial_profiles.txt'), False)
        p = GATProfile(False, fragments = ['1', '7\'', '6\'', '4\'', '5\'', '3', '2\''])
        tg = TypeGenerator(profiles, p, False, False, prefix = 'GS')
        self.assertEqual(2, tg.calculate_orientationless_order())
        self.assertEqual('RED', tg.quality)
        self.assertEqual('GS2.122', tg.calculate_type())
    
