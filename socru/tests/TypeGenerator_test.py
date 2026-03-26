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

    def test_novel_profile_order_not_in_db(self):
        """A fragment order not in the database should yield order 0."""
        profiles = Profiles(os.path.join(data_dir, 'initial_profiles.txt'), False)
        # Completely novel order: 1 3 2 7 6 5 4 (not matching order 1 or 2)
        p = GATProfile(False, fragments=['1', '3', '2', '7', '6', '5', '4'])
        tg = TypeGenerator(profiles, p, False, True, prefix='GS')
        self.assertEqual(0, tg.calculate_orientationless_order())
        self.assertFalse(tg.has_previously_seen)
        # Novel + valid fragments = AMBER
        self.assertEqual('AMBER', tg.quality)

    def test_novel_profile_invalid_fragments(self):
        """Novel order with invalid fragments should be RED quality."""
        profiles = Profiles(os.path.join(data_dir, 'initial_profiles.txt'), False)
        p = GATProfile(False, fragments=['1', '3', '2', '7', '6', '5', '4'])
        tg = TypeGenerator(profiles, p, False, False, prefix='GS')
        self.assertEqual('RED', tg.quality)

    def test_profile_with_unknown_fragments(self):
        """Profile containing '?' unknown fragments should not crash."""
        profiles = Profiles(os.path.join(data_dir, 'initial_profiles.txt'), False)
        p = GATProfile(False, fragments=['1', '?', '3', '5', '4', '6', '7'])
        tg = TypeGenerator(profiles, p, False, True, prefix='GS')
        # Should not raise and should produce some GS type string
        self.assertTrue(tg.gs_type.startswith('GS'))

    def test_empty_profile(self):
        """Empty fragment list should not crash the type generator."""
        profiles = Profiles(os.path.join(data_dir, 'initial_profiles.txt'), False)
        p = GATProfile(False, fragments=[])
        tg = TypeGenerator(profiles, p, False, True, prefix='GS')
        self.assertTrue(tg.gs_type.startswith('GS'))
        self.assertEqual(0, tg.calculate_orientationless_order())

