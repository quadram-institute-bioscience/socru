import os
import sys
import unittest
from io import StringIO

from socru.Schemas import Schemas


class TestSchemas(unittest.TestCase):

    def test_schemas(self):
        s = Schemas(False)
        self.assertTrue( len(s.all_available()) > 4)

    def test_all_available(self):
        s = Schemas(False)
        available = s.all_available()
        self.assertIsInstance(available, list)
        self.assertTrue(len(available) > 0)

    def test_print_all(self):
        s = Schemas(False)
        # Capture stdout
        captured_output = StringIO()
        old_stdout = sys.stdout
        try:
            sys.stdout = captured_output
            s.print_all()
        finally:
            sys.stdout = old_stdout
        output = captured_output.getvalue()
        self.assertTrue(len(output) > 0)

    def test_extended(self):
        s = Schemas(False)
        extended_info = s.extended()
        self.assertIsInstance(extended_info, dict)
        # Should have at least some species in the extended info
        self.assertTrue(len(extended_info) > 0)

        # Check structure of extended info
        for species, info in extended_info.items():
            self.assertIsInstance(info, list)
            self.assertEqual(len(info), 5)  # species, dnaA_fragment_number, dnaa_forward_orientation, dif_fragment_number, reference_genome

    def test_print_extended(self):
        s = Schemas(False)
        # Capture stdout
        captured_output = StringIO()
        old_stdout = sys.stdout
        try:
            sys.stdout = captured_output
            s.print_extended()
        finally:
            sys.stdout = old_stdout
        output = captured_output.getvalue()
        self.assertTrue(len(output) > 0)
        self.assertIn('Species', output)

    def test_database_directory_with_base_directory(self):
        s = Schemas(False)
        available = s.all_available()
        if len(available) > 0:
            species = available[0]
            db_dir = s.database_directory(None, species)
            self.assertIsNotNone(db_dir)
            self.assertTrue(os.path.isdir(db_dir))

    def test_database_directory_with_custom_dir(self):
        s = Schemas(False)
        available = s.all_available()
        if len(available) > 0:
            species = available[0]
            # Using the actual base directory as custom dir for testing
            db_dir = s.database_directory(s.base_directory, species)
            self.assertIsNotNone(db_dir)

    def test_database_directory_invalid_species(self):
        s = Schemas(False)
        # Capture stdout to suppress error message
        captured_output = StringIO()
        old_stdout = sys.stdout
        try:
            sys.stdout = captured_output
            db_dir = s.database_directory(None, 'NonExistentSpecies_12345')
        finally:
            sys.stdout = old_stdout
        self.assertIsNone(db_dir)

