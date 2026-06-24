import os
import shutil
import tempfile
import unittest
from unittest import mock

from socru.DatabaseManager import DEFAULT_DATA_DIR, DatabaseManager


class TestDatabaseManagerBundled(unittest.TestCase):
    """Tests that interact with the real bundled data directory."""

    def test_list_species_returns_bundled(self):
        dm = DatabaseManager(data_dir=tempfile.mkdtemp())
        species = dm.list_species(include_bundled=True)
        self.assertIsInstance(species, list)
        self.assertTrue(len(species) > 4, "Expected multiple bundled species")

    def test_list_species_sorted(self):
        dm = DatabaseManager(data_dir=tempfile.mkdtemp())
        species = dm.list_species()
        self.assertEqual(species, sorted(species))

    def test_get_database_dir_finds_bundled(self):
        dm = DatabaseManager(data_dir=tempfile.mkdtemp())
        species = dm.list_species()
        self.assertTrue(len(species) > 0)
        db_dir = dm.get_database_dir(species[0])
        self.assertIsNotNone(db_dir)
        self.assertTrue(os.path.isdir(db_dir))

    def test_get_database_dir_returns_none_for_unknown(self):
        dm = DatabaseManager(data_dir=tempfile.mkdtemp())
        result = dm.get_database_dir('Nonexistent_species_xyz_99999')
        self.assertIsNone(result)

    def test_database_info_returns_expected_structure(self):
        dm = DatabaseManager(data_dir=tempfile.mkdtemp())
        species = dm.list_species()
        self.assertTrue(len(species) > 0)
        info = dm.database_info(species[0])
        self.assertIsNotNone(info)
        self.assertIn('species', info)
        self.assertIn('path', info)
        self.assertIn('fragment_count', info)
        self.assertIn('has_profile', info)
        self.assertIn('is_bundled', info)
        self.assertEqual(info['species'], species[0])

    def test_database_info_returns_none_for_unknown(self):
        dm = DatabaseManager(data_dir=tempfile.mkdtemp())
        info = dm.database_info('Nonexistent_species_xyz_99999')
        self.assertIsNone(info)

    def test_database_info_known_types_present_when_profile_exists(self):
        dm = DatabaseManager(data_dir=tempfile.mkdtemp())
        species = dm.list_species()
        for sp in species:
            info = dm.database_info(sp)
            if info and info['has_profile']:
                self.assertIn('known_types', info)
                self.assertIsInstance(info['known_types'], int)
                break


class TestDatabaseManagerUserDir(unittest.TestCase):
    """Tests for user-installed databases and install functionality."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.user_data_dir = os.path.join(self.tmpdir, 'user_data')
        self.source_dir = os.path.join(self.tmpdir, 'source_db')
        os.makedirs(self.source_dir)
        # Create fake database files
        with open(os.path.join(self.source_dir, 'fragment1.fa'), 'w') as f:
            f.write('>frag1\nATCG\n')
        with open(os.path.join(self.source_dir, 'fragment2.fa'), 'w') as f:
            f.write('>frag2\nGCTA\n')
        with open(os.path.join(self.source_dir, 'profile.txt'), 'w') as f:
            f.write('# header\n')
            f.write('type1\t1\t2\t3\n')
            f.write('type2\t3\t2\t1\n')

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_install_database_copies_files(self):
        dm = DatabaseManager(data_dir=self.user_data_dir)
        dest = dm.install_database(self.source_dir, species_name='Test_species')
        self.assertTrue(os.path.isdir(dest))
        self.assertTrue(os.path.exists(os.path.join(dest, 'fragment1.fa')))
        self.assertTrue(os.path.exists(os.path.join(dest, 'fragment2.fa')))
        self.assertTrue(os.path.exists(os.path.join(dest, 'profile.txt')))

    def test_install_database_default_name(self):
        dm = DatabaseManager(data_dir=self.user_data_dir)
        dest = dm.install_database(self.source_dir)
        self.assertEqual(os.path.basename(dest), 'source_db')

    def test_get_database_dir_prefers_user_dir(self):
        dm = DatabaseManager(data_dir=self.user_data_dir)
        dm.install_database(self.source_dir, species_name='Test_species')
        db_dir = dm.get_database_dir('Test_species')
        self.assertIsNotNone(db_dir)
        self.assertTrue(db_dir.startswith(self.user_data_dir))

    def test_list_species_includes_user_installed(self):
        dm = DatabaseManager(data_dir=self.user_data_dir)
        dm.install_database(self.source_dir, species_name='Test_species')
        species = dm.list_species(include_bundled=False)
        self.assertIn('Test_species', species)

    def test_list_species_excludes_hidden_dirs(self):
        dm = DatabaseManager(data_dir=self.user_data_dir)
        os.makedirs(os.path.join(self.user_data_dir, '.hidden'))
        dm.install_database(self.source_dir, species_name='Visible_species')
        species = dm.list_species(include_bundled=False)
        self.assertIn('Visible_species', species)
        self.assertNotIn('.hidden', species)

    def test_database_info_for_installed_db(self):
        dm = DatabaseManager(data_dir=self.user_data_dir)
        dm.install_database(self.source_dir, species_name='Test_species')
        info = dm.database_info('Test_species')
        self.assertIsNotNone(info)
        self.assertEqual(info['species'], 'Test_species')
        self.assertEqual(info['fragment_count'], 2)
        self.assertTrue(info['has_profile'])
        self.assertEqual(info['known_types'], 2)


class TestDatabaseManagerEnvVar(unittest.TestCase):
    """Test SOCRU_DATA_DIR environment variable override."""

    def test_env_var_override(self):
        custom_dir = '/tmp/socru_test_custom_data'
        with mock.patch.dict(os.environ, {'SOCRU_DATA_DIR': custom_dir}):
            dm = DatabaseManager()
            self.assertEqual(dm.data_dir, custom_dir)

    def test_explicit_data_dir_takes_precedence_over_env(self):
        custom_dir = '/tmp/socru_test_custom_data'
        explicit_dir = '/tmp/socru_test_explicit'
        with mock.patch.dict(os.environ, {'SOCRU_DATA_DIR': custom_dir}):
            dm = DatabaseManager(data_dir=explicit_dir)
            self.assertEqual(dm.data_dir, explicit_dir)

    def test_default_data_dir_without_env(self):
        with mock.patch.dict(os.environ, {}, clear=True):
            # Remove SOCRU_DATA_DIR if present
            os.environ.pop('SOCRU_DATA_DIR', None)
            dm = DatabaseManager()
            self.assertEqual(dm.data_dir, DEFAULT_DATA_DIR)
