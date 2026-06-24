import filecmp
import os
import shutil
import tempfile
import unittest

from socru.ProfileGenerator import ProfileGenerator

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data', 'profile_generator')


class TestProfileGenerator(unittest.TestCase):

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        # Copy the database directory contents to the temp directory so that
        # output files are written there instead of polluting the source tree.
        self.tmp_db = os.path.join(self.tmpdir, 'database')
        shutil.copytree(os.path.join(data_dir, 'database'), self.tmp_db)

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_profile_generator(self):
        p = ProfileGenerator(
            self.tmp_db, 7,
            os.path.join(data_dir, 'dnaA.fa.gz'),
            os.path.join(data_dir, 'dif.fa.gz'),
            1, 'reffile.fa', False,
        )
        p.write_output_file()

        profile_path = os.path.join(self.tmp_db, 'profile.txt')
        yml_path = os.path.join(self.tmp_db, 'profile.txt.yml')

        self.assertTrue(os.path.exists(profile_path))
        self.assertTrue(
            filecmp.cmp(profile_path, os.path.join(data_dir, 'expected_profile.txt'))
        )
        self.assertTrue(os.path.exists(yml_path))
