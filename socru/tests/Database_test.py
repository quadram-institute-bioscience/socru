import unittest
import os
import shutil
from socru.Database  import Database

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data','database')

class TestDatabase(unittest.TestCase):
   
    def test_database(self):
        d = Database( data_dir)
        self.assertEqual(2, len(d.get_database_files()))
        self.assertEqual(0, len(d.get_database_files_compressed()))
        self.assertTrue(os.path.exists(d.db_prefix+'.nhr'))
        
    def test_database_compressed_inputs(self):
        d = Database( os.path.join(data_dir,'compressed'))
        self.assertEqual(2, len(d.get_database_files_compressed()))
        self.assertEqual(0, len(d.get_database_files()))
        self.assertTrue(os.path.exists(d.db_prefix+'.nhr'))
