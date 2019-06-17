import unittest
import os
import shutil
from socru.Schemas  import Schemas

class TestSchemas(unittest.TestCase):
   
    def test_schemas(self):
        s = Schemas(False)
        self.assertTrue( len(s.all_available()) > 4)

        