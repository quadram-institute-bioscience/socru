import unittest
import os
import shutil
from socru.Barrnap  import Barrnap

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data','barrnap')


class TestBarrnap(unittest.TestCase):
   
    def test_barrnap_normal(self):
        barrnap = Barrnap('xxx', 1, False)
        b = barrnap.read_barrnap_output( os.path.join(data_dir,'barrnap_output.txt'))
        boundries = []
        for bobj in b:
            boundries.append([bobj.start,bobj.end, bobj.direction])
        self.assertEqual(boundries, [[36740, 42014, False], [185730, 190843, False], [698430, 703543, True], [1535716, 1540829, True], [3935714, 3940989, False], [4651666, 4657321, False]])

    def test_barrnap_70s_over_end(self):
        barrnap = Barrnap(os.path.join(data_dir, 'test.fa'), 1, False)
        b = barrnap.read_barrnap_output( os.path.join(data_dir,'70s_over_end.txt'))
        boundries = []
        for bobj in b:
            boundries.append([bobj.start,bobj.end, bobj.direction])
        self.assertEqual(boundries,  [[2805, 5721, True], [9200, 1889, True]])
        
    def test_repeated_5s(self):
        barrnap = Barrnap(os.path.join(data_dir, 'test.fa'), 1, False)
        b = barrnap.read_barrnap_output( os.path.join(data_dir,'repeated_5S.txt'))
        boundries = []
        for bobj in b:
            boundries.append([bobj.start,bobj.end, bobj.direction])
        self.assertEqual(boundries,  [[1254, 6909, True]])
		
    def test_orphan_5s(self):
        barrnap = Barrnap('xxx', 1, False, chromosome_length = 2500000)
        b = barrnap.read_barrnap_output( os.path.join(data_dir,'orphan_5S.txt'))
        boundries = []
        for bobj in b:
            boundries.append([bobj.start,bobj.end, bobj.direction])
        self.assertEqual(boundries, [[804629, 809557, False], [809807, 814769, False], [854058, 859048, False], [1947089, 1952017, True], [2063399, 2068327, True], [2241391, 2246490, True]] )
        