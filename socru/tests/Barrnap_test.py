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
        self.assertEqual(b, [[36740, 42014], [185730, 190843], [698430, 703543], [1535716, 1540829], [3935714, 3940989], [4651666, 4657321]])


    def test_barrnap_70s_over_end(self):
        barrnap = Barrnap(os.path.join(data_dir, 'test.fa'), 1, False)
        b = barrnap.read_barrnap_output( os.path.join(data_dir,'70s_over_end.txt'))
        self.assertEqual(b,  [[2805, 5721], [9200, 1889]])
        
    def test_repeated_5s(self):
        barrnap = Barrnap(os.path.join(data_dir, 'test.fa'), 1, False)
        b = barrnap.read_barrnap_output( os.path.join(data_dir,'repeated_5S.txt'))
        self.assertEqual(b,  [[1254, 6909]])
		
    def test_orphan_5s(self):
        barrnap = Barrnap('xxx', 1, False, chromosome_length = 2500000)
        b = barrnap.read_barrnap_output( os.path.join(data_dir,'orphan_5S.txt'))
        self.assertEqual(b,  [[804629, 809557],[809807, 814769],[854058, 859048],[1947089, 1952017],[2063399, 2068327],[2241391, 2246490]]  )
        