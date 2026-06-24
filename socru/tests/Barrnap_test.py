import os
import unittest

from socru.Barrnap import Barrnap

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data','barrnap')


class TestBarrnap(unittest.TestCase):

    def test_barrnap_normal(self):
        barrnap = Barrnap('xxx', 1, False)
        b = barrnap.read_barrnap_output( os.path.join(data_dir,'barrnap_output.txt'))
        boundries = []
        for bobj in b:
            boundries.append([bobj.start,bobj.end, bobj.direction])
        self.assertEqual(boundries, [[36739, 42014, False], [185729, 190843, False], [698429, 703543, True], [1535715, 1540829, True], [3935713, 3940989, False], [4651665, 4657321, False]])

    def test_barrnap_70s_over_end(self):
        barrnap = Barrnap(os.path.join(data_dir, 'test.fa'), 1, False)
        b = barrnap.read_barrnap_output( os.path.join(data_dir,'70s_over_end.txt'))
        boundries = []
        for bobj in b:
            boundries.append([bobj.start,bobj.end, bobj.direction])
        self.assertEqual(boundries,  [[2804, 5721, True], [9199, 1889, True]])

    def test_repeated_5s(self):
        barrnap = Barrnap(os.path.join(data_dir, 'test.fa'), 1, False)
        b = barrnap.read_barrnap_output( os.path.join(data_dir,'repeated_5S.txt'))
        boundries = []
        for bobj in b:
            boundries.append([bobj.start,bobj.end, bobj.direction])
        self.assertEqual(boundries,  [[1253, 6909, True]])

    def test_orphan_5s(self):
        barrnap = Barrnap('xxx', 1, False, chromosome_length = 2500000)
        b = barrnap.read_barrnap_output( os.path.join(data_dir,'orphan_5S.txt'))
        boundries = []
        for bobj in b:
            boundries.append([bobj.start,bobj.end, bobj.direction])
        self.assertEqual(boundries, [[804628, 809557, False], [809806, 814769, False], [854057, 859048, False], [1947088, 1952017, True], [2063398, 2068327, True], [2241390, 2246490, True]] )
