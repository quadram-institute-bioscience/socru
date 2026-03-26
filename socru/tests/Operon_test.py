import unittest

from socru.Operon import Operon


class TestOperon(unittest.TestCase):

    def test_creation_forward(self):
        o = Operon(100, 5000, True)
        self.assertEqual(o.start, 100)
        self.assertEqual(o.end, 5000)
        self.assertTrue(o.direction)

    def test_creation_reverse(self):
        o = Operon(5000, 100, False)
        self.assertEqual(o.start, 5000)
        self.assertEqual(o.end, 100)
        self.assertFalse(o.direction)

    def test_str_forward(self):
        o = Operon(100, 5000, True)
        self.assertEqual(str(o), "100\t5000\tTrue")

    def test_str_reverse(self):
        o = Operon(100, 5000, False)
        self.assertEqual(str(o), "100\t5000\tFalse")

    def test_str_tab_delimited_fields(self):
        o = Operon(42, 999, True)
        fields = str(o).split("\t")
        self.assertEqual(len(fields), 3)
        self.assertEqual(fields[0], "42")
        self.assertEqual(fields[1], "999")
        self.assertEqual(fields[2], "True")

    def test_attributes_mutable(self):
        o = Operon(100, 200, True)
        o.start = 150
        o.direction = False
        self.assertEqual(o.start, 150)
        self.assertFalse(o.direction)
