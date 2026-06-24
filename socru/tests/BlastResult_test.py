import unittest

from socru.BlastResult import BlastResult


class TestBlastResult(unittest.TestCase):

    def _make_result(self, subject_start=2221, subject_end=4740):
        """Helper to create a BlastResult with sensible defaults."""
        return BlastResult(
            query_name="query7",
            subject="7",
            identity="100.000",
            alignment_length="2520",
            mismatches="0",
            gap_openings="0",
            query_start="1",
            query_end="2520",
            subject_start=str(subject_start),
            subject_end=str(subject_end),
            e_value="0.0",
            bit_score="4654",
        )

    def test_parse_fields_from_strings(self):
        br = self._make_result()
        self.assertEqual(br.query_name, "query7")
        self.assertEqual(br.subject, "7")
        self.assertAlmostEqual(br.identity, 100.0)
        self.assertEqual(br.alignment_length, 2520)
        self.assertEqual(br.mismatches, 0)
        self.assertEqual(br.gap_openings, 0)
        self.assertEqual(br.query_start, 1)
        self.assertEqual(br.query_end, 2520)
        self.assertEqual(br.subject_start, 2221)
        self.assertEqual(br.subject_end, 4740)
        self.assertAlmostEqual(br.e_value, 0.0)
        self.assertAlmostEqual(br.bit_score, 4654.0)

    def test_parse_from_tab_delimited_line(self):
        line = "query7\t7\t100.000\t2520\t0\t0\t1\t2520\t2221\t4740\t0.0\t4654"
        fields = line.split("\t")
        br = BlastResult(*fields)
        self.assertEqual(br.query_name, "query7")
        self.assertEqual(br.alignment_length, 2520)
        self.assertEqual(br.subject_start, 2221)
        self.assertEqual(br.subject_end, 4740)

    def test_is_forward_true_when_start_less_than_end(self):
        br = self._make_result(subject_start=100, subject_end=500)
        self.assertTrue(br.is_forward())

    def test_is_forward_false_when_start_greater_than_end(self):
        br = self._make_result(subject_start=4740, subject_end=2221)
        self.assertFalse(br.is_forward())

    def test_is_forward_true_when_equal(self):
        br = self._make_result(subject_start=100, subject_end=100)
        self.assertTrue(br.is_forward())

    def test_str_produces_tab_delimited_output(self):
        br = self._make_result()
        output = str(br)
        fields = output.split("\t")
        self.assertEqual(len(fields), 12)
        self.assertEqual(fields[0], "query7")
        self.assertEqual(fields[1], "7")
        self.assertEqual(fields[8], "2221")
        self.assertEqual(fields[9], "4740")

    def test_str_roundtrip(self):
        """Parsing a str() output should yield equivalent values."""
        br1 = self._make_result()
        fields = str(br1).split("\t")
        br2 = BlastResult(*fields)
        self.assertEqual(br2.query_name, br1.query_name)
        self.assertEqual(br2.alignment_length, br1.alignment_length)
        self.assertEqual(br2.subject_start, br1.subject_start)
        self.assertEqual(br2.subject_end, br1.subject_end)

    def test_type_conversions(self):
        """Fields passed as strings should be converted to correct types."""
        br = BlastResult("q", "s", "99.5", "100", "2", "1", "10", "110", "200", "300", "1e-50", "200.5")
        self.assertIsInstance(br.identity, float)
        self.assertIsInstance(br.alignment_length, int)
        self.assertIsInstance(br.mismatches, int)
        self.assertIsInstance(br.gap_openings, int)
        self.assertIsInstance(br.query_start, int)
        self.assertIsInstance(br.query_end, int)
        self.assertIsInstance(br.subject_start, int)
        self.assertIsInstance(br.subject_end, int)
        self.assertIsInstance(br.e_value, float)
        self.assertIsInstance(br.bit_score, float)
