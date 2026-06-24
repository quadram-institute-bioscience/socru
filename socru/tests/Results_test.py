import os
import tempfile
import unittest

from socru.Results import Results

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data','results')

class TestResults(unittest.TestCase):

    def test_results(self):
        r = Results(os.path.join(data_dir, 'yersinia'), False)

        valid_patterns  = [str(p)for p in r.filter(7)]
        self.assertEqual(valid_patterns, ["1'\t2\t3\t4\t5\t6\t7",
                                          "1'\t7'\t6'\t5'\t3\t4\t2'",
                                          "1'\t7'\t3\t4\t5\t6\t2'",
                                          "1\t7'\t3\t4\t5\t6\t2'"])

    def test_all_rows_processed(self):
        """Verify that all rows in a results file are processed, not just the first."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as f:
            # Write 3 distinct novel profiles, each with 3 fragments
            f.write("dir1\tGS0.0\t1\t2\t3\n")
            f.write("dir2\tGS0.0\t1\t3\t2\n")
            f.write("dir3\tGS0.0\t1\t2'\t3\n")
            tmpfile = f.name
        try:
            r = Results(tmpfile, False)
            self.assertEqual(len(r.profiles), 3)
        finally:
            os.unlink(tmpfile)


