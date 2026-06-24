import os
import shutil
import tempfile
import unittest

from socru.Fragment import Fragment
from socru.PlotProfile import PlotProfile

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data', 'plot_profile')


class TestPlotProfile(unittest.TestCase):

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_plot_profile(self):
        fragments = []
        fragments.append(Fragment([], sequence="AAAA", number=1, reversed_frag=False, dna_A=False))
        fragments.append(Fragment([], sequence="AAA", number=3, reversed_frag=False, dna_A=False))
        fragments.append(Fragment([], sequence="AAAAAAAAAA", number=2, reversed_frag=False, dna_A=True))
        fragments.append(Fragment([], sequence="A", number=4, reversed_frag=False, dna_A=False))

        outfile = os.path.join(self.tmpdir, 'plot.pdf')
        p = PlotProfile(fragments, outfile, False)
        p.create_plot()
        self.assertTrue(os.path.exists(outfile))

    def test_plot_profile_some_reversed(self):
        fragments = []
        fragments.append(Fragment([], sequence="AAAA", number=1, reversed_frag=True, dna_A=True))
        fragments.append(Fragment([], sequence="AAA", number=3, reversed_frag=False, dna_A=False))
        fragments.append(Fragment([], sequence="AAAAAAAAAA", number=2, reversed_frag=True, dna_A=False))
        fragments.append(Fragment([], sequence="A", number=4, reversed_frag=False, dna_A=False))

        outfile = os.path.join(self.tmpdir, 'plot2.pdf')
        p = PlotProfile(fragments, outfile, False)
        p.create_plot()

        self.assertTrue(os.path.exists(outfile))
