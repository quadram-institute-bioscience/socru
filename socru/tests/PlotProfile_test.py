import unittest
import os
import shutil
import filecmp
from socru.PlotProfile  import PlotProfile
from socru.Fragment import Fragment

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data','plot_profile')

class TestPlotProfile(unittest.TestCase):
   
    def test_plot_profile(self):
        fragments = []
        fragments.append(Fragment([], sequence = "AAAA", number = 1, reversed_frag = False, dna_A = False))
        fragments.append(Fragment([], sequence = "AAA", number = 3, reversed_frag = False, dna_A = False))
        fragments.append(Fragment([], sequence = "AAAAAAAAAA", number = 2, reversed_frag = False, dna_A = True))
        fragments.append(Fragment([], sequence = "A", number = 4, reversed_frag = False, dna_A = False))

        p = PlotProfile(fragments, 'plot.pdf', False)
        p.create_plot()
        self.assertTrue(os.path.exists('plot.pdf'))
        os.remove('plot.pdf')
        
    def test_plot_profile_some_reversed(self):
        fragments = []
        fragments.append(Fragment([], sequence = "AAAA", number = 1, reversed_frag = True, dna_A = True))
        fragments.append(Fragment([], sequence = "AAA", number = 3, reversed_frag = False, dna_A = False))
        fragments.append(Fragment([], sequence = "AAAAAAAAAA", number = 2, reversed_frag = True, dna_A = False))
        fragments.append(Fragment([], sequence = "A", number = 4, reversed_frag = False, dna_A = False))

        p = PlotProfile(fragments, 'plot2.pdf', False)
        p.create_plot()
        
        self.assertTrue(os.path.exists('plot2.pdf'))
        os.remove('plot2.pdf')
    
        

