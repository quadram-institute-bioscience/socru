import unittest
import os
import shutil
from socru.ValidateFragments import ValidateFragments
from socru.Fragment import Fragment

class TestValidateFragments(unittest.TestCase):
   
    def test_term_first_valid(self):
        # --> 1'(Ter) <-- 6' <-- 3 <-- 4 <-- 5(Ori) --> 2'
        fragments = []
        fragments.append(Fragment([], number = 1, operon_forward_start = False, dna_A = False, dif = True))
        fragments.append(Fragment([], number = 6, operon_forward_start = False, dna_A = False, dif = False))
        fragments.append(Fragment([], number = 3, operon_forward_start = False, dna_A = False, dif = False))    
        fragments.append(Fragment([], number = 4, operon_forward_start = False, dna_A = False, dif = False))         
        fragments.append(Fragment([], number = 5, operon_forward_start = True, dna_A = True, dif = False))   
        fragments.append(Fragment([], number = 2, operon_forward_start = True, dna_A = False, dif = False)) 
        
        vf = ValidateFragments(fragments)
        self.assertTrue(vf.validate())
        
    def test_term_first_invalid(self):
        # --> 1'(Ter) <-- 6' <-- 3 --> 4 <-- 5(Ori) --> 2'
        fragments = []
        fragments.append(Fragment([], number = 1, operon_forward_start = False, dna_A = False, dif = True))
        fragments.append(Fragment([], number = 6, operon_forward_start = False, dna_A = False, dif = False))
        fragments.append(Fragment([], number = 3, operon_forward_start = True, dna_A = False, dif = False))    
        fragments.append(Fragment([], number = 4, operon_forward_start = False, dna_A = False, dif = False))         
        fragments.append(Fragment([], number = 5, operon_forward_start = True, dna_A = True, dif = False))   
        fragments.append(Fragment([], number = 2, operon_forward_start = True, dna_A = False, dif = False)) 
        
        vf = ValidateFragments(fragments)
        self.assertFalse(vf.validate())

        