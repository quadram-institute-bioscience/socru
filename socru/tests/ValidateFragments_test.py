import unittest
import os
import shutil
from socru.ValidateFragments import ValidateFragments
from socru.Fragment import Fragment

class TestValidateFragments(unittest.TestCase):
   
    def test_term_first_valid(self):
        # --> 1'(Ter) <-- 6' <-- 3 <-- 4 <-- 5(Ori) --> 2'
        fragments = []
        fragments.append(Fragment([], number = 1, operon_forward_start = True, dna_A = False, dif = True))
        fragments.append(Fragment([], number = 6, operon_forward_start = False, dna_A = False, dif = False))
        fragments.append(Fragment([], number = 3, operon_forward_start = False, dna_A = False, dif = False))    
        fragments.append(Fragment([], number = 4, operon_forward_start = False, dna_A = False, dif = False))         
        fragments.append(Fragment([], number = 5, operon_forward_start = False, dna_A = True, dif = False))   
        fragments.append(Fragment([], number = 2, operon_forward_start = True, dna_A = False, dif = False)) 
        
        vf = ValidateFragments(fragments)
        self.assertTrue(vf.validate())

    def test_valid_ef(self):
        #E. faecium: --> 1(Ter) <-- 2 <-- 3 <-- 4 <-- 5(Ori) --> 6
        fragments = []
        fragments.append(Fragment([], number = 1, operon_forward_start = True, dna_A = False, dif = True))
        fragments.append(Fragment([], number = 2, operon_forward_start = False, dna_A = False, dif = False))
        fragments.append(Fragment([], number = 3, operon_forward_start = False, dna_A = False, dif = False))    
        fragments.append(Fragment([], number = 4, operon_forward_start = False, dna_A = False, dif = False))         
        fragments.append(Fragment([], number = 5, operon_forward_start = False, dna_A = True, dif = False))   
        fragments.append(Fragment([], number = 6, operon_forward_start = True, dna_A = False, dif = False)) 
        vf = ValidateFragments(fragments)
        self.assertTrue(vf.validate())
        
    def test_valid_ec(self):
        #E. coli:   --> 1'(Ter) <-- 4' <-- 5' <-- 2 <-- 3(Ori) --> 6 --> 7
        fragments = []
        fragments.append(Fragment([], number = 1, operon_forward_start = True, dna_A = False, dif = True))
        fragments.append(Fragment([], number = 4, operon_forward_start = False, dna_A = False, dif = False))
        fragments.append(Fragment([], number = 5, operon_forward_start = False, dna_A = False, dif = False))    
        fragments.append(Fragment([], number = 2, operon_forward_start = False, dna_A = False, dif = False))         
        fragments.append(Fragment([], number = 3, operon_forward_start = False, dna_A = True, dif = False))   
        fragments.append(Fragment([], number = 6, operon_forward_start = True, dna_A = False, dif = False)) 
        fragments.append(Fragment([], number = 7, operon_forward_start = True, dna_A = False, dif = False)) 
        vf = ValidateFragments(fragments)
        self.assertTrue(vf.validate())
 
    def test_valid_kp(self):
        #K. pneumoniae:  --> 1(Ter) <-- 8 <-- 7 <-- 6 <-- 5 <-- 4 <-- 3(Ori) --> 2
        fragments = []
        fragments.append(Fragment([], number = 1, operon_forward_start = True, dna_A = False, dif = True))
        fragments.append(Fragment([], number = 8, operon_forward_start = False, dna_A = False, dif = False))
        fragments.append(Fragment([], number = 7, operon_forward_start = False, dna_A = False, dif = False))    
        fragments.append(Fragment([], number = 6, operon_forward_start = False, dna_A = False, dif = False))         
        fragments.append(Fragment([], number = 5, operon_forward_start = False, dna_A = False, dif = False))   
        fragments.append(Fragment([], number = 4, operon_forward_start = False, dna_A = False, dif = False)) 
        fragments.append(Fragment([], number = 3, operon_forward_start = False, dna_A = True, dif = False)) 
        fragments.append(Fragment([], number = 2, operon_forward_start = True, dna_A = False, dif = False)) 
        vf = ValidateFragments(fragments)
        self.assertTrue(vf.validate())
 
    def test_valid_se(self):
        #S. enterica: --> 1'(Ter) <-- 2' <-- 6 <-- 5 <-- 4 <-- 3(Ori) --> 7'
        fragments = []
        fragments.append(Fragment([], number = 1, operon_forward_start = True, dna_A = False, dif = True))
        fragments.append(Fragment([], number = 2, operon_forward_start = False, dna_A = False, dif = False))
        fragments.append(Fragment([], number = 6, operon_forward_start = False, dna_A = False, dif = False))    
        fragments.append(Fragment([], number = 5, operon_forward_start = False, dna_A = False, dif = False))         
        fragments.append(Fragment([], number = 4, operon_forward_start = False, dna_A = False, dif = False))   
        fragments.append(Fragment([], number = 3, operon_forward_start = False, dna_A = True, dif = False)) 
        fragments.append(Fragment([], number = 7, operon_forward_start = True, dna_A = False, dif = False)) 
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
        
    def test_invalid_operons_but_flagged_as_valid_ef(self):
        #E. faecium: Aus0004 <-- 1(Ter) <-- 2 <-- 3 <-- 4 <-- 5(Ori) --> 6
        fragments = []
        fragments.append(Fragment([], number = 1, operon_forward_start = False, dna_A = False, dif = True))
        fragments.append(Fragment([], number = 2, operon_forward_start = False, dna_A = False, dif = False))
        fragments.append(Fragment([], number = 3, operon_forward_start = False, dna_A = False, dif = False))    
        fragments.append(Fragment([], number = 4, operon_forward_start = False, dna_A = False, dif = False))         
        fragments.append(Fragment([], number = 5, operon_forward_start = False, dna_A = True, dif = False))   
        fragments.append(Fragment([], number = 6, operon_forward_start = True, dna_A = False, dif = False)) 
        vf = ValidateFragments(fragments)
        self.assertFalse(vf.validate())
        
    def test_invalid_operons_but_flagged_as_valid_ec(self):
        #E. coli: IAI39  <-- 1'(Ter) <-- 4' <-- 5' <-- 2 <-- 3(Ori) --> 6 --> 7
        fragments = []
        fragments.append(Fragment([], number = 1, operon_forward_start = False, dna_A = False, dif = True))
        fragments.append(Fragment([], number = 4, operon_forward_start = False, dna_A = False, dif = False))
        fragments.append(Fragment([], number = 5, operon_forward_start = False, dna_A = False, dif = False))    
        fragments.append(Fragment([], number = 2, operon_forward_start = False, dna_A = False, dif = False))         
        fragments.append(Fragment([], number = 3, operon_forward_start = False, dna_A = True, dif = False))   
        fragments.append(Fragment([], number = 6, operon_forward_start = True, dna_A = False, dif = False)) 
        fragments.append(Fragment([], number = 7, operon_forward_start = True, dna_A = False, dif = False)) 
        vf = ValidateFragments(fragments)
        self.assertFalse(vf.validate())
 
    def test_invalid_operons_but_flagged_as_valid_kp(self):
        #K. pneumoniae: WCHKP015625 <-- 1(Ter) <-- 8 <-- 7 <-- 6 <-- 5 <-- 4 <-- 3(Ori) --> 2
        fragments = []
        fragments.append(Fragment([], number = 1, operon_forward_start = False, dna_A = False, dif = True))
        fragments.append(Fragment([], number = 8, operon_forward_start = False, dna_A = False, dif = False))
        fragments.append(Fragment([], number = 7, operon_forward_start = False, dna_A = False, dif = False))    
        fragments.append(Fragment([], number = 6, operon_forward_start = False, dna_A = False, dif = False))         
        fragments.append(Fragment([], number = 5, operon_forward_start = False, dna_A = False, dif = False))   
        fragments.append(Fragment([], number = 4, operon_forward_start = False, dna_A = False, dif = False)) 
        fragments.append(Fragment([], number = 3, operon_forward_start = False, dna_A = True, dif = False)) 
        fragments.append(Fragment([], number = 2, operon_forward_start = True, dna_A = False, dif = False)) 
        vf = ValidateFragments(fragments)
        self.assertFalse(vf.validate())
 
    def test_invalid_operons_but_flagged_as_valid_se(self):
        #S. enterica: FORC_007 - <-- 1'(Ter) <-- 2' <-- 6 <-- 5 <-- 4 <-- 3(Ori) --> 7'
        fragments = []
        fragments.append(Fragment([], number = 1, operon_forward_start = False, dna_A = False, dif = True))
        fragments.append(Fragment([], number = 2, operon_forward_start = False, dna_A = False, dif = False))
        fragments.append(Fragment([], number = 6, operon_forward_start = False, dna_A = False, dif = False))    
        fragments.append(Fragment([], number = 5, operon_forward_start = False, dna_A = False, dif = False))         
        fragments.append(Fragment([], number = 4, operon_forward_start = False, dna_A = False, dif = False))   
        fragments.append(Fragment([], number = 3, operon_forward_start = False, dna_A = True, dif = False)) 
        fragments.append(Fragment([], number = 7, operon_forward_start = True, dna_A = False, dif = False)) 
        vf = ValidateFragments(fragments)
        self.assertFalse(vf.validate())
    

        