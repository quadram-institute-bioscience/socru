import os
import sys
from socru.GATProfile  import GATProfile
from socru.Profiles import Profiles
from socru.TypeGenerator import TypeGenerator
from socru.ValidateFragments import ValidateFragments
from socru.Fragment import Fragment

class SocruLookup:
    def __init__(self,options):
        self.fragments = options.fragments
        self.db_dir = options.db_dir
        self.verbose = options.verbose
        
        self.tg = self.type_generator()
        
    def create_fragments(self, input_profile, standard_fragments):
        frags = []
        # get the basic details
        
        reversed_frags = input_profile.orientation_array()
        for i, frag_number in enumerate(input_profile.orientationless_fragments()):
            f = Fragment([], number = frag_number, operon_forward_start = False, dna_A = False, dif = False, reversed_frag = not reversed_frags[i])
            if str(frag_number) == str(input_profile.dnaA_fragment_number):
                f.dna_A = True
            elif str(frag_number) == str(input_profile.dif_fragment_number):
                f.dif = True
            frags.append(f)

        # loop over the standard and set the direction (if inverted, then invert)
        for f in frags:
            for s in standard_fragments:
                if f.number == s.number:
                    f.operon_forward_start = s.operon_forward_start
            if f.reversed_frag:
                f.operon_forward_start = not f.operon_forward_start
                
            if s.dna_A:
                f.operon_forward_start = True
            if s.dif:
                f.operon_forward_start = False
            
        return frags
  
    def find_dnaa(self, frags):
        for i, frag in enumerate(frags):
            if frag.dna_A:
                return i
        
        print("DnaA couldnt be found")
        return -1
        
    def find_dif(self, frags):
        for i, frag in enumerate(frags):
            if frag.dif:
                return i
        
        print("Dif couldnt be found")
        return -1
        
    def reorientate_ori(self, frags):
        dnaa_index = self.find_dnaa(frags)
        if dnaa_index == -1:
            return frags

        reorientated_frags = frags[dnaa_index: len(frags)] + frags[0: dnaa_index]
        return reorientated_frags
        
    # Take reference case 1 2 3 4 5 6 7
    # assign terminus and origin
    # assign direction to each number
    # dnaa points outwards
    # dif points inwards
    # This is the standard
    def create_standard_fragments(self, profile_db):
        reference_profile = profile_db.gats[0]
        
        ref_frags = []
        # get the basic details
        for frag_number in reference_profile.orientationless_fragments():
            f = Fragment([], number = frag_number, operon_forward_start = False, dna_A = False, dif = False)
            if str(frag_number) == str(profile_db.dnaA_fragment_number):
                f.dna_A = True
            elif str(frag_number) == str(profile_db.dif_fragment_number):
                f.dif = True
            ref_frags.append(f)
            
        # reorientate to dnaA
        dna_start_frags = self.reorientate_ori(ref_frags)
        dnaa_index = self.find_dnaa(dna_start_frags)
        dif_index = self.find_dif(dna_start_frags)
        if dnaa_index == -1 or dif_index == -1:
            return []
        
        # forward walk (reverse is already set )
        for i in range(dnaa_index, dif_index):
            if i >= dif_index:
                # we have reached the end
                continue
            else:
                dna_start_frags[i].operon_forward_start = True

        return dna_start_frags 
        
    def validate_profile(self, profile_db, input_profile):
        standard_fragments = self.create_standard_fragments(profile_db)
        frags = self.create_fragments( input_profile, standard_fragments)
        return ValidateFragments(frags).validate()
        
    def type_generator(self):
        profile_db = Profiles(os.path.join(self.db_dir, 'profile.txt'), self.verbose)

        split_fragments = self.fragments.split('-')
        input_profile = GATProfile(self.verbose, fragments = split_fragments, dnaA_fragment_number = profile_db.dnaA_fragment_number, dif_fragment_number = profile_db.dif_fragment_number)
        input_profile.orientate_for_dnaA()
        
        is_profile_valid = self.validate_profile(profile_db, input_profile)
        tg = TypeGenerator(profile_db,input_profile, self.verbose, is_profile_valid)
        return tg
        
    def calc_type(self):
        return self.tg.calculate_type()

    def calc_quality(self):
        return self.tg.quality
        
    def run(self):
        print(str(self.calc_quality())+"\t"+str(self.tg.calculate_type()))
        