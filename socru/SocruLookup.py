import os
import sys
from socru.GATProfile  import GATProfile
from socru.Profiles import Profiles
from socru.TypeGenerator import TypeGenerator

class SocruLookup:
    def __init__(self,options):
        self.fragments = options.fragments
        self.db_dir = options.db_dir
        
    def calc_type(self):  
        profile_db = Profiles(os.path.join(self.db_dir, 'profile.txt'))
        
        split_fragments = self.fragments.split('-')
        input_profile = GATProfile(fragments = split_fragments)
        
        tg = TypeGenerator(profile_db,input_profile)
        return tg.calculate_type()

    def run(self):
        print(calc_type())
        