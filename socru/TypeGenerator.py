from Bio import SeqIO
from socru.GATProfile  import GATProfile
import re
import gzip

class TypeGenerator:
    def __init__(self, profile_db, gat_profile, verbose, prefix = 'GS'):
        self.profile_db = profile_db
        self.gat_profile = gat_profile
        self.prefix = prefix
        self.has_previously_seen = False
        self.verbose = verbose
        
        self.gs_type = self.calculate_type()
        
    def find_order_orientationless(self, orientationless_fragment):
        for db_profile in self.profile_db.gats:
            orientationless_db_profile = GATProfile(self.verbose, fragments = db_profile.orientationless_fragments())
            if orientationless_db_profile.does_the_profile_match(orientationless_fragment):
                return db_profile.order()
        return 0
        
    def calculate_orientationless_order(self):
        orientationless_fragment = GATProfile(self.verbose, fragments = self.gat_profile.orientationless_fragments())
        
        order = self.find_order_orientationless(orientationless_fragment)
        if order > 0:
            return order
        
        # invert it
        inverted_orientationless_fragment = GATProfile(self.verbose, fragments = orientationless_fragment.inverted_orientationless_fragments())
        
        order = self.find_order_orientationless(inverted_orientationless_fragment)
        if order > 0:
            return order
        else:
            return 0

    def calculate_type(self):
        if self.gat_profile.gat_number == 0:
            self.gat_profile.gat_number = self.calculate_orientationless_order()
            
        
        # lookup the gat_profile to get the number
        for db_profile in self.profile_db.gats:
            if db_profile.does_the_profile_match(self.gat_profile):
                self.gat_profile.gat_number = db_profile.gat_number
                continue
                
        self.has_previously_seen = self.previously_seen(self.gat_profile.gat_number, self.gat_profile.orientation_binary())
        return self.create_gs_type(self.gat_profile.gat_number, self.gat_profile.orientation_binary())
    
    def previously_seen(self, gat_number, orientation_binary):
        m = re.search(r'([\d]+)\.([\d]+)$', str(gat_number) )
        
        if m:
            if str(m.group(1)) == str(0):
                return False
            
            if str(m.group(2)) == str(orientation_binary):
                return True
            else:
                return False
        else:
            return False
    
    def create_gs_type(self, gat_number, orientation_binary):
        m = re.search(r'([\d]+)\.([\d]+)$', str(gat_number) )
        
        if m:
            if str(m.group(2)) == str(orientation_binary):
                return self.prefix + str(gat_number) 
            else:
                return self.prefix + str(m.group(1)) + "." + str(orientation_binary)
                
        else:
            return self.prefix + str(gat_number) + "." + str(orientation_binary)
