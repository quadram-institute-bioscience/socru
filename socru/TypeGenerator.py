from Bio import SeqIO
from socru.Fragment import Fragment
from socru.GATProfile  import GATProfile
import re
import gzip

class TypeGenerator:
    def __init__(self, profile_db, gat_profile, prefix = 'GS'):
        self.profile_db = profile_db
        self.gat_profile = gat_profile
        self.prefix = prefix
        self.has_previously_seen = False
        
        self.gs_type = self.calculate_type()
        
    def calculate_orientationless_order(self):
        orientationless_fragment = GATProfile(fragments = self.gat_profile.orientationless_fragments())
        
        for db_profile in self.profile_db.gats:
            orientationless_db_profile = GATProfile(fragments = db_profile.orientationless_fragments())
            if orientationless_db_profile.does_the_profile_match(orientationless_fragment):
                return db_profile.order()
        return 0
        
    def calculate_type(self):
        # lookup the gat_profile to get the number
        for db_profile in self.profile_db.gats:
            if db_profile.does_the_profile_match(self.gat_profile):
                self.gat_profile.gat_number = db_profile.gat_number
                continue
                
        self.has_previously_seen = self.previously_seen(self.gat_profile.gat_number, self.gat_profile.orientation_binary())
        return self.create_gs_type(self.gat_profile.gat_number, self.gat_profile.orientation_binary())
    
    def previously_seen(self, gat_number, orientation_binary):
        m = re.search('([\d]+)\.([\d]+)$', str(gat_number) )
        
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
        m = re.search('([\d]+)\.([\d]+)$', str(gat_number) )
        
        if m:
            if str(m.group(2)) == str(orientation_binary):
                return self.prefix + str(gat_number) 
            else:
                return self.prefix + str(m.group(1)) + "." + str(orientation_binary)
                
        else:
            return self.prefix + str(gat_number) + "." + str(orientation_binary)
