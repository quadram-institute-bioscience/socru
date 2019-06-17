import os
import csv
import re
from socru.GATProfile  import GATProfile

class Results:
    def __init__(self, results_file, verbose):
        self.results_file = results_file
        self.verbose = verbose
        self.profiles = self.create_profiles()
    
    def create_profiles(self):
        seen_profiles = []
        profiles = []
        with open(self.results_file, newline='') as csvfile:
            profile_reader = csv.reader(csvfile, delimiter='\t')
            for row in profile_reader:
                for row in profile_reader:
                    if len(row) > 3:
                        # 1: directory of schema
                        # 2: GS number with GS at the start
                        # 3..N: fragment pattern
                        m = re.match("GS(.+)", row[1])
                        if m:
                            gat_number = m.group(1)
                            fragments = [row[f] for f in range(2, len(row))]
                            unknown = [f for f in fragments if f == '?' or f == "?'"]
                            if len(unknown) > 0:
                                continue
                            
                            g = GATProfile(self.verbose, gat_number = gat_number, fragments = fragments)
                            if str(g) not in seen_profiles:
                                seen_profiles.append(str(g))
                                profiles.append(g)
        return profiles
    
    def filter(self, num_fragments):
        profiles_novel = self.filter_previously_seen_profiles(self.profiles)
        valid_profiles = self.all_fragments_present_once(num_fragments, profiles_novel)
        return valid_profiles
        
    def all_fragments_present_once(self, num_fragments, profiles):
        filtered_profiles = []
        for p in profiles:
            if len(p.fragments) != num_fragments:
                continue
            
            sorted_orientationless = sorted(p.orientationless_fragments())
            out_of_order = [ i+1 for i in range(0,len(sorted_orientationless)) if int(sorted_orientationless[i]) != i+1]
            if len(out_of_order) > 0:
                continue
            else:
                filtered_profiles.append(p)
            
        return filtered_profiles
        

    def filter_previously_seen_profiles(self, profiles):
        filtered_profiles = []
        for p in profiles:
            m = re.match(r'([\d]+)\.([\d]+)', p.gat_number)
            if m:
                if int(m.group(1)) == 0:
                     # not seen before so keep
                     filtered_profiles.append(p)
        return filtered_profiles
  