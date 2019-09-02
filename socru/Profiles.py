import csv
import yaml
import os
from socru.GATProfile import GATProfile


class Profiles:
    def __init__(self, input_file, verbose, metadata_file_suffix = '.yml', default_dnaA_fragment_number = 3, default_dif_fragment_number = 1):
        self.input_file = input_file
        self.metadata_file = self.input_file + metadata_file_suffix
        self.dnaA_fragment_number = default_dnaA_fragment_number
        self.dif_fragment_number = default_dif_fragment_number
        self.dnaa_forward_orientation = False
        self.verbose = verbose
        
        self.gats = self.read_profiles()
        self.read_metadata()
        self.num_fragments = self.expected_num_fragments()
        
    def read_profiles(self):
        with open(self.input_file, newline='') as csvfile:
            profile_reader = csv.reader(csvfile, delimiter='\t')
            # skip the header
            next(profile_reader) 
            profiles = []
            for row in profile_reader:
                if len(row) > 2:
                    fragments = [row[f] for f in range(1, len(row)) if row[f] != '']
                    g = GATProfile(self.verbose, gat_number = row[0], fragments = fragments)
                    g.orientate_for_dnaA()
                    profiles.append(g)
            return profiles
            
    def read_metadata(self):
        if not os.path.exists(self.metadata_file):
            return
        with open(self.metadata_file, 'r') as metadatafh:
            try:
                metadata = yaml.load(metadatafh)
                if 'dnaa_fragment' in metadata:
                    self.dnaA_fragment_number = int(metadata['dnaa_fragment'])
                    self.dnaa_forward_orientation = metadata['dnaa_forward_orientation']
					
                if 'dif_fragment' in metadata:
                    self.dif_fragment_number = int(metadata['dif_fragment'])
					
            except yaml.YAMLError as exc:
                print(exc)
        return
        
    def expected_num_fragments(self):
        if len(self.gats)> 0:
            return len(self.gats[0].fragments)
        else:
            return -1
            
    
    def next_order_number(self):
        return max([p.order() for p in self.gats]) + 1
            