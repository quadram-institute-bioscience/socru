import os
from os import listdir
from os.path import isdir
import pkg_resources
import yaml

class Schemas:
    def __init__(self, verbose):
        self.verbose = verbose
        self.base_directory = str(pkg_resources.resource_filename( __name__, 'data/'))
    
    def all_available(self):
        return [ d for d in listdir(self.base_directory) if isdir(os.path.join(self.base_directory, d))]
        
    def print_all(self):
        for d in sorted(self.all_available()):
            print(d)
            
    def extended(self):
        db_info = {}
        for species in listdir(self.base_directory):
            
            db_dir = os.path.join(self.base_directory, species)
            db_file = os.path.join(self.base_directory, species,'profile.txt.yml')

            if not os.path.exists(db_file):
                continue
                
            with open(db_file, 'r') as metadatafh:
                try:
                    metadata = yaml.load(metadatafh, yaml.SafeLoader)
                    dnaA_fragment_number = 0
                    dnaa_forward_orientation = False
                    dif_fragment_number = 0
                    reference_genome = ''
                    
                    if 'dnaa_fragment' in metadata:
                        dnaA_fragment_number = int(metadata['dnaa_fragment'])
                        dnaa_forward_orientation = metadata['dnaa_forward_orientation']
					
                    if 'dif_fragment' in metadata:
                        dif_fragment_number = int(metadata['dif_fragment'])
                        
                    if 'reference_genome' in metadata:
                        reference_genome = str(metadata['reference_genome'])
                        
                    db_info[species] = [species, dnaA_fragment_number, dnaa_forward_orientation, dif_fragment_number, reference_genome]
                    
					
                except yaml.YAMLError as exc:
                    print(exc)
        return db_info
        
    def print_extended(self):
        extended_details = self.extended()
        print("\t".join(['Species','dnaA fragment No.', 'dnaA forward orientation', 'dif fragment No.', 'Reference']))
        for species in sorted(extended_details):
            print("\t".join( [ str(a) for a in extended_details[species] ]))
        
        
    def database_directory(self, db_dir, species):
        proposed_db_dir = os.path.join(self.base_directory, species)
        if db_dir is not None:
            proposed_db_dir = os.path.join(db_dir, species)

        if isdir(proposed_db_dir):
            return proposed_db_dir
        else:
            print("Cannot access the database directory for the given species")
            return None
            

                
         