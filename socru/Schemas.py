import os
from os import listdir
from os.path import isdir
import pkg_resources

class Schemas:
    def __init__(self):
        self.base_directory = str(pkg_resources.resource_filename( __name__, 'data/'))
    
    def all_available(self):
        return [ d for d in listdir(self.base_directory) if isdir(os.path.join(self.base_directory, d))]
        
    def print_all(self):
        for d in sorted(self.all_available()):
            print(d)
        
        
    def database_directory(self, db_dir, species):
        proposed_db_dir = os.path.join(self.base_directory, species)
        if db_dir is not None:
            proposed_db_dir = os.path.join(db_dir, species)

        if isdir(proposed_db_dir):
            return proposed_db_dir
        else:
            print("Cannot access the database directory for the given species")
            return None
            

                
         