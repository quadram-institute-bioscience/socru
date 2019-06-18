import os
import sys

from socru.ShrinkDatabase import ShrinkDatabase

class SocruShrinkDatabase:
    def __init__(self,options):
        self.blast_results = options.blast_results
        self.input_database = options.input_database
        self.output_database = options.output_database
        self.min_fragment_size = options.min_fragment_size
        self.verbose = options.verbose
        
        if not os.path.exists(self.input_database):
             print(
             "Cannot access the input database you specified, please check again")
             sys.exit(1)
             
        if os.path.exists(self.output_database):
             print(
             "Output directory already exists, please choose a new name")
             sys.exit(1)
        else:
            os.makedirs(self.output_database)
        self.output_database = os.path.abspath(self.output_database)

    def run(self):
        input_dir_size =self.directory_size(self.input_database) 
        s = ShrinkDatabase(self.input_database, self.output_database, self.blast_results, self.min_fragment_size, self.verbose)
        s.shrink_files()
            
        output_dir_size = self.directory_size(self.output_database)
        
        percentage_reduction = round(100 - (output_dir_size*100/input_dir_size), 2)
        mb_reduction = round((input_dir_size - output_dir_size)/1000000, 3)
        stats = ( os.path.basename(self.input_database) + "\t" + str(mb_reduction) + "\t" + str(percentage_reduction))
        
        return stats
        
    def directory_size(self, directory):
        return sum(os.path.getsize(os.path.join(directory ,f)) for f in os.listdir(directory ) if os.path.isfile(os.path.join(directory ,f)))
        
      
        