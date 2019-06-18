from socru.Profiles import Profiles
from socru.Results import Results
from socru.TypeGenerator import TypeGenerator
import shutil

class SocruUpdateProfile:
    def __init__(self,options):
        self.socru_output_filename = options.socru_output_filename
        self.profile_filename = options.profile_filename
        self.output_file = options.output_file
        self.verbose = options.verbose
        
        shutil.copy(self.profile_filename, self.output_file)
        self.profiles = Profiles( self.profile_filename, self.verbose)
        self.results = Results( self.socru_output_filename, self.verbose)
            
    def run(self):
        valid_profiles_not_db_checked = self.results.filter(self.profiles.num_fragments)
        valid_profiles = []
        
        for input_profile in valid_profiles_not_db_checked:
            match_found = False
            for db_profile in self.profiles.gats:
                if db_profile.does_the_profile_match(input_profile):
                    match_found = True
                    continue
                    
            if not match_found:
                valid_profiles.append(input_profile)
        
        with open(self.output_file, "a+") as output_fh:
            for p in valid_profiles:
                tg = TypeGenerator(self.profiles, p, self.verbose, prefix = '')
                # get orientationless versions of p and db
                # check if there is an order number assigned, if not generate the next one.
                order_to_use = tg.calculate_orientationless_order()
                if order_to_use == 0:
                    order_to_use = self.profiles.next_order_number()
                    
                p.gat_number = str(order_to_use) + '.' + str(p.orientation_binary())
                self.profiles.gats.append(p)
                
                novel_type = p.gat_number + "\t" + str(p)
                output_fh.write(novel_type + "\n")
   
    