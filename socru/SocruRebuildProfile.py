from socru.Profiles import Profiles
from socru.Results import Results
from socru.TypeGenerator import TypeGenerator
from socru.SocruUpdateProfile import SocruUpdateProfile
import shutil
import os
import re
from tempfile import mkstemp
from tempfile import mkdtemp

class SocruUpdateProfileOptions:
    def __init__(self, socru_output_filename, profile_filename, output_file, verbose):
        self.socru_output_filename = socru_output_filename
        self.profile_filename = profile_filename
        self.output_file = output_file
        self.verbose = verbose

class SocruRebuildProfile:
    def __init__(self,options):
        self.profile_filename = options.profile_filename
        self.output_file = options.output_file
        self.prefix = options.prefix
        self.verbose = options.verbose
        self.missing_character = '0'
        self.files_to_cleanup = []
        
    def run(self):
        # Get all the results from line 3 onwards
        # Save first 2 lines into output profile first
        # Set the GS number to be zero and save results to a file
        fd_ipf, intermediate_profile_file = mkstemp()
        fd_apf, additional_profiles_file = mkstemp()
        self.files_to_cleanup.append(intermediate_profile_file)
        self.files_to_cleanup.append(additional_profiles_file)
        
        regex_prefix_order = re.compile('^[\d]+\.')
        
        line_count = 1
        with open(self.profile_filename) as profile_fh:
            for line in profile_fh:
            
                if line_count <=2:
                    with open(intermediate_profile_file, '+a') as fh:
                        fh.write(line)
                else:
                    updated_profile = regex_prefix_order.sub( self.missing_character+'.', line)
                    fake_socru_result = "file\t" + self.prefix + updated_profile
                    with open(additional_profiles_file, '+a') as fh:
                        fh.write(fake_socru_result)
                line_count += 1
        
        os.close(fd_ipf)
        os.close(fd_apf)
        
        # incrementally update/add back to the profile
        fd_upf, unsorted_profile_file = mkstemp()
        self.files_to_cleanup.append(unsorted_profile_file)
        
        sup = SocruUpdateProfile(SocruUpdateProfileOptions(additional_profiles_file, intermediate_profile_file,  unsorted_profile_file, self.verbose))
        sup.run()
        
        # sort the profile file in numerical order
        with open(unsorted_profile_file) as unsorted_profile_fh:
            # print out the header
            line = unsorted_profile_fh.readline()
            with open(self.output_file, '+a') as fh:
                fh.write(line)
            
            # read in the rest of the lines (new line already on the end)
            lineList = unsorted_profile_fh.readlines()
            
            def atoi(text):
                return int(text) if text.isdigit() else text

            def sort_key(text):
                return [ atoi(c) for c in re.split(r'(\d+)', text) ]
            
            lineList.sort(key=sort_key)
            with open(self.output_file, '+a') as fh:
                fh.write(''.join(lineList)) 
        
        os.close(fd_upf)
        

        
    def __del__(self):
        for f in self.files_to_cleanup:
            if os.path.exists(f):
                os.remove(f)
        