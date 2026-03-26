"""
Profile database rebuild and renumbering workflow.

This module rebuilds a profile database by systematically renumbering all
GS types. It preserves fragment patterns but assigns sequential order numbers
and recalculates orientation numbers, resulting in a clean, standardized
profile database.

Classes:
    SocruUpdateProfileOptions: Options container for SocruUpdateProfile
    SocruRebuildProfile: Manages profile rebuild workflow
"""

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
    """
    Simple options container for SocruUpdateProfile.
    
    Attributes:
        socru_output_filename (str): Path to Socru results
        profile_filename (str): Path to profile database
        output_file (str): Path to output file
        verbose (bool): Enable verbose output
    """
    def __init__(self, socru_output_filename, profile_filename, output_file, verbose):
        """
        Initialize options.
        
        Args:
            socru_output_filename (str): Socru results path
            profile_filename (str): Profile database path
            output_file (str): Output path
            verbose (bool): Verbose flag
        """
        self.socru_output_filename = socru_output_filename
        self.profile_filename = profile_filename
        self.output_file = output_file
        self.verbose = verbose

class SocruRebuildProfile:
    """
    Rebuild profile database with systematic renumbering.
    
    This class takes an existing profile.txt and completely rebuilds it by:
    1. Preserving header rows
    2. Resetting all GS types to 0.X (marking as "new")
    3. Using SocruUpdateProfile to systematically reassign numbers
    4. Sorting results by GS type number
    
    This is useful for cleaning up databases after manual edits or merging.
    
    Attributes:
        profile_filename (str): Input profile.txt path
        output_file (str): Output profile.txt path
        prefix (str): GS type prefix (e.g., 'GS')
        verbose (bool): Enable verbose output
        missing_character (str): Character for "new" profiles (default '0')
        files_to_cleanup (list): Temporary files to delete
    """
    def __init__(self,options):
        """
        Initialize SocruRebuildProfile.
        
        Args:
            options: Parsed command-line arguments
        """
        self.profile_filename = options.profile_filename
        self.output_file = options.output_file
        self.prefix = options.prefix
        self.verbose = options.verbose
        self.missing_character = '0'
        self.files_to_cleanup = []
        
    def run(self):
        """
        Execute profile rebuild workflow.
        
        Steps:
        1. Split input into header (lines 1-2) and profiles (line 3+)
        2. Mark all profiles as new (0.X)
        3. Use update workflow to reassign systematic numbers
        4. Sort output by GS type number
        """
        # Create temporary files
        fd_ipf, intermediate_profile_file = mkstemp()
        fd_apf, additional_profiles_file = mkstemp()
        self.files_to_cleanup.append(intermediate_profile_file)
        self.files_to_cleanup.append(additional_profiles_file)
        
        # Regex to match order number (e.g., "1." in "1.3")
        regex_prefix_order = re.compile(r'^[\d]+\.')
        
        # Split input file: header vs profiles
        line_count = 1
        with open(self.profile_filename) as profile_fh:
            for line in profile_fh:
            
                if line_count <=2:
                    # Copy header rows to intermediate file
                    with open(intermediate_profile_file, '+a') as fh:
                        fh.write(line)
                else:
                    # Reset GS type to 0.X (mark as new)
                    updated_profile = regex_prefix_order.sub( self.missing_character+'.', line)
                    # Format as fake Socru result for update workflow
                    fake_socru_result = "file\t" + self.prefix + updated_profile
                    with open(additional_profiles_file, '+a') as fh:
                        fh.write(fake_socru_result)
                line_count += 1
        
        os.close(fd_ipf)
        os.close(fd_apf)
        
        # Use update workflow to reassign numbers systematically
        fd_upf, unsorted_profile_file = mkstemp()
        self.files_to_cleanup.append(unsorted_profile_file)
        
        sup = SocruUpdateProfile(SocruUpdateProfileOptions(additional_profiles_file, intermediate_profile_file,  unsorted_profile_file, self.verbose))
        sup.run()
        
        # Sort the profile file in numerical order by GS type
        with open(unsorted_profile_file) as unsorted_profile_fh:
            # Write header to output
            line = unsorted_profile_fh.readline()
            with open(self.output_file, '+a') as fh:
                fh.write(line)
            
            # Read remaining lines for sorting
            lineList = unsorted_profile_fh.readlines()
            
            # Define natural sort key function
            def atoi(text):
                """Convert text to int if numeric."""
                return int(text) if text.isdigit() else text

            def sort_key(text):
                """Natural sort key splitting on numbers."""
                return [ atoi(c) for c in re.split(r'(\d+)', text) ]
            
            # Sort lines naturally (1.1, 1.2, ..., 2.1, 2.2, ...)
            lineList.sort(key=sort_key)
            
            # Write sorted profiles to output
            with open(self.output_file, '+a') as fh:
                fh.write(''.join(lineList)) 
        
        os.close(fd_upf)
        

        
    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.cleanup()
        return False

    def cleanup(self):
        """Clean up temporary files and directories."""
        for f in self.files_to_cleanup:
            if os.path.exists(f):
                os.remove(f)
        self.files_to_cleanup = []

    def __del__(self):
        """Safety net cleanup -- prefer using as context manager."""
        self.cleanup()