"""
Profile database update workflow with novel patterns.

This module updates profile databases by adding validated novel genome structure
patterns from Socru analysis results. It filters novel patterns, assigns them
proper GS type numbers, and appends them to the profile database.

Classes:
    SocruUpdateProfile: Updates profile databases with novel types
"""

import shutil

from socru.Profiles import Profiles
from socru.Results import Results
from socru.TypeGenerator import TypeGenerator


class SocruUpdateProfile:
    """
    Update profile database with validated novel genome structures.

    This class processes Socru output containing novel profiles (GS0.X types),
    validates them, assigns proper GS type numbers, and adds them to the
    profile database. It handles:
    - Filtering to remove duplicates
    - Assigning order numbers (new or matching existing orientationless patterns)
    - Assigning orientation numbers
    - Appending to profile.txt

    Attributes:
        socru_output_filename (str): Path to Socru results with novel profiles
        profile_filename (str): Path to existing profile.txt
        output_file (str): Path to updated profile.txt
        verbose (bool): Enable verbose output
        profiles (Profiles): Loaded profile database
        results (Results): Parsed Socru results
    """
    def __init__(self,options):
        """
        Initialize SocruUpdateProfile and load databases.

        Args:
            options: Parsed command-line arguments
        """
        self.socru_output_filename = options.socru_output_filename
        self.profile_filename = options.profile_filename
        self.output_file = options.output_file
        self.verbose = options.verbose

        # Copy existing profile to output (will append to it)
        shutil.copy(self.profile_filename, self.output_file)

        # Load existing profiles
        self.profiles = Profiles( self.profile_filename, self.verbose)

        # Parse results with novel profiles
        self.results = Results( self.socru_output_filename, self.verbose)

    def run(self):
        """
        Process novel profiles and add them to database.

        Workflow:
        1. Filter results to get valid novel profiles
        2. Remove profiles already in database
        3. For each truly novel profile:
           - Determine order number (orientationless matching)
           - Assign full GS type (order.orientation)
           - Append to output file
        """
        # Get validated novel profiles (complete patterns, correct fragment count)
        valid_profiles_not_db_checked = self.results.filter(self.profiles.num_fragments)
        valid_profiles = []

        # Check against existing database to avoid duplicates
        for input_profile in valid_profiles_not_db_checked:
            match_found = False
            for db_profile in self.profiles.gats:
                if db_profile.does_the_profile_match(input_profile):
                    # Already in database
                    match_found = True
                    continue

            if not match_found:
                # Truly novel profile
                valid_profiles.append(input_profile)

        # Add each novel profile to output file
        with open(self.output_file, "a+") as output_fh:
            for p in valid_profiles:
                # Create type generator to calculate order number
                tg = TypeGenerator(self.profiles, p, self.verbose, True, prefix = '')

                # Try to match orientationless pattern to existing order
                order_to_use = tg.calculate_orientationless_order()
                if order_to_use == 0:
                    # Completely new pattern - assign next order number
                    order_to_use = self.profiles.next_order_number()

                # Assign full GS type: order.orientation
                p.gat_number = str(order_to_use) + '.' + str(p.orientation_binary())

                # Add to in-memory database (for subsequent matches)
                self.profiles.gats.append(p)

                # Write to file
                novel_type = p.gat_number + "\t" + str(p)
                output_fh.write(novel_type + "\n")
