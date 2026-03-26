"""
Novel profile results parsing and validation.

This module parses Socru output files containing novel genome structure profiles
and validates them for inclusion in databases. It filters profiles to ensure
all fragments are present exactly once and that the patterns are truly novel.

Classes:
    Results: Parses and validates novel profile results
"""

from __future__ import annotations

import csv
import re

from socru.GATProfile import GATProfile


class Results:
    """
    Parse and validate novel genome structure profiles from results.

    This class reads Socru output files (especially those containing novel
    profiles), extracts the genome structure patterns, and validates them
    for potential inclusion in a database. It filters out incomplete patterns,
    duplicates, and profiles that don't meet validation criteria.

    Attributes:
        results_file (str): Path to Socru results file
        verbose (bool): Enable verbose output
        profiles (list): List of GATProfile objects extracted from results
    """
    def __init__(self, results_file: str, verbose: bool) -> None:
        """
        Initialize Results parser and extract profiles.

        Args:
            results_file (str): Path to Socru output file
            verbose (bool): Enable verbose output
        """
        self.results_file = results_file
        self.verbose = verbose
        # Parse profiles from results file
        self.profiles = self.create_profiles()

    def create_profiles(self) -> list[GATProfile]:
        """
        Parse Socru results file to extract genome structure profiles.

        Reads tab-delimited results and extracts profile information:
        - Column 1: Database directory
        - Column 2: GS type (e.g., "GS0.1")
        - Columns 3+: Fragment pattern

        Skips profiles with unknown fragments ('?') and removes duplicates.

        Returns:
            list: Unique GATProfile objects from results
        """
        seen_profiles = []
        profiles = []
        with open(self.results_file, newline='') as csvfile:
            profile_reader = csv.reader(csvfile, delimiter='\t')
            for row in profile_reader:
                if len(row) > 3:
                    # Parse results format:
                    # 1: directory of schema
                    # 2: GS number with GS at the start
                    # 3..N: fragment pattern
                    m = re.match("GS(.+)", row[1])
                    if m:
                        gat_number = m.group(1)
                        fragments = [row[f] for f in range(2, len(row))]

                        # Skip profiles with unknown fragments
                        unknown = [f for f in fragments if f == '?' or f == "?'"]
                        if len(unknown) > 0:
                            continue

                        # Create profile and check for duplicates
                        g = GATProfile(self.verbose, gat_number = gat_number, fragments = fragments)
                        if str(g) not in seen_profiles:
                            seen_profiles.append(str(g))
                            profiles.append(g)
        return profiles

    def filter(self, num_fragments: int) -> list[GATProfile]:
        """
        Filter profiles to keep only valid, novel patterns.

        Applies two filters:
        1. Removes profiles seen before (non-zero major order number)
        2. Ensures all fragments 1..N are present exactly once

        Args:
            num_fragments (int): Expected number of fragments in valid profiles

        Returns:
            list: Filtered list of valid novel profiles
        """
        # Filter out previously seen profiles
        profiles_novel = self.filter_previously_seen_profiles(self.profiles)
        # Ensure all fragments present once
        valid_profiles = self.all_fragments_present_once(num_fragments, profiles_novel)
        return valid_profiles

    def all_fragments_present_once(self, num_fragments: int, profiles: list[GATProfile]) -> list[GATProfile]:
        """
        Validate that each profile has all fragments 1..N exactly once.

        Checks:
        1. Profile has correct number of fragments
        2. When orientations removed and sorted, fragments are [1,2,3,...,N]

        Args:
            num_fragments (int): Expected number of fragments
            profiles (list): Profiles to validate

        Returns:
            list: Profiles where all fragments present exactly once
        """
        filtered_profiles = []
        for p in profiles:
            # Check fragment count
            if len(p.fragments) != num_fragments:
                continue

            # Check that fragments 1..N are all present once
            sorted_orientationless = sorted(p.orientationless_fragments())
            out_of_order = [ i+1 for i in range(0,len(sorted_orientationless)) if int(sorted_orientationless[i]) != i+1]
            if len(out_of_order) > 0:
                # Missing or duplicate fragments
                continue
            else:
                filtered_profiles.append(p)

        return filtered_profiles


    def filter_previously_seen_profiles(self, profiles: list[GATProfile]) -> list[GATProfile]:
        """
        Keep only truly novel profiles (major order number is 0).

        GS types like "0.X" indicate novel patterns not in the database.
        GS types like "1.X" indicate known order with different orientation.
        This method keeps only "0.X" patterns.

        Args:
            profiles (list): Profiles to filter

        Returns:
            list: Profiles with order number 0 (novel)
        """
        filtered_profiles = []
        for p in profiles:
            m = re.match(r'([\d]+)\.([\d]+)', p.gat_number)
            if m:
                if int(m.group(1)) == 0:
                     # Order 0 means not seen before - keep it
                     filtered_profiles.append(p)
        return filtered_profiles
