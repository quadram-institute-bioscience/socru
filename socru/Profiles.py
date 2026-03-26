"""
Profile database management and loading.

This module handles loading and managing the profile database that contains
known genome arrangement types (GS types). It reads profile files listing
fragment patterns, loads metadata about special fragments (dnaA, dif), and
provides access to profile information for genome typing.

Classes:
    Profiles: Manages profile database loading and access
"""

import csv
import logging
import yaml
import os
from socru.GATProfile import GATProfile

logger = logging.getLogger(__name__)


class Profiles:
    """
    Manage profile database of known genome arrangement types.
    
    This class loads a profile database from a tab-delimited file, where each
    row represents a known GS type and its fragment pattern. It also reads
    metadata about special fragments like dnaA and dif positions. The profiles
    are used to match query genomes to known types.
    
    Attributes:
        input_file (str): Path to profile database file (tab-delimited)
        metadata_file (str): Path to metadata YAML file
        dnaA_fragment_number (int): Fragment containing dnaA origin
        dif_fragment_number (int): Fragment containing dif terminus  
        dnaa_forward_orientation (bool): Expected orientation of dnaA
        verbose (bool): Enable verbose output
        gats (list): List of GATProfile objects from database
        num_fragments (int): Expected number of fragments in profiles
    """
    def __init__(self, input_file, verbose, metadata_file_suffix = '.yml', default_dnaA_fragment_number = 3, default_dif_fragment_number = 1):
        """
        Initialize Profiles by loading database and metadata.
        
        Args:
            input_file (str): Path to profile database file
            verbose (bool): Enable verbose output
            metadata_file_suffix (str): Suffix for metadata file (default '.yml')
            default_dnaA_fragment_number (int): Default dnaA fragment (default 3)
            default_dif_fragment_number (int): Default dif fragment (default 1)
        """
        self.input_file = input_file
        self.metadata_file = self.input_file + metadata_file_suffix
        self.dnaA_fragment_number = default_dnaA_fragment_number
        self.dif_fragment_number = default_dif_fragment_number
        self.dnaa_forward_orientation = False
        self.verbose = verbose
        
        # Load profiles from database file
        self.gats = self.read_profiles()
        # Load metadata (dnaA, dif positions)
        self.read_metadata()
        # Determine expected number of fragments
        self.num_fragments = self.expected_num_fragments()
        
    def read_profiles(self):
        """
        Read profile database from tab-delimited file.
        
        File format: Each row has GS type in first column, followed by
        fragment pattern (e.g., "1  2'  3  4")
        
        Returns:
            list: List of GATProfile objects from database
        """
        with open(self.input_file, newline='') as csvfile:
            profile_reader = csv.reader(csvfile, delimiter='\t')
            # skip the header line
            next(profile_reader) 
            profiles = []
            for row in profile_reader:
                if len(row) > 2:
                    # Extract fragment pattern from columns after GS type
                    fragments = [row[f] for f in range(1, len(row)) if row[f] != '']
                    # Create GATProfile with GS type number and fragment pattern
                    g = GATProfile(self.verbose, gat_number = row[0], fragments = fragments)
                    # Ensure profile is oriented correctly relative to dnaA
                    g.orientate_for_dnaA()
                    profiles.append(g)
            return profiles
            
    def read_metadata(self):
        """
        Read metadata YAML file for database configuration.
        
        Metadata includes positions of special fragments:
        - dnaa_fragment: Fragment containing origin
        - dnaa_forward_orientation: Expected orientation
        - dif_fragment: Fragment containing terminus
        
        Returns:
            None
        """
        if not os.path.exists(self.metadata_file):
            return
        with open(self.metadata_file, 'r') as metadatafh:
            try:
                metadata = yaml.load(metadatafh, yaml.SafeLoader)
                # Load dnaA fragment information
                if 'dnaa_fragment' in metadata:
                    self.dnaA_fragment_number = int(metadata['dnaa_fragment'])
                    self.dnaa_forward_orientation = metadata['dnaa_forward_orientation']
                
                # Load dif fragment information	
                if 'dif_fragment' in metadata:
                    self.dif_fragment_number = int(metadata['dif_fragment'])
					
            except yaml.YAMLError as exc:
                logger.warning("Failed to parse YAML metadata: %s", exc)
        return
        
    def expected_num_fragments(self):
        """
        Determine expected number of fragments from first profile.
        
        All profiles should have the same number of fragments. This gets
        the count from the first profile to validate query genomes.
        
        Returns:
            int: Expected number of fragments, or -1 if no profiles loaded
        """
        if len(self.gats)> 0:
            return len(self.gats[0].fragments)
        else:
            return -1
            
    def next_order_number(self):
        """
        Get next available GS type order number.
        
        Used when creating new GS types for novel patterns. Returns one
        more than the highest existing order number.
        
        Returns:
            int: Next available order number
        """
        return max([p.order() for p in self.gats]) + 1
            