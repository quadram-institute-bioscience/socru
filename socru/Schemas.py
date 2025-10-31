"""
Species database schema management and lookup.

This module manages the bundled species databases that come with Socru.
It provides functions to list available species, display database metadata,
and locate database directories. Each species has a directory containing
fragment files and a profile database.

Classes:
    Schemas: Manages species database discovery and access
"""

import os
from os import listdir
from os.path import isdir
import pkg_resources
import yaml

class Schemas:
    """
    Manage species database schemas and metadata.
    
    This class provides access to the bundled species databases included
    with Socru. It can list available species, display extended metadata
    (dnaA/dif positions, reference genomes), and locate database directories.
    
    Attributes:
        verbose (bool): Enable verbose output
        base_directory (str): Path to bundled database directory
    """
    def __init__(self, verbose):
        """
        Initialize Schemas with database base directory.
        
        Args:
            verbose (bool): Enable verbose output
        """
        self.verbose = verbose
        # Get path to bundled data directory from package resources
        self.base_directory = str(pkg_resources.resource_filename( __name__, 'data/'))
    
    def all_available(self):
        """
        Get list of all available species databases.
        
        Returns:
            list: Names of available species (directory names)
        """
        return [ d for d in listdir(self.base_directory) if isdir(os.path.join(self.base_directory, d))]
        
    def print_all(self):
        """
        Print simple list of available species names.
        """
        for d in sorted(self.all_available()):
            print(d)
            
    def extended(self):
        """
        Get extended metadata for all species databases.
        
        Reads metadata YAML files to extract information about:
        - dnaA fragment number and orientation
        - dif fragment number
        - Reference genome used
        
        Returns:
            dict: Species name -> [species, dnaA_frag, dnaA_orient, dif_frag, reference]
        """
        db_info = {}
        for species in listdir(self.base_directory):
            
            db_dir = os.path.join(self.base_directory, species)
            db_file = os.path.join(self.base_directory, species,'profile.txt.yml')

            if not os.path.exists(db_file):
                continue
            
            # Read metadata YAML file    
            with open(db_file, 'r') as metadatafh:
                try:
                    metadata = yaml.load(metadatafh, yaml.SafeLoader)
                    dnaA_fragment_number = 0
                    dnaa_forward_orientation = False
                    dif_fragment_number = 0
                    reference_genome = ''
                    
                    # Extract dnaA information
                    if 'dnaa_fragment' in metadata:
                        dnaA_fragment_number = int(metadata['dnaa_fragment'])
                        dnaa_forward_orientation = metadata['dnaa_forward_orientation']
                    
                    # Extract dif information
                    if 'dif_fragment' in metadata:
                        dif_fragment_number = int(metadata['dif_fragment'])
                    
                    # Extract reference genome    
                    if 'reference_genome' in metadata:
                        reference_genome = str(metadata['reference_genome'])
                    
                    # Store database info    
                    db_info[species] = [species, dnaA_fragment_number, dnaa_forward_orientation, dif_fragment_number, reference_genome]
                    
					
                except yaml.YAMLError as exc:
                    print(exc)
        return db_info
        
    def print_extended(self):
        """
        Print detailed table of all species databases with metadata.
        
        Displays tab-delimited table with columns:
        Species, dnaA fragment No., dnaA forward orientation, dif fragment No., Reference
        """
        extended_details = self.extended()
        # Print header
        print("\t".join(['Species','dnaA fragment No.', 'dnaA forward orientation', 'dif fragment No.', 'Reference']))
        # Print each species
        for species in sorted(extended_details):
            print("\t".join( [ str(a) for a in extended_details[species] ]))
        
        
    def database_directory(self, db_dir, species):
        """
        Locate database directory for a species.
        
        First checks custom db_dir if provided, otherwise uses bundled
        databases. Returns path if directory exists.
        
        Args:
            db_dir (str): Custom database directory (optional)
            species (str): Species name
            
        Returns:
            str: Path to database directory, or None if not found
        """
        # Default to bundled database
        proposed_db_dir = os.path.join(self.base_directory, species)
        
        # Use custom directory if provided
        if db_dir is not None:
            proposed_db_dir = os.path.join(db_dir, species)

        # Check if directory exists
        if isdir(proposed_db_dir):
            return proposed_db_dir
        else:
            print("Cannot access the database directory for the given species")
            return None