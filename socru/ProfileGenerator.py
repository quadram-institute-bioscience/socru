"""
Profile database initialization for new species.

This module creates the initial profile database file when setting up a new
species. It generates a default profile (GS1.0 with fragments in order 1-2-3...),
identifies dnaA and dif positions by BLASTing against fragments, and creates
the metadata file with this information.

Classes:
    ProfileGenerator: Creates initial profile and metadata files
"""

import re
import csv
import os
import yaml
from socru.DnaA  import DnaA
from socru.Dif  import Dif

class ProfileGenerator:
    """
    Generate initial profile database and metadata for new species.
    
    This class creates the foundational files for a new species database:
    - profile.txt: Contains the default GS type (1.0 with forward fragments)
    - profile.txt.yml: Metadata including dnaA/dif positions
    
    It BLASTs dnaA and dif sequences against the extracted fragments to
    identify their locations.
    
    Attributes:
        output_directory (str): Database directory where files will be written
        output_filename (str): Full path to profile.txt
        input_file (str): Reference genome used for database creation
        prefix (str): Prefix for GS types (default 'GS')
        num_fragments (int): Number of fragments in genome
        verbose (bool): Enable verbose output
        dnaa_fasta (str): Path to dnaA query sequence
        dif_fasta (str): Path to dif query sequence
        threads (int): Number of CPU threads
        metadata_file_suffix (str): Suffix for metadata file (default '.yml')
    """
    def __init__(self, output_directory, num_fragments, dnaa_fasta, dif_fasta, threads, input_file,verbose, prefix = 'GS', output_filename = 'profile.txt', metadata_file_suffix = '.yml'):
        """
        Initialize ProfileGenerator with output location and parameters.
        
        Args:
            output_directory (str): Database directory
            num_fragments (int): Number of fragments
            dnaa_fasta (str): Path to dnaA sequence
            dif_fasta (str): Path to dif sequence
            threads (int): Number of threads
            input_file (str): Reference genome path
            verbose (bool): Enable verbose output
            prefix (str): GS type prefix (default 'GS')
            output_filename (str): Profile filename (default 'profile.txt')
            metadata_file_suffix (str): Metadata suffix (default '.yml')
        """
        self.output_directory = output_directory
        self.output_filename = os.path.join(self.output_directory, output_filename)
        self.input_file = input_file
        self.prefix = prefix
        self.num_fragments = num_fragments
        self.verbose = verbose
        self.dnaa_fasta = dnaa_fasta
        self.dif_fasta = dif_fasta
        self.threads = threads
        self.metadata_file_suffix = metadata_file_suffix
        
    def header(self):
        """
        Generate header row for profile file.
        
        Returns:
            list: Header with GS type column followed by fragment columns
        """
        header_row = [self.prefix]
        for i in range(self.num_fragments):
            header_row.append('Frag_'+str(i +1)) 
        return header_row
        
    def default_profile(self):
        """
        Generate default profile (GS1.0 with all fragments forward).
        
        The default profile represents the reference genome's structure
        with all fragments in order 1-2-3-... without reversals.
        
        Returns:
            list: Default profile row [1.0, 1, 2, 3, ...]
        """
        default_row = ['1.0']
        for i in range(self.num_fragments):
            default_row.append(str(i + 1)) 
        return default_row
    
    def write_output_file(self):
        """
        Write profile.txt file with header and default profile.
        
        Creates tab-delimited file with:
        - Header row
        - Default GS1.0 profile
        
        Also triggers metadata file generation.
        """
        content = [self.header(), self.default_profile()]
        with open(self.output_filename, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, delimiter='\t')
            writer.writerows(content)
            
        # Generate metadata file
        self.write_metadata_file()
    
    def write_metadata_file(self):
        """
        Create metadata YAML file with dnaA and dif positions.
        
        BLASTs dnaA and dif sequences against fragments to identify their
        locations and orientations. Writes metadata including:
        - dnaa_fragment: Fragment number containing origin
        - dnaa_forward_orientation: Orientation of dnaA
        - dif_fragment: Fragment number containing terminus
        - reference_genome: Source genome path
        """
        # Find dnaA position
        d = DnaA(self.dnaa_fasta, self.output_directory, self.threads, self.verbose)
        # Find dif position
        dif = Dif(self.dif_fasta, self.output_directory, self.threads, self.verbose)
        
        # Prepare metadata
        metadata_file =  self.output_filename + self.metadata_file_suffix
        metadata_content = { 
            'dnaa_fragment': int(d.fragment_with_dnaa), 
            'dif_fragment': int(dif.fragment_with_dif), 
            'dnaa_forward_orientation': d.forward_orientation, 
            'reference_genome':  self.input_file 
        }
        
        # Write YAML file
        with open(metadata_file, 'w') as yaml_file:
            yaml.dump(metadata_content, yaml_file, default_flow_style=False)