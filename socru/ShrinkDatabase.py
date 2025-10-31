"""
Database optimization by extracting high-coverage regions.

This module shrinks fragment databases by analyzing BLAST coverage from multiple
genome analyses and extracting only the most-covered regions. This reduces
database size while maintaining typing accuracy by keeping the conserved
portions of fragments that are consistently found across many genomes.

Classes:
    ShrinkDatabase: Optimizes databases based on coverage analysis
"""

import os
from os import listdir
from os.path import isfile, join
import re
import shutil
import subprocess
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from socru.FilterBlast import FilterBlast
from socru.Fasta import Fasta

class ShrinkDatabase:
    """
    Optimize fragment database by extracting high-coverage regions.
    
    This class analyzes BLAST results from analyzing many genomes to identify
    which parts of each fragment are consistently found. It then creates a
    new database containing only these high-coverage regions, reducing storage
    and search time while maintaining typing accuracy.
    
    Small fragments are kept intact; large fragments are trimmed to only the
    most-covered regions.
    
    Attributes:
        input_database (str): Source database directory
        output_database (str): Destination database directory
        blast_results (str): BLAST results file from multiple genome analyses
        target_bases (int): Minimum bases to retain per fragment
        verbose (bool): Enable verbose output
    """
    def __init__(self,input_database, output_database, blast_results, target_bases, verbose):
        """
        Initialize ShrinkDatabase with paths and parameters.
        
        Args:
            input_database (str): Source database directory
            output_database (str): Destination directory (must not exist)
            blast_results (str): Path to accumulated BLAST results
            target_bases (int): Minimum bases to keep per fragment
            verbose (bool): Enable verbose output
        """
        self.input_database = input_database
        self.output_database = output_database
        self.blast_results = blast_results
        self.target_bases = target_bases
        self.verbose = verbose
    
    def get_database_files(self):
        """
        Get list of uncompressed fragment FASTA files.
        
        Returns:
            list: Paths to plain .fa files
        """
        return [os.path.join(self.input_database,f) for f in listdir(self.input_database) if isfile(join(self.input_database, f)) and re.match(r'[\d]+\.fa$', f)]
        
    def get_database_files_compressed(self):
        """
        Get list of compressed fragment FASTA files.
        
        Returns:
            list: Paths to .fa.gz files
        """
        return [os.path.join(self.input_database,f) for f in listdir(self.input_database) if isfile(join(self.input_database, f)) and re.match(r'[\d]+\.fa.gz$', f)]

    def shrink_files(self):
        """
        Process all fragments: copy small ones, extract regions from large ones.
        
        For each fragment:
        - If smaller than target_bases: copy entire fragment
        - If larger than target_bases: analyze coverage and extract high-coverage regions
        
        Also copies profile.txt and metadata files to output database.
        
        Returns:
            list: Paths to compressed output files
        """
        output_filenames = []
        
        # Copy profile files to destination unchanged
        shutil.copy(os.path.join(self.input_database, 'profile.txt'), self.output_database)
        shutil.copy(os.path.join(self.input_database, 'profile.txt.yml'), self.output_database)
        
        # Get all fragment files
        fasta_file_names_compressed = self.get_database_files_compressed()
        fasta_file_names_uncompressed = self.get_database_files()
        fasta_file_names = fasta_file_names_uncompressed + fasta_file_names_compressed
        
        # Process each fragment
        fasta_obj = [ Fasta(f, self.verbose) for f in fasta_file_names]
        for f in fasta_obj:
            # Analyze coverage for this fragment
            fb = FilterBlast(self.blast_results, 1, 1, self.verbose)
            destination_filename = os.path.join(self.output_database, str(f.fragment_number()) + '.fa')
            
            if len(f.chromosome.seq) < self.target_bases:
                # Fragment is small - copy entire file
                if f.input_file in fasta_file_names_compressed:
                    shutil.copy(f.input_file, destination_filename + '.gz')
                    output_filenames.append(destination_filename + '.gz')
                else:
                    shutil.copy(f.input_file, destination_filename)
                    output_filenames.append(destination_filename)
            else:
                # Fragment is large - extract high-coverage regions
                blocks = fb.identify_regions(f.fragment_number(), self.target_bases)
                sequence = ""
                for b in blocks:
                    # Concatenate high-coverage blocks
                    sequence += f.chromosome.seq[(b[0]):(b[1])]
                
                # Write extracted sequence
                record = [SeqRecord(sequence, str(f.fragment_number()) , '', '')]
                SeqIO.write(record, destination_filename, "fasta")
                output_filenames.append(destination_filename)
                
        # Compress all output files
        return self.compress_files(output_filenames)
            
    def compress_files(self, filenames):
        """
        Compress all uncompressed FASTA files with gzip.
        
        Args:
            filenames (list): Paths to files to compress
            
        Returns:
            list: Paths to compressed files (*.fa.gz)
        """
        compressed_filenames = []
        for filename in filenames:
            # Check if file needs compression
            m = re.search(r"[\d]+\.fa$", filename)
            if m:
                # Compress with gzip
                subprocess.check_output( 'gzip '+ filename,  shell=True)
                compressed_filenames.append(filename + '.gz')
            else:
                # Already compressed
                compressed_filenames.append(filename)
        return compressed_filenames