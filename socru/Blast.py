"""
BLAST wrapper for fragment database searching.

This module provides a wrapper around the NCBI BLAST+ blastn tool for comparing
inter-operon fragments against a database of known fragments. It constructs BLAST
commands, executes them with appropriate parameters, and returns results sorted
by bit score.

Classes:
    Blast: Wrapper for executing BLAST searches
"""

from tempfile import mkstemp
import subprocess
import os
import re

class Blast:
    """
    Wrapper for NCBI BLAST+ blastn command-line tool.
    
    This class manages BLAST database searching for fragment identification.
    It handles compressed input files, constructs optimized BLAST commands,
    and sorts results by bit score for downstream filtering.
    
    Attributes:
        blast_db (str): Path to BLAST database prefix
        evalue (float): E-value threshold for BLAST
        threads (int): Number of CPU threads to use
        word_size (int): BLAST word size parameter
        exec (str): BLAST executable name (default 'blastn')
        task (str): BLAST task type (default 'megablast')
        verbose (bool): Enable verbose output
        files_to_cleanup (list): Temporary files to delete on cleanup
    """
    def __init__(self, blast_db, threads, verbose, word_size = 28, exec = 'blastn', evalue = 0.000001, task = 'megablast'):
        """
        Initialize BLAST wrapper.
        
        Args:
            blast_db (str): Path to BLAST database
            threads (int): Number of CPU threads
            verbose (bool): Enable verbose output
            word_size (int): BLAST word size (default 28)
            exec (str): BLAST executable name (default 'blastn')
            evalue (float): E-value cutoff (default 0.000001)
            task (str): BLAST task algorithm (default 'megablast')
        """
        self.blast_db = blast_db
        self.evalue = evalue
        self.threads = threads
        self.word_size = word_size
        self.exec = exec
        self.task = task
        self.verbose = verbose
        self.files_to_cleanup = []

    def decompress_file_to_tmp(self, input_file):
        """
        Decompress gzipped query file if necessary.
        
        Args:
            input_file (str): Path to query file (possibly gzipped)
            
        Returns:
            str: Path to decompressed file or original if not compressed
        """
        m = re.search(r".gz$", input_file)
        if m:
            # Create temporary file and decompress
            fd, decompressed_input_file = mkstemp()
            self.files_to_cleanup.append(decompressed_input_file)
            cmd = " ".join(['gunzip', '-c', input_file, '>', decompressed_input_file])
            if self.verbose:
                print("Decompress file before blasting:\t" + cmd)
            subprocess.check_output( cmd,  shell=True)
            
            return decompressed_input_file
        else:
            return input_file

    def blast_command(self, query, blast_results):
        """
        Construct the BLAST command line.
        
        Builds a blastn command with tabular output (format 6), sorts results
        by bit score in descending order.
        
        Args:
            query (str): Path to query FASTA file
            blast_results (str): Path where results should be written
            
        Returns:
            str: Complete shell command to run BLAST
        """
        return " ".join([self.exec, '-outfmt', str(6), '-evalue', str(self.evalue), '-db', self.blast_db, '-word_size', str(self.word_size), '-num_threads', str(self.threads), '-task', self.task, '-query', self.decompress_file_to_tmp(query), '|', 'sort', '-k', str(12), '-g', '-r', '>', blast_results])
        
    def run_blast(self, query):
        """
        Execute BLAST search and return path to results.
        
        Args:
            query (str): Path to query FASTA file
            
        Returns:
            str: Path to BLAST results file (tabular format, sorted by bit score)
        """
        # Create temporary results file
        fd, blast_results = mkstemp()
        cmd = self.blast_command(query, blast_results)
        if self.verbose:
            print("Run blastn:\t" + cmd )
        subprocess.check_output( cmd,  shell=True)
        os.close(fd)
        return blast_results

    def __del__(self):
        """
        Clean up temporary decompressed files.
        """
        for f in self.files_to_cleanup:
            if os.path.exists(f):
                os.remove(f)    
                