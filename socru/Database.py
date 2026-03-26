"""
Fragment database management and BLAST database creation.

This module handles loading fragment databases and creating BLAST databases
from fragment FASTA files. It concatenates individual fragment files (both
plain and gzipped), creates a temporary BLAST database, and manages cleanup.

Classes:
    Database: Manages fragment database loading and BLAST database creation
"""

import gzip
import os
from os import listdir
from os.path import isfile, join
import re
from tempfile import mkstemp
from tempfile import mkdtemp
import subprocess
import shutil

class Database:
    """
    Manage fragment databases and BLAST database creation.

    This class locates fragment files in a database directory, concatenates
    them (handling both plain and gzipped formats), and creates a BLAST
    database for searching. The BLAST database is temporary and cleaned
    up automatically.

    Attributes:
        directory (str): Path to database directory containing fragment files
        verbose (bool): Enable verbose output
        concat_fasta (str): Path to concatenated temporary FASTA file
        db_prefix (str): Path prefix for BLAST database files
    """
    def __init__(self,directory, verbose):
        """
        Initialize Database by loading and indexing fragments.

        Args:
            directory (str): Path to database directory
            verbose (bool): Enable verbose output
        """
        self.directory = directory
        self.verbose = verbose
        self._tmpdir = None
        # Concatenate all fragment files into one FASTA
        self.concat_fasta = self.concat_db_files()
        # Create BLAST database from concatenated FASTA
        self.db_prefix = self.make_blastdb(self.concat_fasta)

    def get_database_files(self):
        """
        Get list of plain (uncompressed) fragment FASTA files.

        Finds files matching pattern: [digit]+.fa

        Returns:
            list: Paths to plain FASTA fragment files
        """
        return [os.path.join(self.directory,f) for f in listdir(self.directory) if isfile(join(self.directory, f)) and re.match(r'[\d]+\.fa$', f)]

    def get_database_files_compressed(self):
        """
        Get list of gzipped fragment FASTA files.

        Finds files matching pattern: [digit]+.fa.gz

        Returns:
            list: Paths to gzipped FASTA fragment files
        """
        return [os.path.join(self.directory,f) for f in listdir(self.directory) if isfile(join(self.directory, f)) and re.match(r'[\d]+\.fa.gz$', f)]

    def concat_db_files(self):
        """
        Concatenate all fragment files into a single temporary FASTA.

        Handles a mixture of gzipped and plain FASTA files by decompressing
        gzipped files on the fly and appending all to one file.

        Returns:
            str: Path to concatenated temporary FASTA file
        """
        # Create temporary file for concatenated sequences
        fd, concat_db_fasta = mkstemp()
        os.close(fd)

        with open(concat_db_fasta, 'wb') as out_fh:
            # Decompress and concatenate gzipped files
            for gz_file in self.get_database_files_compressed():
                with gzip.open(gz_file, 'rb') as gz_in:
                    shutil.copyfileobj(gz_in, out_fh)

            # Append plain files
            for plain_file in self.get_database_files():
                with open(plain_file, 'rb') as f_in:
                    shutil.copyfileobj(f_in, out_fh)

        return concat_db_fasta

    def make_blastdb(self, concat_fasta):
        """
        Create BLAST nucleotide database from concatenated FASTA.

        Runs makeblastdb to create a searchable BLAST database in a
        temporary directory.

        Args:
            concat_fasta (str): Path to concatenated FASTA file

        Returns:
            str: Path prefix for BLAST database files
        """
        # Create temporary directory for BLAST database
        tmpdir = mkdtemp()
        self._tmpdir = tmpdir
        output_prefix = os.path.join(tmpdir, 'all')

        # Run makeblastdb command
        cmd = ['makeblastdb', '-in', concat_fasta, '-dbtype', 'nucl', '-out', output_prefix]
        if self.verbose:
            print("Creating blast database:\t" + ' '.join(cmd))
        subprocess.run(cmd, check=True, capture_output=True)
        return output_prefix

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.cleanup()
        return False

    def cleanup(self):
        """Clean up temporary files and directories."""
        # Remove BLAST database directory
        if self._tmpdir is not None and os.path.exists(self._tmpdir):
            shutil.rmtree(self._tmpdir, ignore_errors=True)
            self._tmpdir = None
        # Remove concatenated FASTA file
        if hasattr(self, 'concat_fasta') and self.concat_fasta and os.path.exists(self.concat_fasta):
            os.remove(self.concat_fasta)
            self.concat_fasta = None

    def __del__(self):
        """Safety net cleanup -- prefer using as context manager."""
        self.cleanup()
