"""
BLAST wrapper for fragment database searching.

This module provides a wrapper around the NCBI BLAST+ blastn tool for comparing
inter-operon fragments against a database of known fragments. It constructs BLAST
commands, executes them with appropriate parameters, and returns results sorted
by bit score.

Classes:
    Blast: Wrapper for executing BLAST searches
"""

import gzip
import logging
import os
import re
import shutil
import subprocess
from tempfile import mkstemp

from socru.ToolCheck import MissingToolError

logger = logging.getLogger(__name__)

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
            os.close(fd)
            self.files_to_cleanup.append(decompressed_input_file)

            logger.info("Decompress file before blasting:\tgunzip -c %s > %s", input_file, decompressed_input_file)

            with gzip.open(input_file, 'rb') as gz_in:
                with open(decompressed_input_file, 'wb') as f_out:
                    shutil.copyfileobj(gz_in, f_out)

            return decompressed_input_file
        else:
            return input_file

    def run_blast(self, query):
        """
        Execute BLAST search and return path to results.

        Runs blastn as a subprocess with safe argument list, captures output,
        sorts by bit score (column 12) in descending order, and writes to
        a temporary results file.

        Args:
            query (str): Path to query FASTA file

        Returns:
            str: Path to BLAST results file (tabular format, sorted by bit score)
        """
        # Create temporary results file
        fd, blast_results = mkstemp()
        os.close(fd)

        decompressed_query = self.decompress_file_to_tmp(query)

        cmd = [
            self.exec, '-outfmt', str(6), '-evalue', str(self.evalue),
            '-db', self.blast_db, '-word_size', str(self.word_size),
            '-num_threads', str(self.threads), '-task', self.task,
            '-query', decompressed_query,
        ]
        logger.info("Run blastn:\t%s", ' '.join(cmd))

        try:
            result = subprocess.run(cmd, capture_output=True, text=True)
            result.check_returncode()
        except subprocess.CalledProcessError as e:
            logger.error("blastn failed with exit code %d on %s", e.returncode, query)
            raise
        except FileNotFoundError:
            raise MissingToolError(f"Tool not found: {cmd[0]}")

        # Sort output by bit score (column 12, 0-indexed 11) descending
        lines = result.stdout.strip().split('\n') if result.stdout.strip() else []
        if not lines or (len(lines) == 1 and not lines[0].strip()):
            logger.warning("blastn produced no output for query %s", query)
            lines = []
        sorted_lines = sorted(
            lines,
            key=lambda line: float(line.split('\t')[11]) if len(line.split('\t')) > 11 else 0.0,
            reverse=True,
        )

        with open(blast_results, 'w') as f:
            for line in sorted_lines:
                f.write(line + '\n')

        return blast_results

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
