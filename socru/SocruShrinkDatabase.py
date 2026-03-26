"""
Database shrinking workflow orchestration.

This module provides the workflow wrapper for optimizing fragment databases
based on coverage analysis. It handles directory validation, invokes the
shrinking process, and reports size reduction statistics.

Classes:
    SocruShrinkDatabase: Orchestrates database shrinking workflow
"""

import logging
import os
import sys

from socru.ShrinkDatabase import ShrinkDatabase

logger = logging.getLogger(__name__)

class SocruShrinkDatabase:
    """
    Orchestrate database shrinking workflow with validation and reporting.
    
    This class wraps the ShrinkDatabase functionality with:
    - Input/output directory validation
    - Directory creation
    - Size comparison and statistics reporting
    
    Attributes:
        blast_results (str): Path to accumulated BLAST results
        input_database (str): Source database directory
        output_database (str): Destination database directory
        min_fragment_size (int): Minimum bases to retain per fragment
        verbose (bool): Enable verbose output
    """
    def __init__(self,options):
        """
        Initialize SocruShrinkDatabase and validate directories.
        
        Args:
            options: Parsed command-line arguments
        """
        self.blast_results = options.blast_results
        self.input_database = options.input_database
        self.output_database = options.output_database
        self.min_fragment_size = options.min_fragment_size
        self.verbose = options.verbose
        
        # Validate input database exists
        if not os.path.exists(self.input_database):
             logger.error("Cannot access the input database you specified, please check again")
             sys.exit(1)

        # Validate output database doesn't exist
        if os.path.exists(self.output_database):
             logger.error("Output directory already exists, please choose a new name")
             sys.exit(1)
        else:
            # Create output directory
            os.makedirs(self.output_database)
        
        # Get absolute path for output
        self.output_database = os.path.abspath(self.output_database)

    def run(self):
        """
        Execute database shrinking and report statistics.
        
        Measures input and output database sizes, performs shrinking,
        and returns a tab-delimited statistics string.
        
        Returns:
            str: Tab-delimited statistics: "species\tMB_reduction\tpercent_reduction"
        """
        # Measure input size
        input_dir_size =self.directory_size(self.input_database) 
        
        # Perform shrinking
        s = ShrinkDatabase(self.input_database, self.output_database, self.blast_results, self.min_fragment_size, self.verbose)
        s.shrink_files()
        
        # Measure output size    
        output_dir_size = self.directory_size(self.output_database)
        
        # Calculate reduction statistics
        percentage_reduction = round(100 - (output_dir_size*100/input_dir_size), 2)
        mb_reduction = round((input_dir_size - output_dir_size)/1000000, 3)
        
        # Format statistics string
        stats = ( os.path.basename(self.input_database) + "\t" + str(mb_reduction) + "\t" + str(percentage_reduction))
        
        return stats
        
    def directory_size(self, directory):
        """
        Calculate total size of all files in directory.
        
        Args:
            directory (str): Path to directory
            
        Returns:
            int: Total size in bytes
        """
        return sum(os.path.getsize(os.path.join(directory ,f)) for f in os.listdir(directory ) if os.path.isfile(os.path.join(directory ,f)))