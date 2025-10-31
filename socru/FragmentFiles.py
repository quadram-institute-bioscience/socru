"""
Fragment file management for creating individual FASTA files.

This module manages the creation of individual FASTA files for each inter-operon
fragment. It reorders fragments starting with the largest, assigns numbers,
handles user-specified fragment orders, and writes out FASTA files for each
fragment to be used in BLAST searches.

Classes:
    FragmentFiles: Manages fragment ordering and FASTA file creation
"""

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

import os
import re

class FragmentFiles:
    """
    Create and manage individual FASTA files for genome fragments.
    
    This class takes a list of fragments, reorders them (starting with the
    largest), assigns fragment numbers, and creates individual FASTA files
    for each. Supports custom fragment ordering and reverse complementing.
    
    Attributes:
        fragments (list): List of Fragment objects
        output_directory (str): Directory where FASTA files will be written
        fragment_order (str): Optional custom fragment order (e.g., "1-2-3-4")
        verbose (bool): Enable verbose output
        ordered_fragments (list): Fragments reordered starting with largest
        output_filenames (list): Paths to created FASTA files
    """
    def __init__(self, fragments, output_directory, verbose, fragment_order = None):
        """
        Initialize FragmentFiles with fragments and output location.
        
        Args:
            fragments (list): List of Fragment objects to process
            output_directory (str): Directory for output FASTA files
            verbose (bool): Enable verbose output
            fragment_order (str): Optional custom order string (e.g., "1-2'-3")
        """
        self.fragments = fragments
        self.output_directory = output_directory
        self.fragment_order = fragment_order
        self.verbose = verbose
        # Reorder fragments with largest first and assign numbers
        self.ordered_fragments = self.fragments_with_largest_first()
        self.output_filenames = []
        
    def fragments_with_largest_first(self):
        """
        Reorder fragments starting with the largest and assign numbers.
        
        Rotates the fragment list so the largest fragment is first. This
        provides a consistent starting point for fragment numbering. If a
        custom fragment order is specified, applies it and handles reverse
        complementing (indicated by prime notation).
        
        Returns:
            list: Reordered Fragment objects with numbers assigned
        """
        # Find the largest fragment
        fragment_sizes = [f.num_bases() for f in self.fragments]
        largest_size = max(fragment_sizes)
        
        largest_index = 0
        for f in range(0,len(fragment_sizes)):
            if largest_size == fragment_sizes[f]:
                largest_index = f
                break

        # Rotate list to start with largest fragment
        reordered_fragments = self.fragments[largest_index:len(self.fragments)] + self.fragments[0:largest_index]
        
        if self.fragment_order is None:
            # Assign sequential numbers starting from 1
            for i in range(0,len(reordered_fragments)):
                reordered_fragments[i].number = i+1
        else:
            # Apply custom fragment order
            input_order = self.split_fragment_order()
            for i in range(0,len(reordered_fragments)):
                # Check if fragment should be reverse complemented (ends with ')
                m = re.match(r"([\d]+)'", input_order[i])
                if m:
                    # Reverse complement this fragment
                    reordered_fragments[i].number = m.group(1)
                    reordered_fragments[i].sequence = reordered_fragments[i].sequence.reverse_complement()
                    reordered_fragments[i].reversed = True
                else:
                    # Use fragment as-is
                    reordered_fragments[i].number = input_order[i]
                
        return reordered_fragments
      
    def create_fragment_fastas(self):
        """
        Write each fragment to its own FASTA file.
        
        Creates individual FASTA files in the output directory for each
        fragment. Files are named based on fragment number (e.g., "1.fa").
        Updates fragment objects with their output filenames.
        """
        for f in self.ordered_fragments:
            # Create SeqRecord for BioPython
            record = [SeqRecord(f.sequence, str(f) , '', '')]
            
            # Generate output filename and path
            outname = os.path.join(self.output_directory, f.output_filename())
            self.output_filenames.append(outname)
            f.output_filename = outname
            
            # Write FASTA file
            SeqIO.write(record, outname, "fasta")
    
    def split_fragment_order(self):
        """
        Parse custom fragment order string into list.
        
        Splits a dash-delimited string like "1-2-3'" into individual
        fragment identifiers.
        
        Returns:
            list: Individual fragment identifiers
        """
        return self.fragment_order.split('-')
        