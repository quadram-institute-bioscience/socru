"""
FASTA file handling and fragment coordinate calculation for Socru.

This module handles reading genome FASTA files (including gzipped versions),
identifying the largest contig (chromosome), calculating inter-operon fragment
coordinates, and extracting fragment sequences. It supports both circular and
linear chromosomes.

Classes:
    Fasta: Main class for FASTA parsing and fragment extraction
"""

from Bio import SeqIO
from socru.Fragment import Fragment
import logging
import os
import re
import sys
import gzip

logger = logging.getLogger(__name__)

class Fasta:
    """
    Handle FASTA file reading and fragment sequence extraction.
    
    This class reads genome assemblies in FASTA format, identifies the largest
    contig (assumed to be the main chromosome), and extracts inter-operon
    fragment sequences based on rRNA boundary positions.
    
    Attributes:
        input_file (str): Path to input FASTA file (may be gzipped)
        is_circular (bool): Whether to treat chromosome as circular
        verbose (bool): Enable verbose output
        chromosome (SeqRecord): BioPython sequence record for the chromosome
    """
    def __init__(self, input_file,verbose, is_circular=True):
        """
        Initialize Fasta parser.
        
        Args:
            input_file (str): Path to input FASTA file
            verbose (bool): Enable verbose output
            is_circular (bool): Assume circular chromosome (default True)
        """
        self.input_file = input_file
        self.is_circular = is_circular
        self.verbose = verbose

        # Validate input file exists
        if not os.path.exists(input_file):
            raise FileNotFoundError(f"Input file not found: {input_file}")

        # Validate input file is not empty
        if os.path.getsize(input_file) == 0:
            raise ValueError(f"Input file is empty: {input_file}")

        self.chromosome = self.get_chromosome_from_fasta()
    
    def get_chromosome_from_fasta(self):
        """
        Extract the largest contig from FASTA file as the chromosome.
        
        Handles both plain and gzipped FASTA files. Selects the largest
        sequence by length, assuming this is the main chromosome.
        
        Returns:
            SeqRecord: BioPython sequence record of the largest contig
        """
        # find largest contig and ignore the rest
        largest_contig = None
		
        # Check if file is gzipped based on extension
        m = re.search(".gz$", self.input_file)
        if m:
            # Handle gzipped FASTA
            with gzip.open(self.input_file, "rt") as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    largest_contig = self.largest_contig_check( largest_contig, record)
        else:
            # Handle plain FASTA
            for record in SeqIO.parse(self.input_file, "fasta"):
                largest_contig = self.largest_contig_check( largest_contig, record)

        if largest_contig is None:
            raise ValueError(f"No valid sequences found in {self.input_file}")

        if len(largest_contig.seq) < 100000:
            logger.warning(
                "Largest contig in %s is only %d bp, which is unlikely to be a "
                "complete chromosome",
                self.input_file, len(largest_contig.seq),
            )

        return largest_contig
        
    def fragment_number(self):
        """
        Get fragment number from chromosome ID.
        
        Returns:
            int: Fragment number parsed from chromosome identifier
        """
        return int(self.chromosome.id)
        
    def largest_contig_check(self, largest_contig, record):
        """
        Compare current record to largest seen so far, update if bigger.
        
        Args:
            largest_contig (SeqRecord): Current largest contig (or None)
            record (SeqRecord): New record to compare
            
        Returns:
            SeqRecord: The larger of the two contigs
        """
        sequence_length  = len(record.seq)
        if largest_contig is None \
           or sequence_length > len(largest_contig.seq):
           largest_contig = record
        return largest_contig
        
    def calc_fragment_coords(self, boundries):
        """
        Calculate coordinates of inter-operon fragments from rRNA boundaries.
        
        This method determines the genomic regions between rRNA operons and
        creates Fragment objects with appropriate coordinates. Handles both
        circular and linear chromosomes, including edge cases where fragments
        span the origin.
        
        Args:
            boundries (list): List of operon boundary objects with start/end positions
            
        Returns:
            list: List of Fragment objects with coordinate ranges
        """
        genome_length = len(self.chromosome.seq)
        fragments = []
        if not boundries:
            return []
        
        start_coord = 0
        end_coord = genome_length
        
        # Check if any fragment wraps around chromosome end (start > end)
        end_frag = [ b for b in boundries if b.start > b.end ]
        if len(end_frag) > 0:
            # Fragment wrapping indicates linear or improperly assembled chromosome
            self.is_circular = False
            del boundries[-1]
            start_coord = end_frag[0].end
            end_coord = end_frag[0].start
            
        if self.is_circular:
            # For circular chromosomes, create fragment spanning origin
            # from last operon end -> chromosome end -> 0 -> first operon start
            f = Fragment([[boundries[-1].end, genome_length], [0, boundries[0].start]], operon_forward_start = boundries[-1].direction, operon_forward_end = boundries[0].direction)
            fragments.append(f)
        else:
            # For linear chromosomes, create separate fragments at each end
            f = Fragment([[boundries[-1].end, end_coord]], operon_forward_start = boundries[-1].direction)
            fragments.append(f)
            f = Fragment([ [start_coord, boundries[0].start] ], operon_forward_end = boundries[0].direction )
            fragments.append(f)
        
        # Create fragments for all regions between consecutive operons
        for i in range(0,len(boundries)-1):
            fragments.append(Fragment([[boundries[i].end, boundries[i+1].start]], operon_forward_start = boundries[i].direction, operon_forward_end = boundries[i+1].direction))
        return fragments

    def populate_fragments_from_chromosome(self, fragments, max_bases_from_ends):
        """
        Extract actual DNA sequences for all fragments from chromosome.
        
        Fills in the sequence attribute for each Fragment object. Optionally
        subsamples long fragments by taking only sequence from ends.
        
        Args:
            fragments (list): List of Fragment objects with coordinates
            max_bases_from_ends (int): If set, only extract this many bases
                                      from each end of long fragments
        """        
        # Extract sequence for each coordinate range in each fragment
        for f in fragments:
            for c in f.coords:
               f.sequence += self.chromosome.seq[(c[0]):(c[1])]
        
        # Optionally subsample long fragments from ends only       
        for f in fragments:  
           if max_bases_from_ends is not None and max_bases_from_ends > 0 and (len(f.sequence) > 2*max_bases_from_ends):
           	   # Take max_bases_from_ends from start and end, join with NNN
               filtered_seq = f.sequence[0:max_bases_from_ends] +"NNN" + f.sequence[len(f.sequence) - max_bases_from_ends:len(f.sequence)]
               f.sequence =filtered_seq
        return
       
