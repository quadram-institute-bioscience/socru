"""
Dif terminus identification in fragment databases.

This module identifies which fragment contains the dif site (replication terminus)
by BLASTing a dif sequence against fragment databases. The dif site marks the
region where replication forks meet, which helps validate genome structure and
operon arrangements.

Classes:
    Dif: Identifies dif-containing fragments
"""

from socru.Database import Database
from socru.Blast import Blast
from socru.FilterBlast import FilterBlast

class Dif:
    """
    Identify fragment containing the dif replication terminus.
    
    This class searches for the dif site in a set of genome fragments by
    BLASTing a known dif sequence. The dif site marks the chromosomal
    terminus of replication (opposite the origin), where the two replication
    forks meet. This is used to validate that rRNA operons are correctly
    oriented relative to replication direction.
    
    Attributes:
        dif_fasta (str): Path to dif query FASTA file
        directory_of_fasta_files (str): Directory containing fragment files
        threads (int): Number of CPU threads for BLAST
        min_bit_score (int): Minimum BLAST bit score threshold
        min_alignment_length (int): Minimum alignment length threshold
        verbose (bool): Enable verbose output
        fragment_with_dif (int): Fragment number containing dif
        dif_orientation (bool): Orientation of dif (forward/reverse)
    """
    def __init__(self, dif_fasta, directory_of_fasta_files, threads, verbose, min_bit_score = 30, min_alignment_length = 25):
        """
        Initialize Dif finder and search for dif site.
        
        Args:
            dif_fasta (str): Path to dif query sequence
            directory_of_fasta_files (str): Directory with fragment files
            threads (int): Number of threads for BLAST
            verbose (bool): Enable verbose output
            min_bit_score (int): Minimum bit score (default 30, lower than dnaA)
            min_alignment_length (int): Minimum alignment length (default 25, dif is shorter)
        """
        self.dif_fasta = dif_fasta
        self.directory_of_fasta_files = directory_of_fasta_files
        self.threads = threads 
        self.min_bit_score = min_bit_score
        self.min_alignment_length = min_alignment_length
        self.verbose = verbose
        
        # Default values
        self.fragment_with_dif = 0
        self.dif_orientation = False
        
        # Search for dif
        self.find_dif()
    
    def find_dif(self):
        """
        BLAST dif sequence against fragments to find terminus location.
        
        Creates a BLAST database from fragments, searches with dif query,
        and identifies which fragment contains the best match. Dif is much
        shorter than dnaA (~28bp), so lower thresholds are used. If not found,
        defaults to fragment 1 (usually the largest, often containing terminus).
        """
        if self.verbose:
            print("Finding fragment FASTA containing dif (near terminus of replication)")
        
        # Create BLAST database from fragments
        blastdb = Database(self.directory_of_fasta_files, self.verbose)
        
        # Run BLAST search with sensitive parameters for short sequence
        blast = Blast(blastdb.db_prefix, self.threads, self.verbose, task = 'blastn',  word_size = 11)
        blast_results = blast.run_blast(self.dif_fasta)
        
        # Filter results and get top hit
        fb = FilterBlast(blast_results, self.min_bit_score, self.min_alignment_length, self.verbose)
        top_result = fb.return_top_result()
        
        if top_result is None:
            # No match found - default to fragment 1
            # Fragment 1 is the largest, and usually contains the terminus
            self.forward_orientation = None
            self.fragment_with_dif = 1
            print("Dif cannot be found in\t" + self.directory_of_fasta_files)
        else:
            # Record fragment containing dif
            self.fragment_with_dif  = top_result.subject
            
            # Determine orientation
            if top_result.is_forward():
                self.forward_orientation = True
            else:
                self.forward_orientation = False
                
            if self.verbose:
                print("Found dif on fragment:\t" + str(self.fragment_with_dif))
                print("dif in forward orientation:\t" + str(self.forward_orientation))
