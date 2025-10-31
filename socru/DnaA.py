"""
DnaA origin identification in fragment databases.

This module identifies which fragment contains the dnaA gene (replication origin)
by BLASTing a dnaA sequence against fragment databases. This information is used
to orient genome structures consistently and validate fragment arrangements.

Classes:
    DnaA: Identifies dnaA-containing fragments
"""

from socru.Database import Database
from socru.Blast import Blast
from socru.FilterBlast import FilterBlast

class DnaA:
    """
    Identify fragment containing the dnaA replication origin.
    
    This class searches for the dnaA gene in a set of genome fragments by
    BLASTing a known dnaA sequence against the fragments. The dnaA gene marks
    the chromosomal origin of replication, which is used as a reference point
    for standardizing genome structure representations.
    
    Attributes:
        dnaa_fasta (str): Path to dnaA query FASTA file
        directory_of_fasta_files (str): Directory containing fragment files
        threads (int): Number of CPU threads for BLAST
        min_bit_score (int): Minimum BLAST bit score threshold
        min_alignment_length (int): Minimum alignment length threshold
        verbose (bool): Enable verbose output
        fragment_with_dnaa (int): Fragment number containing dnaA
        dnaa_orientation (bool): Orientation of dnaA (forward/reverse)
    """
    def __init__(self, dnaa_fasta, directory_of_fasta_files, threads, verbose, min_bit_score = 100, min_alignment_length = 1000):
        """
        Initialize DnaA finder and search for dnaA.
        
        Args:
            dnaa_fasta (str): Path to dnaA query sequence
            directory_of_fasta_files (str): Directory with fragment files
            threads (int): Number of threads for BLAST
            verbose (bool): Enable verbose output
            min_bit_score (int): Minimum bit score (default 100)
            min_alignment_length (int): Minimum alignment length (default 1000)
        """
        self.dnaa_fasta = dnaa_fasta
        self.directory_of_fasta_files = directory_of_fasta_files
        self.threads = threads 
        self.min_bit_score = min_bit_score
        self.min_alignment_length = min_alignment_length
        self.verbose = verbose
        
        # Default values
        self.fragment_with_dnaa = 1
        self.dnaa_orientation = False
        
        # Search for dnaA
        self.find_dnaa()
    
    def find_dnaa(self):
        """
        BLAST dnaA sequence against fragments to find origin location.
        
        Creates a BLAST database from fragments, searches with dnaA query,
        and identifies which fragment contains the best match. Also determines
        the orientation of dnaA (forward or reverse strand).
        """
        if self.verbose:
            print("Finding fragment FASTA containing dnaA (origin of replication)")
        
        # Create BLAST database from fragments
        blastdb = Database(self.directory_of_fasta_files, self.verbose)
        
        # Run BLAST search with sensitive parameters
        blast = Blast(blastdb.db_prefix, self.threads, self.verbose, task = 'blastn',  word_size = 11)
        blast_results = blast.run_blast(self.dnaa_fasta)
        
        # Filter results and get top hit
        fb = FilterBlast(blast_results, self.min_bit_score, self.min_alignment_length, self.verbose)
        top_result = fb.return_top_result()
        
        if top_result is None:
            # No match found - use defaults
            self.forward_orientation = None
            self.fragment_with_dnaa = 1
        else:
            # Record fragment containing dnaA
            self.fragment_with_dnaa  = top_result.subject
            
            # Determine orientation
            if top_result.is_forward():
                self.forward_orientation = True
            else:
                self.forward_orientation = False
                
            if self.verbose:
                print("Found dnaA on fragment:\t" + str(self.fragment_with_dnaa))
                print("dnaA in forward orientation:\t" + str(self.forward_orientation))
