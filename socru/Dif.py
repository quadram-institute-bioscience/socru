from socru.Database import Database
from socru.Blast import Blast
from socru.FilterBlast import FilterBlast

class Dif:
    def __init__(self, dif_fasta, directory_of_fasta_files, threads, verbose, min_bit_score = 30, min_alignment_length = 25):
        self.dif_fasta = dif_fasta
        self.directory_of_fasta_files = directory_of_fasta_files
        self.threads = threads 
        self.min_bit_score = min_bit_score
        self.min_alignment_length = min_alignment_length
        self.verbose = verbose
        
        self.fragment_with_dif = 1
        self.dif_orientation = False
        self.find_dif()
    
    def find_dif(self):
        if self.verbose:
            print("Finding fragment FASTA containing dif (near terminus of replication)")
        blastdb = Database(self.directory_of_fasta_files, self.verbose)
        blast = Blast(blastdb.db_prefix, self.threads, self.verbose, task = 'blastn',  word_size = 11)
        blast_results = blast.run_blast(self.dif_fasta)
        fb = FilterBlast(blast_results, self.min_bit_score, self.min_alignment_length, self.verbose)
        top_result = fb.return_top_result()
        
        if top_result is None:
            self.forward_orientation = None
            self.fragment_with_dif = 1
        else:
            self.fragment_with_dif  = top_result.subject
            if top_result.is_forward():
                self.forward_orientation = True
            else:
                self.forward_orientation = False
                
            if self.verbose:
                print("Found dif on fragment:\t" + str(self.fragment_with_dif))
                print("dif in forward orientation:\t" + str(self.forward_orientation))
            
