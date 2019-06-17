from socru.Database import Database
from socru.Blast import Blast
from socru.FilterBlast import FilterBlast

class DnaA:
    def __init__(self, dnaa_fasta, directory_of_fasta_files, threads, verbose, min_bit_score = 100, min_alignment_length = 1000):
        self.dnaa_fasta = dnaa_fasta
        self.directory_of_fasta_files = directory_of_fasta_files
        self.threads = threads 
        self.min_bit_score = min_bit_score
        self.min_alignment_length = min_alignment_length
        self.verbose = verbose
        
        self.fragment_with_dnaa = 1
        self.dnaa_orientation = False
        self.find_dnaa()
    
    def find_dnaa(self):
        blastdb = Database(self.directory_of_fasta_files, self.verbose)
        blast = Blast(blastdb.db_prefix, self.threads, self.verbose, task = 'blastn',  word_size = 11)
        blast_results = blast.run_blast(self.dnaa_fasta)
        fb = FilterBlast(blast_results, self.min_bit_score, self.min_alignment_length, self.verbose)
        top_result = fb.return_top_result()
        
        if top_result is None:
            self.forward_orientation = None
            self.fragment_with_dnaa = 1
        else:
            self.fragment_with_dnaa  = top_result.subject
            if top_result.is_forward():
                self.forward_orientation = True
            else:
                self.forward_orientation = False
            
