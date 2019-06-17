from Bio import SeqIO
from socru.Fragment import Fragment
import re
import gzip

class Fasta:
    def __init__(self, input_file,verbose, is_circular=True):
        self.input_file = input_file
        self.is_circular = is_circular
        self.verbose = verbose
        self.chromosome = self.get_chromosome_from_fasta()
    
    def get_chromosome_from_fasta(self):
        # find largest contig and ignore the rest
        largest_contig = None
		
        m = re.search(".gz$", self.input_file)
        if m:
            with gzip.open(self.input_file, "rt") as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    largest_contig = self.largest_contig_check( largest_contig, record)
        else:
            for record in SeqIO.parse(self.input_file, "fasta"):
                largest_contig = self.largest_contig_check( largest_contig, record)

        return largest_contig
        
    def fragment_number(self):
        return int(self.chromosome.id)
        
    def largest_contig_check(self, largest_contig, record):
        sequence_length  = len(record.seq)
        if largest_contig is None \
           or sequence_length > len(largest_contig.seq):
           largest_contig = record
        return largest_contig
        
    def calc_fragment_coords(self, boundries):
        genome_length = len(self.chromosome.seq)
        fragments = []
        
        start_coord = 0
        end_coord = genome_length
        
        # check to see if one of the fragments goes over the end
        end_frag = [ b for b in boundries if b[0] > b[1] ]
        if len(end_frag) > 0:
            self.is_circular = False
            del boundries[-1]
            start_coord = end_frag[0][1]
            end_coord = end_frag[0][0]
            
        if self.is_circular:
            # first - we assume its circular
            f = Fragment([[boundries[-1][1],genome_length], [0, boundries[0][0]]])
            fragments.append(f)
        else:
            f = Fragment([[boundries[-1][1],end_coord]])
            fragments.append(f)
            f = Fragment([ [start_coord, boundries[0][0]] ])
            fragments.append(f)
        
        for i in range(0,len(boundries)-1):
            fragments.append(Fragment([[boundries[i][1], boundries[i+1][0]]]))
        return fragments

    def populate_fragments_from_chromosome(self, fragments, max_bases_from_ends):        
        for f in fragments:
            for c in f.coords:
               f.sequence += self.chromosome.seq[(c[0]):(c[1])]
               
        for f in fragments:  
           if max_bases_from_ends is not None and max_bases_from_ends > 0 and (len(f.sequence) > 2*max_bases_from_ends):
           	   #subsample from start and end
               filtered_seq = f.sequence[0:max_bases_from_ends] +"NNN" + f.sequence[len(f.sequence) - max_bases_from_ends:len(f.sequence)]
               f.sequence =filtered_seq
        return
       
