from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

import os
import re

class FragmentFiles:
    def __init__(self, fragments, output_directory, fragment_order = None):
        self.fragments = fragments
        self.output_directory = output_directory
        self.fragment_order = fragment_order
        self.ordered_fragments = self.fragments_with_largest_first()
        self.output_filenames = []
        
    def fragments_with_largest_first(self):
        fragment_sizes = [f.num_bases() for f in self.fragments]
        largest_size = max(fragment_sizes)
        
        largest_index = 0
        for f in range(0,len(fragment_sizes)):
            if largest_size == fragment_sizes[f]:
                largest_index = f
                break

        reordered_fragments = self.fragments[largest_index:len(self.fragments)] + self.fragments[0:largest_index]
        
        if self.fragment_order is None:
            for i in range(0,len(reordered_fragments)):
                reordered_fragments[i].number = i+1
        else:
            input_order = self.split_fragment_order()
            for i in range(0,len(reordered_fragments)):
                # reverse complement if prime
                m = re.match(r"([\d]+)'", input_order[i])
                if m:
                    reordered_fragments[i].number = m.group(1)
                    reordered_fragments[i].sequence = reordered_fragments[i].sequence.reverse_complement()
                    
                else:
                    reordered_fragments[i].number = input_order[i]
                
        return reordered_fragments
      
    def create_fragment_fastas(self):
        for f in self.ordered_fragments:
            record = [SeqRecord(f.sequence, str(f) , '', '')]
            outname = os.path.join(self.output_directory, f.output_filename())
            self.output_filenames.append(outname)
            SeqIO.write(record, outname, "fasta")
    
    def split_fragment_order(self):
        return self.fragment_order.split('-')
        