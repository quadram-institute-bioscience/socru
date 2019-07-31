
class Fragment:
    def __init__(self, coords, sequence = "", number = 0, reversed_frag = False, dna_A = False, operon_forward_start = True, operon_forward_end = True):
        self.coords = coords
        self.sequence = sequence
        self.number = number
        self.reversed_frag = reversed_frag
        self.dna_A = dna_A
        
        self.operon_forward_start = operon_forward_start
        self.operon_forward_end = operon_forward_end
        
    def num_bases(self):
        return len(self.sequence)
        
    def output_filename(self):
        return str(self.number) + '.fa'
    
    def __str__(self):
        seqname = str(self.number)+ " " +"__".join([str(i[0]) + "_" + str(i[1]) for i in self.coords])
        return seqname
        