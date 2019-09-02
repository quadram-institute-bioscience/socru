
class Fragment:
    def __init__(self, coords, sequence = "", number = 0, reversed_frag = False, dna_A = False, operon_forward_start = True, operon_forward_end = True, dif = False):
        self.coords = coords
        self.sequence = sequence
        self.number = number
        self.reversed_frag = reversed_frag
        self.dna_A = dna_A
        self.dif = dif
        
        self.operon_forward_start = operon_forward_start
        self.operon_forward_end = operon_forward_end
        
    def num_bases(self):
        return len(self.sequence)
        
    def output_filename(self):
        return str(self.number) + '.fa'
        
    # <-- 3(Ori) 
    def operon_direction_str(self):
        output_str = ''
        if self.operon_forward_start:
            output_str += '--> '
        else:
            output_str += '<-- '

        if self.reversed_frag:
            output_str += str(self.number)+"'"
        else:
            output_str += str(self.number)
            
        if self.dna_A:
            output_str += '(Ori)'
			
        if self.dif:
            output_str += '(Ter)'
            
        return output_str
        
    def __str__(self):
        seqname = str(self.number)+ " " +"__".join([str(i[0]) + "_" + str(i[1]) for i in self.coords])
        return seqname
        