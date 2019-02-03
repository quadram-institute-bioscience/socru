
class BlastResult:
    def __init__(self, query_name, subject, identity, alignment_length, mismatches, gap_openings, query_start, query_end, subject_start, subject_end, e_value, bit_score):
        self.query_name = query_name
        self.subject = subject
        self.identity = float(identity)
        self.alignment_length = int(alignment_length)
        self.mismatches =  int(mismatches)
        self.gap_openings =  int(gap_openings)
        self.query_start =  int(query_start)
        self.query_end =  int(query_end)
        self.subject_start =  int(subject_start)
        self.subject_end =  int(subject_end)
        self.e_value = float(e_value)
        self.bit_score = float(bit_score)
        
        #query7	7	100.000	2520	0	0	1	2520	2221	4740	0.04654
        #reverse
        #query7	7	100.000	2520	0	0	1	2520	4740	2221	0.04654
    
    def is_forward(self):
        if self.subject_start > self.subject_end:
            return False
        else:
            return True
			
    def __str__(self):
       return "\t".join([str(self.query_name ), str(self.subject ), str(self.identity ), 
                   str(self.alignment_length ), str(self.mismatches ), str(self.gap_openings ), 
                   str(self.query_start ), str(self.query_end ), str(self.subject_start ), 
                   str(self.subject_end ), str(self.e_value ), str(self.bit_score)])
            