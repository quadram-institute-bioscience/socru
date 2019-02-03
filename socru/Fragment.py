
class Fragment:
    def __init__(self, coords):
        self.coords = coords
        self.sequence = ""
        self.number = 0
        
    def num_bases(self):
        return len(self.sequence)
        
    def output_filename(self):
        return str(self.number) + '.fa'
    
    def __str__(self):
        seqname = str(self.number)+ " " +"__".join([str(i[0]) + "_" + str(i[1]) for i in self.coords])
        return seqname
        