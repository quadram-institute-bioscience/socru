import re
class Operon:
    
    def __init__(self, start, end, direction ):
        self.start = start
        self.end = end
        self.direction = direction
        
    def __str__(self):
        return str(self.start) + "\t" + str(self.end) + "\t" + str(self.direction)
      