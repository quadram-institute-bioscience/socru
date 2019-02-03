import re
class GATProfile:
    
    def __init__(self, gat_number = 0, fragments = [], orientation_number = 0, dnaA_fragment_number = 3 ):
        self.gat_number = gat_number
        self.fragments = fragments
        self.orientation_number = orientation_number
        self.dnaA_fragment_number = dnaA_fragment_number
        
    def order(self):
        m = re.search('([\d]+)\.[\d]+$', str(self.gat_number) )
        
        if m:
            return int(m.group(1)) 
        else:
            return 0
        
    def fragment_str(self):
        return "\t".join(self.fragments)
        
    def invert_fragments(self):
        inverted = []
        for f in self.fragments:
            m = re.match("([\d]+)'", f)
            if m:
                inverted.append(m.group(1))
            else:
                inverted.append(f + '\'')

        inverted.reverse()
        return self.reorientate_list_to_start_with_one(inverted)
        
    def reorientate_list_to_start_with_one(self, raw_fragments):
        # find the index where fragment 1 or 1' is
        chop_index = 0
        for i in range(0,len(raw_fragments)):
            if raw_fragments[i] == '1' or raw_fragments[i] == "1'":
                chop_index = i
                continue
        
        reorientated = raw_fragments
        if chop_index > 0:
            reorientated = raw_fragments[chop_index:len(raw_fragments)] + raw_fragments[0:chop_index]
        return reorientated
        
    def inverted_fragment_str(self):
        return "\t".join(self.invert_fragments())
        
    def is_profile_in_correct_orientation(self):
        for f in self.fragments:
            if f == str(self.dnaA_fragment_number):
                return True
            elif f == str(self.dnaA_fragment_number) + "'":
                return False
        # dnaA hasnt been found so really should raise an exception
        return True
        
    def orientate_for_dnaA(self):
        if not self.is_profile_in_correct_orientation():
            self.fragments  = self.invert_fragments()
        
    def does_the_profile_match(self, gat_query):
        if gat_query.fragment_str() == self.fragment_str():
            return True
        else:
            return False
    
    def __str__(self):
        return self.fragment_str()
        
    def orientationless_fragments(self):
        orientationless = []
        for f in self.fragments:
            m = re.match("([\d]+)'", f)
            if m:
                orientationless.append(m.group(1))
            else:
                orientationless.append(f)
        return orientationless
		
        # 1' - 00000001
        # 2' - 00000010
        # 3' - 00000100
        # 4' - 00001000
        # 5' - 00010000
        # 6' - 00100000
        # 7' - 01000000
        # 8' - 10000000
    def orientation_binary(self):
        # If a fragment is inverted it gets a 1, otherwise zero. The fragments are then shifted by the fragment number and added.
        total = 0
        for f in self.fragments:
            m = re.match("([\d]+)'", f)
            if m:
                bits_to_shift = int(m.group(1)) -1
                total += 1 << bits_to_shift
        return total 
        
    
		