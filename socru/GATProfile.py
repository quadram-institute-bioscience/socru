import re
class GATProfile:
    
    def __init__(self, verbose, gat_number = 0, fragments = [], orientation_number = 0, dnaA_fragment_number = 3, dif_fragment_number = 1 ):
        self.gat_number = gat_number
        self.fragments = fragments
        self.verbose = verbose
        self.orientation_number = orientation_number
        self.dnaA_fragment_number = dnaA_fragment_number
        self.dif_fragment_number = dif_fragment_number
        
    def order(self):
        m = re.search(r'([\d]+)\.[\d]+$', str(self.gat_number) )
        
        if m:
            return int(m.group(1)) 
        else:
            return 0
        
    def fragment_str(self):
        return "\t".join(self.fragments)
        
    def invert_fragments(self):
        inverted = []
        for f in self.fragments:
            m = re.match(r"([\d]+)'", f)
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
        
    def reorder_fragment_objects_based_on_fragment_name_array(self, fragment_objects):
        reordered_fragment_objects = []
        # inefficient but small searches
        # need to check for unmatched
        for frag_name in self.fragments:
            frag_name_orientationless = frag_name
            frag_name_reversed = False
            m = re.match(r"([\d]+)'", frag_name)
            if m:
                frag_name_orientationless = m.group(1)
                frag_name_reversed = True

            for frag_obj in fragment_objects:
                if str(frag_obj.number) == str(frag_name_orientationless):
                    if frag_name_reversed:
                        frag_obj.reversed_frag = True
                    else:
                        frag_obj.reversed_frag = False
                    reordered_fragment_objects.append(frag_obj)
                    
        return reordered_fragment_objects
        
        
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

    def remove_orientation(self, fragments):
        orientationless = []
        for f in fragments:
            m = re.match(r"([\d]+)'", f)
            if m:
                orientationless.append(m.group(1))
            else:
                orientationless.append(f)
        return orientationless
        
    def orientationless_fragments(self):
        return self.remove_orientation(self.fragments)
        
    def inverted_orientationless_fragments(self):
        inverted = self.invert_fragments()
        return self.remove_orientation(inverted)
		
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
            m = re.match(r"([\d]+)'", f)
            if m:
                bits_to_shift = int(m.group(1)) -1
                total += 1 << bits_to_shift
        return total 
        
    
		