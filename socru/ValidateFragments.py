class ValidateFragments:
    def __init__(self, fragments):
        self.fragments = fragments
    
    def find_dnaa(self, frags):
        for i, frag in enumerate(frags):
            if frag.dna_A:
                return i
        
        print("DnaA couldnt be found")
        return -1
        
    def find_dif(self, frags):
        for i, frag in enumerate(frags):
            if frag.dif:
                return i
        
        print("Dif couldnt be found")
        return -1
        
    def reorientate_ori(self):
        dnaa_index = self.find_dnaa(self.fragments)
        if dnaa_index == -1:
            return self.fragments

        reorientated_frags = self.fragments[dnaa_index: len(self.fragments)] + self.fragments[0: dnaa_index]
        return reorientated_frags
        
    def validate(self):        
        # reorinate so that ori is at index zero
        reorientated_frags = self.reorientate_ori()
        
        dnaa_index = self.find_dnaa(reorientated_frags)
        dif_index = self.find_dif(reorientated_frags)

        if dnaa_index == -1 or dif_index == -1:
            return False
            
        # Walk from the origin in both directions until the terminus

        # forward walk
        for i in range(dnaa_index + 1, len(self.fragments)):
            #  direction should be -->
            # if the operon is not forward then its invalid
            if not reorientated_frags[i].operon_forward_start:
                return False
            if i >= dif_index:
                # we have reached the end
                break
       
        # reverse walk
        for i in range( len(self.fragments) -1, dif_index, -1):
            #  direction should be <--
            if reorientated_frags[i].operon_forward_start:
                return False    
            if i <= dif_index:
                # we have reached the end
                break

        return True
        
        