"""
Fragment arrangement validation relative to origin and terminus.

This module validates that genome fragments are arranged correctly relative
to the replication origin (dnaA) and terminus (dif). It checks that rRNA
operons are oriented as expected: transcribing away from the origin toward
the terminus on both replichores.

Classes:
    ValidateFragments: Validates fragment and operon orientation
"""

class ValidateFragments:
    """
    Validate genome fragment arrangement and operon orientations.
    
    This class checks that fragments are arranged correctly in relation to
    the chromosome's replication origin and terminus. In bacteria, rRNA
    operons should be transcribed in the direction of replication fork
    movement (away from origin, toward terminus).
    
    Attributes:
        fragments (list): List of Fragment objects to validate
    """
    def __init__(self, fragments):
        """
        Initialize validator with fragments to check.
        
        Args:
            fragments (list): List of Fragment objects
        """
        self.fragments = fragments
    
    def find_dnaa(self, frags):
        """
        Find the fragment containing the dnaA origin marker.
        
        Args:
            frags (list): List of Fragment objects
            
        Returns:
            int: Index of dnaA-containing fragment, or -1 if not found
        """
        for i, frag in enumerate(frags):
            if frag.dna_A:
                return i
        
        print("DnaA couldnt be found")
        return -1
        
    def find_dif(self, frags):
        """
        Find the fragment containing the dif terminus marker.
        
        Args:
            frags (list): List of Fragment objects
            
        Returns:
            int: Index of dif-containing fragment, or -1 if not found
        """
        for i, frag in enumerate(frags):
            if frag.dif:
                return i
        
        print("Dif couldnt be found")
        return -1
        
    def reorientate_ori(self):
        """
        Rotate fragment list to start with the origin (dnaA).
        
        This provides a consistent reference point for validation,
        starting from the replication origin.
        
        Returns:
            list: Fragments reordered to start with dnaA
        """
        dnaa_index = self.find_dnaa(self.fragments)
        if dnaa_index == -1:
            return self.fragments

        # Rotate to put dnaA at position 0
        reorientated_frags = self.fragments[dnaa_index: len(self.fragments)] + self.fragments[0: dnaa_index]
        return reorientated_frags
        
    def validate(self):
        """
        Validate that fragments and operons are correctly oriented.
        
        Checks that:
        1. dnaA and dif markers are present
        2. Operons from origin to terminus are forward-oriented (-->)
        3. Operons from terminus back to origin are reverse-oriented (<--)
        
        This validates the biological expectation that rRNA operons are
        transcribed in the direction of replication fork movement.
        
        Returns:
            bool: True if arrangement is valid, False otherwise
        """        
        # Reorient so that origin (dnaA) is at index zero
        reorientated_frags = self.reorientate_ori()
        
        # Find positions of origin and terminus
        dnaa_index = self.find_dnaa(reorientated_frags)
        dif_index = self.find_dif(reorientated_frags)

        # Must have both markers
        if dnaa_index == -1 or dif_index == -1:
            return False
            
        # Walk from the origin in both directions until the terminus
        # checking operon orientations

        # Forward walk (right replichore): origin -> terminus
        # All operons should be forward (-->)
        for i in range(dnaa_index + 1, len(self.fragments)):
            # If the operon is not forward then it's invalid
            if not reorientated_frags[i].operon_forward_start:
                return False
            if i >= dif_index:
                # We have reached the terminus
                break
       
        # Reverse walk (left replichore): going backwards from end to terminus
        # All operons should be reverse (<--)
        for i in range( len(self.fragments) -1, dif_index, -1):
            # If the operon is forward (should be reverse) then invalid
            if reorientated_frags[i].operon_forward_start:
                return False    
            if i <= dif_index:
                # We have reached the terminus
                break

        # All checks passed
        return True
        
        