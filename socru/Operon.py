"""
Operon representation for ribosomal RNA operons.

This module defines a simple data class to represent ribosomal operons identified
by Barrnap. Each operon has genomic coordinates and a direction indicating whether
the rRNA genes are transcribed on the forward or reverse strand.

Classes:
    Operon: Represents a single ribosomal operon
"""

from __future__ import annotations

import re

class Operon:
    """
    Represents a ribosomal RNA operon in a bacterial genome.
    
    An operon is a cluster of rRNA genes (16S, 23S, 5S) transcribed together.
    The boundaries of operons define the inter-operon fragments that Socru analyzes.
    
    Attributes:
        start (int): Genomic start coordinate of the operon
        end (int): Genomic end coordinate of the operon
        direction (bool): True if operon is on forward strand, False if reverse
    """
    
    def __init__(self, start: int, end: int, direction: bool) -> None:
        """
        Initialize an Operon object.
        
        Args:
            start (int): Start coordinate in genome
            end (int): End coordinate in genome
            direction (bool): True for forward strand, False for reverse
        """
        self.start = start
        self.end = end
        # Forward - boolean True
        self.direction = direction
        
    def __str__(self) -> str:
        """
        String representation for debugging/logging.
        
        Returns:
            str: Tab-delimited string of start, end, and direction
        """
        return str(self.start) + "\t" + str(self.end) + "\t" + str(self.direction)
      