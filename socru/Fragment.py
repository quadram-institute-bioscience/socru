"""
Fragment representation for inter-operon genomic regions.

This module defines the Fragment class which represents a genomic region between
rRNA operons. Fragments store their coordinates, sequences, orientation, and
special markers (dnaA origin, dif terminator) used for genome structure typing.

Classes:
    Fragment: Represents an inter-operon genomic fragment
"""

from __future__ import annotations

from typing import Any


class Fragment:
    """
    Represents a genomic fragment between rRNA operons.

    Fragments are the key units analyzed in Socru. Each fragment represents
    the DNA sequence between two adjacent rRNA operons. Fragments can span
    multiple coordinate ranges (e.g., wrapping around chromosome origin) and
    can be in forward or reverse orientation.

    Attributes:
        coords (list): List of [start, end] coordinate pairs
        sequence (str): DNA sequence of the fragment
        number (int): Fragment identifier from database matching
        reversed_frag (bool): True if fragment is reverse complemented
        dna_A (bool): True if this fragment contains the dnaA origin
        dif (bool): True if this fragment contains the dif terminator
        operon_forward_start (bool): Orientation of operon at fragment start
        operon_forward_end (bool): Orientation of operon at fragment end
    """
    def __init__(self, coords: list[list[int]], sequence: Any = "", number: int = 0, reversed_frag: bool = False, dna_A: bool = False, operon_forward_start: bool = True, operon_forward_end: bool = True, dif: bool = False) -> None:
        """
        Initialize a Fragment object.

        Args:
            coords (list): List of [start, end] coordinate pairs
            sequence (str): DNA sequence (default empty)
            number (int): Fragment database identifier (default 0)
            reversed_frag (bool): Fragment is reversed (default False)
            dna_A (bool): Contains dnaA origin marker (default False)
            operon_forward_start (bool): Start operon forward strand (default True)
            operon_forward_end (bool): End operon forward strand (default True)
            dif (bool): Contains dif terminator marker (default False)
        """
        self.coords = coords
        self.sequence = sequence
        self.number = number
        self.reversed_frag = reversed_frag
        self.dna_A = dna_A
        self.dif = dif

        self.operon_forward_start = operon_forward_start
        self.operon_forward_end = operon_forward_end

    def num_bases(self) -> int:
        """
        Get the length of the fragment sequence.

        Returns:
            int: Number of bases in the sequence
        """
        return len(self.sequence)

    def output_filename(self) -> str:
        """
        Generate output filename for this fragment.

        Returns:
            str: Filename in format "{number}.fa"
        """
        return str(self.number) + '.fa'

    def operon_direction_str(self) -> str:
        """
        Create a string representation of fragment with operon orientations.

        Generates a human-readable string showing:
        - Start operon direction (-->) or (<--)
        - Fragment number (with ' if reversed)
        - Special markers: (Ori) for dnaA, (Ter) for dif

        Example: "--> 3'(Ori)" means forward start operon, reversed fragment 3 with origin

        Returns:
            str: Formatted direction string
        """
        output_str = ''
        # Show start operon direction
        if self.operon_forward_start:
            output_str += '--> '
        else:
            output_str += '<-- '

        # Show fragment number with orientation
        if self.reversed_frag:
            output_str += str(self.number)+"'"
        else:
            output_str += str(self.number)

        # Add special position markers
        if self.dna_A:
            output_str += '(Ori)'

        if self.dif:
            output_str += '(Ter)'

        return output_str

    def __str__(self) -> str:
        """
        String representation for debugging/logging.

        Returns:
            str: Fragment number followed by coordinates
        """
        seqname = str(self.number)+ " " +"__".join([str(i[0]) + "_" + str(i[1]) for i in self.coords])
        return seqname
