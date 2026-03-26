"""
BLAST result representation and parsing.

This module defines the BlastResult class which represents a single BLAST hit
in tabular format (outfmt 6). It stores all standard BLAST fields and provides
methods for determining hit orientation (forward/reverse).

Classes:
    BlastResult: Represents a single BLAST alignment result
"""

from __future__ import annotations

from typing import Union


class BlastResult:
    """
    Represents a single BLAST hit from tabular output.
    
    Stores all standard BLAST tabular output fields (format 6) including
    sequence identifiers, alignment statistics, and coordinates. Provides
    methods to determine hit orientation.
    
    Attributes:
        query_name (str): Name of query sequence
        subject (str): Name of subject (database) sequence
        identity (float): Percent identity of alignment
        alignment_length (int): Length of alignment
        mismatches (int): Number of mismatches
        gap_openings (int): Number of gap openings
        query_start (int): Start position in query
        query_end (int): End position in query
        subject_start (int): Start position in subject
        subject_end (int): End position in subject
        e_value (float): E-value (expectation value)
        bit_score (float): Bit score
    """
    def __init__(self, query_name: str, subject: str, identity: Union[str, float], alignment_length: Union[str, int], mismatches: Union[str, int], gap_openings: Union[str, int], query_start: Union[str, int], query_end: Union[str, int], subject_start: Union[str, int], subject_end: Union[str, int], e_value: Union[str, float], bit_score: Union[str, float]) -> None:
        """
        Initialize BlastResult from BLAST tabular output fields.
        
        Args:
            query_name (str): Query sequence name
            subject (str): Subject sequence name
            identity: Percent identity (converted to float)
            alignment_length: Alignment length (converted to int)
            mismatches: Number of mismatches (converted to int)
            gap_openings: Number of gap openings (converted to int)
            query_start: Query start position (converted to int)
            query_end: Query end position (converted to int)
            subject_start: Subject start position (converted to int)
            subject_end: Subject end position (converted to int)
            e_value: E-value (converted to float)
            bit_score: Bit score (converted to float)
        """
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
        
        # Example forward hit:
        # query7	7	100.000	2520	0	0	1	2520	2221	4740	0.04654
        # Example reverse hit (subject_start > subject_end):
        # query7	7	100.000	2520	0	0	1	2520	4740	2221	0.04654
    
    def is_forward(self) -> bool:
        """
        Determine if the hit is in forward or reverse orientation.
        
        In BLAST output, if subject_start > subject_end, the hit is on the
        reverse strand. This is important for determining fragment orientation.
        
        Returns:
            bool: True if forward orientation, False if reverse
        """
        if self.subject_start > self.subject_end:
            return False
        else:
            return True
			
    def __str__(self) -> str:
        """
        Format as tab-delimited string matching BLAST tabular output.
        
        Returns:
            str: Tab-separated values of all BLAST fields
        """
        return "\t".join([str(self.query_name ), str(self.subject ), str(self.identity ), 
                   str(self.alignment_length ), str(self.mismatches ), str(self.gap_openings ), 
                   str(self.query_start ), str(self.query_end ), str(self.subject_start ), 
                   str(self.subject_end ), str(self.e_value ), str(self.bit_score)])
            