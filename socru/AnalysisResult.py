"""
Structured analysis result representation for Socru genome typing.

This module defines dataclasses for representing the complete output of a Socru
analysis run, including per-fragment BLAST details, operon positions, quality
control flags, and confidence scoring. Results can be serialized to JSON for
downstream consumption by visualization and reporting tools.

Classes:
    FragmentResult: Per-fragment analysis metrics
    OperonResult: Per-operon position and direction
    QCFlag: A quality control flag or warning
    AnalysisResult: Complete result for one genome analysis
"""

import json
from dataclasses import asdict, dataclass, field
from typing import List, Optional


@dataclass
class FragmentResult:
    """Per-fragment analysis result with BLAST alignment details.

    Attributes:
        number: Fragment identifier from database matching
        reversed: Whether fragment is in reverse orientation
        is_dnaA: Whether this fragment contains the replication origin
        is_dif: Whether this fragment contains the replication terminus
        length: Length of the fragment sequence in base pairs
        coords: List of (start, end) coordinate tuples
        blast_identity: Percent identity of best BLAST hit
        blast_alignment_length: Alignment length of best BLAST hit
        blast_bit_score: Bit score of best BLAST hit
        blast_e_value: E-value of best BLAST hit
        blast_mismatches: Number of mismatches in best BLAST hit
        blast_subject: Subject sequence name from best BLAST hit
    """

    number: int
    reversed: bool
    is_dnaA: bool
    is_dif: bool
    length: int
    coords: list
    blast_identity: Optional[float] = None
    blast_alignment_length: Optional[int] = None
    blast_bit_score: Optional[float] = None
    blast_e_value: Optional[float] = None
    blast_mismatches: Optional[int] = None
    blast_subject: Optional[str] = None


@dataclass
class OperonResult:
    """Per-operon position and direction result.

    Attributes:
        start: Genomic start coordinate
        end: Genomic end coordinate
        direction: Operon transcription direction ('forward' or 'reverse')
    """

    start: int
    end: int
    direction: str


@dataclass
class QCFlag:
    """A quality control flag or warning raised during analysis.

    Attributes:
        code: Short machine-readable flag code (e.g. 'LOW_IDENTITY')
        severity: Either 'warning' or 'error'
        message: Human-readable description of the issue
        details: Optional additional context
    """

    code: str
    severity: str
    message: str
    details: Optional[str] = None


@dataclass
class AnalysisResult:
    """Complete structured result for one genome analysis.

    Aggregates all information produced during a single Socru run: the genome
    structure type, quality assessment, per-fragment BLAST metrics, operon
    positions, confidence score, and any quality control flags.

    Attributes:
        genome_file: Path to the input genome FASTA file
        genome_length: Length of the chromosome in base pairs
        is_circular: Whether the chromosome was treated as circular
        num_operons: Number of rRNA operons identified
        gs_type: Genome structure type string (e.g. 'GS1.3')
        quality: Traffic-light quality label (GREEN, AMBER, RED)
        is_novel: Whether this is a previously unseen profile
        fragment_pattern: Tab-delimited fragment order string
        orientation_binary: Binary encoding of fragment orientations
        confidence_score: Confidence score from 0 to 100
        fragments: Per-fragment analysis results
        operons: Per-operon position results
        qc_flags: Quality control flags raised during analysis
        validation_passed: Whether fragment arrangement passed validation
        operon_direction_string: Human-readable operon direction summary
    """

    genome_file: str
    genome_length: int
    is_circular: bool
    num_operons: int
    gs_type: str
    quality: str
    is_novel: bool
    fragment_pattern: str
    orientation_binary: int
    confidence_score: float
    fragments: List[FragmentResult] = field(default_factory=list)
    operons: List[OperonResult] = field(default_factory=list)
    qc_flags: List[QCFlag] = field(default_factory=list)
    validation_passed: bool = True
    operon_direction_string: str = ''
    novelty_assessment: Optional[dict] = None

    def to_dict(self) -> dict:
        """Convert to a plain dictionary suitable for JSON serialization.

        Returns:
            dict: Nested dictionary representation of the result.
        """
        return asdict(self)

    def to_json(self, indent: int = 2) -> str:
        """Serialize the result to a JSON string.

        Args:
            indent: Number of spaces for JSON indentation.

        Returns:
            str: JSON string representation of the result.
        """
        return json.dumps(self.to_dict(), indent=indent)
