"""Configuration dataclasses for Socru analysis and database creation.

These dataclasses provide a clean, typed interface for configuring Socru
programmatically as a library, without requiring fabrication of argparse
Namespace objects. Backward compatibility with argparse is preserved via
the ``from_options()`` class methods.

Classes:
    SocruConfig: Configuration for running Socru analysis.
    SocruCreateConfig: Configuration for creating a new species database.
"""

from dataclasses import dataclass, field
from typing import List, Optional


@dataclass
class SocruConfig:
    """Configuration for running Socru analysis.

    Can be constructed directly for library use, or from an argparse
    Namespace via ``from_options()``.

    Attributes:
        input_files: List of input FASTA file paths to analyze.
        db_dir: Base directory for species databases.  ``None`` uses the
            bundled default.
        species: Species name used to locate the database when *db_dir*
            is not given.
        output_file: Path to main tab-separated output file.  ``None``
            prints to stdout.
        output_json: Path to JSON structured results output.
        output_svg: Path to SVG circular genome diagram.
        output_html: Path to self-contained HTML report.
        output_plot_file: Path to PDF plot of genome structure.
        output_operon_directions_file: Path to operon direction results.
        novel_profiles: Path to file for novel profile patterns.
        new_fragments: Path to file for unmatched novel fragments.
        top_blast_hits: Path to file for detailed BLAST hit output.
        min_bit_score: Minimum BLAST bit score for fragment matching.
        min_alignment_length: Minimum BLAST alignment length.
        max_bases_from_ends: If non-zero, only use this many bases from
            each fragment end for matching.
        threads: Number of threads for parallel processing.
        not_circular: When ``True``, do not assume the chromosome is
            circular.
        verbose: Enable verbose logging output.
    """

    input_files: List[str] = field(default_factory=list)
    db_dir: Optional[str] = None
    species: str = ""
    output_file: Optional[str] = None
    output_json: Optional[str] = None
    output_svg: Optional[str] = None
    output_html: Optional[str] = None
    output_plot_file: str = "genome_structure.pdf"
    output_operon_directions_file: str = "operon_directions.txt"
    novel_profiles: str = "profile.txt.novel"
    new_fragments: str = "novel_fragments.fa"
    top_blast_hits: Optional[str] = None
    min_bit_score: int = 100
    min_alignment_length: int = 100
    max_bases_from_ends: Optional[int] = None
    threads: int = 1
    not_circular: bool = False
    verbose: bool = False
    output_dir: Optional[str] = None

    @classmethod
    def from_options(cls, options: object) -> "SocruConfig":
        """Create a config from an argparse Namespace (backward compatible).

        Args:
            options: An argparse ``Namespace`` or any object whose
                attributes match the expected command-line option names.

        Returns:
            A fully populated ``SocruConfig`` instance.
        """
        return cls(
            input_files=getattr(options, "input_files", []),
            db_dir=getattr(options, "db_dir", None),
            species=getattr(options, "species", ""),
            output_file=getattr(options, "output_file", None),
            output_json=getattr(options, "output_json", None),
            output_svg=getattr(options, "output_svg", None),
            output_html=getattr(options, "output_html", None),
            output_plot_file=getattr(options, "output_plot_file", "genome_structure.pdf"),
            output_operon_directions_file=getattr(
                options, "output_operon_directions_file", "operon_directions.txt"
            ),
            novel_profiles=getattr(options, "novel_profiles", "profile.txt.novel"),
            new_fragments=getattr(options, "new_fragments", "novel_fragments.fa"),
            top_blast_hits=getattr(options, "top_blast_hits", None),
            min_bit_score=getattr(options, "min_bit_score", 100),
            min_alignment_length=getattr(options, "min_alignment_length", 100),
            max_bases_from_ends=getattr(options, "max_bases_from_ends", None),
            threads=getattr(options, "threads", 1),
            not_circular=getattr(options, "not_circular", False),
            verbose=getattr(options, "verbose", False),
            output_dir=getattr(options, "output_dir", None),
        )


@dataclass
class SocruCreateConfig:
    """Configuration for creating a new species database.

    Can be constructed directly for library use, or from an argparse
    Namespace via ``from_options()``.

    Attributes:
        input_file: Path to reference genome FASTA.
        output_directory: Directory where the new database will be written.
        fragment_order: Optional custom fragment numbering (e.g.
            ``"1-2-3-4-5-6-7"``).
        threads: Number of CPU threads for parallel processing.
        dnaa_fasta: Path to dnaA query sequence.  ``None`` uses the
            bundled default.
        dif_fasta: Path to dif query sequence.  ``None`` uses the
            bundled default.
        verbose: Enable verbose logging output.
        max_bases_from_ends: Optional fragment end trimming.
    """

    input_file: str = ""
    output_directory: str = ""
    fragment_order: Optional[str] = None
    threads: int = 1
    dnaa_fasta: Optional[str] = None
    dif_fasta: Optional[str] = None
    verbose: bool = False
    max_bases_from_ends: Optional[int] = None

    @classmethod
    def from_options(cls, options: object) -> "SocruCreateConfig":
        """Create a config from an argparse Namespace (backward compatible).

        Args:
            options: An argparse ``Namespace`` or any object whose
                attributes match the expected command-line option names.

        Returns:
            A fully populated ``SocruCreateConfig`` instance.
        """
        return cls(
            input_file=getattr(options, "input_file", ""),
            output_directory=getattr(options, "output_directory", ""),
            fragment_order=getattr(options, "fragment_order", None),
            threads=getattr(options, "threads", 1),
            dnaa_fasta=getattr(options, "dnaa_fasta", None),
            dif_fasta=getattr(options, "dif_fasta", None),
            verbose=getattr(options, "verbose", False),
            max_bases_from_ends=getattr(options, "max_bases_from_ends", None),
        )
