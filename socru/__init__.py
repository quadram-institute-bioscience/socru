"""
socru - Typing of genome-level order and orientation around ribosomal operons in bacteria.

Usage as a library:
    from socru import Socru, SocruConfig, AnalysisResult

    config = SocruConfig(
        input_files=['genome.fasta'],
        species='Escherichia_coli',
    )
    with Socru(config) as s:
        s.run()
"""
from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("socru")
except PackageNotFoundError:
    __version__ = "unknown"

# Core analysis
# Data models
from socru.AnalysisResult import AnalysisResult, FragmentResult, OperonResult, QCFlag
from socru.BatchStats import BatchStats
from socru.BlastResult import BlastResult

# Analysis modules
from socru.ConfidenceScore import calculate_confidence
from socru.DatabaseManager import DatabaseManager
from socru.Fragment import Fragment
from socru.GATProfile import GATProfile

# Reporting
from socru.HtmlReport import HtmlReport
from socru.NoveltyDetector import NoveltyAssessment, assess_novelty
from socru.Operon import Operon
from socru.QCFlags import generate_qc_flags
from socru.RearrangementDistance import pairwise_distance_matrix, rearrangement_distance
from socru.Socru import Socru
from socru.SocruConfig import SocruConfig, SocruCreateConfig
from socru.SocruCreate import SocruCreate
from socru.SvgConfidenceHeatmap import generate_confidence_heatmap_svg
from socru.SvgCoveragePileup import generate_coverage_pileup_svg
from socru.SvgFragmentQuality import generate_fragment_quality_svg

# Visualization
from socru.SvgGenomePlot import generate_genome_svg, save_genome_svg
from socru.SvgSynteny import generate_synteny_svg
from socru.SvgTypeDistribution import generate_type_distribution_svg

__all__ = [
    # Core
    'Socru', 'SocruCreate', 'SocruConfig', 'SocruCreateConfig',
    # Models
    'AnalysisResult', 'FragmentResult', 'OperonResult', 'QCFlag',
    'Fragment', 'BlastResult', 'Operon', 'GATProfile',
    'NoveltyAssessment',
    # Analysis
    'calculate_confidence', 'generate_qc_flags',
    'assess_novelty', 'rearrangement_distance', 'pairwise_distance_matrix',
    'BatchStats', 'DatabaseManager',
    # Visualization
    'generate_genome_svg', 'save_genome_svg',
    'generate_synteny_svg', 'generate_fragment_quality_svg',
    'generate_type_distribution_svg', 'generate_confidence_heatmap_svg',
    'generate_coverage_pileup_svg',
    # Reporting
    'HtmlReport',
    # Meta
    '__version__',
]
