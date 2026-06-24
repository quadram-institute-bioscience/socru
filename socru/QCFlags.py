"""
Quality control flag generation for Socru analysis results.

This module inspects an AnalysisResult and generates QCFlag objects that
highlight potential issues with the genome structure typing. Flags are
used by downstream reporting and visualization to draw attention to
results that may need manual review.

Functions:
    generate_qc_flags: Generate QC flags from an analysis result.
"""

from typing import List

from socru.AnalysisResult import QCFlag

# Thresholds
IDENTITY_THRESHOLD = 95.0
COVERAGE_THRESHOLD = 0.80
CONFIDENCE_THRESHOLD = 50.0
SMALL_CHROMOSOME_BP = 500000


def generate_qc_flags(
    analysis_result,
    expected_fragment_count: int,
) -> List[QCFlag]:
    """Generate quality control flags for a completed analysis result.

    Inspects the analysis result and returns a list of QCFlag objects for
    any conditions that warrant attention.

    Flags generated:
        LOW_IDENTITY: Any fragment BLAST identity below 95%.
        SHORT_ALIGNMENT: Any fragment alignment covers <80% of its length.
        MISSING_FRAGMENT: One or more fragments returned '?'.
        UNEXPECTED_OPERON_COUNT: Operon count differs from expected.
        NOVEL_PATTERN: Fragment order not found in database.
        NOVEL_ORIENTATION: Fragment order known but orientation is new.
        INVALID_ARRANGEMENT: Operon orientations violate replichore rules.
        LOW_CONFIDENCE: Confidence score below 50.
        SMALL_CHROMOSOME: Chromosome shorter than 500,000 bp.

    Args:
        analysis_result: An AnalysisResult object.
        expected_fragment_count: Expected number of fragments from the
            profile database.

    Returns:
        List of QCFlag objects.
    """
    flags: List[QCFlag] = []

    # LOW_IDENTITY: any fragment with BLAST identity below threshold
    for frag in analysis_result.fragments:
        if frag.blast_identity is not None and frag.blast_identity < IDENTITY_THRESHOLD:
            flags.append(QCFlag(
                code='LOW_IDENTITY',
                severity='warning',
                message=f'Fragment {frag.number} has BLAST identity {frag.blast_identity:.1f}% (threshold {IDENTITY_THRESHOLD}%)',
                details=f'subject={frag.blast_subject}',
            ))

    # SHORT_ALIGNMENT: alignment covers less than 80% of fragment length
    for frag in analysis_result.fragments:
        if (frag.blast_alignment_length is not None
                and frag.length > 0
                and (frag.blast_alignment_length / frag.length) < COVERAGE_THRESHOLD):
            coverage_pct = (frag.blast_alignment_length / frag.length) * 100
            flags.append(QCFlag(
                code='SHORT_ALIGNMENT',
                severity='warning',
                message=f'Fragment {frag.number} alignment covers {coverage_pct:.1f}% of fragment length (threshold {COVERAGE_THRESHOLD * 100:.0f}%)',
            ))

    # MISSING_FRAGMENT: fragments that could not be matched
    missing = [f for f in analysis_result.fragments if str(f.number) == '?']
    if missing:
        flags.append(QCFlag(
            code='MISSING_FRAGMENT',
            severity='error',
            message=f'{len(missing)} fragment(s) could not be matched to the database',
        ))

    # UNEXPECTED_OPERON_COUNT: number of operons differs from expected
    # Expected operons = expected_fragment_count (for circular) or expected_fragment_count - 1
    if expected_fragment_count > 0:
        expected_operons = expected_fragment_count
        if analysis_result.num_operons != expected_operons:
            flags.append(QCFlag(
                code='UNEXPECTED_OPERON_COUNT',
                severity='warning',
                message=f'Found {analysis_result.num_operons} operons, expected {expected_operons}',
            ))

    # NOVEL_PATTERN: fragment order not in database
    if analysis_result.is_novel and analysis_result.quality in ('RED', 'AMBER'):
        # Distinguish between novel order and novel orientation
        # If quality is AMBER, the order is known but orientation differs
        if analysis_result.quality == 'AMBER':
            flags.append(QCFlag(
                code='NOVEL_ORIENTATION',
                severity='warning',
                message='Fragment order is known but this orientation combination is new',
            ))
        else:
            flags.append(QCFlag(
                code='NOVEL_PATTERN',
                severity='warning',
                message='Fragment order not found in the profile database',
            ))

    # INVALID_ARRANGEMENT: validation failed
    if not analysis_result.validation_passed:
        flags.append(QCFlag(
            code='INVALID_ARRANGEMENT',
            severity='error',
            message='Operon orientations violate replichore rules',
        ))

    # LOW_CONFIDENCE
    if analysis_result.confidence_score < CONFIDENCE_THRESHOLD:
        flags.append(QCFlag(
            code='LOW_CONFIDENCE',
            severity='warning',
            message=f'Confidence score {analysis_result.confidence_score:.1f} is below threshold {CONFIDENCE_THRESHOLD}',
        ))

    # SMALL_CHROMOSOME
    if analysis_result.genome_length < SMALL_CHROMOSOME_BP:
        flags.append(QCFlag(
            code='SMALL_CHROMOSOME',
            severity='warning',
            message=f'Chromosome length {analysis_result.genome_length:,} bp is below {SMALL_CHROMOSOME_BP:,} bp',
        ))

    return flags
