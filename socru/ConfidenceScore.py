"""
Confidence scoring for Socru genome structure analysis results.

This module computes a composite confidence score (0-100) that summarizes
how reliable a genome structure typing result is. The score combines BLAST
alignment quality, profile match level, and validation status.

Functions:
    calculate_confidence: Compute a 0-100 confidence score for an analysis.
"""

from typing import List

# Weights for each component of the confidence score
WEIGHT_IDENTITY = 0.30
WEIGHT_COVERAGE = 0.25
WEIGHT_PROFILE_MATCH = 0.25
WEIGHT_VALIDATION = 0.20

# Penalty subtracted for each unknown ("?") fragment
UNKNOWN_FRAGMENT_PENALTY = 10


def _mean_blast_identity(fragment_results: List) -> float:
    """Compute mean BLAST identity across fragments with hits.

    Args:
        fragment_results: List of FragmentResult objects.

    Returns:
        Mean identity as a 0-100 value. Returns 0.0 if no fragments have hits.
    """
    identities = [
        f.blast_identity for f in fragment_results
        if f.blast_identity is not None
    ]
    if not identities:
        return 0.0
    return sum(identities) / len(identities)


def _mean_coverage_ratio(fragment_results: List) -> float:
    """Compute mean alignment-length / fragment-length ratio.

    Args:
        fragment_results: List of FragmentResult objects.

    Returns:
        Mean coverage ratio scaled to 0-100. Returns 0.0 if no data available.
    """
    ratios = []
    for f in fragment_results:
        if f.blast_alignment_length is not None and f.length > 0:
            ratio = f.blast_alignment_length / f.length
            # Cap at 1.0 (alignment can exceed fragment with extensions)
            ratios.append(min(ratio, 1.0))
    if not ratios:
        return 0.0
    return (sum(ratios) / len(ratios)) * 100.0


def _profile_match_score(quality: str, is_novel: bool) -> float:
    """Score based on how well the profile matches known types.

    Args:
        quality: Quality label (GREEN, AMBER, RED).
        is_novel: Whether the profile is novel (not in database).

    Returns:
        Score from 0-100 reflecting profile match level.
    """
    if quality == 'GREEN':
        return 100.0
    elif quality == 'AMBER' and not is_novel:
        # Known order, different orientation
        return 60.0
    elif quality == 'AMBER' and is_novel:
        # Valid arrangement but novel pattern
        return 60.0
    elif quality == 'RED' and not is_novel:
        # Invalid arrangement but recognized order
        return 20.0
    else:
        # RED and novel or unknown fragments
        return 0.0


def _count_unknown_fragments(fragment_results: List) -> int:
    """Count fragments that could not be matched (marked as '?').

    Args:
        fragment_results: List of FragmentResult objects.

    Returns:
        Number of unknown fragments.
    """
    return sum(
        1 for f in fragment_results
        if f.blast_identity is None and str(f.number) == '?'
    )


def calculate_confidence(
    fragment_results: List,
    quality: str,
    validation_passed: bool,
    is_novel: bool = False,
) -> float:
    """Compute a composite confidence score for a genome structure result.

    The score is a weighted combination of:
    - Mean fragment BLAST identity (30%)
    - Mean fragment coverage ratio (25%)
    - Profile match level (25%)
    - Validation status (20%)

    A penalty of 10 points per unknown ("?") fragment is subtracted.

    Args:
        fragment_results: List of FragmentResult objects from the analysis.
        quality: Quality label (GREEN, AMBER, RED).
        validation_passed: Whether fragment arrangement passed validation.
        is_novel: Whether the profile is novel.

    Returns:
        Confidence score clamped to [0.0, 100.0].
    """
    identity_score = _mean_blast_identity(fragment_results)
    coverage_score = _mean_coverage_ratio(fragment_results)
    profile_score = _profile_match_score(quality, is_novel)
    validation_score = 100.0 if validation_passed else 0.0

    raw_score = (
        WEIGHT_IDENTITY * identity_score
        + WEIGHT_COVERAGE * coverage_score
        + WEIGHT_PROFILE_MATCH * profile_score
        + WEIGHT_VALIDATION * validation_score
    )

    # Apply penalty for unknown fragments
    num_unknown = _count_unknown_fragments(fragment_results)
    raw_score -= num_unknown * UNKNOWN_FRAGMENT_PENALTY

    return max(0.0, min(100.0, round(raw_score, 1)))
