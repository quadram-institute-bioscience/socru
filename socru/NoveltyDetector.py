"""
Novelty detection for genome structure profiles.

When a genome receives a novel GS type (not found in the profile database),
this module assesses whether the novel arrangement is likely a genuine
biological rearrangement or a possible assembly artifact.  The assessment
combines rearrangement distance to the nearest known type, assembly
confidence scores, and per-fragment BLAST identity values.

Classes:
    NoveltyAssessment: Dataclass holding the assessment result

Functions:
    assess_novelty: Perform the full novelty assessment
"""

from dataclasses import dataclass, field
from typing import List, Optional

from socru.RearrangementDistance import rearrangement_distance


@dataclass
class NoveltyAssessment:
    """Assessment of a novel genome structure.

    Attributes:
        nearest_known_type: Closest GS type identifier in the database.
        edit_distance: Minimum fragment swaps/inversions to reach the
            nearest known type.
        is_likely_real: True when the assessment is ``'likely_real'``.
        assessment: One of ``'likely_real'``, ``'possibly_artifactual'``,
            or ``'uncertain'``.
        reasoning: Human-readable explanation of the assessment.
        fragment_differences: List of individual difference descriptions.
    """

    nearest_known_type: str
    edit_distance: int
    is_likely_real: bool
    assessment: str
    reasoning: str
    fragment_differences: List[str] = field(default_factory=list)


def _find_nearest_profile(
    query_fragments: List[str],
    known_profiles: List[List[str]],
) -> tuple:
    """Find the known profile closest to the query.

    Args:
        query_fragments: Query fragment list.
        known_profiles: List of known fragment lists from the profile DB.

    Returns:
        Tuple of ``(best_index, best_distance_dict)`` for the closest
        known profile, or ``(None, None)`` when *known_profiles* is empty.
    """
    best_index: Optional[int] = None
    best_dist: Optional[dict] = None
    best_total: int = 999_999

    for idx, known in enumerate(known_profiles):
        dist = rearrangement_distance(query_fragments, known)
        total = dist["distance"]
        if total < best_total:
            best_total = total
            best_dist = dist
            best_index = idx

    return best_index, best_dist


def assess_novelty(
    query_fragments: List[str],
    known_profiles: List[List[str]],
    confidence_score: float,
    blast_identities: List[float],
) -> NoveltyAssessment:
    """Assess whether a novel profile is likely real or an artifact.

    The function locates the nearest known profile, computes the
    rearrangement distance, and combines that with assembly confidence
    and BLAST identity metrics to classify the novelty.

    Args:
        query_fragments: Fragment strings for the query genome, e.g.
            ``["1", "2'", "3", "4"]``.
        known_profiles: List of fragment-string lists from the profile
            database.
        confidence_score: Assembly confidence value in the range 0--100.
        blast_identities: Per-fragment BLAST identity percentages (0--100).

    Returns:
        A :class:`NoveltyAssessment` with the classification result.
    """
    # Edge case: no known profiles to compare against
    if not known_profiles:
        return NoveltyAssessment(
            nearest_known_type="none",
            edit_distance=len(query_fragments),
            is_likely_real=False,
            assessment="uncertain",
            reasoning="No known profiles available for comparison.",
            fragment_differences=[],
        )

    best_index, best_dist = _find_nearest_profile(query_fragments, known_profiles)
    edit_distance: int = best_dist["distance"]  # type: ignore[index]
    details: List[str] = best_dist["details"]  # type: ignore[index]

    nearest_label = f"profile_{best_index}"

    # --- Classification logic ---
    min_blast = min(blast_identities) if blast_identities else 0.0

    if confidence_score > 80 and min_blast > 95:
        assessment = "likely_real"
        is_likely_real = True
        reasoning = (
            f"High confidence ({confidence_score:.1f}) and all BLAST identities "
            f"above 95% (min {min_blast:.1f}%). Edit distance to nearest known "
            f"type is {edit_distance}."
        )
    elif confidence_score < 50 or min_blast < 90:
        assessment = "possibly_artifactual"
        is_likely_real = False
        reasons = []
        if confidence_score < 50:
            reasons.append(f"low confidence ({confidence_score:.1f})")
        if min_blast < 90:
            reasons.append(f"low BLAST identity (min {min_blast:.1f}%)")
        reasoning = (
            f"Possibly artifactual due to {' and '.join(reasons)}. "
            f"Edit distance to nearest known type is {edit_distance}."
        )
    else:
        assessment = "uncertain"
        is_likely_real = False
        reasoning = (
            f"Confidence ({confidence_score:.1f}) and BLAST identities "
            f"(min {min_blast:.1f}%) are in the ambiguous range. "
            f"Edit distance to nearest known type is {edit_distance}. "
            f"Manual review recommended."
        )

    return NoveltyAssessment(
        nearest_known_type=nearest_label,
        edit_distance=edit_distance,
        is_likely_real=is_likely_real,
        assessment=assessment,
        reasoning=reasoning,
        fragment_differences=details,
    )
