"""
Rearrangement distance computation between genome structure profiles.

This module computes the minimum number of operations (inversions and
translocations) needed to transform one GS type fragment order into another.
Inspired by the GRIMM reversal distance but simplified for operon-level
fragments used in socru.

Functions:
    rearrangement_distance: Compute edit distance between two profiles
    pairwise_distance_matrix: Compute all pairwise distances between named profiles
"""

import re
from typing import Dict, List, Tuple


def _strip_orientation(fragment: str) -> str:
    """Remove orientation marker (prime) from a fragment identifier.

    Args:
        fragment: Fragment string like "2" or "2'".

    Returns:
        Base fragment number as string, e.g. "2".
    """
    m = re.match(r"([\d]+)'", fragment)
    if m:
        return m.group(1)
    return fragment


def _is_inverted(fragment: str) -> bool:
    """Check whether a fragment has the inverted (prime) orientation.

    Args:
        fragment: Fragment string like "2" or "2'".

    Returns:
        True if the fragment carries a prime marker.
    """
    return fragment.endswith("'")


def rearrangement_distance(
    profile_a: List[str], profile_b: List[str]
) -> Dict[str, object]:
    """Compute minimum edit distance between two fragment profiles.

    The distance is decomposed into:
    - **inversions**: fragments at the same position but with different
      orientation (one forward, one reversed).
    - **translocations**: fragments that appear at different positions in
      the two profiles.
    - **missing/extra**: fragments present in one profile but absent from
      the other (represented by "?" in either profile).

    Args:
        profile_a: List of fragment strings, e.g. ``["1", "2'", "3"]``.
        profile_b: List of fragment strings, e.g. ``["1", "3'", "2"]``.

    Returns:
        Dictionary with keys:
            - ``distance`` (int): total operations
            - ``inversions`` (int): orientation flip count
            - ``translocations`` (int): position change count
            - ``missing`` (int): fragments present in only one profile
            - ``details`` (list of str): human-readable difference descriptions
    """
    inversions = 0
    translocations = 0
    missing_count = 0
    details: List[str] = []

    # Build position maps keyed by base fragment number
    base_a = [_strip_orientation(f) for f in profile_a]
    base_b = [_strip_orientation(f) for f in profile_b]

    pos_a = {base: idx for idx, base in enumerate(base_a)}
    pos_b = {base: idx for idx, base in enumerate(base_b)}

    all_bases = sorted(set(base_a) | set(base_b), key=lambda x: int(x) if x.isdigit() else 0)

    for base in all_bases:
        # Handle missing / unknown fragments
        in_a = base in pos_a and profile_a[pos_a[base]] != "?"
        in_b = base in pos_b and profile_b[pos_b[base]] != "?"

        # Also treat "?" entries in the profile list itself
        a_is_missing = not in_a or base == "?"
        b_is_missing = not in_b or base == "?"

        if a_is_missing and not b_is_missing:
            missing_count += 1
            details.append(f"Fragment {base}: missing from profile A")
            continue
        if b_is_missing and not a_is_missing:
            missing_count += 1
            details.append(f"Fragment {base}: missing from profile B")
            continue
        if a_is_missing and b_is_missing:
            continue

        idx_a = pos_a[base]
        idx_b = pos_b[base]
        frag_a = profile_a[idx_a]
        frag_b = profile_b[idx_b]

        # Check positional difference (translocation)
        if idx_a != idx_b:
            translocations += 1
            details.append(
                f"Fragment {base}: translocated from position {idx_a} to {idx_b}"
            )

        # Check orientation difference (inversion)
        if _is_inverted(frag_a) != _is_inverted(frag_b):
            inversions += 1
            details.append(
                f"Fragment {base}: orientation flipped "
                f"({'reverse' if _is_inverted(frag_a) else 'forward'} -> "
                f"{'reverse' if _is_inverted(frag_b) else 'forward'})"
            )

    # Handle "?" entries that appear literally in the fragment lists
    for profile_label, profile in [("A", profile_a), ("B", profile_b)]:
        for frag in profile:
            if frag == "?":
                missing_count += 1
                details.append(f"Unknown fragment '?' in profile {profile_label}")

    distance = inversions + translocations + missing_count
    return {
        "distance": distance,
        "inversions": inversions,
        "translocations": translocations,
        "missing": missing_count,
        "details": details,
    }


def pairwise_distance_matrix(
    profiles_dict: Dict[str, List[str]],
) -> Dict[Tuple[str, str], Dict[str, object]]:
    """Compute all pairwise distances between named profiles.

    Args:
        profiles_dict: Mapping of GS type name to fragment list, e.g.
            ``{"GS1.0": ["1", "2", "3"], "GS2.0": ["1", "3", "2"]}``.

    Returns:
        Dictionary keyed by ``(type_a, type_b)`` tuples with distance
        dictionaries as values.  Both orderings are included (the matrix
        is symmetric).
    """
    names = sorted(profiles_dict.keys())
    matrix: Dict[Tuple[str, str], Dict[str, object]] = {}

    for i, name_a in enumerate(names):
        for j, name_b in enumerate(names):
            if i <= j:
                dist = rearrangement_distance(
                    profiles_dict[name_a], profiles_dict[name_b]
                )
                matrix[(name_a, name_b)] = dist
                matrix[(name_b, name_a)] = dist

    return matrix
