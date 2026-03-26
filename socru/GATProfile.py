"""
Genome Arrangement Type (GAT) profile representation and manipulation.

This module defines the GATProfile class which represents a genome structure
profile as an ordered list of inter-operon fragments with their orientations.
It handles profile reorientation relative to dnaA (origin), profile matching,
and various transformations needed for genome structure typing.

Classes:
    GATProfile: Represents and manipulates genome arrangement profiles
"""

from __future__ import annotations

import re
from typing import Any, Union


class GATProfile:
    """
    Represents a genome arrangement (GAT) profile.

    A GAT profile is an ordered list of fragment numbers with orientations
    (indicated by prime ' for reverse). Profiles can be compared, inverted,
    and reoriented. The profile must start from the dnaA origin fragment
    for standardization.

    Attributes:
        gat_number: GS type identifier (e.g., "1.1")
        fragments (list): Ordered list of fragment identifiers (e.g., ["1", "2'", "3"])
        verbose (bool): Enable verbose output
        orientation_number (int): Binary encoding of fragment orientations
        dnaA_fragment_number (int): Fragment containing dnaA origin (default 3)
        dif_fragment_number (int): Fragment containing dif terminator (default 1)
    """

    def __init__(self, verbose: bool, gat_number: Union[str, int] = 0, fragments: list[str] | None = None, orientation_number: int = 0, dnaA_fragment_number: int = 3, dif_fragment_number: int = 1) -> None:
        """
        Initialize a GATProfile.

        Args:
            verbose (bool): Enable verbose output
            gat_number: GS type number (default 0)
            fragments (list): List of fragment identifiers (default empty)
            orientation_number (int): Binary orientation encoding (default 0)
            dnaA_fragment_number (int): Origin fragment ID (default 3)
            dif_fragment_number (int): Terminus fragment ID (default 1)
        """
        self.gat_number = gat_number
        self.fragments = fragments if fragments is not None else []
        self.verbose = verbose
        self.orientation_number = orientation_number
        self.dnaA_fragment_number = dnaA_fragment_number
        self.dif_fragment_number = dif_fragment_number

    def order(self) -> int:
        """
        Extract the major order number from GS type (e.g., 1 from "1.2").

        Returns:
            int: Major order number, or 0 if not found
        """
        m = re.search(r'([\d]+)\.[\d]+$', str(self.gat_number) )

        if m:
            return int(m.group(1))
        else:
            return 0

    def fragment_str(self) -> str:
        """
        Format fragments as tab-delimited string.

        Returns:
            str: Tab-separated fragment identifiers
        """
        return "\t".join(self.fragments)

    def invert_fragments(self) -> list[str]:
        """
        Invert the fragment profile (reverse order and orientations).

        Reverses the fragment list and flips each fragment's orientation
        (removes or adds prime notation). Then reorients to start with
        fragment 1 or 1' for standardization.

        Returns:
            list: Inverted and reoriented fragment list
        """
        inverted = []
        for f in self.fragments:
            # Check if fragment is reversed (has prime)
            m = re.match(r"([\d]+)'", f)
            if m:
                # Remove prime (now forward)
                inverted.append(m.group(1))
            else:
                # Add prime (now reverse)
                inverted.append(f + '\'')

        # Reverse the order
        inverted.reverse()
        # Reorient to start with fragment 1
        return self.reorientate_list_to_start_with_one(inverted)

    def reorientate_list_to_start_with_one(self, raw_fragments: list[str]) -> list[str]:
        """
        Rotate fragment list to start with fragment 1 or 1'.

        This provides a consistent starting point for fragment lists,
        regardless of the original circular order.

        Args:
            raw_fragments (list): Fragment list in any order

        Returns:
            list: Rotated list starting with fragment 1 or 1'
        """
        # find the index where fragment 1 or 1' is
        chop_index = 0
        for i in range(0,len(raw_fragments)):
            if raw_fragments[i] == '1' or raw_fragments[i] == "1'":
                chop_index = i
                continue

        # Rotate if fragment 1 is not at the start
        reorientated = raw_fragments
        if chop_index > 0:
            reorientated = raw_fragments[chop_index:len(raw_fragments)] + raw_fragments[0:chop_index]
        return reorientated

    def inverted_fragment_str(self) -> str:
        """
        Get inverted profile as tab-delimited string.

        Returns:
            str: Tab-separated inverted fragment identifiers
        """
        return "\t".join(self.invert_fragments())

    def is_profile_in_correct_orientation(self) -> bool:
        """
        Check if profile starts in correct orientation relative to dnaA.

        Correct orientation means dnaA fragment is forward (not reversed).
        This ensures consistent profile representation.

        Returns:
            bool: True if correctly oriented, False if needs inversion
        """
        for f in self.fragments:
            if f == str(self.dnaA_fragment_number):
                # DnaA is forward - correct orientation
                return True
            elif f == str(self.dnaA_fragment_number) + "'":
                # DnaA is reversed - needs inversion
                return False
        # dnaA hasnt been found - assume correct (should perhaps raise exception)
        return True

    def reorder_fragment_objects_based_on_fragment_name_array(self, fragment_objects: list[Any]) -> list[Any]:
        """
        Reorder Fragment objects to match this profile's fragment order.

        Takes a list of Fragment objects and reorders them to match the
        fragment order in this profile. Also sets the reversed_frag flag
        based on prime notation.

        Args:
            fragment_objects (list): List of Fragment objects

        Returns:
            list: Reordered Fragment objects matching profile order
        """
        reordered_fragment_objects = []

        # For each fragment name in profile, find corresponding object
        for frag_name in self.fragments:
            frag_name_orientationless = frag_name
            frag_name_reversed = False

            # Check for prime notation (reversed)
            m = re.match(r"([\d]+)'", frag_name)
            if m:
                frag_name_orientationless = m.group(1)
                frag_name_reversed = True

            # Find matching fragment object
            for frag_obj in fragment_objects:
                if str(frag_obj.number) == str(frag_name_orientationless):
                    # Set orientation flag
                    if frag_name_reversed:
                        frag_obj.reversed_frag = True
                    else:
                        frag_obj.reversed_frag = False
                    reordered_fragment_objects.append(frag_obj)

        return reordered_fragment_objects

    def orientate_for_dnaA(self) -> None:
        """
        Ensure profile is oriented correctly relative to dnaA origin.

        If profile is inverted (dnaA reversed), inverts the entire profile
        so dnaA is forward. This standardizes profile representation.
        """
        if not self.is_profile_in_correct_orientation():
            self.fragments  = self.invert_fragments()

    def does_the_profile_match(self, gat_query: GATProfile) -> bool:
        """
        Check if another profile matches this one exactly.

        Args:
            gat_query (GATProfile): Profile to compare

        Returns:
            bool: True if profiles match, False otherwise
        """
        if gat_query.fragment_str() == self.fragment_str():
            return True
        else:
            return False

    def __str__(self) -> str:
        """
        String representation as tab-delimited fragments.

        Returns:
            str: Fragment string representation
        """
        return self.fragment_str()

    def remove_orientation(self, fragments: list[str]) -> list[str]:
        """
        Remove orientation indicators (primes) from fragment list.

        Converts fragments like ["1", "2'", "3"] to ["1", "2", "3"],
        ignoring strand orientation.

        Args:
            fragments (list): Fragment list with orientations

        Returns:
            list: Fragment list without orientation indicators
        """
        orientationless = []
        for f in fragments:
            m = re.match(r"([\d]+)'", f)
            if m:
                # Remove prime
                orientationless.append(m.group(1))
            else:
                orientationless.append(f)
        return orientationless

    def orientationless_fragments(self) -> list[str]:
        """
        Get current fragments without orientation information.

        Returns:
            list: Fragment numbers without orientation
        """
        return self.remove_orientation(self.fragments)

    def inverted_orientationless_fragments(self) -> list[str]:
        """
        Get inverted fragments without orientation information.

        Returns:
            list: Inverted fragment numbers without orientation
        """
        inverted = self.invert_fragments()
        return self.remove_orientation(inverted)

    def orientation_binary(self) -> int:
        """
        Encode fragment orientations as binary number.

        Each fragment position gets a bit: 1 if reversed (prime), 0 if forward.
        Fragment N uses bit position N-1. This provides a compact numerical
        representation of the orientation pattern.

        Example:
            Fragments ["1", "2'", "3"] -> binary 010 -> decimal 2

        Bit positions:
            1' - 00000001 (bit 0)
            2' - 00000010 (bit 1)
            3' - 00000100 (bit 2)
            etc.

        Returns:
            int: Binary encoding of orientations
        """
        total = 0
        for f in self.fragments:
            m = re.match(r"([\d]+)'", f)
            if m:
                # Fragment is reversed - set its bit
                bits_to_shift = int(m.group(1)) -1
                total += 1 << bits_to_shift
        return total

    def orientation_array(self) -> list[bool]:
        """
        Get fragment orientations as boolean array.

        Returns:
            list: Boolean list where True=forward, False=reverse
        """
        orientation_bools = []
        for f in self.fragments:
            m = re.match(r"([\d]+)'", f)
            if m:
                # Fragment is reversed
                orientation_bools.append(False)
            else:
                # Fragment is forward
                orientation_bools.append(True)
        return orientation_bools
