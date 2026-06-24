"""
GS type generation from fragment patterns.

This module determines the Genome Structure (GS) type for a query genome by
matching its fragment pattern against a profile database. It handles exact
matches, orientation-independent matching, and assigns quality scores based
on how well the pattern matches known types and whether fragment arrangement
is valid.

Classes:
    TypeGenerator: Generates GS type assignments from fragment patterns
"""

from __future__ import annotations

import re
from typing import Union

from socru.GATProfile import GATProfile
from socru.Profiles import Profiles


class TypeGenerator:
    """
    Determine GS type from fragment pattern.

    This class matches a query genome's fragment pattern against a profile
    database to assign a GS (Genome Structure) type. It handles:
    - Exact matches to known types
    - Orientation-independent matching (ignoring fragment reversals)
    - Novel pattern detection
    - Quality scoring (GREEN/AMBER/RED) based on match quality and validation

    Quality levels:
    - GREEN: Known type, correct orientation, valid fragment arrangement
    - AMBER: Known order but orientation differs, or valid but novel pattern
    - RED: Invalid fragment arrangement or completely novel pattern

    Attributes:
        profile_db (Profiles): Database of known GS types
        gat_profile (GATProfile): Query genome's fragment pattern
        prefix (str): Prefix for GS type names (default 'GS')
        has_previously_seen (bool): True if exact match found in database
        verbose (bool): Enable verbose output
        is_frag_valid (bool): Whether fragment arrangement is biologically valid
        quality (str): Quality score (RED/AMBER/GREEN)
        gs_type (str): Assigned GS type string
    """
    def __init__(self, profile_db: Profiles, gat_profile: GATProfile, verbose: bool, is_frag_valid: bool, prefix: str = 'GS') -> None:
        """
        Initialize TypeGenerator with database and query profile.

        Args:
            profile_db (Profiles): Profile database to match against
            gat_profile (GATProfile): Query genome's fragment pattern
            verbose (bool): Enable verbose output
            is_frag_valid (bool): Whether fragment arrangement is valid
            prefix (str): Prefix for GS type names (default 'GS')
        """
        self.profile_db = profile_db
        self.gat_profile = gat_profile
        self.prefix = prefix
        self.has_previously_seen = False
        self.verbose = verbose
        self.is_frag_valid = is_frag_valid

        # Start with RED quality (invalid)
        self.quality = 'RED'
        if self.is_frag_valid:
            # Upgrade to AMBER if fragments are validly arranged
            self.quality = 'AMBER'

        # Calculate the GS type
        self.gs_type = self.calculate_type()

    def find_order_orientationless(self, orientationless_fragment: GATProfile) -> int:
        """
        Find matching order number ignoring fragment orientations.

        Searches database for a profile that matches the fragment order,
        regardless of whether fragments are forward or reversed.

        Args:
            orientationless_fragment (GATProfile): Query without orientation info

        Returns:
            int: Order number if match found, 0 otherwise
        """
        for db_profile in self.profile_db.gats:
            # Remove orientation from database profile for comparison
            orientationless_db_profile = GATProfile(self.verbose, fragments = db_profile.orientationless_fragments())
            if orientationless_db_profile.does_the_profile_match(orientationless_fragment):
                return db_profile.order()
        return 0

    def calculate_orientationless_order(self) -> int:
        """
        Calculate order number by matching without considering orientations.

        Tries to match the query's fragment order (ignoring reversals) to
        the database. Also tries the inverted order in case the genome is
        sequenced in the opposite direction.

        Returns:
            int: Order number if match found, 0 for completely novel pattern
        """
        # Try matching with orientations removed
        orientationless_fragment = GATProfile(self.verbose, fragments = self.gat_profile.orientationless_fragments())

        order = self.find_order_orientationless(orientationless_fragment)
        if order > 0:
            return order

        # Try inverted order
        inverted_orientationless_fragment = GATProfile(self.verbose, fragments = orientationless_fragment.inverted_orientationless_fragments())

        order = self.find_order_orientationless(inverted_orientationless_fragment)
        if order > 0:
            return order
        else:
            return 0

    def calculate_type(self) -> str:
        """
        Calculate the GS type for the query genome.

        This is the main method that:
        1. Determines the order number (fragment arrangement pattern)
        2. Checks for exact match in database (including orientations)
        3. Determines if pattern has been seen before
        4. Creates the final GS type string (e.g., "GS1.3")

        Returns:
            str: GS type string
        """
        # First, determine order number if not already set
        if self.gat_profile.gat_number == 0:
            self.gat_profile.gat_number = self.calculate_orientationless_order()

        # Check for exact match in database (including orientations)
        for db_profile in self.profile_db.gats:
            if db_profile.does_the_profile_match(self.gat_profile):
                self.gat_profile.gat_number = db_profile.gat_number
                continue

        # Check if exact pattern has been seen before
        self.has_previously_seen = self.previously_seen(self.gat_profile.gat_number, self.gat_profile.orientation_binary())

        # Create final GS type string
        return self.create_gs_type(self.gat_profile.gat_number, self.gat_profile.orientation_binary())

    def previously_seen(self, gat_number: Union[str, int], orientation_binary: int) -> bool:
        """
        Check if exact pattern (order + orientation) is in database.

        GS types are formatted as "order.orientation" (e.g., "1.3").
        This checks if the orientation number matches the database entry.
        If it matches and fragments are valid, upgrades quality to GREEN.

        Args:
            gat_number: GS type number from database
            orientation_binary (int): Binary encoding of fragment orientations

        Returns:
            bool: True if exact match found, False otherwise
        """
        m = re.search(r'([\d]+)\.([\d]+)$', str(gat_number) )

        if m:
            if str(m.group(1)) == str(0):
                # Order 0 means completely novel
                return False

            if str(m.group(2)) == str(orientation_binary):
                # Orientation matches - exact match found!
                # Upgrade to GREEN quality if fragments are valid
                if self.is_frag_valid:
                    self.quality = 'GREEN'

                return True
            else:
                # Order matches but orientation differs
                return False
        else:
            return False

    def create_gs_type(self, gat_number: Union[str, int], orientation_binary: int) -> str:
        """
        Create GS type string from order and orientation numbers.

        Format: "GS{order}.{orientation}"
        Examples: "GS1.3", "GS2.5", "GS0.0" (unknown)

        Args:
            gat_number: Order number or full GS type from database
            orientation_binary (int): Binary encoding of orientations

        Returns:
            str: Formatted GS type string
        """
        m = re.search(r'([\d]+)\.([\d]+)$', str(gat_number) )

        if m:
            if str(m.group(2)) == str(orientation_binary):
                # Orientation matches database
                return self.prefix + str(gat_number)
            else:
                # Same order, different orientation
                return self.prefix + str(m.group(1)) + "." + str(orientation_binary)

        else:
            # New format needed
            return self.prefix + str(gat_number) + "." + str(orientation_binary)
