"""
BLAST result filtering and coverage analysis.

This module filters BLAST results based on quality thresholds (bit score, alignment
length) and performs coverage analysis to identify high-confidence regions in
fragment alignments. It supports selecting the best hit and analyzing fragment
coverage depth across genomic regions.

Classes:
    FilterBlast: Filters BLAST results and analyzes alignment coverage
"""

from __future__ import annotations

import csv
import logging
from typing import Optional, Union

import numpy
import numpy.typing as npt

from socru.BlastResult import BlastResult

logger = logging.getLogger(__name__)


class FilterBlast:
    """
    Filter and analyze BLAST search results.

    This class reads BLAST tabular output, applies quality filters, and can
    perform coverage analysis on fragment hits. It identifies the top hit
    for each query and can calculate coverage depth across subject sequences.

    Attributes:
        results_file (str): Path to BLAST results file (tabular format)
        min_bit_score (int): Minimum bit score threshold
        min_alignment_length (int): Minimum alignment length threshold
        verbose (bool): Enable verbose output
        results (list): List of BlastResult objects from input file
    """
    def __init__(self, results_file: str, min_bit_score: int, min_alignment_length: int, verbose: bool) -> None:
        """
        Initialize FilterBlast with results file and thresholds.

        Args:
            results_file (str): Path to BLAST tabular output
            min_bit_score (int): Minimum bit score for filtering
            min_alignment_length (int): Minimum alignment length for filtering
            verbose (bool): Enable verbose output
        """
        self.results_file = results_file
        self.min_bit_score = min_bit_score
        self.min_alignment_length =  min_alignment_length
        self.verbose = verbose
        self.results = self.readin_results()

    def readin_results(self) -> list[BlastResult]:
        """
        Parse BLAST tabular output file into BlastResult objects.

        Reads tab-delimited BLAST output (format 6) and creates BlastResult
        objects for each hit.

        Returns:
            list: List of BlastResult objects
        """
        results = []
        with open(self.results_file , newline='') as blastfile:
            blastreader = csv.reader(blastfile, delimiter='\t')
            for line_num, row in enumerate(blastreader, start=1):
                if not row or (len(row) == 1 and not row[0].strip()):
                    # Skip blank lines
                    continue
                if len(row) < 12:
                    logger.warning(
                        "Skipping malformed BLAST line %d in %s: expected 12 columns, got %d",
                        line_num, self.results_file, len(row),
                    )
                    continue
                try:
                    # Create BlastResult from 12 standard BLAST columns
                    results.append(BlastResult(row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7], row[8], row[9], row[10], row[11]))
                except (ValueError, IndexError) as e:
                    logger.warning(
                        "Skipping unparseable BLAST line %d in %s: %s",
                        line_num, self.results_file, e,
                    )
        return results

    def filter_results(self) -> list[BlastResult]:
        """
        Apply quality filters to BLAST results.

        Filters results based on minimum bit score and alignment length
        thresholds to remove low-quality hits.

        Returns:
            list: Filtered list of BlastResult objects
        """
        filtered = [r for r in self.results if r.bit_score > self.min_bit_score and r.alignment_length > self.min_alignment_length ]
        return filtered

    def return_top_result(self) -> Optional[BlastResult]:
        """
        Get the best BLAST hit after filtering.

        Assumes results are sorted by bit score (from BLAST command piped
        through sort). Returns the first result after filtering.

        Returns:
            BlastResult: Top scoring result, or None if no results pass filters
        """
        results = self.filter_results()
        if len(results) > 0:
            return results[0]
        else:
            return None

    def max_coord_from_blast_results(self, blast_result_objs: list[BlastResult]) -> int:
        """
        Find the maximum coordinate across all BLAST hits.

        Used to determine the size needed for coverage arrays.

        Args:
            blast_result_objs (list): List of BlastResult objects

        Returns:
            int: Maximum coordinate value
        """
        start = [int(r.subject_start) for r in blast_result_objs]
        end = [int(r.subject_end) for r in blast_result_objs]
        return max(start + end)

    def pileup_fragment(self, fragment_number: Union[str, int]) -> npt.NDArray[numpy.floating]:
        """
        Create coverage pileup array for a specific fragment.

        Counts how many times each base position in the fragment is covered
        by BLAST alignments. This is similar to read depth in sequencing.

        Args:
            fragment_number: Fragment identifier to analyze

        Returns:
            numpy.array: Coverage depth at each position
        """
        # Get all hits to this fragment
        fragment_results = [r for r in self.results if str(r.subject) == str(fragment_number)]
        max_coord = self.max_coord_from_blast_results(fragment_results)

        # Initialize coverage array
        fragment_pileup = numpy.zeros(max_coord)

        # Add coverage for each alignment
        for r in fragment_results:
            start = int(r.subject_start)
            end = int(r.subject_end)

            # Handle reverse strand hits (start > end)
            if start > end:
                start = int(r.subject_end)
                end = int(r.subject_start)

            # Increment coverage for aligned region
            for i in range(start, end):
                fragment_pileup[i-1] += 1
        return fragment_pileup

    def num_bases_with_at_least_this_coverage(self, target_coverage: int, pileup: npt.NDArray[numpy.floating]) -> int:
        """
        Count bases meeting minimum coverage threshold.

        Args:
            target_coverage (int): Minimum coverage depth required
            pileup (numpy.array): Coverage depth array

        Returns:
            int: Number of bases with at least target_coverage depth
        """
        if target_coverage == 0:
            return 0
        base_count = sum([1 for p in pileup if p >= target_coverage])
        return base_count

    def calc_coverage_threshold(self, max_coverage: float, pileup: npt.NDArray[numpy.floating], target_bases: int) -> int:
        """
        Calculate coverage threshold to achieve target number of bases.

        Finds the highest coverage threshold where at least target_bases
        positions meet or exceed that coverage level.

        Args:
            max_coverage (float): Maximum coverage in pileup
            pileup (numpy.array): Coverage depth array
            target_bases (int): Minimum number of bases required

        Returns:
            int: Coverage threshold that achieves target_bases coverage
        """
        coverage_threshold = 0
        # Try coverage levels from max down to 1
        for c in range(int(max_coverage), 0,-1 ):
            base_count = self.num_bases_with_at_least_this_coverage(c, pileup)
            if base_count >= target_bases:
                coverage_threshold = c
                return coverage_threshold
        return coverage_threshold

    def identify_regions(self, fragment_number: Union[str, int], target_bases: int) -> list[list[int]]:
        """
        Identify high-coverage regions in a fragment.

        Finds contiguous blocks where coverage meets the calculated threshold
        needed to cover at least target_bases positions.

        Args:
            fragment_number: Fragment identifier to analyze
            target_bases (int): Minimum number of bases to cover

        Returns:
            list: List of [start, end] coordinate pairs for high-coverage blocks
        """
        # Get coverage pileup
        pileup = self.pileup_fragment(fragment_number)
        max_coverage = numpy.max(pileup)

        # Calculate coverage threshold
        coverage_threshold = self.calc_coverage_threshold(max_coverage, pileup, target_bases)

        # Identify contiguous regions above threshold
        blocks = []
        in_block = False
        block_start = 0
        block_end = 0

        for b in range(0, len(pileup)):
            if pileup[b] >= coverage_threshold and in_block == False:
                # Start of a high-coverage block
                in_block = True
                block_start = b+1
            elif in_block == True and pileup[b] < coverage_threshold:
                # End of a high-coverage block
                in_block = False
                block_end = b+1
                blocks.append([block_start, block_end])

        # Handle case where block extends to end
        if in_block:
            block_end = len(pileup) + 1
            blocks.append([block_start, block_end])
        return blocks
