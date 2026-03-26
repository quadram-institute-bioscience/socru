"""
Barrnap integration for rRNA operon identification.

This module wraps the Barrnap tool to identify ribosomal RNA (16S, 23S, 5S) genes
in bacterial genomes. It parses Barrnap's GFF3 output to determine rRNA operon
boundaries and orientations, which define the inter-operon fragments analyzed by Socru.

Classes:
    Barrnap: Manages Barrnap execution and output parsing
"""

import csv
import gzip
import logging
import re
import os
import shutil
import subprocess
from tempfile import mkstemp
from socru.Fasta import Fasta
from socru.Operon import Operon

logger = logging.getLogger(__name__)

class Barrnap:
    """
    Interface to Barrnap tool for rRNA operon detection.

    This class constructs Barrnap commands, executes them, and parses the output
    to identify complete ribosomal operons (70S ribosomes) in bacterial genomes.
    It handles both compressed and uncompressed input files, filters overlapping
    predictions, and determines operon boundaries and orientations.

    Attributes:
        input_file (str): Path to input genome FASTA file
        threads (int): Number of CPU threads for Barrnap
        overlap_margin (int): Distance threshold for filtering overlaps (default 1500bp)
        len_70s (int): Expected size of 70S ribosome operon (default 8000bp)
        chromosome_length (int): Length of chromosome if known
        verbose (bool): Enable verbose output
        files_to_cleanup (list): Temporary files to delete on cleanup
    """
    def __init__(self,input_file, threads, verbose, overlap_margin = 1500, len_70s = 8000, chromosome_length = 0):
        """
        Initialize Barrnap wrapper.

        Args:
            input_file (str): Path to genome FASTA file
            threads (int): Number of threads for parallel processing
            verbose (bool): Enable verbose output
            overlap_margin (int): Distance for filtering overlapping hits (default 1500)
            len_70s (int): Expected 70S operon size in bases (default 8000)
            chromosome_length (int): Chromosome length if pre-known (default 0)
        """
        self.input_file = input_file
        self.threads = threads
        self.overlap_margin = overlap_margin
        self.len_70s = len_70s
        self.chromosome_length = chromosome_length
        self.verbose = verbose
        self.files_to_cleanup = []

    def decompress_to_file(self):
        """
        Decompress gzipped input if necessary.

        Returns:
            str: Path to decompressed file, or original path if not gzipped
        """
        m = re.search(".gz$", self.input_file)
        if m:
            # Create temporary file for decompressed content
            fd, decompressed_input_file = mkstemp()
            os.close(fd)
            self.files_to_cleanup.append(decompressed_input_file)
            # Decompress using Python gzip module
            with gzip.open(self.input_file, 'rb') as gz_in:
                with open(decompressed_input_file, 'wb') as f_out:
                    shutil.copyfileobj(gz_in, f_out)
            return decompressed_input_file
        else:
            return self.input_file

    def read_barrnap_output(self,barrnap_outputfile):
        """
        Read and process Barrnap output file.

        Args:
            barrnap_outputfile (str): Path to Barrnap GFF3 output

        Returns:
            list: List of Operon objects representing rRNA operon boundaries
        """
        coords = self.parse_barrnap_output(barrnap_outputfile)
        return self.find_boundries(coords)

    def parse_barrnap_output(self,barrnap_outputfile):
        """
        Parse Barrnap GFF3 output file to extract rRNA gene coordinates.

        Reads tab-delimited GFF3 format and extracts:
        - Start and end coordinates
        - rRNA type (5S, 16S, 23S)
        - Strand orientation (+/-)

        Args:
            barrnap_outputfile (str): Path to Barrnap output file

        Returns:
            list: List of tuples (start, end, rRNA_size, strand)
        """
        with open(barrnap_outputfile, newline='') as csvfile:
            bnreader = csv.reader(csvfile, delimiter='\t')
            coords = []
            for row in bnreader:
                # Skip GFF header line
                if row[0] == '##gff-version 3':
                    continue

                # Extract rRNA type from Name attribute (e.g., "Name=16S_rRNA")
                m = re.search(r"Name=([\d]+)S_rRNA", row[8])
                if m:
                    # Store (start, end, rRNA_type, strand)
                    # GFF uses 1-based inclusive coordinates; convert start to
                    # 0-based so that Python slicing works correctly. The end
                    # coordinate needs no adjustment because GFF's 1-based
                    # inclusive end equals Python's 0-based exclusive end.
                    coords.append((int(row[3]) - 1, int(row[4]), int(m.group(1)), row[6] ))
            return coords

    def filter_out_close_start_coords(self, coords):
        """
        Filter out overlapping start coordinates by proximity.

        When multiple rRNA predictions have start coordinates within overlap_margin
        of each other, only keeps the first one. This removes redundant predictions.

        Args:
            coords (list): List of start coordinates

        Returns:
            list: Filtered list with overlaps removed
        """
        # Check for coords that are close together - split results.
        for i in range(len(coords)-1):
            if coords[i] < 0:
                continue
            if coords[i] + self.overlap_margin >= coords[i+1]:
                # Mark next coord for removal (take the first of the pair)
                coords[i + 1] = -1
        # Filter out marked coordinates (negative values)
        filtered_coords = [s for s in coords if s >= 0]
        return filtered_coords

    def filter_out_close_end_coords(self, coords):
        """
        Filter out overlapping end coordinates by proximity.

        When multiple rRNA predictions have end coordinates within overlap_margin
        of each other, only keeps the last one. This removes redundant predictions.

        Args:
            coords (list): List of end coordinates

        Returns:
            list: Filtered list with overlaps removed
        """
        for i in range(len(coords)-1):
            if coords[i] < 0:
                continue
            if coords[i] + self.overlap_margin >= coords[i+1]:
                # Mark current coord for removal (take the last of the pair)
                coords[i] = -1
        # Filter out marked coordinates (negative values)
        filtered_coords = [s for s in coords if s >= 0]
        return filtered_coords

    def five_or_23s(self,coords):
        """
        Determine whether to use 5S or 23S rRNA for operon end detection.

        5S rRNA is not present in every genome, so check if it exists.
        If 5S is found, use it; otherwise fall back to 23S.

        Args:
            coords (list): List of rRNA coordinate tuples

        Returns:
            int: Either 5 or 23 indicating which rRNA type to use
        """
        fives = [c for c in coords if c[2] == 5]
        if len(fives) > 0:
            return 5
        else:
            return 23

    def find_boundries(self, coords):
        """
        Identify complete ribosomal operon boundaries from rRNA gene coordinates.

        This complex method matches start and end rRNA genes to form complete
        operons. It:
        1. Identifies operon starts (16S+ or 5S-/23S-)
        2. Identifies operon ends (16S- or 5S+/23S+)
        3. Pairs starts with ends based on expected 70S size
        4. Handles operons that wrap around chromosome origin
        5. Determines operon orientation (forward/reverse)

        Args:
            coords (list): List of (start, end, rRNA_type, strand) tuples

        Returns:
            list: List of Operon objects with start, end, and direction
        """
        boundries = []

        # Track whether each coordinate is on forward strand
        coord_forward = {}
        starting_coords = []
        ending_coords = []

        # Determine which rRNA marks operon end (5S or 23S)
        variable_s = self.five_or_23s(coords)

        # Classify each rRNA gene as operon start or end based on type and strand
        for c in coords:
            if (c[2] == 16 and c[3] == '+'):
                # 16S on + strand marks operon start (forward direction)
                coord_forward[c[0]] = True
                starting_coords.append(c[0])
            elif (c[2] == variable_s and c[3] == '-'):
                # 5S/23S on - strand marks operon start (reverse direction)
                coord_forward[c[0]] = False
                starting_coords.append(c[0])
            elif (c[2] == 16 and c[3] == '-'):
                # 16S on - strand marks operon end (reverse direction)
                coord_forward[c[0]] = False
                ending_coords.append(c[1])
            elif (c[2] == variable_s and c[3] == '+'):
                # 5S/23S on + strand marks operon end (forward direction)
                coord_forward[c[0]] = True
                ending_coords.append(c[1])

        # Filter out overlapping/redundant coordinates
        starting_coords = self.filter_out_close_start_coords(starting_coords)
        ending_coords = self.filter_out_close_end_coords(ending_coords)

        # Match starts with ends to form complete operons
        for start_index in range(len(starting_coords)):
            start = starting_coords[start_index]
            if start < 0:
                continue
            for end_index in range(len(ending_coords)):
                end = ending_coords[end_index]
                if end < 0:
                    continue
                # Check if distance matches expected 70S size
                if end - start < self.len_70s and end - start > 0:

                    # Determine operon orientation
                    direction = True
                    if start in coord_forward:
                        direction = coord_forward[start]
                    elif end in coord_forward:
                        direction = coord_forward[end]

                    boundries.append(Operon(start, end, direction))
                    # Mark these coordinates as used
                    ending_coords[end_index] = -1
                    starting_coords[start_index] = -1
                    continue

        # Handle operons that wrap around the chromosome end/start
        remaining_start_coords = [s for s in starting_coords if s >= 0]
        remaining_end_coords = [e for e in ending_coords if e >= 0]
        if len(remaining_start_coords) > 0 and len(remaining_end_coords) > 0:
            # Get chromosome length if not already known
            chromosome_length = self.chromosome_length
            if self.chromosome_length <= 0:
                chromosome_length = len(Fasta(self.input_file, self.verbose).chromosome)

            # Check for operons spanning origin (start near end, end near beginning)
            for start_index in range(len(remaining_start_coords)):
                start = remaining_start_coords[start_index]
                for end_index in range(len(remaining_end_coords)):
                    end = remaining_end_coords[end_index]
                    if end < 0:
                        continue
                    # Check if operon wraps: start near end + end near beginning = ~70S size
                    if (chromosome_length - start) < self.len_70s and end < self.len_70s and (chromosome_length - start) + end  < self.len_70s:
                        # Determine direction
                        direction = True
                        if start in coord_forward:
                            direction = coord_forward[start]
                        elif end in coord_forward:
                            direction = coord_forward[end]
                        boundries.append(Operon(start, end, direction))

                        remaining_end_coords[end_index] = -1
                        remaining_start_coords[start_index] = -1
                        continue

        return boundries

    def construct_barrnap_command(self, barrnap_outputfile):
        """
        Build the command list and run Barrnap.

        Handles decompression of gzipped input files using Python gzip module.
        Runs barrnap as a subprocess with safe argument list (no shell=True).

        Args:
            barrnap_outputfile (str): Path where Barrnap output should be written
        """
        decompressed_file = self.decompress_to_file()
        cmd = ['barrnap', '--quiet', '--threads', str(self.threads), decompressed_file]
        logger.info("Run barrnap:\t%s > %s", ' '.join(cmd), barrnap_outputfile)
        with open(barrnap_outputfile, 'w') as out_fh:
            subprocess.run(cmd, stdout=out_fh, check=True)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.cleanup()
        return False

    def cleanup(self):
        """Clean up temporary files and directories."""
        for f in self.files_to_cleanup:
            if os.path.exists(f):
                os.remove(f)
        self.files_to_cleanup = []

    def __del__(self):
        """Safety net cleanup -- prefer using as context manager."""
        self.cleanup()
