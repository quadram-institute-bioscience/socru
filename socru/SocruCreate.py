"""
New species database creation workflow.

This module orchestrates the complete workflow for creating a new Socru species
database from a reference genome. It extracts inter-operon fragments, numbers
them, identifies dnaA/dif positions, and generates the initial profile database.

Classes:
    SocruCreate: Manages database creation workflow
"""

import os
import sys
import pkg_resources
from tempfile import mkstemp


from socru.Fasta import Fasta
from socru.FragmentFiles import FragmentFiles
from socru.Barrnap  import Barrnap
from socru.ProfileGenerator import ProfileGenerator

class SocruCreate:
    """
    Create a new Socru species database from a reference genome.

    This class coordinates the complete database creation process:
    1. Identify rRNA operons with Barrnap
    2. Extract inter-operon fragments
    3. Number fragments (largest first)
    4. Write fragment FASTA files
    5. Identify dnaA and dif positions
    6. Generate initial profile.txt with GS1.0

    Attributes:
        input_file (str): Path to reference genome FASTA
        output_directory (str): Directory for new database
        fragment_order (str): Optional custom fragment numbering
        threads (int): Number of CPU threads
        dnaa_fasta (str): Path to dnaA query sequence
        dif_fasta (str): Path to dif query sequence
        verbose (bool): Enable verbose output
        max_bases_from_ends (int): Optional fragment end trimming
        files_to_cleanup (list): Temporary files to delete
    """
    def __init__(self,options):
        """
        Initialize SocruCreate with command-line options.

        Args:
            options: Parsed command-line arguments
        """
        self.input_file = options.input_file
        self.output_directory = options.output_directory
        self.fragment_order = options.fragment_order
        self.threads = options.threads
        self.dnaa_fasta = options.dnaa_fasta
        self.dif_fasta = options.dif_fasta
        self.verbose = options.verbose
        self.max_bases_from_ends = options.max_bases_from_ends
        self.files_to_cleanup = []

        # Validate output directory doesn't exist
        if os.path.exists(self.output_directory):
             print(
             "The output directory already exists, "
             "please choose another name: "
             + self.output_directory)
             sys.exit(1)

        else:
            # Create output directory
            os.makedirs(self.output_directory)

        # Use bundled dnaA and dif sequences if not provided
        if self.dnaa_fasta is None:
            self.dnaa_fasta = str(pkg_resources.resource_filename( __name__, 'data/dnaA.fa.gz'))
        if self.dif_fasta is None:
            self.dif_fasta = str(pkg_resources.resource_filename( __name__, 'data/dif.fa.gz'))

    def run(self):
        """
        Execute complete database creation workflow.

        Steps:
        1. Run Barrnap to find rRNA operons
        2. Calculate fragment coordinates
        3. Extract fragment sequences
        4. Write fragment FASTA files
        5. Generate profile.txt and metadata
        """
        # Step 1: Run Barrnap to identify rRNA operons
        fd, barrnap_outputfile = mkstemp()
        os.close(fd)
        self.files_to_cleanup.append(barrnap_outputfile)
        b = Barrnap(self.input_file, self.threads, self.verbose)
        b.construct_barrnap_command(barrnap_outputfile)

        # Parse Barrnap output to get operon boundaries
        boundries = b.read_barrnap_output(barrnap_outputfile)

        # Step 2: Calculate fragment coordinates from operons
        f = Fasta(self.input_file, self.verbose)
        fragments = f.calc_fragment_coords( boundries)

        # Step 3: Extract fragment sequences
        f.populate_fragments_from_chromosome(fragments, self.max_bases_from_ends)

        # Step 4: Number fragments and write FASTA files
        ff = FragmentFiles(fragments, self.output_directory, self.verbose, fragment_order =  self.fragment_order)
        ff.create_fragment_fastas()

        # Step 5: Create default profile.txt file with dnaA/dif positions
        default_profile = ProfileGenerator(self.output_directory, len(ff.ordered_fragments), self.dnaa_fasta, self.dif_fasta, self.threads, self.input_file, self.verbose )
        default_profile.write_output_file()

    def __del__(self):
        """
        Clean up temporary files when object is destroyed.
        """
        for f in self.files_to_cleanup:
            if os.path.exists(f):
                os.remove(f)
