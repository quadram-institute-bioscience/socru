"""
Main Socru analysis module for genome structure typing.

This module orchestrates the complete Socru workflow for analyzing bacterial genomes
to determine their genome structure (GS) types based on the arrangement of fragments
around ribosomal operons. It coordinates rRNA identification, fragment extraction,
BLAST searching, and profile matching.

Classes:
    Socru: Main class that runs the complete analysis pipeline
"""

import logging
import os
import sys
import time

import subprocess
from tempfile import mkstemp
from tempfile import mkdtemp
import pkg_resources
import shutil

from socru.Fasta import Fasta
from socru.FragmentFiles import FragmentFiles
from socru.Barrnap  import Barrnap
from socru.Database import Database
from socru.Blast import Blast
from socru.FilterBlast import FilterBlast
from socru.Profiles import Profiles
from socru.GATProfile import GATProfile
from socru.TypeGenerator import TypeGenerator
from socru.Schemas import Schemas
from socru.PlotProfile import PlotProfile
from socru.ValidateFragments import ValidateFragments

class Socru:
    """
    Main class for running Socru genome structure analysis.
    
    This class manages the complete workflow: loading databases, finding rRNA operons,
    extracting inter-operon fragments, BLASTing against the fragment database, matching
    to known profiles, and generating output including plots and novel fragment reports.
    
    Attributes:
        input_files (list): List of input FASTA files to analyze
        min_bit_score (int): Minimum BLAST bit score for fragment matching
        min_alignment_length (int): Minimum BLAST alignment length
        threads (int): Number of threads for parallel processing
        output_file (str): Path to main output file
        novel_profiles (str): Path to file for novel profile patterns
        new_fragments (str): Path to file for unmatched fragments
        max_bases_from_ends (int): If set, only use this many bases from fragment ends
        top_blast_hits (str): Path to file for detailed BLAST results
        output_plot_file (str): Path to PDF plot of genome structure
        output_operon_directions_file (str): Path to operon direction results
        verbose (bool): Enable verbose output
        db_dir (str): Path to species database directory
        is_circular (bool): Whether chromosomes are assumed circular
        dirs_to_cleanup (list): Temporary directories to remove on cleanup
        top_results (list): Collection of top BLAST hits for all fragments
    """
    def __init__(self,options):
        """
        Initialize Socru with command-line options.
        
        Args:
            options: Parsed command-line arguments containing all configuration
        """
        self.input_files = options.input_files
        
        # BLAST filtering parameters
        self.min_bit_score = options.min_bit_score
        self.min_alignment_length = options.min_alignment_length
        self.threads = options.threads
        
        # Output file paths
        self.output_file = options.output_file
        self.novel_profiles = options.novel_profiles
        self.new_fragments = options.new_fragments
        self.max_bases_from_ends = options.max_bases_from_ends
        self.top_blast_hits = options.top_blast_hits
        self.output_plot_file = options.output_plot_file
        self.output_operon_directions_file = options.output_operon_directions_file
        
        self.verbose = options.verbose
        self.dirs_to_cleanup = []
        self.top_results = []
        
        # Locate and validate the species database
        self.db_dir =  Schemas(self.verbose).database_directory(options.db_dir, options.species)
        if self.db_dir is None:
             print(
             "Cannot access the database you specified, please check again")
             sys.exit(1)
        
        # Set chromosome circularity assumption
        if options.not_circular:
            self.is_circular = False
        else:
            self.is_circular = True
            
    def run(self):
        """
        Execute the main Socru analysis workflow for all input files.
        
        This method:
        1. Loads the profile database and fragment database
        2. Iterates through each input genome file
        3. Runs complete analysis on each file
        4. Outputs results to file or stdout
        5. Optionally saves detailed BLAST hits
        """
        # load the profiles
        if self.verbose:
            print("Loading profiles:\t" + os.path.join(self.db_dir, 'profile.txt'))
        p = Profiles(os.path.join(self.db_dir, 'profile.txt'), self.verbose)
        
        if self.verbose:
            print("Loading database:\t" + self.db_dir)
        d = Database(self.db_dir, self.verbose)
        
        # Process each input genome file
        for i in self.input_files:
            if self.verbose:
                print("Beginning analysis of input file:\t" + i)
            output_type = self.run_analysis(i, p, d)
            self.output_results(i, output_type)
        
        # Write all top BLAST hits to file if requested    
        if self.top_blast_hits is not None:
            with open(self.top_blast_hits, "a+") as output_fh:
                for h in self.top_results:
                    output_fh.write(str(h)+"\n")
            
    def output_results(self, input_file, profile_type):
        """
        Output analysis results for a single genome.
        
        Args:
            input_file (str): Name of input genome file
            profile_type (str): Determined GS type and quality (e.g., "GREEN\tGS1.1")
        """
        # Default to RED quality with GS0.0 if no valid type found
        if profile_type == '':
            profile_type = "RED\tGS0.0"
        
        # Write to stdout or to specified output file
        if self.output_file is None:
            print(input_file + "\t" + profile_type)
        else:
            with open(self.output_file, "a+") as output_fh:
                output_fh.write(input_file + "\t" + profile_type + "\n") 
            
    def find_rrna_boundries(self, input_file):
        """
        Identify rRNA operon boundaries in the genome using Barrnap.
        
        Args:
            input_file (str): Path to input genome FASTA file
            
        Returns:
            list: List of boundary objects marking rRNA operon locations
        """
        # run the fasta through barrnap
        fd, barrnap_outputfile = mkstemp()
        b = Barrnap(input_file, self.threads, self.verbose)
        cmd = b.construct_barrnap_command(barrnap_outputfile)
        if self.verbose:
            print("Finding rRNA boundries:\t" + cmd)
        subprocess.check_output(cmd, shell=True)

        # Parse barrnap output to extract boundary coordinates
        boundries = b.read_barrnap_output(barrnap_outputfile)
        if self.verbose:
            print("Boundries:")
            for b in boundries:
                print(b)
        os.close(fd)
        return boundries
    
    def populate_fragments_from_chromosome(self, input_file, boundries):
        """
        Extract inter-operon fragment sequences from the chromosome.
        
        Args:
            input_file (str): Path to input genome FASTA file
            boundries (list): rRNA operon boundary locations
            
        Returns:
            list: Fragment objects containing coordinates and sequences
        """
        f = Fasta(input_file, self.verbose, is_circular = self.is_circular)
        # Calculate fragment coordinates based on operon boundaries
        fragments = f.calc_fragment_coords( boundries)
        # Extract actual sequences for each fragment
        f.populate_fragments_from_chromosome(fragments, self.max_bases_from_ends)
        return fragments
        
    def write_novel_profile_to_file(self, tg, type_output_string):
        """
        Save novel (previously unseen) genome structure profiles to file.
        
        Args:
            tg (TypeGenerator): Type generator containing novelty information
            type_output_string (str): Full profile description string
        """
        if not tg.has_previously_seen:
            with open(self.novel_profiles, "a+") as output_fh:
                output_fh.write(self.db_dir + "\t" + type_output_string + "\n")
    
    # refactor
    def run_analysis(self, input_file, p, d):
        """
        Run complete analysis workflow for a single genome.
        
        This is the core analysis method that:
        1. Finds rRNA operons
        2. Extracts inter-operon fragments
        3. BLASTs fragments against database
        4. Matches fragment pattern to known profiles
        5. Validates fragment ordering
        6. Generates plots and saves novel data
        
        Args:
            input_file (str): Path to input genome FASTA
            p (Profiles): Loaded profile database
            d (Database): Fragment BLAST database
            
        Returns:
            str: Type output string (quality + GS type + fragment pattern)
        """
        # Step 1: Find rRNA operon boundaries
        boundries = self.find_rrna_boundries(input_file)
        if not boundries:
            return ''
        
        # Step 2: Extract fragment sequences
        fragments = self.populate_fragments_from_chromosome(input_file, boundries)
        
        # Step 3: Create temporary FASTA files for each fragment
        tmpdir = mkdtemp()
        self.dirs_to_cleanup.append(tmpdir)
        ff = FragmentFiles(fragments, tmpdir, self.verbose)
        ff.create_fragment_fastas()
        
        # Step 4: BLAST each fragment against the database
        blast = Blast(d.db_prefix, self.threads, self.verbose)
        
        # Build the fragment profile pattern
        gat_profile = GATProfile(self.verbose, fragments = [])
        for current_fragment in ff.ordered_fragments:
            fasta_file = current_fragment.output_filename
            
            # BLAST the fragment and filter results
            blast_results = blast.run_blast(fasta_file)
            fb = FilterBlast(blast_results, self.min_bit_score, self.min_alignment_length, self.verbose)
            top_result = fb.return_top_result()
            
            if top_result is None:
                # No match found - mark as unknown and save for manual review
                gat_profile.fragments.append('?')
                current_fragment.number = '?'
                
                with open(fasta_file, "r") as fasta_file_fh:
                    with open(self.new_fragments, "a+") as newfrag_fh:
                        newfrag_fh.write(fasta_file_fh.read())
                continue
            else:
                self.top_results.append(top_result) 
                
                # Record fragment number from BLAST hit
                current_fragment.number = str(top_result.subject)
                
                # Mark special fragments (dnaA and dif)
                if str(p.dnaA_fragment_number) == str(current_fragment.number):
                    current_fragment.dna_A = True
					
                if str(p.dif_fragment_number) == str(current_fragment.number):
                    current_fragment.dif = True
                
                # Add fragment to profile with orientation indicator
                if top_result.is_forward():
                    gat_profile.fragments.append( str(top_result.subject))
                else:
                    # Mark reversed fragments with prime (')
                    current_fragment.reversed_frag = True
                    gat_profile.fragments.append( str(top_result.subject)+ '\'')
        
        # Step 5: Orient profile starting from dnaA and reorder fragments
        gat_profile.orientate_for_dnaA()
        reordered_frag_objs = gat_profile.reorder_fragment_objects_based_on_fragment_name_array( ff.ordered_fragments )
        
        # Step 6: Validate that fragments are ordered correctly
        validate_fragments = ValidateFragments(ff.ordered_fragments)
        is_frag_valid = validate_fragments.validate()
        
        # Build operon direction string for output
        operon_directions_str = " ".join([current_fragment.operon_direction_str() for current_fragment in ff.ordered_fragments])
        if is_frag_valid:
            operon_directions_str = "Valid\t" + operon_directions_str
        else:
            operon_directions_str = "Invalid\t" + operon_directions_str
        
        if self.verbose:
            print("Operon directions:\t" +  operon_directions_str)
        self.output_operon_direction(input_file, operon_directions_str)
        
        # Step 7: Match profile pattern to known GS type
        tg = TypeGenerator(p, gat_profile, self.verbose, is_frag_valid)
        type_output_string  =  tg.quality + "\t" + tg.calculate_type() + "\t" + str(gat_profile)
        self.write_novel_profile_to_file(tg, type_output_string)
        
        # Step 8: Generate genome structure plot for high-quality results
        if tg.quality == 'GREEN':
            pp = PlotProfile(reordered_frag_objs, self.output_plot_file, self.verbose)
            pp.create_plot()
        
        return type_output_string
        
    def output_operon_direction(self, input_file, operon_directions):
        """
        Write operon direction information to output file.
        
        Args:
            input_file (str): Name of input genome file
            operon_directions (str): String describing operon orientations
        """
        with open(self.output_operon_directions_file, "a+") as output_fh:
            output_fh.write(input_file + "\t" + operon_directions + "\n")

    def __del__(self):
        """
        Cleanup temporary directories when object is destroyed.
        """
        for f in self.dirs_to_cleanup:
            if os.path.exists(f):
                shutil.rmtree(f)
        