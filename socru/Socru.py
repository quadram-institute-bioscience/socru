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

class Socru:
    def __init__(self,options):
        self.input_files = options.input_files
        
        self.min_bit_score = options.min_bit_score
        self.min_alignment_length = options.min_alignment_length
        self.threads = options.threads
        self.output_file = options.output_file
        self.novel_profiles = options.novel_profiles
        self.new_fragments = options.new_fragments
        self.max_bases_from_ends = options.max_bases_from_ends
        self.top_blast_hits = options.top_blast_hits
        self.output_plot_file = options.output_plot_file
        self.dirs_to_cleanup = []
        self.top_results = []
        
        self.db_dir =  Schemas().database_directory(options.db_dir, options.species)
        if self.db_dir is None:
             print(
             "Cannot access the database you specified, please check again")
             sys.exit(1)
        
        
        if options.not_circular:
            self.is_circular = False
        else:
            self.is_circular = True
            
    def run(self):
        # load the profiles
        p = Profiles(os.path.join(self.db_dir, 'profile.txt'))
        d = Database(self.db_dir)
        
        for i in self.input_files:
            output_type = self.run_analysis(i, p, d)
            self.output_results(i, output_type)
            
        if self.top_blast_hits is not None:
            with open(self.top_blast_hits, "a+") as output_fh:
                for h in self.top_results:
                    output_fh.write(str(h)+"\n")
            
    def output_results(self, input_file, profile_type):
        if self.output_file is None:
            print(input_file + "\t" + profile_type)
        else:
            with open(self.output_file, "a+") as output_fh:
                output_fh.write(input_file + "\t" + profile_type + "\n") 
            
    def find_rrna_boundries(self, input_file):
        # run the fasta through barrnap
        fd, barrnap_outputfile = mkstemp()
        b = Barrnap(input_file, self.threads)
        subprocess.check_output(
           b.construct_barrnap_command(barrnap_outputfile), 
           shell=True)

        boundries = b.read_barrnap_output(barrnap_outputfile)
        os.close(fd)
        return boundries
    
    # refactor
    def run_analysis(self, input_file, p, d):
        boundries = self.find_rrna_boundries(input_file)
        
        f = Fasta(input_file, is_circular = self.is_circular)
        fragments = f.calc_fragment_coords( boundries)
        f.populate_fragments_from_chromosome(fragments, self.max_bases_from_ends)
        
        tmpdir = mkdtemp()
        self.dirs_to_cleanup.append(tmpdir)

        ff = FragmentFiles(fragments, tmpdir)
        ff.create_fragment_fastas()
        
         # take each fasta file and blast it against the database
        blast = Blast(d.db_prefix, self.threads)
        
        gat_profile = GATProfile(fragments = [])
        for current_fragment in ff.ordered_fragments:
            fasta_file = current_fragment.output_filename
            
            blast_results = blast.run_blast(fasta_file)
            fb = FilterBlast(blast_results, self.min_bit_score, self.min_alignment_length)
            top_result = fb.return_top_result()
            if top_result is None:
                gat_profile.fragments.append('?')
                current_fragment.number = '?'
                
                with open(fasta_file, "r") as fasta_file_fh:
                    with open(self.new_fragments, "a+") as newfrag_fh:
                        newfrag_fh.write(fasta_file_fh.read())
                continue
            else:
                self.top_results.append(top_result) 
                
                current_fragment.number = str(top_result.subject)
                if str(p.dnaA_fragment_number) == str(current_fragment.number):
                    current_fragment.dna_A = True
                
                if top_result.is_forward():
                    gat_profile.fragments.append( str(top_result.subject))
                else:
                    current_fragment.reversed_frag = True
                    gat_profile.fragments.append( str(top_result.subject)+ '\'')
        
        gat_profile.orientate_for_dnaA()
        reordered_frag_obj = gat_profile.reorder_fragment_objects_based_on_fragment_name_array( ff.ordered_fragments )
        pp = PlotProfile(ff.ordered_fragments, self.output_plot_file)
        pp.create_plot()
        
        # lookup the gat_profile to get the number
        tg = TypeGenerator(p, gat_profile)
        type_output_string  =  tg.calculate_type() + "\t" + str(gat_profile)
        if not tg.has_previously_seen:
            with open(self.novel_profiles, "a+") as output_fh:
                output_fh.write(self.db_dir + "\t" + type_output_string + "\n")
        
        return type_output_string

    def __del__(self):
        for f in self.dirs_to_cleanup:
            if os.path.exists(f):
                shutil.rmtree(f)
        