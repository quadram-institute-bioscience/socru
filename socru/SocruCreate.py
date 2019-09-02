
import os
import sys
import time
import pkg_resources
import subprocess
from tempfile import mkstemp


from socru.Fasta import Fasta
from socru.FragmentFiles import FragmentFiles
from socru.Barrnap  import Barrnap
from socru.ProfileGenerator import ProfileGenerator

class SocruCreate:
    def __init__(self,options):
        self.input_file = options.input_file
        self.output_directory = options.output_directory
        self.fragment_order = options.fragment_order
        self.threads = options.threads
        self.dnaa_fasta = options.dnaa_fasta
        self.dif_fasta = options.dif_fasta
        self.verbose = options.verbose
        self.max_bases_from_ends = options.max_bases_from_ends
        self.files_to_cleanup = []
        
        if os.path.exists(self.output_directory):
             print(
             "The output directory already exists, "
             "please choose another name: "
             + self.output_directory)
             sys.exit(1)
         
        else:
            os.makedirs(self.output_directory)
            
        if self.dnaa_fasta is None:
            self.dnaa_fasta = str(pkg_resources.resource_filename( __name__, 'data/dnaA.fa.gz'))
        if self.dif_fasta is None:
            self.dif_fasta = str(pkg_resources.resource_filename( __name__, 'data/dif.fa.gz'))
        
    def run(self):
        # run the fasta through barrnap
        fd, barrnap_outputfile = mkstemp()
        self.files_to_cleanup.append(barrnap_outputfile)
        b = Barrnap(self.input_file, self.threads, self.verbose)
        subprocess.check_output(
           b.construct_barrnap_command(barrnap_outputfile), 
           shell=True)

        boundries = b.read_barrnap_output(barrnap_outputfile)
        
        f = Fasta(self.input_file, self.verbose)
        fragments = f.calc_fragment_coords( boundries)
        f.populate_fragments_from_chromosome(fragments, self.max_bases_from_ends)
        
        ff = FragmentFiles(fragments, self.output_directory, self.verbose, fragment_order =  self.fragment_order)
        ff.create_fragment_fastas()
        
        # create a default profile.txt file
        default_profile = ProfileGenerator(self.output_directory, len(ff.ordered_fragments), self.dnaa_fasta, self.dif_fasta, self.threads, self.input_file, self.verbose )
        default_profile.write_output_file()
        
    def __del__(self):
        for f in self.files_to_cleanup:
            if os.path.exists(f):
                os.remove(f)
            
        