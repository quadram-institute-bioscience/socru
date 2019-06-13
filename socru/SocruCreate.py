
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
        
    def run(self):
        # run the fasta through barrnap
        fd, barrnap_outputfile = mkstemp()
        self.files_to_cleanup.append(barrnap_outputfile)
        b = Barrnap(self.input_file, self.threads)
        subprocess.check_output(
           b.construct_barrnap_command(barrnap_outputfile), 
           shell=True)

        boundries = b.read_barrnap_output(barrnap_outputfile)
        
        f = Fasta(self.input_file)
        fragments = f.calc_fragment_coords( boundries)
        f.populate_fragments_from_chromosome(fragments, self.max_bases_from_ends)
        
        ff = FragmentFiles(fragments, self.output_directory, fragment_order =  self.fragment_order)
        ff.create_fragment_fastas()
        
        # create a default profile.txt file
        default_profile = ProfileGenerator(self.output_directory, len(ff.ordered_fragments), self.dnaa_fasta, self.threads)
        default_profile.write_output_file()
        
    def __del__(self):
        for f in self.files_to_cleanup:
            if os.path.exists(f):
                os.remove(f)
            
        