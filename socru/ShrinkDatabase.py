import os
from os import listdir
from os.path import isfile, join
import re
import shutil
import subprocess
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from socru.FilterBlast import FilterBlast
from socru.Fasta import Fasta

class ShrinkDatabase:
    def __init__(self,input_database, output_database, blast_results, target_bases):
        self.input_database = input_database
        self.output_database = output_database
        self.blast_results = blast_results
        self.target_bases = target_bases
    
    # read in the files in the directory starting with number and ending in fa
    def get_database_files(self):
        return [os.path.join(self.input_database,f) for f in listdir(self.input_database) if isfile(join(self.input_database, f)) and re.match("[\d]+\.fa$", f)]
        
    def get_database_files_compressed(self):
        return [os.path.join(self.input_database,f) for f in listdir(self.input_database) if isfile(join(self.input_database, f)) and re.match("[\d]+\.fa.gz$", f)]

    def shrink_files(self):
        output_filenames = []
        # copy profile files to destination
        shutil.copy(os.path.join(self.input_database, 'profile.txt'), self.output_database)
        shutil.copy(os.path.join(self.input_database, 'profile.txt.yml'), self.output_database)
        
        fasta_file_names_compressed = self.get_database_files_compressed()
        fasta_file_names_uncompressed = self.get_database_files()
        fasta_file_names = fasta_file_names_uncompressed + fasta_file_names_compressed
        
        fasta_obj = [ Fasta(f) for f in fasta_file_names]
        for f in fasta_obj:
            fb = FilterBlast(self.blast_results, 1, 1)
            destination_filename = os.path.join(self.output_database, str(f.fragment_number()) + '.fa')
            if len(f.chromosome.seq) < self.target_bases:
                if f.input_file in fasta_file_names_compressed:
                    shutil.copy(f.input_file, destination_filename + '.gz')
                    output_filenames.append(destination_filename + '.gz')
                else:
                    shutil.copy(f.input_file, destination_filename)
                    output_filenames.append(destination_filename)
            else:
                blocks = fb.identify_regions(f.fragment_number(), self.target_bases)
                sequence = ""
                for b in blocks:
                    sequence += f.chromosome.seq[(b[0]):(b[1])]
                record = [SeqRecord(sequence, str(f.fragment_number()) , '', '')]
                SeqIO.write(record, destination_filename, "fasta")
                output_filenames.append(destination_filename)
                
        return self.compress_files(output_filenames)
            
    def compress_files(self, filenames):
        compressed_filenames = []
        for filename in filenames:
            m = re.search("[\d]+\.fa$", filename)
            if m:
                subprocess.check_output( 'gzip '+ filename,  shell=True)
                compressed_filenames.append(filename + '.gz')
            else:
                compressed_filenames.append(filename)
        return compressed_filenames
           
                
