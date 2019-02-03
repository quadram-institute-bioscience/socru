import re
import csv
import os
import yaml
from socru.DnaA  import DnaA

class ProfileGenerator:
    def __init__(self, output_directory, num_fragments, dnaa_fasta, threads, prefix = 'GS', output_filename = 'profile.txt', metadata_file_suffix = '.yml',):
        self.output_directory = output_directory
        self.output_filename = os.path.join(self.output_directory, output_filename)
        self.prefix = prefix
        self.num_fragments = num_fragments
        self.dnaa_fasta = dnaa_fasta
        self.threads = threads
        self.metadata_file_suffix = metadata_file_suffix
        
    def header(self):
        header_row = [self.prefix]
        for i in range(self.num_fragments):
            header_row.append('Frag_'+str(i +1)) 
        return header_row
        
    def default_profile(self):
        default_row = ['1.0']
        for i in range(self.num_fragments):
            default_row.append(str(i + 1)) 
        return default_row
    
    def write_output_file(self):
        content = [self.header(), self.default_profile()]
        with open(self.output_filename, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, delimiter='\t')
            writer.writerows(content)
            
        self.write_metadata_file()
    
    def write_metadata_file(self):
        d = DnaA(self.dnaa_fasta, self.output_directory, self.threads)
        metadata_file =  self.output_filename + self.metadata_file_suffix
        metadata_content = { 'dnaa_fragment': int(d.fragment_with_dnaa), 'dnaa_forward_orientation': d.forward_orientation }
        
        with open(metadata_file, 'w') as yaml_file:
            yaml.dump(metadata_content, yaml_file, default_flow_style=False)
        
        