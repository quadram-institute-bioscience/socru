import csv
import re
import os
from tempfile import mkstemp
from socru.Fasta import Fasta

class Barrnap:
    def __init__(self,input_file, threads, verbose, overlap_margin = 1500, len_70s = 8000, chromosome_length = 0):
        self.input_file = input_file
        self.threads = threads
        self.overlap_margin = overlap_margin
        self.len_70s = len_70s
        self.chromosome_length = chromosome_length
        self.verbose = verbose
        self.files_to_cleanup = []
        
    def decompress_to_file(self):
        m = re.search(".gz$", self.input_file)
        if m:
            fd, decompressed_input_file = mkstemp()
            self.files_to_cleanup.append(decompressed_input_file)
            return decompressed_input_file
        else:
            return self.input_file
        
    def read_barrnap_output(self,barrnap_outputfile):
        coords = self.parse_barrnap_output(barrnap_outputfile)
        return self.find_boundries(coords)
    
    def parse_barrnap_output(self,barrnap_outputfile):
        with open(barrnap_outputfile, newline='') as csvfile:
            bnreader = csv.reader(csvfile, delimiter='\t')
            coords = []
            for row in bnreader:
                if row[0] == '##gff-version 3':
                    continue
                    
                m = re.search(r"Name=([\d]+)S_rRNA", row[8])
                if m:
                    coords.append((int(row[3]), int(row[4]), int(m.group(1)), row[6] ))
            return coords
            
    # Check for coords that are close together - split results.    
    def filter_out_close_start_coords(self, coords):
        for i in range(len(coords)-1):
            if coords[i] < 0:
                continue
            if coords[i] + self.overlap_margin >= coords[i+1]:
                #take the biggest range
                coords[i + 1] = -1
        filtered_coords = [s for s in coords if s >= 0]   
        return filtered_coords
        
    def filter_out_close_end_coords(self, coords):
        for i in range(len(coords)-1):
            if coords[i] < 0:
                continue
            if coords[i] + self.overlap_margin >= coords[i+1]:
                #take the biggest range
                coords[i] = -1
        filtered_coords = [s for s in coords if s >= 0]   
        return filtered_coords
            
    # 5S isnt in every genome, so check its there, otherwise use 23S
    def five_or_23s(self,coords):
        fives = [c for c in coords if c[2] == 5]
        if len(fives) > 0:
            return 5
        else:
            return 23
            
    # todo: refactor
    def find_boundries(self, coords):
        boundries = []
        
        starting_coords = []
        ending_coords = []
        variable_s = self.five_or_23s(coords)
        for c in coords:
            if (c[2] == 16 and c[3] == '+') or (c[2] == variable_s and c[3] == '-'):
                # start of ribo
                starting_coords.append(c[0])  
            elif (c[2] == 16 and c[3] == '-') or (c[2] == variable_s and c[3] == '+'):
                # end of ribo
                ending_coords.append(c[1])
                    
        starting_coords = self.filter_out_close_start_coords(starting_coords)
        ending_coords = self.filter_out_close_end_coords(ending_coords)
           
        for start_index in range(len(starting_coords)):
            start = starting_coords[start_index]
            if start < 0:
                continue
            for end_index in range(len(ending_coords)):
                end = ending_coords[end_index]
                if end < 0:
                    continue
                if end - start < self.len_70s and end - start > 0:
                    boundries.append([start, end])
                    ending_coords[end_index] = -1
                    starting_coords[start_index] = -1
                    continue
        
        # check for 70S that goes over the end of the genome and for errors
        remaining_start_coords = [s for s in starting_coords if s >= 0]   
        remaining_end_coords = [e for e in ending_coords if e >= 0]    
        if len(remaining_start_coords) > 0 and len(remaining_end_coords) > 0:
            chromosome_length = self.chromosome_length
            if self.chromosome_length <= 0:
                chromosome_length = len(Fasta(self.input_file, self.verbose).chromosome)
            
            for start_index in range(len(remaining_start_coords)):
                start = remaining_start_coords[start_index]
                for end_index in range(len(remaining_end_coords)):
                    end = remaining_end_coords[end_index]
                    if end < 0:
                        continue
                    if (chromosome_length - start) < self.len_70s and end < self.len_70s and (chromosome_length - start) + end  < self.len_70s:
                        boundries.append([start, end])
                        remaining_end_coords[end_index] = -1
                        remaining_start_coords[start_index] = -1
                        continue
        
        return boundries
            
    def construct_barrnap_command(self, barrnap_outputfile):
        decompressed_file = self.decompress_to_file()
        if decompressed_file == self.input_file:
            return " ".join(['barrnap', '--quiet', '--threads', str(self.threads),  self.input_file, '>', barrnap_outputfile])
        else:
            return " ".join(['gunzip', '-c', self.input_file, '>', decompressed_file, '&&','barrnap', '--quiet', '--threads', str(self.threads),  decompressed_file, '>', barrnap_outputfile])
        
    def __del__(self):
        for f in self.files_to_cleanup:
            if os.path.exists(f):
                os.remove(f)
                