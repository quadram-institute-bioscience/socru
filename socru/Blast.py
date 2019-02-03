from tempfile import mkstemp
import subprocess
import os
import re

class Blast:
    def __init__(self, blast_db, threads, word_size = 28, exec = 'blastn', evalue = 0.000001, task = 'megablast'):
        self.blast_db = blast_db
        self.evalue = evalue
        self.threads = threads
        self.word_size = word_size
        self.exec = exec
        self.task = task
        self.files_to_cleanup = []

    def decompress_file_to_tmp(self, input_file):
        m = re.search(".gz$", input_file)
        if m:
            fd, decompressed_input_file = mkstemp()
            self.files_to_cleanup.append(decompressed_input_file)
            cmd = " ".join(['gunzip', '-c', input_file, '>', decompressed_input_file])
            subprocess.check_output( cmd,  shell=True)
            
            return decompressed_input_file
        else:
            return input_file

    def blast_command(self, query, blast_results):
        return " ".join([self.exec, '-outfmt', str(6), '-evalue', str(self.evalue), '-db', self.blast_db, '-word_size', str(self.word_size), '-num_threads', str(self.threads), '-task', self.task, '-query', self.decompress_file_to_tmp(query), '|', 'sort', '-k', str(12), '-g', '-r', '>', blast_results])
        
    def run_blast(self, query):
        fd, blast_results = mkstemp()
        subprocess.check_output( self.blast_command(query, blast_results),  shell=True)
        return blast_results


    def __del__(self):
        for f in self.files_to_cleanup:
            if os.path.exists(f):
                os.remove(f)    
                