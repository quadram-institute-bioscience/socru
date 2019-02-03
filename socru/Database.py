import os
from os import listdir
from os.path import isfile, join
import re
from tempfile import mkstemp
from tempfile import mkdtemp
import subprocess
import shutil

class Database:
    def __init__(self,directory):
        self.directory = directory
        self.concat_fasta = self.concat_db_files()
        self.db_prefix = self.make_blastdb(self.concat_fasta)
    
    # read in the files in the directory starting with number and ending in fa
    def get_database_files(self):
        return [os.path.join(self.directory,f) for f in listdir(self.directory) if isfile(join(self.directory, f)) and re.match("[\d]+\.fa$", f)]
        
    def get_database_files_compressed(self):
        return [os.path.join(self.directory,f) for f in listdir(self.directory) if isfile(join(self.directory, f)) and re.match("[\d]+\.fa.gz$", f)]
        
    # concat into a temp file. It can take a mixture of gz and non gzip files
    def concat_db_files(self):
        fd, concat_db_fasta = mkstemp()
        
        if len(self.get_database_files_compressed()) > 0:
            cmd = " ".join(["gunzip", '-c'] + self.get_database_files_compressed() + [ '>' + concat_db_fasta])
            subprocess.check_output( cmd, shell=True)
            
        if len(self.get_database_files()) > 0:
            cmd = " ".join(["cat"] + self.get_database_files() + [ '>>' + concat_db_fasta])
            subprocess.check_output( cmd, shell=True)
        return concat_db_fasta
        
     # make a blast database
    def make_blastdb(self, concat_fasta):
        tmpdir = mkdtemp()
        output_prefix = os.path.join(tmpdir, 'all')
        cmd = " ".join(['makeblastdb', '-in', concat_fasta, '-dbtype', 'nucl',  '-out', output_prefix])
        subprocess.check_output(cmd, shell=True)
        return output_prefix
    
    def __del__(self):
        if os.path.exists(self.db_prefix):
            shutil.rmtree(self.db_prefix)
        if os.path.exists(self.concat_fasta):
            os.remove(self.concat_fasta)
            