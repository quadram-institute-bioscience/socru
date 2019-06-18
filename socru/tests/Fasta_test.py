import unittest
import os
import shutil
from socru.Fasta  import Fasta

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data','fasta')

class TestFasta(unittest.TestCase):
   
    def test_fasta_get_largest_contig(self):
        f = Fasta(os.path.join(data_dir,'get_largest_contig.fa'), False)
        largest_contig_record = f.get_chromosome_from_fasta()
        
        self.assertEqual(len(largest_contig_record.seq), 60)
        
    def test_calc_fragment_coords(self):
        f = Fasta(os.path.join(data_dir,'calc_fragment_coords.fa'), False)
        boundries = [[45,55],[90,110],[150,180]]
        fragments = f.calc_fragment_coords(boundries)
        
        coords = [f.coords for f in fragments]
        self.assertEqual(coords, [[[180, 200], [0, 45]], [[55, 90]],  [[110,150]]])
        
        f.populate_fragments_from_chromosome(fragments, None)
        
    def test_calc_fragment_coords_gz(self):
        f = Fasta(os.path.join(data_dir,'calc_fragment_coords.fa.gz'), False)
        boundries = [[45,55],[90,110],[150,180]]
        fragments = f.calc_fragment_coords(boundries)
        
        coords = [f.coords for f in fragments]
        self.assertEqual(coords, [[[180, 200], [0, 45]], [[55, 90]],  [[110,150]]])
        
    def test_populate_fragments_from_chromosome(self):
        f = Fasta(os.path.join(data_dir,'calc_fragment_coords.fa'), False)
        fragments = f.calc_fragment_coords([[45,55],[90,110],[150,180]])
        sequences = [str(f.sequence) for f in fragments]
        self.assertEqual(sequences, ["","",""])
        
        f.populate_fragments_from_chromosome(fragments, None)
        sequences = [str(f.sequence) for f in fragments]
        self.assertEqual(sequences, [
            "TTTTTTTTTTTTTTTTTTTTAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
            "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC",
            "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG"])
			
    def test_chop_from_ends(self):
        f = Fasta(os.path.join(data_dir,'calc_fragment_coords.fa'), False)
        fragments = f.calc_fragment_coords([[45,55],[90,110],[150,180]])
        sequences = [str(f.sequence) for f in fragments]
        
        f.populate_fragments_from_chromosome(fragments, 5)
        sequences = [str(f.sequence) for f in fragments]
        self.assertEqual(sequences, ['TTTTTNNNAAAAA', 'CCCCCNNNCCCCC', 'GGGGGNNNGGGGG'])	
        
        