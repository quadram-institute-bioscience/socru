import unittest
import os
import shutil
import filecmp
from socru.SocruCreate  import SocruCreate

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data','create')


class TestOptions:
    def __init__(self, output_directory, input_file, threads,dnaa_fasta, max_bases_from_ends, fragment_order = None):
        self.output_directory = output_directory
        self.input_file = input_file
        self.fragment_order = fragment_order
        self.threads = threads
        self.max_bases_from_ends = max_bases_from_ends
        self.dnaa_fasta = dnaa_fasta
        self.verbose = False

class TestSocruCreate(unittest.TestCase):

    def test_compressed_input(self):
        if os.path.exists('test_output'):
            shutil.rmtree('test_output')
            
        g = SocruCreate(TestOptions('test_output', os.path.join(data_dir, 'test.fa.gz'), 1, os.path.join(data_dir, 'dnaA.fa.gz'), None, fragment_order = None))
        g.run()
        self.assertTrue(os.path.exists('test_output/1.fa'))
        self.assertTrue(os.path.exists('test_output/2.fa'))
        self.assertTrue(os.path.exists('test_output/3.fa'))
        self.assertTrue(os.path.exists('test_output/4.fa'))
        self.assertFalse(os.path.exists('test_output/5.fa'))
        
        self.assertTrue(filecmp.cmp('test_output/1.fa', os.path.join(data_dir, 'expected_create', '1.fa')))
        self.assertTrue(filecmp.cmp('test_output/2.fa', os.path.join(data_dir, 'expected_create', '2.fa')))
        self.assertTrue(filecmp.cmp('test_output/3.fa', os.path.join(data_dir, 'expected_create', '3.fa')))
        self.assertTrue(filecmp.cmp('test_output/4.fa', os.path.join(data_dir, 'expected_create', '4.fa')))
        self.assertTrue(os.path.exists('test_output/profile.txt'))
        
        if os.path.exists('test_output'):
            shutil.rmtree('test_output')

    def test_socru_create(self):
        if os.path.exists('test_output'):
            shutil.rmtree('test_output')
            
        g = SocruCreate(TestOptions('test_output', os.path.join(data_dir, 'test.fa'), 1, os.path.join(data_dir, 'dnaA.fa.gz'),None, fragment_order = None))
        g.run()
        self.assertTrue(os.path.exists('test_output/1.fa'))
        self.assertTrue(os.path.exists('test_output/2.fa'))
        self.assertTrue(os.path.exists('test_output/3.fa'))
        self.assertTrue(os.path.exists('test_output/4.fa'))
        self.assertFalse(os.path.exists('test_output/5.fa'))
        
        self.assertTrue(filecmp.cmp('test_output/1.fa', os.path.join(data_dir, 'expected_create', '1.fa')))
        self.assertTrue(filecmp.cmp('test_output/2.fa', os.path.join(data_dir, 'expected_create', '2.fa')))
        self.assertTrue(filecmp.cmp('test_output/3.fa', os.path.join(data_dir, 'expected_create', '3.fa')))
        self.assertTrue(filecmp.cmp('test_output/4.fa', os.path.join(data_dir, 'expected_create', '4.fa')))
        
        self.assertTrue(os.path.exists('test_output/profile.txt'))
        
        if os.path.exists('test_output'):
            shutil.rmtree('test_output')
            
    def test_socru_create_set_fragment_order(self):
        if os.path.exists('test_output_reorder'):
            shutil.rmtree('test_output_reorder')  

        g = SocruCreate(TestOptions('test_output_reorder', os.path.join(data_dir, 'test.fa'), 1,os.path.join(data_dir, 'dnaA.fa.gz'),None, fragment_order = "1-4-3'-2"))
        g.run()
        
        self.assertTrue(filecmp.cmp('test_output_reorder/1.fa', os.path.join(data_dir, 'expected_fragment_order', '1.fa')))
        self.assertTrue(filecmp.cmp('test_output_reorder/2.fa', os.path.join(data_dir, 'expected_fragment_order', '2.fa')))
        self.assertTrue(filecmp.cmp('test_output_reorder/3.fa', os.path.join(data_dir, 'expected_fragment_order', '3.fa')))
        self.assertTrue(filecmp.cmp('test_output_reorder/4.fa', os.path.join(data_dir, 'expected_fragment_order', '4.fa')))

        if os.path.exists('test_output_reorder'):
            shutil.rmtree('test_output_reorder')
        
    