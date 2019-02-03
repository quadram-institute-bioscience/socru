import csv         
import numpy
import math
from socru.BlastResult  import BlastResult

class FilterBlast:
    def __init__(self, results_file, min_bit_score, min_alignment_length):
        self.results_file = results_file
        self.min_bit_score = min_bit_score
        self.min_alignment_length =  min_alignment_length
        self.results = self.readin_results()
    
    def readin_results(self):
        results = []
        with open(self.results_file , newline='') as blastfile:
            blastreader = csv.reader(blastfile, delimiter='\t')
            for row in blastreader:
                results.append(BlastResult(row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7], row[8], row[9], row[10], row[11]))
        return results
                
    def filter_results(self):
        filtered = [r for r in self.results if r.bit_score > self.min_bit_score and r.alignment_length > self.min_alignment_length ]
        return filtered
        
    def return_top_result(self):
        # assume they are sorted by bit score and return the first
        results = self.filter_results()
        if len(results) > 0:
            return (self.filter_results())[0]
        else:
            return None
            
    def max_coord_from_blast_results(self, blast_result_objs):
        start = [int(r.subject_start) for r in blast_result_objs]
        end = [int(r.subject_end) for r in blast_result_objs]
        return max(start + end)
    
    def pileup_fragment(self, fragment_number):
        fragment_results = [r for r in self.results if str(r.subject) == str(fragment_number)]
        max_coord = self.max_coord_from_blast_results(fragment_results)
        fragment_pileup = numpy.zeros(max_coord)
        
        for r in fragment_results:
            start = int(r.subject_start)
            end = int(r.subject_end)
            if start > end:
                start = int(r.subject_end)
                end = int(r.subject_start)   
            
            for i in range(start, end):
                fragment_pileup[i-1] += 1
        return fragment_pileup
        
    def num_bases_with_at_least_this_coverage(self, target_coverage, pileup):
        if target_coverage == 0:
            return 0
        base_count = sum([1 for p in pileup if p >= target_coverage])
        return base_count
        
        
    def calc_coverage_threshold(self, max_coverage, pileup, target_bases):
        coverage_threshold = 0
        for c in range(int(max_coverage), 0,-1 ):
            base_count = self.num_bases_with_at_least_this_coverage(c, pileup)
            if base_count >= target_bases:
                coverage_threshold = c
                return coverage_threshold
        return coverage_threshold
        
    def identify_regions(self, fragment_number, target_bases):
        pileup = self.pileup_fragment(fragment_number)
        max_coverage = numpy.max(pileup)
        coverage_threshold = self.calc_coverage_threshold(max_coverage, pileup, target_bases)
                
        blocks = []
        in_block = False
        block_start = 0
        block_end = 0
        for b in range(0, len(pileup)):
            if pileup[b] >= coverage_threshold and in_block == False:
                in_block = True
                block_start = b+1
            elif in_block == True and pileup[b] < coverage_threshold:
                in_block = False
                block_end = b+1
                blocks.append([block_start, block_end])
        
        if in_block:
            block_end = len(pileup) + 1 
            blocks.append([block_start, block_end])
        return blocks
                