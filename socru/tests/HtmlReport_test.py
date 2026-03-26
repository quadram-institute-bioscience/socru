import unittest
from socru.HtmlReport import HtmlReport


class TestHtmlReport(unittest.TestCase):
    def test_basic_report_generation(self):
        results = [{
            'genome_file': 'test.fasta',
            'gs_type': 'GS1.0',
            'quality': 'GREEN',
            'confidence_score': 95.0,
            'fragment_pattern': '1 2 3 4 5',
            'num_operons': 5,
            'is_novel': False,
            'qc_flags': [],
            'fragments': [
                {'number': 1, 'reversed': False, 'blast_identity': 99.5,
                 'blast_alignment_length': 50000, 'blast_bit_score': 92000,
                 'length': 50000}
            ],
            'genome_length': 4800000,
            'validation_passed': True,
        }]
        report = HtmlReport(results, species='Escherichia_coli')
        html = report.generate()
        self.assertIn('<!DOCTYPE html>', html)
        self.assertIn('Socru', html)
        self.assertIn('GS1.0', html)
        self.assertIn('Escherichia_coli', html)
        self.assertIn('GREEN', html)
        self.assertIn('test.fasta', html)

    def test_multiple_results(self):
        results = [
            {'genome_file': f'sample_{i}.fasta', 'gs_type': f'GS{i}.0',
             'quality': q, 'confidence_score': 90 - i * 10,
             'fragment_pattern': '1 2 3', 'num_operons': 3,
             'is_novel': i > 2, 'qc_flags': [], 'fragments': [],
             'genome_length': 3000000, 'validation_passed': True}
            for i, q in enumerate(['GREEN', 'GREEN', 'AMBER', 'RED'], 1)
        ]
        report = HtmlReport(results, species='Test_species')
        html = report.generate()
        self.assertIn('4', html)  # total count
        for r in results:
            self.assertIn(r['genome_file'], html)

    def test_qc_flags_displayed(self):
        results = [{
            'genome_file': 'test.fasta', 'gs_type': 'GS1.0',
            'quality': 'AMBER', 'confidence_score': 65.0,
            'fragment_pattern': '1 2 ? 4', 'num_operons': 4,
            'is_novel': False,
            'qc_flags': [{'code': 'MISSING_FRAGMENT', 'severity': 'warning',
                          'message': 'Fragment 3 not identified'}],
            'fragments': [], 'genome_length': 3000000,
            'validation_passed': True,
        }]
        report = HtmlReport(results, species='Test_species')
        html = report.generate()
        self.assertIn('MISSING_FRAGMENT', html)

    def test_save_to_file(self):
        import tempfile
        import os

        results = [{
            'genome_file': 'test.fasta', 'gs_type': 'GS1.0',
            'quality': 'GREEN', 'confidence_score': 95.0,
            'fragment_pattern': '1 2 3', 'num_operons': 3,
            'is_novel': False, 'qc_flags': [], 'fragments': [],
            'genome_length': 3000000, 'validation_passed': True,
        }]
        report = HtmlReport(results, species='Test')
        fd, path = tempfile.mkstemp(suffix='.html')
        os.close(fd)
        try:
            report.save(path)
            self.assertTrue(os.path.exists(path))
            with open(path) as f:
                content = f.read()
            self.assertIn('<!DOCTYPE html>', content)
        finally:
            os.unlink(path)

    def test_empty_results(self):
        report = HtmlReport([], species='Empty_species')
        html = report.generate()
        self.assertIn('<!DOCTYPE html>', html)
        self.assertIn('0', html)  # total count should be 0
        self.assertIn('Empty_species', html)

    def test_html_escaping(self):
        results = [{
            'genome_file': 'test<b>bold</b>.fasta',
            'gs_type': 'GS1.0', 'quality': 'GREEN',
            'confidence_score': 95.0, 'fragment_pattern': '1 2 3',
            'num_operons': 3, 'is_novel': False, 'qc_flags': [],
            'fragments': [], 'genome_length': 3000000,
            'validation_passed': True,
        }]
        report = HtmlReport(results, species='<em>Test</em>')
        html = report.generate()
        # Raw HTML tags must not appear unescaped
        self.assertNotIn('<b>bold</b>', html)
        self.assertIn('&lt;b&gt;bold&lt;/b&gt;', html)
        self.assertNotIn('<em>Test</em>', html)
        self.assertIn('&lt;em&gt;Test&lt;/em&gt;', html)

    def test_fragment_visualization(self):
        results = [{
            'genome_file': 'test.fasta', 'gs_type': 'GS1.0',
            'quality': 'GREEN', 'confidence_score': 95.0,
            'fragment_pattern': "1 2' ? 4", 'num_operons': 4,
            'is_novel': False, 'qc_flags': [], 'fragments': [],
            'genome_length': 3000000, 'validation_passed': True,
        }]
        report = HtmlReport(results, species='Test')
        html = report.generate()
        self.assertIn('frag-fwd', html)
        self.assertIn('frag-rev', html)
        self.assertIn('frag-unknown', html)

    def test_sortable_table_js(self):
        results = [{
            'genome_file': 'test.fasta', 'gs_type': 'GS1.0',
            'quality': 'GREEN', 'confidence_score': 95.0,
            'fragment_pattern': '1 2 3', 'num_operons': 3,
            'is_novel': False, 'qc_flags': [], 'fragments': [],
            'genome_length': 3000000, 'validation_passed': True,
        }]
        report = HtmlReport(results, species='Test')
        html = report.generate()
        self.assertIn('<script>', html)
        self.assertIn('sort', html.lower())

    def test_qc_summary_with_flags(self):
        results = [
            {'genome_file': 'a.fasta', 'gs_type': 'GS1.0', 'quality': 'AMBER',
             'confidence_score': 60.0, 'fragment_pattern': '1 ? 3',
             'num_operons': 3, 'is_novel': False,
             'qc_flags': [{'code': 'MISSING_FRAGMENT', 'severity': 'warning',
                           'message': 'frag missing'}],
             'fragments': [], 'genome_length': 3000000, 'validation_passed': True},
            {'genome_file': 'b.fasta', 'gs_type': 'GS2.0', 'quality': 'AMBER',
             'confidence_score': 55.0, 'fragment_pattern': '1 ? 3',
             'num_operons': 3, 'is_novel': False,
             'qc_flags': [{'code': 'MISSING_FRAGMENT', 'severity': 'warning',
                           'message': 'frag missing'}],
             'fragments': [], 'genome_length': 3000000, 'validation_passed': True},
        ]
        report = HtmlReport(results, species='Test')
        html = report.generate()
        self.assertIn('QC Summary', html)
        self.assertIn('MISSING_FRAGMENT', html)
        # Count should be 2
        self.assertIn('2', html)

    def test_tool_version_in_report(self):
        results = [{
            'genome_file': 'test.fasta', 'gs_type': 'GS1.0',
            'quality': 'GREEN', 'confidence_score': 95.0,
            'fragment_pattern': '1 2 3', 'num_operons': 3,
            'is_novel': False, 'qc_flags': [], 'fragments': [],
            'genome_length': 3000000, 'validation_passed': True,
        }]
        report = HtmlReport(results, species='Test', tool_version='3.0.0')
        html = report.generate()
        self.assertIn('3.0.0', html)


if __name__ == '__main__':
    unittest.main()
