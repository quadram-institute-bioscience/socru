import unittest


class TestPublicAPI(unittest.TestCase):
    def test_version_available(self):
        import socru
        self.assertIsInstance(socru.__version__, str)

    def test_core_imports(self):
        pass

    def test_model_imports(self):
        pass

    def test_analysis_imports(self):
        pass

    def test_viz_imports(self):
        pass

    def test_report_imports(self):
        pass

    def test_all_exports(self):
        import socru
        self.assertIn('Socru', socru.__all__)
        self.assertIn('AnalysisResult', socru.__all__)
        self.assertIn('generate_genome_svg', socru.__all__)
