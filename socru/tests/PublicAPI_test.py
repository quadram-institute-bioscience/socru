import unittest


class TestPublicAPI(unittest.TestCase):
    def test_version_available(self):
        import socru
        self.assertIsInstance(socru.__version__, str)

    def test_core_imports(self):
        from socru import Socru, SocruConfig, SocruCreate, SocruCreateConfig

    def test_model_imports(self):
        from socru import AnalysisResult, FragmentResult, OperonResult, QCFlag
        from socru import Fragment, BlastResult, Operon, GATProfile

    def test_analysis_imports(self):
        from socru import calculate_confidence, generate_qc_flags
        from socru import assess_novelty, rearrangement_distance
        from socru import BatchStats, DatabaseManager

    def test_viz_imports(self):
        from socru import generate_genome_svg, generate_synteny_svg
        from socru import generate_fragment_quality_svg
        from socru import generate_type_distribution_svg

    def test_report_imports(self):
        from socru import HtmlReport

    def test_all_exports(self):
        import socru
        self.assertIn('Socru', socru.__all__)
        self.assertIn('AnalysisResult', socru.__all__)
        self.assertIn('generate_genome_svg', socru.__all__)
