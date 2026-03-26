"""Integration tests for novelty detection, batch stats, and batch visualizations.

These tests verify that:
- NoveltyDetector fields appear in AnalysisResult.to_dict() when novel
- BatchStats can consume AnalysisResult.to_dict() output
- Batch visualization functions accept real AnalysisResult data
"""

import json
import os
import tempfile

from socru.AnalysisResult import AnalysisResult, FragmentResult, OperonResult
from socru.BatchStats import BatchStats
from socru.NoveltyDetector import assess_novelty
from socru.SvgConfidenceHeatmap import generate_confidence_heatmap_svg
from socru.SvgFragmentQuality import generate_fragment_quality_svg
from socru.SvgSynteny import generate_synteny_svg
from socru.SvgTypeDistribution import generate_type_distribution_svg


def _make_fragment(number=1, reversed_=False, identity=99.5, length=50000):
    """Helper to create a FragmentResult with plausible values."""
    return FragmentResult(
        number=number,
        reversed=reversed_,
        is_dnaA=(number == 3),
        is_dif=(number == 1),
        length=length,
        coords=[(0, length)],
        blast_identity=identity,
        blast_alignment_length=length - 100,
        blast_bit_score=8000.0,
        blast_e_value=0.0,
        blast_mismatches=5,
        blast_subject=str(number),
    )


def _make_analysis_result(
    genome_file="assembly_1.fasta",
    gs_type="GS1.0",
    quality="GREEN",
    is_novel=False,
    confidence_score=95.0,
    fragment_numbers=None,
    novelty_assessment=None,
):
    """Helper to build an AnalysisResult with reasonable defaults."""
    if fragment_numbers is None:
        fragment_numbers = [1, 2, 3, 4]
    fragments = [_make_fragment(number=n) for n in fragment_numbers]
    operons = [
        OperonResult(start=i * 100000, end=i * 100000 + 5000, direction="forward")
        for i in range(len(fragment_numbers))
    ]
    pattern = " ".join(str(n) for n in fragment_numbers)
    result = AnalysisResult(
        genome_file=genome_file,
        genome_length=5000000,
        is_circular=True,
        num_operons=len(fragment_numbers),
        gs_type=gs_type,
        quality=quality,
        is_novel=is_novel,
        fragment_pattern=pattern,
        orientation_binary=0,
        confidence_score=confidence_score,
        fragments=fragments,
        operons=operons,
        validation_passed=True,
        operon_direction_string="Valid\tforward forward forward forward",
        novelty_assessment=novelty_assessment,
    )
    result.qc_flags = []
    return result


# ---------------------------------------------------------------------------
# TASK 1: Novelty assessment appears in AnalysisResult.to_dict()
# ---------------------------------------------------------------------------


class TestNoveltyInAnalysisResult:
    """Verify novelty_assessment field round-trips through to_dict()."""

    def test_novelty_assessment_none_by_default(self):
        result = _make_analysis_result()
        d = result.to_dict()
        assert "novelty_assessment" in d
        assert d["novelty_assessment"] is None

    def test_novelty_assessment_present_when_novel(self):
        from dataclasses import asdict

        assessment = assess_novelty(
            query_fragments=["1", "3'", "2", "4"],
            known_profiles=[["1", "2", "3", "4"]],
            confidence_score=85.0,
            blast_identities=[99.0, 98.5, 97.0, 99.2],
        )
        assessment_dict = asdict(assessment)

        result = _make_analysis_result(
            is_novel=True,
            quality="AMBER",
            novelty_assessment=assessment_dict,
        )
        d = result.to_dict()

        assert d["novelty_assessment"] is not None
        assert "nearest_known_type" in d["novelty_assessment"]
        assert "edit_distance" in d["novelty_assessment"]
        assert "assessment" in d["novelty_assessment"]
        assert "reasoning" in d["novelty_assessment"]
        assert "fragment_differences" in d["novelty_assessment"]
        assert isinstance(d["novelty_assessment"]["is_likely_real"], bool)

    def test_novelty_assessment_json_serializable(self):
        from dataclasses import asdict

        assessment = assess_novelty(
            query_fragments=["1", "2'", "3"],
            known_profiles=[["1", "2", "3"]],
            confidence_score=90.0,
            blast_identities=[99.0, 98.0, 97.5],
        )
        result = _make_analysis_result(
            is_novel=True,
            novelty_assessment=asdict(assessment),
        )
        # Should not raise
        json_str = json.dumps(result.to_dict())
        parsed = json.loads(json_str)
        assert parsed["novelty_assessment"]["edit_distance"] >= 0


# ---------------------------------------------------------------------------
# TASK 2: BatchStats consumes AnalysisResult.to_dict() output
# ---------------------------------------------------------------------------


class TestBatchStatsWithAnalysisResult:
    """Verify BatchStats works with real AnalysisResult.to_dict() data."""

    def _make_batch(self):
        results = [
            _make_analysis_result(
                genome_file="asm1.fa", gs_type="GS1.0", quality="GREEN",
                confidence_score=95.0,
            ),
            _make_analysis_result(
                genome_file="asm2.fa", gs_type="GS1.0", quality="GREEN",
                confidence_score=90.0,
            ),
            _make_analysis_result(
                genome_file="asm3.fa", gs_type="GS2.1", quality="AMBER",
                confidence_score=70.0, is_novel=True,
            ),
        ]
        return [r.to_dict() for r in results]

    def test_type_distribution(self):
        dicts = self._make_batch()
        stats = BatchStats(dicts)
        td = stats.type_distribution()
        assert td["GS1.0"] == 2
        assert td["GS2.1"] == 1

    def test_quality_summary(self):
        dicts = self._make_batch()
        stats = BatchStats(dicts)
        qs = stats.quality_summary()
        assert qs["GREEN"] == 2
        assert qs["AMBER"] == 1
        assert qs["RED"] == 0

    def test_mean_confidence(self):
        dicts = self._make_batch()
        stats = BatchStats(dicts)
        mean = stats.mean_confidence()
        assert mean is not None
        assert abs(mean - 85.0) < 0.01

    def test_flag_summary_empty_flags(self):
        dicts = self._make_batch()
        stats = BatchStats(dicts)
        fs = stats.flag_summary()
        # Our test results have no QC flags
        assert isinstance(fs, dict)

    def test_outlier_assemblies(self):
        dicts = self._make_batch()
        stats = BatchStats(dicts)
        outliers = stats.outlier_assemblies()
        assert isinstance(outliers, list)


# ---------------------------------------------------------------------------
# TASK 3: Batch visualization functions accept real AnalysisResult data
# ---------------------------------------------------------------------------


def _build_assemblies_from_results():
    """Build assembly dicts suitable for SVG generators from AnalysisResults."""
    results = [
        _make_analysis_result(
            genome_file="sample_A.fasta", gs_type="GS1.0", quality="GREEN",
            fragment_numbers=[1, 2, 3, 4],
        ),
        _make_analysis_result(
            genome_file="sample_B.fasta", gs_type="GS2.1", quality="AMBER",
            fragment_numbers=[1, 3, 2, 4], is_novel=True,
        ),
    ]
    assemblies = []
    for r in results:
        rd = r.to_dict()
        frags = []
        for f in rd["fragments"]:
            fnum = f["number"]
            try:
                fnum = int(fnum)
            except (ValueError, TypeError):
                fnum = 0
            frags.append({
                "number": fnum,
                "reversed": f["reversed"],
                "length": f["length"],
                "blast_identity": f.get("blast_identity"),
                "blast_alignment_length": f.get("blast_alignment_length"),
                "is_dnaA": f.get("is_dnaA", False),
                "is_dif": f.get("is_dif", False),
            })
        assemblies.append({
            "name": os.path.basename(rd["genome_file"]),
            "gs_type": rd["gs_type"],
            "quality": rd["quality"],
            "fragments": frags,
        })
    return assemblies


class TestBatchVisualizationsWithRealData:
    """Verify SVG generators produce valid output from AnalysisResult data."""

    def test_type_distribution_svg(self):
        results = [
            _make_analysis_result(gs_type="GS1.0"),
            _make_analysis_result(gs_type="GS1.0"),
            _make_analysis_result(gs_type="GS2.1"),
        ]
        dicts = [r.to_dict() for r in results]
        stats = BatchStats(dicts)
        svg = generate_type_distribution_svg(stats.type_distribution())
        assert svg.startswith("<svg")
        assert "GS1.0" in svg
        assert "</svg>" in svg

    def test_confidence_heatmap_svg(self):
        assemblies = _build_assemblies_from_results()
        svg = generate_confidence_heatmap_svg(assemblies)
        assert svg.startswith("<svg")
        assert "Fragment Identity Heatmap" in svg
        assert "</svg>" in svg

    def test_synteny_svg(self):
        assemblies = _build_assemblies_from_results()
        svg = generate_synteny_svg(assemblies)
        assert svg.startswith("<svg")
        assert "</svg>" in svg

    def test_fragment_quality_svg(self):
        assemblies = _build_assemblies_from_results()
        for asm in assemblies:
            svg = generate_fragment_quality_svg(
                asm["fragments"], genome_name=asm["name"],
            )
            assert svg.startswith("<svg")
            assert "</svg>" in svg

    def test_write_batch_outputs_to_directory(self):
        """End-to-end: write all batch outputs to a temp directory."""
        assemblies = _build_assemblies_from_results()
        results = [
            _make_analysis_result(genome_file="a.fa", gs_type="GS1.0"),
            _make_analysis_result(genome_file="b.fa", gs_type="GS2.1"),
        ]
        dicts = [r.to_dict() for r in results]
        stats = BatchStats(dicts)

        with tempfile.TemporaryDirectory() as tmpdir:
            # Type distribution
            path = os.path.join(tmpdir, "type_distribution.svg")
            with open(path, "w") as fh:
                fh.write(generate_type_distribution_svg(stats.type_distribution()))
            assert os.path.exists(path)

            # Confidence heatmap
            path = os.path.join(tmpdir, "confidence_heatmap.svg")
            with open(path, "w") as fh:
                fh.write(generate_confidence_heatmap_svg(assemblies))
            assert os.path.exists(path)

            # Synteny
            path = os.path.join(tmpdir, "synteny.svg")
            with open(path, "w") as fh:
                fh.write(generate_synteny_svg(assemblies))
            assert os.path.exists(path)

            # Fragment quality per assembly
            for asm in assemblies:
                name = asm["name"].replace("/", "_").replace(" ", "_")
                path = os.path.join(tmpdir, f"fragment_quality_{name}.svg")
                with open(path, "w") as fh:
                    fh.write(generate_fragment_quality_svg(
                        asm["fragments"], genome_name=asm["name"],
                    ))
                assert os.path.exists(path)

            # Batch stats JSON
            path = os.path.join(tmpdir, "batch_stats.json")
            batch_stats_dict = {
                "type_distribution": stats.type_distribution(),
                "quality_summary": stats.quality_summary(),
                "mean_confidence": stats.mean_confidence(),
                "flag_summary": stats.flag_summary(),
                "outlier_assemblies": stats.outlier_assemblies(),
                "total_assemblies": len(dicts),
            }
            with open(path, "w") as fh:
                fh.write(json.dumps(batch_stats_dict, indent=2))
            assert os.path.exists(path)

            with open(path) as fh:
                loaded = json.load(fh)
            assert "type_distribution" in loaded
            assert loaded["total_assemblies"] == 2
