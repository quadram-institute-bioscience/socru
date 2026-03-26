# socru - Genome Structural Typing

## What this project does
socru classifies bacterial genome structures by analyzing the order and orientation of chromosomal fragments between ribosomal RNA operons. It assigns standardized GS type identifiers.

## Architecture
- **Core pipeline**: Barrnap (rRNA detection) -> Fragment extraction -> BLAST identification -> Type assignment -> Validation
- **Key modules**: Socru.py (orchestrator), GATProfile.py (core data model), TypeGenerator.py (type assignment)
- **New modules**: AnalysisResult, ConfidenceScore, QCFlags, NoveltyDetector, RearrangementDistance, BatchStats
- **Visualization**: SvgGenomePlot, SvgSynteny, SvgFragmentQuality, SvgTypeDistribution, SvgConfidenceHeatmap, SvgCoveragePileup
- **Reporting**: HtmlReport (self-contained HTML), JSON output

## Testing
- Run: `python -m pytest socru/tests/ -p no:asyncio -v`
- 330+ tests, all should pass when barrnap + blast are installed
- Tests requiring barrnap: SocruCreate_test, Socru_test, SocruRebuild, SocruUpdate

## External dependencies
- barrnap (rRNA prediction, requires Perl + HMMer)
- BLAST+ (blastn, makeblastdb)
- Install via conda: `conda install -c bioconda barrnap blast`

## Key conventions
- Config via SocruConfig dataclasses (or argparse namespace for backward compat)
- Context managers for resource cleanup (Socru, Database, Blast, Barrnap)
- All subprocess calls use argument lists, never shell=True
- SVG generation is pure Python (no external deps)
- Species databases in socru/data/ (433 species, ~238MB)
