# socru

Genome structural typing via ribosomal operon arrangement in bacteria.

## Quick reference

- **Language**: Python 3.9+
- **License**: GPL-3.0
- **External deps**: barrnap (rRNA detection), BLAST+ (fragment identification)
- **Install deps**: `conda install -c bioconda barrnap blast`
- **Install package**: `pip install -e ".[dev]"`
- **Run tests**: `python -m pytest socru/tests/ -p no:asyncio -v`
- **Lint**: `ruff check socru/`
- **CI matrix**: Python 3.9, 3.10, 3.11, 3.12 on ubuntu-latest

## Architecture

**Core pipeline**: Barrnap (rRNA detection) -> Fragment extraction -> BLAST identification -> Type assignment -> Validation

**Module groups**:
- Core pipeline: Socru, SocruCreate, SocruConfig, Barrnap, Blast, Database
- Data models: AnalysisResult, Fragment, BlastResult, Operon, GATProfile
- Analysis: ConfidenceScore, QCFlags, NoveltyDetector, RearrangementDistance, BatchStats
- Visualization (pure Python SVG): SvgGenomePlot, SvgSynteny, SvgFragmentQuality, SvgTypeDistribution, SvgConfidenceHeatmap, SvgCoveragePileup
- Reporting: HtmlReport (self-contained HTML)

## Module layout

```
socru/
  Socru.py                # Main pipeline orchestrator
  SocruCreate.py          # Database creation pipeline
  SocruConfig.py          # Config dataclasses (SocruConfig, SocruCreateConfig)
  AnalysisResult.py       # Structured result model (AnalysisResult, FragmentResult, OperonResult, QCFlag)
  Barrnap.py              # rRNA operon detection via barrnap
  Blast.py                # BLAST+ wrapper (blastn, makeblastdb)
  BlastResult.py          # BLAST hit data model
  Database.py             # Species database access
  DatabaseManager.py      # Bundled + user database path resolution
  Fragment.py             # Chromosome fragment model
  FragmentFiles.py        # Fragment file I/O
  FilterBlast.py          # BLAST result filtering
  Operon.py               # Ribosomal operon model
  GATProfile.py           # Core GS-type profile data model
  TypeGenerator.py        # GS type assignment logic
  Profiles.py             # Profile collection management
  ProfileGenerator.py     # Profile generation from fragments
  Schemas.py              # Database schema definitions
  ValidateFragments.py    # Fragment validation checks
  ConfidenceScore.py      # Confidence scoring for type assignments
  QCFlags.py              # Quality control flag generation
  NoveltyDetector.py      # Novel arrangement detection
  RearrangementDistance.py # Distance metrics between arrangements
  BatchStats.py           # Batch analysis statistics
  DnaA.py                 # DnaA origin detection
  Dif.py                  # Dif site detection
  Fasta.py                # FASTA file utilities
  Results.py              # Legacy results formatting
  PlotProfile.py          # Legacy matplotlib plotting
  ShrinkDatabase.py       # Database size reduction
  Svg*.py                 # Six SVG visualization modules (pure Python)
  HtmlReport.py           # Self-contained HTML report generation
  cli.py                  # CLI entry points (socru, socru_create, etc.)
```

## Key conventions

- Context managers for all resource cleanup (Socru, Database, Blast, Barrnap)
- No `shell=True` -- all subprocess calls use argument lists
- Type hints on all core modules (py.typed marker present)
- Logging via `logging` module, not print()
- SocruConfig dataclass or argparse namespace for configuration
- SVG generation is pure Python, no external deps
- Tests use unittest, run with pytest

## Testing

- ~337 tests across 35 test files
- Integration tests require barrnap + BLAST: SocruCreate_test, Socru_test, Integration_test
- Unit-only filter: `pytest -k "not SocruCreate_test and not Socru_test and not Integration_test"`
- Coverage: `pytest --cov=socru --cov-report=term-missing socru/tests/`

## Species databases

- 433 bundled in socru/data/ (~238MB)
- User databases: ~/.socru/data/ or SOCRU_DATA_DIR env var
- DatabaseManager class handles both bundled and user paths

## CLI entry points

- `socru` -- main typing tool
- `socru_create` -- create new species database
- `socru_rebuild_profile` -- rebuild profile database
- `socru_shrink_database` -- reduce database size
- `socru_species` -- list available species
- `socru_update_profile` -- update profile database
