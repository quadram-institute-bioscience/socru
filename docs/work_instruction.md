# Developer Guide

## Project Structure

```
socru/
  socru/                   # Main package
    __init__.py            # Public API exports
    cli.py                 # Entry points for all CLI commands
    Socru.py               # Main analysis runner (context manager)
    SocruConfig.py         # Config dataclasses (SocruConfig, SocruCreateConfig)
    SocruCreate.py         # Database creation
    AnalysisResult.py      # Result dataclasses (AnalysisResult, FragmentResult, etc.)
    ConfidenceScore.py     # Confidence scoring (calculate_confidence)
    QCFlags.py             # QC flag generation (generate_qc_flags)
    NoveltyDetector.py     # Novel type assessment (assess_novelty)
    RearrangementDistance.py  # Edit distance between profiles
    BatchStats.py          # Batch aggregate statistics
    DatabaseManager.py     # Database discovery and installation
    HtmlReport.py          # Self-contained HTML report generator
    SvgGenomePlot.py       # Circular genome SVG diagram
    SvgSynteny.py          # Synteny comparison SVG
    SvgFragmentQuality.py  # Fragment quality bar chart SVG
    SvgTypeDistribution.py # Type distribution chart SVG
    SvgConfidenceHeatmap.py  # Confidence heatmap SVG
    SvgCoveragePileup.py   # Coverage pileup SVG
    Barrnap.py, Blast.py, Database.py  # External tool wrappers
    Fragment.py, Operon.py, BlastResult.py, GATProfile.py  # Internal data models
    Profiles.py, ProfileGenerator.py, TypeGenerator.py     # Profile matching
    data/                  # Bundled species databases
    tests/                 # Test suite
  scripts/                 # Legacy CLI scripts (entry points now in cli.py)
  docs/                    # This documentation
  pyproject.toml           # Build config, dependencies, entry points
  Dockerfile               # Container build
```

## Running Tests

```bash
python -m pytest socru/tests/ -p no:asyncio
```

For coverage:

```bash
python -m pytest socru/tests/ -p no:asyncio --cov=socru --cov-report=term-missing
```

Install dev dependencies first: `pip install -e ".[dev]"`

## Adding a New Species Database

1. Obtain a complete reference genome FASTA for the species.
2. Run `socru_create output_dir reference.fasta`.
3. Validate by running `socru output_dir -d . test_genomes/*.fasta`.
4. Place the database directory under `socru/data/Species_name/`.
5. Run tests to verify it is discoverable.

## Adding a New Visualization

1. Create `socru/SvgNewViz.py` following the pattern in `SvgGenomePlot.py`:
   - Pure-string SVG generation, no external SVG libraries.
   - Export a `generate_*_svg(...)` function returning an SVG string.
2. Add the import and function name to `socru/__init__.py` and `__all__`.
3. Wire it into `Socru.py` if it should run automatically during analysis.
4. Add tests in `socru/tests/`.

## Code Conventions

- **Context managers**: `Socru` uses `__enter__`/`__exit__` for temp file cleanup. Always use `with Socru(config) as s:`.
- **No `shell=True`**: All subprocess calls use list-form arguments.
- **Type hints**: All public functions and methods have type annotations.
- **Logging**: Use `logging.getLogger(__name__)`, not print statements.
- **Dataclasses**: Prefer `@dataclass` for structured data (see `SocruConfig`, `AnalysisResult`).
- **Linting**: Ruff with `target-version = "py39"`, line length 120. Run `ruff check socru/`.

## CI/CD

GitHub Actions runs on push and PR:
- Lint with ruff
- Test with pytest across Python 3.9-3.12
- Build Docker image

The Dockerfile uses `condaforge/miniforge3` as the base, installs barrnap and BLAST+ via mamba, then pip-installs socru.
