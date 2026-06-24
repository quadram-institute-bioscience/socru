# socru

Typing of genome-level order and orientation around ribosomal operons in bacteria.

## Quick Start

```bash
pip install git+https://github.com/quadram-institute-bioscience/socru
socru_species                              # list available species
socru Escherichia_coli genome.fasta        # analyze a genome
```

## Installation

Requires **Python 3.9+** and two external tools: **barrnap** and **BLAST+**.

### conda (recommended)

```bash
conda install -c conda-forge -c bioconda socru
```

This installs socru, barrnap, and BLAST+ together.

### pip

```bash
pip install git+https://github.com/quadram-institute-bioscience/socru
```

You must install barrnap and BLAST+ separately (e.g. `conda install -c bioconda barrnap blast`).

### Docker

```bash
docker pull quadraminstitute/socru
docker run --rm -v $(pwd):/data quadraminstitute/socru Escherichia_coli /data/genome.fasta
```

## Usage

### Basic analysis

```bash
# Single genome
socru Escherichia_coli genome.fasta

# Multiple genomes
socru Salmonella_enterica *.fasta -o results.txt

# Use 4 threads
socru Escherichia_coli genome.fasta -t 4
```

### Output formats

```bash
# JSON with confidence scores, QC flags, and BLAST details
socru Escherichia_coli genome.fasta --output_json results.json

# SVG circular genome diagram
socru Escherichia_coli genome.fasta --output_svg genome.svg

# Self-contained HTML report
socru Escherichia_coli genome.fasta --output_html report.html

# Batch output directory with SVGs and stats
socru Escherichia_coli *.fasta --output_dir results/

# All formats at once
socru Escherichia_coli genome.fasta \
  -o results.txt \
  --output_json results.json \
  --output_svg genome.svg \
  --output_html report.html
```

### Other commands

```bash
# List available species databases
socru_species

# Create a custom database from a reference genome
socru_create reference.fasta

# Update a profile database with new patterns
socru_update_profile results.txt profile.txt -o updated_profile.txt

# Shrink a database by removing redundant BLAST entries
socru_shrink_database blast_results.txt input_db/ output_db/
```

## Output Formats

| Format | Flag | Content |
|--------|------|---------|
| Text | default (stdout or `-o`) | Tab-delimited, one line per genome: filename, GS type, quality, fragments |
| JSON | `--output_json FILE` | Structured results with confidence scores, QC flags, per-fragment BLAST details |
| SVG | `--output_svg FILE` | Publication-quality circular genome diagram with fragment arrangement |
| HTML | `--output_html FILE` | Self-contained report with embedded SVG and summary tables |

## Library Usage

```python
from socru import Socru, SocruConfig

config = SocruConfig(
    input_files=["genome.fasta"],
    species="Escherichia_coli",
    output_json="results.json",
    output_svg="genome.svg",
    threads=4,
)

with Socru(config) as s:
    s.run()
```

The `SocruConfig` dataclass accepts all the same parameters as the CLI. See `socru/SocruConfig.py` for the full list.

For structured access to results, use `AnalysisResult`:

```python
from socru import AnalysisResult

# AnalysisResult fields include:
#   gs_type, quality, confidence_score, is_novel,
#   fragments (list of FragmentResult), operons (list of OperonResult),
#   qc_flags (list of QCFlag), novelty_assessment
```

## New Features

- **Confidence scoring** -- numeric 0-100 confidence for each GS type assignment
- **QC flags** -- machine-readable warnings (e.g. `LOW_IDENTITY`, `SHORT_ALIGNMENT`) with severity levels
- **Novelty detection** -- identifies previously unseen fragment arrangements with detailed assessment
- **SVG genome plots** -- circular diagrams showing fragment order, orientation, and operon positions
- **Additional SVG visualizations** -- synteny plots, fragment quality, type distribution, confidence heatmaps, coverage pileups
- **HTML reports** -- self-contained single-file reports with embedded visualizations
- **JSON output** -- structured results with full BLAST details for programmatic consumption
- **Batch statistics** -- aggregate stats across multiple genomes via `--output_dir`
- **Rearrangement distance** -- quantify structural differences between genome arrangements
- **SocruConfig dataclass** -- clean typed configuration for library usage without argparse

## CLI Reference

All options for the `socru` command:

| Option | Short | Default | Description |
|--------|-------|---------|-------------|
| `species` | | (required) | Species database name |
| `input_files` | | (required) | One or more FASTA files |
| `--output_file` | `-o` | stdout | Tab-delimited results output |
| `--output_json` | `-j` | | JSON structured results |
| `--output_svg` | `-s` | | SVG circular genome diagram |
| `--output_html` | | | Self-contained HTML report |
| `--output_dir` | | | Directory for batch visualizations and stats |
| `--output_plot_file` | `-p` | `genome_structure.pdf` | PDF genome structure plot |
| `--threads` | `-t` | `1` | Number of threads |
| `--novel_profiles` | `-n` | `profile.txt.novel` | File for novel profiles |
| `--new_fragments` | `-f` | `novel_fragments.fa` | File for novel fragments |
| `--top_blast_hits` | `-b` | | File for top BLAST hits |
| `--output_operon_directions_file` | `-r` | `operon_directions.txt` | Operon direction output |
| `--not_circular` | `-c` | `false` | Treat chromosome as linear |
| `--min_bit_score` | | `100` | Minimum BLAST bit score |
| `--min_alignment_length` | | `100` | Minimum BLAST alignment length |
| `--max_bases_from_ends` | `-m` | | Use only N bases from fragment ends |
| `--db_dir` | `-d` | bundled | Custom database directory |
| `--verbose` | `-v` | `false` | Verbose logging |
| `--debug` | | `false` | Profiling output |
| `--version` | | | Print version and exit |

## Dependencies

**Python packages** (installed automatically):

- biopython >= 1.68
- PyYAML
- numpy
- matplotlib

**External tools** (must be on PATH):

- [barrnap](https://github.com/tseemann/barrnap) -- rRNA operon prediction
- [BLAST+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html) -- blastn and makeblastdb

## Testing

```bash
python -m pytest socru/tests/ -v
```

330+ tests. Some tests require barrnap and BLAST+ to be installed.

## Citation

> **socru: typing of genome-level order and orientation around ribosomal operons in bacteria**
> Andrew J. Page, Emma V. Ainsworth, Gemma C. Langridge
> *Microbial Genomics* (2020)
> https://doi.org/10.1099/mgen.0.000396

## License

GPLv3. See [LICENSE](LICENSE).
