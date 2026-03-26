# Socru
![GitHub Actions](https://img.shields.io/github/actions/workflow/status/quadram-institute-bioscience/socru/ci.yml?branch=master)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-brightgreen.svg)](https://github.com/quadram-institute-bioscience/socru/blob/master/LICENSE)
[![codecov](https://codecov.io/gh/andrewjpage/socru/branch/master/graph/badge.svg)](https://codecov.io/gh/andrewjpage/socru)
[![Docker Pulls](https://img.shields.io/docker/pulls/quadraminstitute/socru.svg)](https://hub.docker.com/r/quadraminstitute/socru)

## Contents
  * [Introduction](#introduction)
  * [Quick Start](#quick-start)
  * [Installation](#installation)
  * [Documentation](#documentation)
  * [Usage Examples](#usage-examples)
  * [Output Formats](#output-formats)
  * [CLI Options](#cli-options)
  * [Testing](#testing)
  * [License](#license)
  * [Feedback/Issues](#feedbackissues)
  * [Citation](#citation)

## Introduction
Socru allows you to easily identify and communicate the order and orientation of complete genomes around ribosomal operons. These large scale structural variants have real impacts on the phenotype of the organism, and with the advent of long read sequencing, we can now start to delve into the mechanisms at work.

## Quick Start

```bash
# Install via conda
conda install -c conda-forge -c bioconda socru

# List available species
socru_species

# Analyze a genome
socru Escherichia_coli genome.fasta
```

## Documentation

Comprehensive documentation is available in the [docs/](docs/) directory:

- **[Installation Guide](docs/installation.md)** - Detailed installation instructions for all platforms
- **[Tutorial](docs/tutorial.md)** - Step-by-step tutorial for new users
- **[User Guide](docs/user_guide.md)** - Complete usage guide with all features
- **[Work Instruction](docs/work_instruction.md)** - Standard operating procedures for analysis workflows
- **[API Reference](docs/api_reference.md)** - Python API documentation for developers

## Installation

**Requires Python 3.9+**, barrnap, and BLAST+.

### Conda (Recommended)

[![Anaconda-Server Badge](https://anaconda.org/bioconda/socru/badges/latest_release_date.svg)](https://anaconda.org/bioconda/socru)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/socru/badges/platforms.svg)](https://anaconda.org/bioconda/socru)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/socru/badges/downloads.svg)](https://anaconda.org/bioconda/socru)

```bash
conda install -c conda-forge -c bioconda socru
```

### Docker

```bash
docker pull quadraminstitute/socru
docker run --rm -it -v /path/to/data:/data quadraminstitute/socru socru <species> <genome.fasta>
```

### pip

```bash
pip3 install git+https://github.com/quadram-institute-bioscience/socru
```

**Note:** pip installation requires barrnap and BLAST+ to be installed separately.

For detailed installation instructions, see the [Installation Guide](docs/installation.md).

## Usage Examples

### Basic Analysis

```bash
# List available species databases
socru_species

# Analyze a single genome
socru Escherichia_coli genome.fasta

# Analyze multiple genomes
socru Salmonella_enterica *.fasta -o results.txt

# Use multiple threads
socru Klebsiella_pneumoniae genome.fasta -t 4
```

### Advanced Usage

```bash
# Save novel patterns and fragments
socru Escherichia_coli genome.fasta \
  -n novel_profiles.txt \
  -f novel_fragments.fa

# Export structured JSON results
socru Escherichia_coli genome.fasta --output_json results.json

# Generate an SVG circular genome diagram
socru Escherichia_coli genome.fasta --output_svg genome.svg

# Generate a self-contained HTML report
socru Escherichia_coli genome.fasta --output_html report.html

# Create a custom database
socru_create my_species_db reference_genome.fasta

# Update database with validated patterns
socru_update_profile results.txt profile.txt -o updated_profile.txt

# Optimize database size
socru_shrink_database blast_results.txt input_db/ output_db/
```

For complete usage information, see the [User Guide](docs/user_guide.md) and [Tutorial](docs/tutorial.md).

## Output Formats

Socru supports multiple output formats:

| Format | Flag | Description |
|--------|------|-------------|
| **Tab-delimited text** | (default) | Backward-compatible one-line-per-genome results to stdout or `-o` file |
| **JSON** | `--output_json results.json` | Structured output for programmatic consumption, including confidence scores, QC flags, and per-fragment BLAST details |
| **SVG** | `--output_svg genome.svg` | Publication-quality circular genome diagram showing fragment arrangement and operon orientations |
| **HTML** | `--output_html report.html` | Self-contained, shareable report with interactive summary tables |

All formats can be used simultaneously in a single run.

## CLI Options

| Option | Short | Description |
|--------|-------|-------------|
| `--output_file` | `-o` | Output filename (default: stdout) |
| `--output_json` | `-j` | JSON structured results output |
| `--output_svg` | `-s` | SVG circular genome diagram |
| `--output_html` | | Self-contained HTML report |
| `--output_plot_file` | `-p` | Genome structure PDF plot |
| `--threads` | `-t` | Number of threads (default: 1) |
| `--novel_profiles` | `-n` | File for novel profiles |
| `--new_fragments` | `-f` | File for novel fragments |
| `--top_blast_hits` | `-b` | File for top BLAST hits |
| `--not_circular` | `-c` | Assume chromosome is not circular |
| `--min_bit_score` | | Minimum BLAST bit score (default: 100) |
| `--min_alignment_length` | | Minimum alignment length (default: 100) |
| `--verbose` | `-v` | Enable verbose logging output |
| `--debug` | | Enable profiling output |

## Testing

Run the test suite with pytest:

```bash
# Run all tests
python3 -m pytest socru/tests/ -v

# Run with coverage
python3 -m pytest socru/tests/ --cov=socru --cov-report=term-missing
```

Current test coverage: **93%**

## Output Files

Socru generates several output files:

- **results.txt** - Tab-delimited results with GS types and fragment patterns
- **operon_directions.txt** - Direction of rRNA operons and origin location
- **genome_structure.pdf** - Visual representation of genome structure
- **profile.txt.novel** - Novel patterns requiring validation
- **novel_fragments.fa** - Unclassified fragments for investigation

## Command Reference

For detailed command-line options and usage, see:
- **socru** - Main analysis command ([User Guide](docs/user_guide.md#socru---analyze-genomes))
- **socru_species** - List available databases ([User Guide](docs/user_guide.md#socru_species---list-available-databases))
- **socru_create** - Create custom databases ([User Guide](docs/user_guide.md#socru_create---create-new-databases))
- **socru_update_profile** - Update profiles ([User Guide](docs/user_guide.md#socru_update_profile---update-database-profiles))
- **socru_shrink_database** - Optimize databases ([User Guide](docs/user_guide.md#socru_shrink_database---optimize-databases))

## Resources

- **Memory**: ~250MB for typical 5Mb bacterial genome
- **Time**: ~20 seconds per genome (single thread)
- **Threads**: Optimal with 2-4 threads

## License
Socru is free software, licensed under [GPLv3](LICENSE).

## Feedback/Issues
Please report issues or provide feedback via [GitHub Issues](https://github.com/quadram-institute-bioscience/socru/issues).

Contributions welcome! Submit improvements via [pull requests](https://github.com/quadram-institute-bioscience/socru/pulls).

## Etymology
[socru](https://www.focloir.ie/en/dictionary/ei/arrangement) (sock-roo) is the Irish (Gaeilge) word for "arrangement".

## Citation
**socru: typing of genome-level order and orientation around ribosomal operons in bacteria**
Andrew J. Page, Emma V. Ainsworth, Gemma C. Langridge
*Microbial Genomics* (2020)
https://doi.org/10.1099/mgen.0.000396

## Authors
See [AUTHORS](AUTHORS) and [CONTRIBUTORS](CONTRIBUTORS) files.
