# Installation Guide

This guide covers the different ways to install Socru on your system.

## System Requirements

Socru is designed to run on Linux and macOS. It will **not** run on Windows.

**Python Requirements:**
- Python >= 3.6

**External Dependencies:**
- barrnap
- NCBI BLAST+ (blastn, makeblastdb)

## Installation Methods

### Conda (Recommended)

The easiest way to install Socru is through Conda with Python 3.

First, install [Miniconda with Python 3](https://conda.io/en/latest/miniconda.html), then run:

```bash
conda install -c conda-forge -c bioconda socru
```

This will automatically install all dependencies including Python packages and external tools.

### Docker

Docker provides an isolated environment for running Socru without installing dependencies on your host system.

1. Install [Docker](https://www.docker.com/)

2. Pull the Socru container:
```bash
docker pull quadraminstitute/socru
```

3. Run Socru with Docker:
```bash
docker run --rm -it -v /path/to/data:/data quadraminstitute/socru socru [options]
```

Replace `/path/to/data` with the actual path to your data files.

### pip (Manual Installation)

If you prefer manual installation, you can use pip to install Socru and its Python dependencies:

```bash
pip3 install git+https://github.com/quadram-institute-bioscience/socru
```

**Note:** This will install Python dependencies but you must manually install:
- barrnap
- NCBI BLAST+ (blastn and makeblastdb commands)

Ensure these tools are available in your system's PATH.

## Verifying Installation

After installation, verify that Socru is working correctly:

```bash
# Check version
socru --version

# List available species databases
socru_species
```

If you see the version number and a list of species, the installation was successful.

## Troubleshooting

### Command not found

If you get a "command not found" error:
- **Conda:** Ensure your conda environment is activated
- **pip:** Check that `~/.local/bin` is in your PATH
- **Docker:** Ensure Docker is running and the container was pulled successfully

### BLAST not found

If you get errors about missing BLAST tools:
- Install NCBI BLAST+ using your system's package manager
- For Ubuntu/Debian: `sudo apt-get install ncbi-blast+`
- For macOS with Homebrew: `brew install blast`

### barrnap not found

If you get errors about barrnap:
- Install barrnap using your system's package manager or conda
- For Ubuntu/Debian: `sudo apt-get install barrnap`
- For conda: `conda install -c bioconda barrnap`
