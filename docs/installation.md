# Installation

## Requirements

- Python 3.9+
- Linux or macOS (not Windows)
- External tools: BLAST+ (blastn, makeblastdb), barrnap

## Install Socru

### pip (from GitHub)

```bash
pip install git+https://github.com/quadram-institute-bioscience/socru
```

This installs Python dependencies (biopython, numpy, matplotlib, PyYAML) but not the external tools. Install those separately (see below).

### conda (bioconda)

```bash
conda install -c conda-forge -c bioconda socru
```

This installs everything, including BLAST+ and barrnap.

### Docker

```bash
docker pull quadraminstitute/socru
docker run --rm -v /path/to/data:/data quadraminstitute/socru Escherichia_coli /data/genome.fasta
```

The Docker image is based on miniforge3 with Python 3.11, BLAST+ 2.15, and barrnap 0.9.

## Install External Dependencies

If you used pip, install BLAST+ and barrnap manually.

### BLAST+

| Platform | Command |
|---|---|
| Ubuntu/Debian | `sudo apt-get install ncbi-blast+` |
| macOS (Homebrew) | `brew install blast` |
| conda | `conda install -c bioconda blast` |

### barrnap

| Platform | Command |
|---|---|
| Ubuntu/Debian | `sudo apt-get install barrnap` |
| conda | `conda install -c bioconda barrnap` |

Both `blastn` and `barrnap` must be on your `PATH`.

## Verify Installation

```bash
socru --version
socru_species
```

The first command prints the version. The second lists all bundled species databases. If both succeed, the installation is working.

## Troubleshooting

| Problem | Fix |
|---|---|
| `socru: command not found` | Activate your conda env, or ensure `~/.local/bin` is on `PATH` |
| `blastn: command not found` | Install BLAST+ (see above) |
| `barrnap: command not found` | Install barrnap (see above) |
| Docker permission errors | Ensure Docker daemon is running; use `sudo` if needed |
