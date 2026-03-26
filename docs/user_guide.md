# User Guide

## Basic Usage

### Single genome

```bash
socru Escherichia_coli genome.fasta
```

### Multiple genomes

```bash
socru Salmonella_enterica *.fasta -t 4 -o results.txt
```

Gzipped FASTA files are accepted. Use `-t` to set threads (2-4 is optimal).

## Output Formats

Socru can produce several output types simultaneously.

| Flag | File | Description |
|---|---|---|
| (default) | stdout or `-o FILE` | Tab-delimited: filename, GS type, fragment order |
| `-j FILE` | JSON | Structured results with BLAST details, QC flags, confidence |
| `-s FILE` | SVG | Circular genome diagram per assembly |
| `--output_html FILE` | HTML | Self-contained report with sortable table and QC summary |
| `-p FILE` | PDF | Genome structure plot (default: `genome_structure.pdf`) |
| `-r FILE` | text | Operon directions (default: `operon_directions.txt`) |
| `--output_dir DIR` | directory | Batch output: per-assembly SVGs and stats JSON |

### Tab-delimited output

```
genome.fasta    GS1.0    1    2    3    4    5    6    7
```

Columns: input filename, GS type identifier, then one column per fragment. A prime suffix (e.g. `3'`) means the fragment is in reverse orientation.

### Operon directions

```
genome.fasta    --> 1 <-- 2' <-- 4 <-- 5 <-- 3(Ori) --> 7'
```

Arrows show rRNA operon transcription direction. `(Ori)` marks the fragment containing the replication origin.

## Understanding Results

### GS Types

- **GS1.0, GS1.1, ...**: Known, validated genome structure types in the database.
- **GS0.X**: Novel patterns not yet in the database. Require validation.

### Quality Tiers

| Tier | Meaning |
|---|---|
| GREEN | Known type, all fragments matched, validation passed |
| AMBER | Known fragment order but novel orientation, or novel pattern with good BLAST hits |
| RED | Missing fragments, failed validation, or unrecognized arrangement |

### Confidence Scores

A 0-100 composite score combining:
- Mean fragment BLAST identity (30%)
- Mean alignment coverage ratio (25%)
- Profile match level (25%)
- Validation status (20%)

Penalty: -10 points per unmatched fragment.

### QC Flags

Flags are generated automatically and included in JSON and HTML output.

| Flag Code | Severity | Trigger |
|---|---|---|
| `LOW_IDENTITY` | warning | Fragment BLAST identity < 95% |
| `SHORT_ALIGNMENT` | warning | Alignment covers < 80% of fragment |
| `MISSING_FRAGMENT` | error | Fragment returned "?" |
| `UNEXPECTED_OPERON_COUNT` | warning | Operon count differs from expected |
| `NOVEL_PATTERN` | warning | Fragment order not in database |
| `NOVEL_ORIENTATION` | warning | Order known, orientation is new |
| `INVALID_ARRANGEMENT` | error | Operon orientations violate replichore rules |
| `LOW_CONFIDENCE` | warning | Score below 50 |
| `SMALL_CHROMOSOME` | warning | Chromosome < 500,000 bp |

## Species Databases

### List available databases

```bash
socru_species
```

### Use a custom database directory

```bash
socru MySpecies genome.fasta -d /path/to/my_database
```

### Create a new species database

```bash
socru_create my_new_db reference_genome.fasta
```

This produces fragment FASTA files, `profile.txt`, and `profile.txt.yml`.

### DatabaseManager and SOCRU_DATA_DIR

User-installed databases are stored in `~/.socru/data/` by default. Override with the `SOCRU_DATA_DIR` environment variable. The search order is:

1. User data directory (`SOCRU_DATA_DIR` or `~/.socru/data/`)
2. Bundled package data

Install a custom database so it is auto-discovered:

```python
from socru import DatabaseManager
dm = DatabaseManager()
dm.install_database('/path/to/my_db', 'My_species')
print(dm.list_species())
```

## Common Workflows

### Surveillance batch

```bash
socru Salmonella_enterica batch/*.fasta \
    -o results.txt -j results.json --output_html report.html \
    --output_dir batch_viz/ -t 4
```

Review `report.html` in a browser. Check the QC summary for flagged assemblies.

### Single genome investigation

```bash
socru Escherichia_coli isolate.fasta -s isolate.svg -j isolate.json
```

Open the SVG for a circular diagram. Inspect JSON for per-fragment BLAST metrics.

### Novel type assessment

When a result shows GS0.X, check assembly quality first. Use the JSON output to inspect confidence scores and QC flags. If confidence > 80 and all BLAST identities > 95%, the novel type is likely real. Update the profile database after validation:

```bash
socru_update_profile results.txt /path/to/profile.txt -o updated_profile.txt
```
