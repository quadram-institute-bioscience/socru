# Tutorial

This walkthrough covers the main socru workflows. You need socru installed and a complete bacterial genome FASTA file.

## 1. Analyze a Single Genome

```bash
socru Escherichia_coli genome.fasta
```

Output:

```
genome.fasta    GS1.0    1    2    3    4    5    6    7
```

This tells you the genome has structure type GS1.0 with seven fragments all in forward orientation.

## 2. Interpret the Output

- **GS1.0**: a known genome structure type.
- **GS0.X**: a novel type that needs validation.
- **Fragment numbers**: sequential identifiers. A prime (`3'`) means reversed.
- **Operon directions** file (`operon_directions.txt`):

```
genome.fasta    --> 1 <-- 2 <-- 3 <-- 4(Ori) --> 5 --> 6 --> 7
```

`(Ori)` marks the fragment with the replication origin. Arrows show rRNA operon transcription direction.

## 3. Generate an HTML Report

```bash
socru Escherichia_coli genome.fasta --output_html report.html
```

Open `report.html` in a browser. The report contains summary cards (quality distribution, mean confidence), a sortable results table with expandable per-assembly details, and a QC flag summary.

## 4. Run Batch Analysis with Visualizations

```bash
socru Salmonella_enterica genomes/*.fasta \
    -o results.txt \
    -j results.json \
    --output_html batch_report.html \
    --output_dir batch_output/ \
    -t 4
```

This produces:
- `results.txt`: tab-delimited results for all genomes
- `results.json`: structured JSON with per-fragment BLAST stats and QC flags
- `batch_report.html`: interactive HTML report
- `batch_output/`: per-assembly SVG diagrams and batch statistics JSON

Review `results.txt` for a quick summary:

```bash
cut -f2 results.txt | sort | uniq -c | sort -rn
```

## 5. Create a Custom Species Database

If your species is not bundled:

```bash
socru_create my_species_db reference.fasta
```

This creates fragment files (`1.fa`, `2.fa`, ...), `profile.txt`, and `profile.txt.yml` inside `my_species_db/`.

Use the new database:

```bash
socru my_species_db -d . test_genome.fasta
```

Install it for auto-discovery:

```bash
export SOCRU_DATA_DIR=~/.socru/data
mkdir -p "$SOCRU_DATA_DIR"
cp -r my_species_db "$SOCRU_DATA_DIR/My_species"
```

## 6. Use as a Python Library

```python
from socru import Socru, SocruConfig, AnalysisResult

config = SocruConfig(
    input_files=['genome.fasta'],
    species='Escherichia_coli',
    output_json='result.json',
    output_svg='genome.svg',
    threads=2,
)
with Socru(config) as s:
    s.run()
```

For batch statistics:

```python
import json
from socru import BatchStats

with open('results.json') as f:
    results = json.load(f)

stats = BatchStats(results if isinstance(results, list) else [results])
print(stats.type_distribution())
print(stats.quality_summary())
print(stats.mean_confidence())
```

For novelty assessment:

```python
from socru import assess_novelty

assessment = assess_novelty(
    query_fragments=["1", "2'", "3", "4"],
    known_profiles=[["1", "2", "3", "4"], ["1", "3", "2", "4"]],
    confidence_score=85.0,
    blast_identities=[99.1, 97.5, 98.8, 99.3],
)
print(assessment.assessment, assessment.reasoning)
```

## Next Steps

- [User Guide](user_guide.md) for all CLI options and output format details
- [API Reference](api_reference.md) for the complete Python library API
