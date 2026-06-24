# API Reference

All public exports from `socru.__init__`. Import with `from socru import <name>`.

## Core

### Socru

Main analysis runner. Used as a context manager.

```python
Socru(config: SocruConfig)
```

- `config`: A `SocruConfig` instance (or argparse Namespace via `SocruConfig.from_options()`).
- Call `run()` inside a `with` block.

```python
with Socru(config) as s:
    s.run()
```

### SocruConfig

Dataclass configuring an analysis run.

| Field | Type | Default | Description |
|---|---|---|---|
| `input_files` | `List[str]` | `[]` | FASTA paths to analyze |
| `species` | `str` | `""` | Species name for database lookup |
| `db_dir` | `Optional[str]` | `None` | Custom database directory |
| `output_file` | `Optional[str]` | `None` | Tab-delimited output path (None = stdout) |
| `output_json` | `Optional[str]` | `None` | JSON output path |
| `output_svg` | `Optional[str]` | `None` | SVG diagram path |
| `output_html` | `Optional[str]` | `None` | HTML report path |
| `output_dir` | `Optional[str]` | `None` | Batch output directory |
| `output_plot_file` | `str` | `"genome_structure.pdf"` | PDF plot path |
| `output_operon_directions_file` | `str` | `"operon_directions.txt"` | Operon directions path |
| `novel_profiles` | `str` | `"profile.txt.novel"` | Novel profiles output |
| `new_fragments` | `str` | `"novel_fragments.fa"` | Novel fragments FASTA |
| `top_blast_hits` | `Optional[str]` | `None` | BLAST hits output |
| `threads` | `int` | `1` | Thread count |
| `min_bit_score` | `int` | `100` | Min BLAST bit score |
| `min_alignment_length` | `int` | `100` | Min alignment length |
| `max_bases_from_ends` | `Optional[int]` | `None` | Fragment end trimming |
| `not_circular` | `bool` | `False` | Treat chromosome as linear |
| `verbose` | `bool` | `False` | Enable verbose logging |

Class method: `SocruConfig.from_options(namespace)` converts an argparse Namespace.

### SocruCreate

Create a new species database from a reference genome.

```python
SocruCreate(config: SocruCreateConfig).run()
```

### SocruCreateConfig

Dataclass configuring database creation.

| Field | Type | Default | Description |
|---|---|---|---|
| `input_file` | `str` | `""` | Reference genome FASTA |
| `output_directory` | `str` | `""` | Output database directory |
| `fragment_order` | `Optional[str]` | `None` | Custom numbering (e.g. `"1-2-3-4"`) |
| `threads` | `int` | `1` | Thread count |
| `dnaa_fasta` | `Optional[str]` | `None` | Custom dnaA FASTA |
| `dif_fasta` | `Optional[str]` | `None` | Custom dif FASTA |
| `verbose` | `bool` | `False` | Verbose logging |
| `max_bases_from_ends` | `Optional[int]` | `None` | Fragment end trimming |

## Data Models

### AnalysisResult

Complete result for one genome. Key fields:

| Field | Type | Description |
|---|---|---|
| `genome_file` | `str` | Input FASTA path |
| `gs_type` | `str` | GS type (e.g. `"GS1.0"`) |
| `quality` | `str` | `"GREEN"`, `"AMBER"`, or `"RED"` |
| `is_novel` | `bool` | True if profile not in database |
| `confidence_score` | `float` | 0-100 composite score |
| `fragment_pattern` | `str` | Tab-delimited fragment string |
| `fragments` | `List[FragmentResult]` | Per-fragment details |
| `operons` | `List[OperonResult]` | Per-operon positions |
| `qc_flags` | `List[QCFlag]` | Quality control flags |
| `validation_passed` | `bool` | Replichore rule check |
| `novelty_assessment` | `Optional[dict]` | Novelty assessment if novel |

Methods: `to_dict() -> dict`, `to_json(indent=2) -> str`

### FragmentResult

Per-fragment BLAST alignment metrics.

| Field | Type | Description |
|---|---|---|
| `number` | `int` | Fragment ID |
| `reversed` | `bool` | Reverse orientation |
| `is_dnaA` / `is_dif` | `bool` | Contains origin / terminus |
| `length` | `int` | Fragment length (bp) |
| `blast_identity` | `Optional[float]` | Percent identity |
| `blast_alignment_length` | `Optional[int]` | Alignment length |
| `blast_bit_score` | `Optional[float]` | Bit score |

### OperonResult

Per-operon position: `start: int`, `end: int`, `direction: str` (`"forward"` or `"reverse"`).

### QCFlag

Quality control flag: `code: str`, `severity: str` (`"warning"` or `"error"`), `message: str`, `details: Optional[str]`.

### Fragment, BlastResult, Operon, GATProfile

Lower-level data structures used internally. `Fragment` holds coordinates and sequence. `BlastResult` holds parsed BLAST outfmt 6 fields. `Operon` represents an rRNA operon. `GATProfile` holds the genome arrangement type profile.

### NoveltyAssessment

Result of `assess_novelty()`. Fields: `nearest_known_type: str`, `edit_distance: int`, `is_likely_real: bool`, `assessment: str` (`"likely_real"`, `"possibly_artifactual"`, `"uncertain"`), `reasoning: str`, `fragment_differences: List[str]`.

## Analysis Functions

### calculate_confidence

```python
calculate_confidence(fragment_results, quality, validation_passed, is_novel=False) -> float
```

Returns a 0-100 confidence score. Weights: identity 30%, coverage 25%, profile match 25%, validation 20%. Penalty: -10 per unknown fragment.

```python
score = calculate_confidence(result.fragments, "GREEN", True)
```

### generate_qc_flags

```python
generate_qc_flags(analysis_result, expected_fragment_count) -> List[QCFlag]
```

Inspects an `AnalysisResult` and returns flags for issues needing attention. See the QC flags table in the [User Guide](user_guide.md).

### assess_novelty

```python
assess_novelty(query_fragments, known_profiles, confidence_score, blast_identities) -> NoveltyAssessment
```

Classifies a novel profile as `"likely_real"`, `"possibly_artifactual"`, or `"uncertain"` based on rearrangement distance, confidence, and BLAST identities.

### rearrangement_distance

```python
rearrangement_distance(profile_a: List[str], profile_b: List[str]) -> dict
```

Returns `{"distance": int, "inversions": int, "translocations": int, "missing": int, "details": List[str]}`.

### pairwise_distance_matrix

```python
pairwise_distance_matrix(profiles_dict: Dict[str, List[str]]) -> Dict[Tuple[str, str], dict]
```

Computes all pairwise rearrangement distances. Returns a symmetric matrix keyed by `(type_a, type_b)` tuples.

### BatchStats

```python
stats = BatchStats(results: List[dict])
```

| Method | Returns |
|---|---|
| `type_distribution()` | `{gs_type: count}` sorted by frequency |
| `quality_summary()` | `{"GREEN": n, "AMBER": n, "RED": n}` |
| `mean_confidence()` | `Optional[float]` |
| `flag_summary()` | `{flag_code: count}` |
| `outlier_assemblies(z_threshold=2.0)` | `List[str]` of outlier genome filenames |

### DatabaseManager

```python
dm = DatabaseManager(data_dir=None)
```

| Method | Returns |
|---|---|
| `get_database_dir(species)` | `Optional[str]` path to database |
| `list_species(include_bundled=True)` | `List[str]` sorted species names |
| `install_database(source_dir, species_name=None)` | `str` installed path |
| `database_info(species)` | `Optional[dict]` with path, fragment_count, known_types |

## Visualization Functions

All SVG generators return SVG strings. No external libraries required.

| Function | Description |
|---|---|
| `generate_genome_svg(fragments, operons, genome_length, gs_type, quality, ...)` | Circular genome diagram |
| `save_genome_svg(filename, **kwargs)` | Write `generate_genome_svg` output to file |
| `generate_synteny_svg(...)` | Synteny comparison between two profiles |
| `generate_fragment_quality_svg(...)` | Per-fragment BLAST quality bar chart |
| `generate_type_distribution_svg(...)` | GS type frequency distribution chart |
| `generate_confidence_heatmap_svg(...)` | Confidence score heatmap |
| `generate_coverage_pileup_svg(...)` | Alignment coverage pileup |

## Reporting

### HtmlReport

```python
report = HtmlReport(results: List[dict], species: str, tool_version="2.2.4")
html_string = report.generate()
report.save("report.html")
```

Produces a self-contained HTML file with inline CSS/JS: summary cards, sortable results table with expandable per-assembly panels, fragment pattern visualization, and QC flag summary. Works in any modern browser with no external dependencies.
