# API Reference

This document provides reference information for developers working with Socru's Python modules.

## Core Modules

### Socru

Main analysis class for processing genome assemblies.

**Module:** `socru.Socru`

**Class:** `Socru(options)`

**Parameters:**
- `options`: Options object containing:
  - `species`: Species name for database lookup
  - `input_files`: List of input FASTA files
  - `db_dir`: Optional custom database directory
  - `threads`: Number of threads to use
  - `output_file`: Output filename (or None for stdout)
  - `novel_profiles`: Filename for novel profiles
  - `new_fragments`: Filename for novel fragments
  - `top_blast_hits`: Optional BLAST results file
  - `max_bases_from_ends`: Optional fragment trimming
  - `not_circular`: Boolean for linear chromosomes
  - `min_bit_score`: Minimum BLAST bit score
  - `min_alignment_length`: Minimum alignment length
  - `verbose`: Boolean for verbose output

**Methods:**
- `run()`: Execute analysis and generate results

**Example:**
```python
from socru.Socru import Socru

class Options:
    def __init__(self):
        self.species = 'Escherichia_coli'
        self.input_files = ['genome.fasta']
        self.db_dir = None
        self.threads = 1
        self.output_file = None
        self.verbose = False
        # ... other options

opts = Options()
s = Socru(opts)
s.run()
```

### Schemas

Database schema management and species listing.

**Module:** `socru.Schemas`

**Class:** `Schemas(verbose)`

**Parameters:**
- `verbose`: Boolean for verbose output

**Methods:**

`all_available() -> List[str]`
- Returns list of all available species databases
- Returns: List of species names

`print_all()`
- Prints all available species to stdout

`extended() -> Dict`
- Returns extended database information
- Returns: Dictionary with species metadata including:
  - Species name
  - dnaA fragment number
  - dnaA forward orientation
  - dif fragment number  
  - Reference genome

`print_extended()`
- Prints extended database information in tabular format

`database_directory(db_dir, species) -> str`
- Returns path to species database
- Parameters:
  - `db_dir`: Optional custom base directory
  - `species`: Species name
- Returns: Full path to database or None if not found

**Example:**
```python
from socru.Schemas import Schemas

s = Schemas(verbose=False)
species_list = s.all_available()
print(f"Available species: {len(species_list)}")

db_path = s.database_directory(None, 'Escherichia_coli')
print(f"Database path: {db_path}")
```

### Database

BLAST database management.

**Module:** `socru.Database`

**Class:** `Database(db_dir, verbose)`

**Parameters:**
- `db_dir`: Directory containing fragment FASTA files
- `verbose`: Boolean for verbose output

**Attributes:**
- `db_prefix`: Path to BLAST database
- `concat_fasta`: Temporary concatenated FASTA file
- `files_to_cleanup`: List of temporary files

**Methods:**
- `make_blastdb(concat_fasta)`: Creates BLAST database
- `cleanup()`: Removes temporary files

### Blast

BLAST search execution and parsing.

**Module:** `socru.Blast`

**Class:** `Blast(blast_db, query_fasta, threads, verbose, max_bases_from_ends, min_bit_score, min_alignment_length)`

**Parameters:**
- `blast_db`: Path to BLAST database
- `query_fasta`: Query FASTA file
- `threads`: Number of threads
- `verbose`: Boolean for verbose output
- `max_bases_from_ends`: Optional fragment trimming
- `min_bit_score`: Minimum bit score filter
- `min_alignment_length`: Minimum alignment length filter

**Methods:**
- `run_blast()`: Executes BLAST search
- Returns: Path to BLAST results file

### Fasta

FASTA file operations.

**Module:** `socru.Fasta`

**Class:** `Fasta(input_file, verbose)`

**Parameters:**
- `input_file`: Path to FASTA file (can be gzipped)
- `verbose`: Boolean for verbose output

**Methods:**

`get_largest_contig() -> SeqRecord`
- Returns largest contig from FASTA file
- Returns: BioPython SeqRecord object

`calc_fragment_coords(barrnap_file) -> List[Fragment]`
- Calculates fragment coordinates from barrnap output
- Parameters:
  - `barrnap_file`: Path to barrnap GFF output
- Returns: List of Fragment objects

`calc_sequence_length() -> int`
- Returns total sequence length
- Returns: Integer length

**Example:**
```python
from socru.Fasta import Fasta

fasta = Fasta('genome.fasta', verbose=False)
largest = fasta.get_largest_contig()
print(f"Largest contig: {largest.id}, length: {len(largest.seq)}")
```

### Fragment

Fragment data representation.

**Module:** `socru.Fragment`

**Class:** `Fragment(start, end, sequence)`

**Parameters:**
- `start`: Start coordinate
- `end`: End coordinate
- `sequence`: BioPython Seq object

**Attributes:**
- `start`: Fragment start position
- `end`: Fragment end position
- `sequence`: Fragment sequence
- `number`: Fragment number (set later)
- `reversed`: Boolean indicating if reversed
- `output_filename`: Output file path

**Methods:**
- `num_bases() -> int`: Returns fragment length
- `output_filename() -> str`: Returns output filename

### Barrnap

rRNA prediction using barrnap.

**Module:** `socru.Barrnap`

**Class:** `Barrnap(input_file, threads, verbose)`

**Parameters:**
- `input_file`: Input FASTA file
- `threads`: Number of threads
- `verbose`: Boolean for verbose output

**Methods:**
- `run()`: Executes barrnap
- Returns: Path to GFF output file

### DnaA and Dif

Origin and terminus identification.

**Module:** `socru.DnaA`, `socru.Dif`

**Class:** `DnaA(input_file, threads, verbose, dnaa_fasta)` / `Dif(input_file, threads, verbose, dif_fasta)`

**Parameters:**
- `input_file`: Input FASTA file
- `threads`: Number of threads
- `verbose`: Boolean for verbose output
- `dnaa_fasta`/`dif_fasta`: Reference database file

**Methods:**
- `run()`: Executes BLAST search against dnaA/dif
- Returns: Path to BLAST results

### Profiles

Profile management and matching.

**Module:** `socru.Profiles`

**Class:** `Profiles(profile_filename, verbose)`

**Parameters:**
- `profile_filename`: Path to profile.txt file
- `verbose`: Boolean for verbose output

**Methods:**

`get_profile(profile_string) -> Tuple[str, str]`
- Looks up profile string
- Parameters:
  - `profile_string`: Fragment pattern string
- Returns: Tuple of (GS_identifier, profile_string)

### Results

Result formatting and output.

**Module:** `socru.Results`

**Class:** `Results(filename, gat_profile, verbose)`

**Parameters:**
- `filename`: Input filename
- `gat_profile`: GATProfile object with results
- `verbose`: Boolean for verbose output

**Methods:**
- `construct_results_string() -> str`: Formats results as tab-delimited string

### ProfileGenerator

Profile generation for database creation.

**Module:** `socru.ProfileGenerator`

**Class:** `ProfileGenerator(fragment_order, verbose)`

**Parameters:**
- `fragment_order`: Fragment order string
- `verbose`: Boolean for verbose output

**Methods:**
- `construct_profile() -> str`: Generates profile string from fragment order

### TypeGenerator

GS type assignment.

**Module:** `socru.TypeGenerator`

**Class:** `TypeGenerator(profile_filename, profile_string, output_file, verbose)`

**Parameters:**
- `profile_filename`: Path to profile file
- `profile_string`: Profile string to match
- `output_file`: Output file path
- `verbose`: Boolean for verbose output

**Methods:**
- `print_type()`: Determines and prints GS type

## Utility Modules

### SocruCreate

Database creation from reference genome.

**Module:** `socru.SocruCreate`

**Class:** `SocruCreate(options)`

**Parameters:**
- `options`: Options object with output_directory, input_file, etc.

**Methods:**
- `run()`: Creates database

### SocruUpdateProfile

Profile update with novel patterns.

**Module:** `socru.SocruUpdateProfile`

**Class:** `SocruUpdateProfile(options)`

**Parameters:**
- `options`: Options object

**Methods:**
- `run()`: Updates profile file

### SocruShrinkDatabase

Database optimization.

**Module:** `socru.SocruShrinkDatabase`

**Class:** `SocruShrinkDatabase(options)`

**Parameters:**
- `options`: Options object

**Methods:**
- `run()`: Creates optimized database

## Data Structures

### GATProfile

Genome Arrangement Type profile result.

**Module:** `socru.GATProfile`

**Attributes:**
- `fragments`: List of fragment identifiers
- `types`: GS type information
- `orientation`: Fragment orientations

### BlastResult

BLAST result representation.

**Module:** `socru.BlastResult`

**Attributes:**
- `qseqid`: Query sequence ID
- `sseqid`: Subject sequence ID
- `pident`: Percent identity
- `length`: Alignment length
- `mismatch`: Number of mismatches
- `gapopen`: Number of gap openings
- `qstart`: Query start position
- `qend`: Query end position
- `sstart`: Subject start position
- `send`: Subject end position
- `evalue`: E-value
- `bitscore`: Bit score

## Error Handling

All modules raise standard Python exceptions:
- `FileNotFoundError`: When input files don't exist
- `ValueError`: For invalid parameters
- `RuntimeError`: For execution errors

## Testing

The test suite is located in `socru/tests/`. Run tests with:

```bash
python3 -m unittest discover -s socru/tests/ -p '*_test.py'
```

## Contributing

When contributing new modules:

1. Follow existing code style
2. Add comprehensive docstrings
3. Include unit tests
4. Update this API reference
5. Ensure test coverage > 80%

## Version Information

Check version programmatically:

```python
import pkg_resources
version = pkg_resources.get_distribution('socru').version
print(f"Socru version: {version}")
```
