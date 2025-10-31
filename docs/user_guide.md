# User Guide

This comprehensive guide covers how to use Socru to analyze bacterial genome arrangements.

## Overview

Socru identifies and communicates the order and orientation of complete bacterial genomes around ribosomal operons. These large-scale structural variants can impact organism phenotype, and with long-read sequencing, we can now investigate these mechanisms.

## Quick Start

For a quick introduction, see the [Tutorial](tutorial.md).

## Core Concepts

### Genome Structure (GS) Types

Socru assigns a **GS identifier** to each genome arrangement pattern. The format is:
- **GS1.0, GS1.1, etc.**: Known, validated patterns
- **GS0.X**: Novel patterns that need validation

### Fragments

Complete genomes are divided into fragments between ribosomal operons. Each fragment is numbered sequentially (1, 2, 3, ...).

### Orientation

Fragments can be in forward or reverse orientation:
- **No suffix**: Forward orientation (same as reference)
- **Prime (')**: Reverse orientation (e.g., `3'` means fragment 3 is reversed)

## Main Commands

### socru - Analyze Genomes

The main command for analyzing complete genome assemblies.

**Basic Usage:**
```bash
socru <species> <assembly.fasta>
```

**Example:**
```bash
socru Escherichia_coli K12.fasta
```

**Output:**
```
K12.fasta	GS1.0	1	2	3	4	5	6	7
```

This shows:
- Input filename
- GS identifier (GS1.0)
- Fragment order and orientation (1 through 7, all forward)

#### Command Options

**Required Arguments:**
- `species`: Species name (use `socru_species` to see available databases)
- `input_files`: One or more FASTA files (can be gzipped)

**Optional Arguments:**

`--db_dir <path>` or `-d <path>`
- Custom database location
- Default: Uses bundled databases

`--threads <n>` or `-t <n>`
- Number of threads to use
- Default: 1
- Recommended: 4 or fewer (diminishing returns after this)

`--output_file <file>` or `-o <file>`
- Output filename
- Default: Print to stdout
- Note: Appends to file if it exists

`--novel_profiles <file>` or `-n <file>`
- File for novel profile patterns
- Default: `profile.txt.novel`

`--new_fragments <file>` or `-f <file>`
- File for unclassified fragments
- Default: `novel_fragments.fa`

`--top_blast_hits <file>` or `-b <file>`
- Output BLAST results (outfmt 6)
- Useful for `socru_shrink_database`

`--not_circular` or `-c`
- Treat chromosome as linear (not circular)
- Use for incomplete assemblies or linear chromosomes

`--min_bit_score <value>`
- Minimum BLAST bit score
- Default: 100

`--min_alignment_length <value>`
- Minimum BLAST alignment length
- Default: 100

`--verbose` or `-v`
- Enable debug output

`--version`
- Show version and exit

#### Understanding Output

**Primary Output (Tab-delimited):**
```
filename.fna	GS1.0	1	2	3	4	5
```

Columns:
1. Input filename
2. GS identifier
3+ Fragment order and orientation

**Operon Directions File (`operon_directions.txt`):**
```
Ty2_1.fa	--> 1 <-- 2' <-- 4 <-- 5 <-- 3(Ori) --> 7'
```

Shows:
- Direction of rRNA operons (arrows)
- Fragment with origin of replication: `(Ori)`
- Fragment orientations

**Genome Structure Visualization:**
- File: `genome_structure.pdf`
- Visual representation of circular genome structure

### socru_species - List Available Databases

Lists all species databases bundled with Socru.

**Usage:**
```bash
socru_species
```

**Example Output:**
```
Acinetobacter_baumannii
Enterobacter_cloacae
Enterococcus_faecium
Klebsiella_pneumoniae
Salmonella_enterica
Staphylococcus_aureus
```

Copy and paste the species name to use with `socru`.

### socru_create - Create New Databases

Create a custom database from a complete reference genome.

**Usage:**
```bash
socru_create [options] <output_directory> <assembly.fasta>
```

**Required Arguments:**
- `output_directory`: Directory for new database (must not exist)
- `assembly.fasta`: Complete reference genome (can be gzipped)

**Optional Arguments:**

`--threads <n>` or `-t <n>`
- Number of threads
- Default: 1

`--fragment_order <order>` or `-f <order>`
- Force specific fragment numbering
- Example: `1-2-3-4-5-6-7`
- Not recommended for general use

`--dnaa_fasta <file>` or `-d <file>`
- Custom dnaA FASTA file
- Default: Uses bundled dnaA database

`--verbose` or `-v`
- Enable debug output

**Output:**
- Fragment FASTA files (1.fa, 2.fa, 3.fa, ...)
- `profile.txt`: Pattern definitions
- `profile.txt.yml`: Metadata (dnaA location, orientation)

### socru_update_profile - Update Database Profiles

Add validated novel patterns to a database.

**Usage:**
```bash
socru_update_profile [options] <socru_output.txt> <profile.txt>
```

**Arguments:**
- `socru_output.txt`: Output from socru containing novel patterns
- `profile.txt`: Current profile file from database

**Optional Arguments:**

`--output_file <file>` or `-o <file>`
- Output filename
- Default: `updated_profile.txt`

**Workflow:**
1. Run socru on multiple genomes
2. Review novel patterns (GS0.X)
3. Validate they're real (not assembly errors)
4. Run socru_update_profile
5. Replace database profile.txt with updated version

### socru_shrink_database - Optimize Databases

Reduce database size by keeping only conserved regions.

**Usage:**
```bash
socru_shrink_database [options] <blast_results> <input_db> <output_db>
```

**Arguments:**
- `blast_results`: BLAST results from `socru -b`
- `input_db`: Directory of database to shrink
- `output_db`: Output directory (must not exist)

**Optional Arguments:**

`--min_fragment_size <bases>`
- Minimum fragment size in bases
- Default: 100000
- Smaller fragments copied in full

**Typical Workflow:**
1. Run socru with `-b` option on multiple assemblies
2. Run socru_shrink_database
3. Test new database
4. Can reduce storage by 80%+

## Quality Considerations

### Complete vs. Draft Assemblies

**Socru requires complete assemblies:**
- One contig for chromosome
- Can resolve large repeats (rRNA regions)
- Long-read sequencing (PacBio, Nanopore) or hybrid assemblies

**Do not use:**
- Short-read draft assemblies
- Scaffolded assemblies (unless truly complete)
- Assemblies with gaps in rRNA regions

### Identifying Assembly Errors

Novel patterns may indicate:
1. **Genuine biological variation** (exciting!)
2. **Assembly errors** (unfortunately common)

**Red flags:**
- Missing fragments
- Extra fragments (duplications)
- Origin and terminus adjacent (unbalanced replichores)
- Biologically improbable arrangements

**Example of likely assembly error:**
```
Salmonella_enterica.fa	GS1.0	1	2	7	4	5	6	3
```
Origin on fragment 3, adjacent to fragment 1 (circular) = improbable

### Validation Steps

For novel patterns:
1. Check assembly quality metrics
2. Look for chimeric contigs
3. Verify with independent methods if possible
4. Compare with closely related genomes
5. Consider sequencing depth and coverage

## Best Practices

1. **Use appropriate databases**: Match species to database
2. **Validate novel patterns**: Don't automatically trust GS0.X patterns
3. **Check quality**: Review warnings in output
4. **Use sufficient threads**: 2-4 threads optimal
5. **Keep databases updated**: Submit novel validated patterns upstream
6. **Document findings**: Note any unusual patterns

## Common Issues

### Missing Fragments

Small fragments may not be detected by assemblers. Check:
- Fragment size in database
- Assembly parameters
- Sequencing coverage

### Extra Fragments

Usually indicates:
- Assembly error
- Contamination
- True biological duplication (rare)

### Wrong Fragment Count

Single rRNA segment detected when multiple expected:
- Likely scaffolded short-read assembly
- Cannot resolve rRNA repeats
- Not a true complete assembly

## Getting Help

- Report issues: [GitHub Issues](https://github.com/quadram-institute-bioscience/socru/issues)
- Pull requests welcome for improvements
- Share new databases with the community
