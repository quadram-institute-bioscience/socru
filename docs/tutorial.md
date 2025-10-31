# Tutorial: Getting Started with Socru

This tutorial will walk you through using Socru to analyze bacterial genome structures.

## Prerequisites

- Socru installed (see [Installation Guide](installation.md))
- A complete bacterial genome assembly in FASTA format
- Basic command-line knowledge

## Step 1: Check Installation

Verify Socru is installed correctly:

```bash
# Check version
socru --version

# List available species
socru_species
```

You should see:
- Version number (e.g., `2.2.4`)
- List of available species databases

## Step 2: List Available Species Databases

See what species databases are available:

```bash
socru_species
```

Example output:
```
Acinetobacter_baumannii
Campylobacter_jejuni
Enterobacter_cloacae
Enterococcus_faecium
Escherichia_coli
Klebsiella_pneumoniae
Salmonella_enterica
Staphylococcus_aureus
```

## Step 3: Analyze a Genome

Let's analyze an *Escherichia coli* genome:

```bash
socru Escherichia_coli genome.fasta
```

### Example Output

```
genome.fasta	GS1.0	1	2	3	4	5	6	7
```

This tells us:
- **genome.fasta**: Input file
- **GS1.0**: Genome Structure type 1.0 (known pattern)
- **1 2 3 4 5 6 7**: Seven fragments in forward orientation

### With Multiple Genomes

Process multiple genomes at once:

```bash
socru Escherichia_coli genome1.fasta genome2.fasta genome3.fasta
```

Or using wildcards:

```bash
socru Escherichia_coli *.fasta
```

### Compressed Input

Socru accepts gzipped files:

```bash
socru Escherichia_coli genome.fasta.gz
```

## Step 4: Save Results to File

By default, output goes to stdout. Save to a file:

```bash
socru Escherichia_coli genome.fasta -o results.txt
```

Process a batch and save all results:

```bash
socru Escherichia_coli *.fasta -o all_results.txt
```

## Step 5: Understanding Advanced Output

Socru creates several output files:

### Operon Directions File

**File:** `operon_directions.txt`

Shows direction of rRNA operons and origin location:

```
genome.fasta	--> 1 <-- 2 <-- 3 <-- 4(Ori) --> 5 --> 6 --> 7
```

- `-->`: 16S comes first in operon
- `<--`: 16S comes last in operon  
- `(Ori)`: Fragment containing origin of replication
- `'`: Prime indicates reversed fragment

### Genome Structure Visualization

**File:** `genome_structure.pdf`

Visual representation of the circular genome structure showing:
- Fragment arrangement
- Operon directions
- Origin of replication

## Step 6: Working with Novel Patterns

When Socru encounters an unknown pattern, it reports a novel GS:

```
novel_genome.fasta	GS0.1	1	2	4'	3	5	6	7
```

**GS0.1** indicates a novel pattern that needs validation.

### Novel Pattern Files

**File:** `profile.txt.novel`

Contains novel patterns found:
```
0.1	1-2-4'-3-5-6-7
```

**File:** `novel_fragments.fa`

Contains unclassifiable fragments for investigation.

### Validating Novel Patterns

Before accepting a novel pattern:

1. **Check assembly quality**
   - Is it a true complete assembly?
   - Any gaps or errors?

2. **Biological plausibility**
   - Does the arrangement make sense?
   - Is origin/terminus placement reasonable?

3. **Compare with related genomes**
   - Is this pattern seen in related strains?

## Step 7: Using Multiple Threads

Speed up analysis with multiple threads:

```bash
socru Escherichia_coli genome.fasta -t 4
```

**Note:** Diminishing returns after 4 threads.

## Step 8: Working with BLAST Results

Save BLAST hits for database optimization:

```bash
socru Escherichia_coli *.fasta -b blast_results.txt
```

Use these results to shrink databases (see Step 10).

## Step 9: Creating a Custom Database

If your species isn't in the bundled databases:

```bash
socru_create my_species_db reference_genome.fasta
```

This creates:
- Fragment FASTA files (1.fa, 2.fa, ...)
- `profile.txt`: Initial pattern (GS1.0)
- `profile.txt.yml`: Metadata

**Compress fragments for production use:**
```bash
cd my_species_db
gzip *.fa
```

## Step 10: Optimizing a Database

Reduce database size using conserved regions:

### Step 10.1: Collect BLAST Results

Run Socru on many genomes with `-b` option:

```bash
socru Escherichia_coli genome*.fasta -b blast_hits.txt
```

### Step 10.2: Shrink Database

```bash
socru_shrink_database blast_hits.txt original_db/ optimized_db/
```

This can reduce storage by 80%+ while maintaining accuracy.

### Step 10.3: Test Optimized Database

```bash
socru Escherichia_coli -d optimized_db/ genome.fasta
```

Verify results match original database.

## Example Workflows

### Workflow 1: Routine Analysis

Analyzing a set of genomes:

```bash
# List species
socru_species

# Run analysis
socru Salmonella_enterica genomes/*.fasta -o results.txt -t 4

# Review results
cat results.txt
less operon_directions.txt
```

### Workflow 2: Creating and Populating a Database

Setting up a new species database:

```bash
# Create initial database
socru_create species_db reference.fasta

# Analyze test genomes
socru species_db -d . test_genomes/*.fasta -o test_results.txt

# Update with validated novel patterns
socru_update_profile test_results.txt species_db/profile.txt -o updated_profile.txt

# Replace old profile
cp updated_profile.txt species_db/profile.txt

# Compress fragment files
cd species_db
gzip *.fa
cd ..
```

### Workflow 3: Database Optimization

Optimizing an existing database:

```bash
# Collect BLAST results from many genomes
socru Klebsiella_pneumoniae large_dataset/*.fasta -b all_hits.txt -t 4

# Create optimized database
socru_shrink_database all_hits.txt original_kp_db/ optimized_kp_db/

# Test optimized database
socru Klebsiella_pneumoniae -d optimized_kp_db/ test.fasta

# If good, use optimized version
```

## Tips and Tricks

1. **Use tab completion**: Most shells support tab completion for filenames
2. **Batch processing**: Use shell loops for large datasets
3. **Check log files**: Use `-v` verbose mode for debugging
4. **Quality first**: Only use complete, high-quality assemblies
5. **Contribute back**: Share new databases and novel patterns with the community

## Troubleshooting

### "Cannot find species database"

Check spelling and use exact name from `socru_species`:
```bash
socru_species | grep -i salmonella
```

### "blastn: command not found"

Install NCBI BLAST+:
```bash
# Ubuntu/Debian
sudo apt-get install ncbi-blast+

# Or use conda
conda install -c bioconda blast
```

### Results look wrong

Check:
1. Is it a complete assembly?
2. Are you using the correct species database?
3. Is the assembly quality high?

### Novel patterns everywhere

Common causes:
- Using draft assemblies (not complete)
- Wrong species database
- Poor assembly quality
- Genuine novel biology (exciting!)

## Next Steps

Now that you're familiar with the basics:

1. Read the [User Guide](user_guide.md) for detailed information
2. Check the [API Reference](api_reference.md) if developing scripts
3. Review quality considerations in the User Guide
4. Join the community and share your findings

## Getting Help

- **Documentation**: Read this guide and the User Guide
- **Issues**: Report bugs on [GitHub Issues](https://github.com/quadram-institute-bioscience/socru/issues)
- **Questions**: Open a discussion on GitHub
- **Contributions**: Pull requests welcome!
