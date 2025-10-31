# Work Instruction: Socru Analysis Workflow

**Document Version:** 1.0  
**Last Updated:** 2025-10-31  
**Purpose:** Standard operating procedure for analyzing bacterial genome arrangements using Socru

## Scope

This work instruction applies to the analysis of complete bacterial genome assemblies to determine genome arrangement types using the Socru bioinformatics tool.

## Prerequisites

### Required Skills
- Basic command-line proficiency
- Understanding of bacterial genome structure
- Familiarity with FASTA file formats

### Required Software
- Socru (version 2.0 or later)
- BLAST+ tools (automatically installed with conda/docker)
- barrnap (automatically installed with conda/docker)

### Required Data
- Complete bacterial genome assemblies in FASTA format
- Appropriate species database for analysis

## Safety and Quality Considerations

### Data Quality Requirements
- **Assembly completeness**: Must be complete assemblies (single contig for chromosome)
- **Assembly method**: Long-read (PacBio/Nanopore) or hybrid assemblies preferred
- **Quality metrics**: Check N50, L50, and assembly graph for circularity
- **Contamination check**: Screen for contamination before analysis

### Data Integrity
- Maintain original data files (do not modify)
- Use version control for analysis scripts
- Document all analysis steps
- Preserve intermediate output files

## Standard Operating Procedure

### Procedure 1: Initial Setup and Verification

**Objective:** Ensure Socru is correctly installed and operational

**Steps:**

1.1. Verify Socru installation
```bash
socru --version
```
Expected: Version number displayed (e.g., 2.2.4)

1.2. List available species databases
```bash
socru_species
```
Expected: List of species names

1.3. Verify required input files exist
```bash
ls -lh *.fasta
```

1.4. Create analysis directory structure
```bash
mkdir -p analysis_output
mkdir -p quality_reports
mkdir -p logs
```

**Quality Check:** All commands complete without errors

---

### Procedure 2: Pre-Analysis Quality Assessment

**Objective:** Verify assembly quality before Socru analysis

**Steps:**

2.1. Check assembly completeness
- Verify single contig for chromosome
- Confirm circular topology if expected

2.2. Document assembly metadata
```bash
# Record in analysis_log.txt
echo "Sample: sample_id" >> logs/analysis_log.txt
echo "Assembly date: $(date)" >> logs/analysis_log.txt
echo "Assembly method: [method]" >> logs/analysis_log.txt
```

2.3. Verify FASTA format
```bash
head -1 *.fasta
```
Expected: Lines starting with ">"

**Quality Check:** Assemblies meet completeness criteria

---

### Procedure 3: Single Genome Analysis

**Objective:** Analyze a single genome to determine arrangement type

**Steps:**

3.1. Identify appropriate species database
```bash
socru_species | grep -i "genus_species"
```

3.2. Run Socru analysis
```bash
socru [Species_name] genome.fasta \
    -o analysis_output/results.txt \
    -n analysis_output/novel_profiles.txt \
    -f analysis_output/novel_fragments.fa \
    -t 4 \
    -v 2>&1 | tee logs/analysis.log
```

3.3. Verify output files created
```bash
ls -lh analysis_output/
```
Expected files:
- results.txt
- operon_directions.txt
- genome_structure.pdf
- profile.txt.novel (if novel patterns found)
- novel_fragments.fa (if unclassified fragments found)

3.4. Review primary results
```bash
cat analysis_output/results.txt
```

**Quality Check:** Output files exist and contain expected data

---

### Procedure 4: Batch Analysis

**Objective:** Process multiple genomes in a single analysis run

**Steps:**

4.1. Organize input files
```bash
# All FASTA files in input_genomes/
ls input_genomes/*.fasta | wc -l
```

4.2. Run batch analysis
```bash
socru [Species_name] input_genomes/*.fasta \
    -o analysis_output/batch_results.txt \
    -n analysis_output/batch_novel_profiles.txt \
    -f analysis_output/batch_novel_fragments.fa \
    -b analysis_output/blast_results.txt \
    -t 4 \
    2>&1 | tee logs/batch_analysis.log
```

4.3. Verify all genomes processed
```bash
# Count input files
ls input_genomes/*.fasta | wc -l

# Count result lines
wc -l analysis_output/batch_results.txt
```

4.4. Extract summary statistics
```bash
# Count by GS type
cut -f2 analysis_output/batch_results.txt | sort | uniq -c
```

**Quality Check:** Number of results matches number of input files

---

### Procedure 5: Results Interpretation

**Objective:** Interpret and validate Socru results

**Steps:**

5.1. Identify result categories

For each result line:
```
filename.fasta	GS[X].[Y]	fragment_pattern
```

- **GS1.0+**: Known, validated pattern - ACCEPT
- **GS0.X**: Novel pattern - REQUIRES VALIDATION

5.2. Review known patterns
```bash
# Extract known patterns (GS1.0 and higher)
grep -E "GS[1-9]" analysis_output/results.txt > analysis_output/known_patterns.txt
```

5.3. Review novel patterns
```bash
# Extract novel patterns (GS0.X)
grep "GS0" analysis_output/results.txt > analysis_output/novel_patterns.txt
```

5.4. Check operon directions
```bash
cat operon_directions.txt
```

Verify:
- Arrows point outward from origin (Ori)
- Origin and terminus are not adjacent (unless expected)
- Pattern is biologically plausible

5.5. Review visualizations
```bash
# Open PDF files
ls -lh genome_structure.pdf
```

**Quality Check:** All patterns are interpreted and validated

---

### Procedure 6: Novel Pattern Validation

**Objective:** Validate or reject novel genome arrangement patterns

**Steps:**

6.1. List novel patterns
```bash
cat analysis_output/novel_patterns.txt
```

6.2. For each novel pattern, check:

**Assembly Quality:**
- Review assembly metrics
- Check for chimeric contigs
- Verify coverage depth
- Look for structural variants in reads

**Biological Plausibility:**
- Compare with related genomes
- Check origin/terminus placement
- Verify fragment count matches expected
- Review for duplications or deletions

6.3. Document validation decision
```bash
# Create validation report
echo "Pattern: GS0.X" >> quality_reports/validation.txt
echo "Decision: [ACCEPT/REJECT]" >> quality_reports/validation.txt
echo "Reason: [explanation]" >> quality_reports/validation.txt
echo "---" >> quality_reports/validation.txt
```

6.4. For accepted novel patterns, update database (see Procedure 8)

**Quality Check:** All novel patterns reviewed and documented

---

### Procedure 7: Quality Control Review

**Objective:** Final quality check of all results

**Steps:**

7.1. Verify completeness
- [ ] All input genomes analyzed
- [ ] All output files generated
- [ ] All novel patterns reviewed
- [ ] Documentation complete

7.2. Check for common errors

**Missing fragments:**
```bash
# Check for incomplete fragment sets
grep "?" analysis_output/results.txt
```

**Extra fragments:**
```bash
# Look for duplicated fragments
awk '{print NF-2}' analysis_output/results.txt | sort | uniq -c
```

7.3. Review log files for warnings
```bash
grep -i "warning\|error" logs/*.log
```

7.4. Generate summary report
```bash
cat > quality_reports/summary.txt <<EOF
Analysis Date: $(date)
Total Genomes: $(wc -l < analysis_output/results.txt)
Known Patterns: $(grep -c "GS[1-9]" analysis_output/results.txt)
Novel Patterns: $(grep -c "GS0" analysis_output/results.txt)
Failed: $(grep -c "?" analysis_output/results.txt)
EOF
```

**Quality Check:** No unresolved errors or warnings

---

### Procedure 8: Database Update (Advanced)

**Objective:** Update species database with validated novel patterns

**Prerequisites:**
- Novel patterns validated (Procedure 6)
- Database write access
- Backup of original database

**Steps:**

8.1. Backup current database
```bash
cp -r /path/to/species_database /path/to/species_database.backup
```

8.2. Update profile file
```bash
socru_update_profile \
    analysis_output/batch_results.txt \
    /path/to/species_database/profile.txt \
    -o updated_profile.txt
```

8.3. Review updated profile
```bash
diff /path/to/species_database/profile.txt updated_profile.txt
```

8.4. Replace profile file (if changes are correct)
```bash
cp updated_profile.txt /path/to/species_database/profile.txt
```

8.5. Test updated database
```bash
socru [Species_name] -d /path/to/species_database test_genome.fasta
```

**Quality Check:** Updated database produces expected results

---

### Procedure 9: Database Optimization (Advanced)

**Objective:** Reduce database size while maintaining accuracy

**Prerequisites:**
- BLAST results from multiple genomes (50+ recommended)
- Database write access
- Backup of original database

**Steps:**

9.1. Collect BLAST results
```bash
# Already done if -b flag used in Procedure 4
# Otherwise, re-run analysis with -b flag
```

9.2. Create optimized database
```bash
socru_shrink_database \
    analysis_output/blast_results.txt \
    /path/to/original_database \
    /path/to/optimized_database \
    --min_fragment_size 100000
```

9.3. Compare database sizes
```bash
du -sh /path/to/original_database
du -sh /path/to/optimized_database
```

9.4. Validate optimized database
```bash
# Test on known genomes
socru [Species_name] -d /path/to/optimized_database test_genomes/*.fasta
```

9.5. Compare results
```bash
diff original_results.txt optimized_results.txt
```
Expected: Identical results with 80%+ size reduction

**Quality Check:** Optimized database produces identical results

---

## Documentation Requirements

### Required Documentation

For each analysis, document:

1. **Input Data**
   - File names and locations
   - Assembly methods
   - Quality metrics

2. **Analysis Parameters**
   - Socru version
   - Species database used
   - Command-line options

3. **Results**
   - GS types identified
   - Novel patterns (if any)
   - Validation decisions

4. **Quality Issues**
   - Warnings or errors
   - Failed analyses
   - Resolution steps

### Documentation Template

Create a file `analysis_report.txt`:

```
SOCRU ANALYSIS REPORT
=====================

Project: [Project Name]
Date: [Date]
Analyst: [Name]
Socru Version: [Version]

INPUT DATA
----------
Species: [Species Name]
Number of Genomes: [N]
Assembly Method: [Method]
Source: [Source/Study]

ANALYSIS PARAMETERS
-------------------
Database: [Species Database]
Threads: [N]
Additional Options: [Options]

RESULTS SUMMARY
---------------
Known Patterns: [N]
Novel Patterns: [N]
Failed: [N]

GS Type Distribution:
- GS1.0: [N]
- GS1.1: [N]
...

NOVEL PATTERNS
--------------
[Pattern ID]: [Status - Accepted/Rejected]
Reason: [Explanation]

QUALITY ISSUES
--------------
[List any issues and resolutions]

CONCLUSIONS
-----------
[Summary of findings]

NEXT STEPS
----------
[Any follow-up actions needed]
```

---

## Troubleshooting

### Common Issues and Solutions

| Issue | Possible Cause | Solution |
|-------|---------------|----------|
| "Cannot find species database" | Typo in species name | Check spelling with `socru_species` |
| "blastn: not found" | BLAST not installed | Install NCBI BLAST+ |
| All genomes show GS0.0 | Wrong species database | Verify species match |
| Missing fragments | Draft assembly | Use only complete assemblies |
| Extra fragments | Assembly error | Validate assembly quality |
| No output files | Write permission issue | Check directory permissions |

### When to Escalate

Contact bioinformatics support if:
- Repeated analysis failures
- Unexpected results across many samples
- Database corruption suspected
- Novel biological patterns identified

---

## Revision History

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | 2025-10-31 | Documentation Team | Initial release |

---

## References

- Socru User Guide: `docs/user_guide.md`
- Socru Tutorial: `docs/tutorial.md`
- Socru GitHub: https://github.com/quadram-institute-bioscience/socru

---

## Appendix A: Quick Reference Commands

```bash
# List species
socru_species

# Single genome analysis
socru Species_name genome.fasta -o results.txt -t 4

# Batch analysis with BLAST output
socru Species_name *.fasta -o batch_results.txt -b blast_hits.txt -t 4

# Create database
socru_create output_db reference.fasta

# Update profile
socru_update_profile results.txt profile.txt -o updated_profile.txt

# Shrink database
socru_shrink_database blast_hits.txt input_db output_db
```

---

## Appendix B: File Format Specifications

### Input FASTA Format
```
>sequence_identifier optional_description
ATCGTAGCTAGCTAGCTAGC...
```

### Output Results Format (Tab-delimited)
```
filename	GS_type	fragment1	fragment2	...	fragmentN
```

### Profile Format
```
GS_number	fragment_pattern
```

Example:
```
1.0	1-2-3-4-5-6-7
1.1	1-2-3'-4-5-6-7
```

---

**END OF WORK INSTRUCTION**
