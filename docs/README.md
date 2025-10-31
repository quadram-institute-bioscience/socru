# Socru Documentation

Welcome to the Socru documentation. This directory contains comprehensive guides and references for using Socru.

## Documentation Structure

### Getting Started
1. **[Installation Guide](installation.md)** - How to install Socru on your system
2. **[Tutorial](tutorial.md)** - Step-by-step tutorial for new users
3. **[Quick Start](#quick-start)** - Get running in 5 minutes

### Usage Documentation
4. **[User Guide](user_guide.md)** - Comprehensive guide to all features
5. **[Work Instruction](work_instruction.md)** - Standard operating procedure for analysis workflows

### Developer Documentation
6. **[API Reference](api_reference.md)** - Python API documentation for developers

## Quick Start

The fastest way to get started with Socru:

### 1. Install via Conda (Recommended)

```bash
conda install -c conda-forge -c bioconda socru
```

### 2. List Available Species

```bash
socru_species
```

### 3. Analyze a Genome

```bash
socru Escherichia_coli genome.fasta
```

That's it! See the [Tutorial](tutorial.md) for more detailed examples.

## What is Socru?

Socru (**So**me **Cr**azy **U**nambiguous) is a tool for analyzing bacterial genome arrangements, specifically the order and orientation of fragments around ribosomal operons in complete bacterial genomes.

### Key Features

- **Rapid Analysis**: Process complete genomes in seconds
- **Multi-genome Support**: Batch process multiple assemblies
- **Extensible Databases**: Create custom databases for any species
- **Visualization**: Generate PDF visualizations of genome structure
- **Novel Pattern Detection**: Identify and validate new arrangements

### Use Cases

- **Comparative Genomics**: Compare genome arrangements across strains
- **Quality Control**: Validate complete genome assemblies
- **Evolutionary Studies**: Track large-scale structural variations
- **Database Curation**: Build and maintain species-specific databases

## Documentation Overview

### For New Users

Start here if you're new to Socru:

1. Read the [Installation Guide](installation.md) to set up Socru
2. Follow the [Tutorial](tutorial.md) with example data
3. Refer to the [User Guide](user_guide.md) for detailed usage

### For Laboratory Technicians

If you're performing routine analysis:

1. Follow the [Work Instruction](work_instruction.md) for standardized workflows
2. Refer to troubleshooting sections for common issues
3. Document results according to provided templates

### For Developers

If you're integrating Socru into pipelines or extending functionality:

1. Review the [API Reference](api_reference.md)
2. Check the test suite in `socru/tests/`
3. Follow coding standards and contribution guidelines

## System Requirements

- **Operating Systems**: Linux, macOS (not Windows)
- **Python**: 3.6 or higher
- **External Tools**: BLAST+, barrnap (installed automatically with conda)
- **Memory**: Varies by genome size (typically < 1GB)
- **Storage**: Depends on database size

## Getting Help

### Documentation

- **Installation issues**: See [Installation Guide](installation.md)
- **Usage questions**: Check [User Guide](user_guide.md)
- **Step-by-step help**: Follow [Tutorial](tutorial.md)
- **Workflow procedures**: Consult [Work Instruction](work_instruction.md)

### Community Support

- **Bug Reports**: [GitHub Issues](https://github.com/quadram-institute-bioscience/socru/issues)
- **Feature Requests**: [GitHub Issues](https://github.com/quadram-institute-bioscience/socru/issues)
- **Questions**: Open a GitHub Discussion
- **Contributions**: Submit Pull Requests

## Citation

If you use Socru in your research, please cite:

```
[Citation information to be added]
```

## License

Socru is free software, licensed under [GPLv3](../LICENSE).

## Contributing

We welcome contributions! Ways to contribute:

1. **Report Bugs**: Open an issue with details
2. **Suggest Features**: Describe new functionality
3. **Submit Code**: Create pull requests
4. **Add Databases**: Share new species databases
5. **Improve Documentation**: Fix errors or add examples

See contribution guidelines in the main repository.

## Version Information

- **Current Version**: 2.2.4
- **Last Updated**: 2025-10-31
- **Documentation Version**: 1.0

## Changelog

See [CHANGELOG](../CHANGELOG) in the main repository for version history.

## Authors

See [AUTHORS](../AUTHORS) and [CONTRIBUTORS](../CONTRIBUTORS) files.

## Additional Resources

### Related Tools

- **Barrnap**: rRNA prediction (used internally)
- **BLAST+**: Sequence alignment (used internally)
- **Circlator**: For dnaA database generation

### Publications

- [Related publications to be added]

### External Links

- **GitHub Repository**: https://github.com/quadram-institute-bioscience/socru
- **Docker Hub**: https://hub.docker.com/r/quadraminstitute/socru
- **Conda Package**: https://anaconda.org/bioconda/socru

## FAQ

### General Questions

**Q: What does "complete genome" mean?**  
A: A genome assembled into a single contig representing the chromosome, typically from long-read sequencing.

**Q: Can I use draft assemblies?**  
A: No. Socru requires complete assemblies to resolve rRNA operons, which are repetitive regions.

**Q: What species are supported?**  
A: Run `socru_species` to see bundled databases. You can create custom databases for any species.

**Q: How do I cite Socru?**  
A: Citation information is available in the main README and publications section.

### Technical Questions

**Q: How much memory does Socru need?**  
A: Typically less than 1GB for most bacterial genomes.

**Q: Can I run Socru in parallel?**  
A: Yes, use the `-t` option to specify threads. 2-4 threads is optimal.

**Q: What output formats are available?**  
A: Tab-delimited text (results), GFF (operon directions), PDF (visualization), FASTA (novel fragments).

**Q: Can I automate Socru in a pipeline?**  
A: Yes. Socru is designed for pipeline integration. See [API Reference](api_reference.md).

### Troubleshooting

**Q: I get "Cannot find species database"**  
A: Check spelling with `socru_species` and use exact name.

**Q: All my results are GS0.0**  
A: You may be using the wrong species database or draft assemblies.

**Q: I found a novel pattern. Is it real?**  
A: Validate assembly quality first. See [User Guide](user_guide.md) quality section.

## Navigation

- **[← Back to Main README](../README.md)**
- **[Installation Guide →](installation.md)**
- **[Tutorial →](tutorial.md)**
- **[User Guide →](user_guide.md)**

---

*Documentation maintained by the Socru development team*
