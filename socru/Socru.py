"""
Main Socru analysis module for genome structure typing.

This module orchestrates the complete Socru workflow for analyzing bacterial genomes
to determine their genome structure (GS) types based on the arrangement of fragments
around ribosomal operons. It coordinates rRNA identification, fragment extraction,
BLAST searching, and profile matching.

Classes:
    Socru: Main class that runs the complete analysis pipeline
"""

import json
import logging
import os

logger = logging.getLogger(__name__)
import shutil
from tempfile import mkdtemp, mkstemp

from socru.AnalysisResult import AnalysisResult, FragmentResult, OperonResult
from socru.Barrnap import Barrnap
from socru.BatchStats import BatchStats
from socru.Blast import Blast
from socru.ConfidenceScore import calculate_confidence
from socru.Database import Database
from socru.Fasta import Fasta
from socru.FilterBlast import FilterBlast
from socru.FragmentFiles import FragmentFiles
from socru.GATProfile import GATProfile
from socru.HtmlReport import HtmlReport
from socru.NoveltyDetector import assess_novelty
from socru.PlotProfile import PlotProfile
from socru.Profiles import Profiles
from socru.QCFlags import generate_qc_flags
from socru.Schemas import Schemas
from socru.SocruConfig import SocruConfig
from socru.ToolCheck import check_all_tools
from socru.SvgConfidenceHeatmap import generate_confidence_heatmap_svg
from socru.SvgFragmentQuality import generate_fragment_quality_svg
from socru.SvgSynteny import generate_synteny_svg
from socru.SvgTypeDistribution import generate_type_distribution_svg
from socru.TypeGenerator import TypeGenerator
from socru.ValidateFragments import ValidateFragments


class Socru:
    """
    Main class for running Socru genome structure analysis.

    This class manages the complete workflow: loading databases, finding rRNA operons,
    extracting inter-operon fragments, BLASTing against the fragment database, matching
    to known profiles, and generating output including plots and novel fragment reports.

    Attributes:
        input_files (list): List of input FASTA files to analyze
        min_bit_score (int): Minimum BLAST bit score for fragment matching
        min_alignment_length (int): Minimum BLAST alignment length
        threads (int): Number of threads for parallel processing
        output_file (str): Path to main output file
        novel_profiles (str): Path to file for novel profile patterns
        new_fragments (str): Path to file for unmatched fragments
        max_bases_from_ends (int): If set, only use this many bases from fragment ends
        top_blast_hits (str): Path to file for detailed BLAST results
        output_plot_file (str): Path to PDF plot of genome structure
        output_operon_directions_file (str): Path to operon direction results
        verbose (bool): Enable verbose output
        db_dir (str): Path to species database directory
        is_circular (bool): Whether chromosomes are assumed circular
        dirs_to_cleanup (list): Temporary directories to remove on cleanup
        top_results (list): Collection of top BLAST hits for all fragments
    """
    def __init__(self, options_or_config):
        """
        Initialize Socru with a SocruConfig or legacy argparse options.

        Args:
            options_or_config: A ``SocruConfig`` instance for library use,
                or an argparse ``Namespace`` for backward compatibility.
        """
        # Initialize cleanup list early so __del__ never fails
        self.dirs_to_cleanup = []

        check_all_tools()

        if isinstance(options_or_config, SocruConfig):
            config = options_or_config
        else:
            config = SocruConfig.from_options(options_or_config)

        self.input_files = config.input_files

        # BLAST filtering parameters
        self.min_bit_score = config.min_bit_score
        self.min_alignment_length = config.min_alignment_length
        self.threads = config.threads

        # Output file paths
        self.output_file = config.output_file
        self.novel_profiles = config.novel_profiles
        self.new_fragments = config.new_fragments
        self.max_bases_from_ends = config.max_bases_from_ends
        self.top_blast_hits = config.top_blast_hits
        self.output_plot_file = config.output_plot_file
        self.output_operon_directions_file = config.output_operon_directions_file

        self.verbose = config.verbose
        self.output_svg = config.output_svg
        self.output_json = config.output_json
        self.top_results = []
        self.analysis_results = []
        self.html_results = []
        self.output_html = config.output_html
        self.output_dir = config.output_dir

        # Locate and validate the species database
        self.db_dir = Schemas(self.verbose).database_directory(config.db_dir, config.species)
        if self.db_dir is None:
             raise FileNotFoundError(
                 "Cannot access the database for species '{}'. "
                 "Please check the species name and database directory.".format(config.species))

        # Set chromosome circularity assumption
        self.is_circular = not config.not_circular

    def run(self):
        """
        Execute the main Socru analysis workflow for all input files.

        This method:
        1. Loads the profile database and fragment database
        2. Iterates through each input genome file
        3. Runs complete analysis on each file
        4. Outputs results to file or stdout
        5. Optionally saves detailed BLAST hits
        """
        # load the profiles
        logger.info("Loading profiles:\t%s", os.path.join(self.db_dir, 'profile.txt'))
        p = Profiles(os.path.join(self.db_dir, 'profile.txt'), self.verbose)

        logger.info("Loading database:\t%s", self.db_dir)
        d = Database(self.db_dir, self.verbose)

        # Process each input genome file
        for idx, i in enumerate(self.input_files):
            if len(self.input_files) > 1:
                logger.info("Processing %d/%d: %s", idx + 1, len(self.input_files), i)
            else:
                logger.info("Beginning analysis of input file:\t%s", i)
            output_type, analysis_result = self.run_analysis(i, p, d)
            self.output_results(i, output_type)
            if analysis_result is not None:
                self.analysis_results.append(analysis_result)
            self._collect_html_result(i, output_type)

        # Write all top BLAST hits to file if requested
        if self.top_blast_hits is not None:
            with open(self.top_blast_hits, "a+") as output_fh:
                for h in self.top_results:
                    output_fh.write(str(h)+"\n")

        # Compute batch statistics when multiple results are available
        batch_stats_dict = None
        if len(self.analysis_results) > 1:
            result_dicts = [r.to_dict() for r in self.analysis_results]
            stats = BatchStats(result_dicts)
            batch_stats_dict = {
                "type_distribution": stats.type_distribution(),
                "quality_summary": stats.quality_summary(),
                "mean_confidence": stats.mean_confidence(),
                "flag_summary": stats.flag_summary(),
                "outlier_assemblies": stats.outlier_assemblies(),
                "total_assemblies": len(result_dicts),
            }

        # Write JSON output if requested
        if self.output_json is not None:
            result_dicts = [r.to_dict() for r in self.analysis_results]
            output_payload = {
                "results": result_dicts,
            }
            if batch_stats_dict is not None:
                output_payload["batch_stats"] = batch_stats_dict
            with open(self.output_json, 'w') as json_fh:
                json_fh.write(json.dumps(output_payload, indent=2))

        # Print summary table for batch mode (multiple files)
        if len(self.input_files) > 1:
            self._print_summary_table()

        # Generate batch visualizations if output_dir is set and multiple results
        if self.output_dir is not None and len(self.analysis_results) > 1:
            self._generate_batch_outputs(batch_stats_dict)

        # Generate HTML report if requested
        if self.output_html is not None:
            species = os.path.basename(self.db_dir) if self.db_dir is not None else "Unknown"
            # Use AnalysisResult data if available, otherwise fall back to parsed data
            html_data = self.html_results
            if self.analysis_results:
                html_data = [r.to_dict() for r in self.analysis_results]
            report = HtmlReport(html_data, species=species)
            report.save(self.output_html)
            logger.info("HTML report written to:\t%s", self.output_html)

    def _collect_html_result(self, input_file, profile_type):
        """Collect structured result data for the HTML report (fallback when AnalysisResult not available).

        Args:
            input_file (str): Name of input genome file.
            profile_type (str): Tab-separated quality, GS type, and fragment pattern.
        """
        if self.output_html is None:
            return

        parts = profile_type.split("\t") if profile_type else []
        quality = parts[0] if len(parts) > 0 else "RED"
        gs_type = parts[1] if len(parts) > 1 else "GS0.0"
        fragment_pattern = parts[2] if len(parts) > 2 else ""

        # Build fragment list from top_results collected during analysis
        fragments = []
        for tr in self.top_results:
            fragments.append({
                "number": str(getattr(tr, "subject", "?")),
                "reversed": not tr.is_forward() if hasattr(tr, "is_forward") else False,
                "blast_identity": getattr(tr, "percentage_identity", 0),
                "blast_alignment_length": getattr(tr, "alignment_length", 0),
                "blast_bit_score": getattr(tr, "bit_score", 0),
                "length": getattr(tr, "alignment_length", 0),
            })

        num_operons = len(fragment_pattern.split()) if fragment_pattern else 0
        is_novel = quality != "GREEN"
        qc_flags = []
        if "?" in fragment_pattern:
            qc_flags.append({
                "code": "MISSING_FRAGMENT",
                "severity": "warning",
                "message": "One or more fragments could not be identified",
            })
        if quality == "RED":
            qc_flags.append({
                "code": "LOW_QUALITY",
                "severity": "error",
                "message": "Result did not pass quality thresholds",
            })

        self.html_results.append({
            "genome_file": input_file,
            "gs_type": gs_type,
            "quality": quality,
            "confidence_score": 100.0 if quality == "GREEN" else (50.0 if quality == "AMBER" else 0.0),
            "fragment_pattern": fragment_pattern,
            "num_operons": num_operons,
            "is_novel": is_novel,
            "qc_flags": qc_flags,
            "fragments": fragments,
            "genome_length": 0,
            "validation_passed": quality == "GREEN",
        })

    def output_results(self, input_file, profile_type):
        """
        Output analysis results for a single genome.

        Args:
            input_file (str): Name of input genome file
            profile_type (str): Determined GS type and quality (e.g., "GREEN\tGS1.1")
        """
        # Default to RED quality with GS0.0 if no valid type found
        if profile_type == '':
            profile_type = "RED\tGS0.0"

        # Write to stdout or to specified output file
        if self.output_file is None:
            print(input_file + "\t" + profile_type)
        else:
            with open(self.output_file, "a+") as output_fh:
                output_fh.write(input_file + "\t" + profile_type + "\n")

    def find_rrna_boundries(self, input_file):
        """
        Identify rRNA operon boundaries in the genome using Barrnap.

        Args:
            input_file (str): Path to input genome FASTA file

        Returns:
            list: List of boundary objects marking rRNA operon locations
        """
        # run the fasta through barrnap
        fd, barrnap_outputfile = mkstemp()
        os.close(fd)
        b = Barrnap(input_file, self.threads, self.verbose)
        b.construct_barrnap_command(barrnap_outputfile)

        # Parse barrnap output to extract boundary coordinates
        boundries = b.read_barrnap_output(barrnap_outputfile)
        if logger.isEnabledFor(logging.DEBUG):
            logger.debug("Boundries:")
            for b in boundries:
                logger.debug(b)
        return boundries

    def populate_fragments_from_chromosome(self, input_file, boundries):
        """
        Extract inter-operon fragment sequences from the chromosome.

        Args:
            input_file (str): Path to input genome FASTA file
            boundries (list): rRNA operon boundary locations

        Returns:
            list: Fragment objects containing coordinates and sequences
        """
        f = Fasta(input_file, self.verbose, is_circular = self.is_circular)
        # Calculate fragment coordinates based on operon boundaries
        fragments = f.calc_fragment_coords( boundries)
        # Extract actual sequences for each fragment
        f.populate_fragments_from_chromosome(fragments, self.max_bases_from_ends)
        return fragments

    def write_novel_profile_to_file(self, tg, type_output_string):
        """
        Save novel (previously unseen) genome structure profiles to file.

        Args:
            tg (TypeGenerator): Type generator containing novelty information
            type_output_string (str): Full profile description string
        """
        if not tg.has_previously_seen:
            with open(self.novel_profiles, "a+") as output_fh:
                output_fh.write(self.db_dir + "\t" + type_output_string + "\n")

    # refactor
    def run_analysis(self, input_file, p, d):
        """
        Run complete analysis workflow for a single genome.

        This is the core analysis method that:
        1. Finds rRNA operons
        2. Extracts inter-operon fragments
        3. BLASTs fragments against database
        4. Matches fragment pattern to known profiles
        5. Validates fragment ordering
        6. Generates plots and saves novel data
        7. Builds a structured AnalysisResult

        Args:
            input_file (str): Path to input genome FASTA
            p (Profiles): Loaded profile database
            d (Database): Fragment BLAST database

        Returns:
            tuple: (type_output_string, AnalysisResult or None)
        """
        # Step 1: Find rRNA operon boundaries
        boundries = self.find_rrna_boundries(input_file)
        if not boundries:
            return '', None

        # Collect operon results for structured output
        operon_results = [
            OperonResult(
                start=op.start,
                end=op.end,
                direction='forward' if op.direction else 'reverse',
            )
            for op in boundries
        ]

        # Step 2: Extract fragment sequences
        fragments = self.populate_fragments_from_chromosome(input_file, boundries)

        # Determine genome length from the Fasta object
        fasta_obj = Fasta(input_file, self.verbose, is_circular=self.is_circular)
        genome_length = len(fasta_obj.chromosome.seq)

        # Step 3: Create temporary FASTA files for each fragment
        tmpdir = mkdtemp()
        self.dirs_to_cleanup.append(tmpdir)
        ff = FragmentFiles(fragments, tmpdir, self.verbose)
        ff.create_fragment_fastas()

        # Step 4: BLAST each fragment against the database
        blast = Blast(d.db_prefix, self.threads, self.verbose)

        # Build the fragment profile pattern, collecting BLAST results per fragment
        gat_profile = GATProfile(self.verbose, fragments = [])
        blast_results_map = {}  # fragment index -> BlastResult
        for frag_idx, current_fragment in enumerate(ff.ordered_fragments):
            fasta_file = current_fragment.output_filename

            # BLAST the fragment and filter results
            blast_results = blast.run_blast(fasta_file)
            fb = FilterBlast(blast_results, self.min_bit_score, self.min_alignment_length, self.verbose)
            top_result = fb.return_top_result()

            if top_result is None:
                # No match found - mark as unknown and save for manual review
                gat_profile.fragments.append('?')
                current_fragment.number = '?'
                blast_results_map[frag_idx] = None

                with open(fasta_file, "r") as fasta_file_fh:
                    with open(self.new_fragments, "a+") as newfrag_fh:
                        newfrag_fh.write(fasta_file_fh.read())
                continue
            else:
                self.top_results.append(top_result)
                blast_results_map[frag_idx] = top_result

                # Record fragment number from BLAST hit
                current_fragment.number = str(top_result.subject)

                # Mark special fragments (dnaA and dif)
                if str(p.dnaA_fragment_number) == str(current_fragment.number):
                    current_fragment.dna_A = True

                if str(p.dif_fragment_number) == str(current_fragment.number):
                    current_fragment.dif = True

                # Add fragment to profile with orientation indicator
                if top_result.is_forward():
                    gat_profile.fragments.append( str(top_result.subject))
                else:
                    # Mark reversed fragments with prime (')
                    current_fragment.reversed_frag = True
                    gat_profile.fragments.append( str(top_result.subject)+ '\'')

        # Step 5: Orient profile starting from dnaA and reorder fragments
        gat_profile.orientate_for_dnaA()
        reordered_frag_objs = gat_profile.reorder_fragment_objects_based_on_fragment_name_array( ff.ordered_fragments )

        # Step 6: Validate that fragments are ordered correctly
        validate_fragments = ValidateFragments(ff.ordered_fragments, genome_name=input_file)
        is_frag_valid = validate_fragments.validate()

        # Build operon direction string for output
        operon_directions_str = " ".join([current_fragment.operon_direction_str() for current_fragment in ff.ordered_fragments])
        if is_frag_valid:
            operon_directions_str = "Valid\t" + operon_directions_str
        else:
            operon_directions_str = "Invalid\t" + operon_directions_str

        logger.info("Operon directions:\t%s", operon_directions_str)
        self.output_operon_direction(input_file, operon_directions_str)

        # Step 7: Match profile pattern to known GS type
        tg = TypeGenerator(p, gat_profile, self.verbose, is_frag_valid)
        gs_type = tg.calculate_type()
        type_output_string  =  tg.quality + "\t" + gs_type + "\t" + str(gat_profile)
        self.write_novel_profile_to_file(tg, type_output_string)

        # Step 8: Generate genome structure plot for high-quality results
        if tg.quality == 'GREEN':
            pp = PlotProfile(reordered_frag_objs, self.output_plot_file, self.verbose)
            pp.create_plot()

        # Step 9: Generate SVG circular genome diagram if requested
        if self.output_svg is not None:
            pp_svg = PlotProfile(reordered_frag_objs, self.output_plot_file, self.verbose)
            gs_label = "GS" + str(tg.calculate_type())
            pp_svg.create_svg(
                self.output_svg,
                gs_type=gs_label,
                quality=tg.quality,
                genome_name=os.path.basename(input_file),
            )

        # Step 10: Build structured AnalysisResult
        fragment_results = []
        for frag_idx, current_fragment in enumerate(ff.ordered_fragments):
            br = blast_results_map.get(frag_idx)
            fr = FragmentResult(
                number=current_fragment.number,
                reversed=current_fragment.reversed_frag,
                is_dnaA=current_fragment.dna_A,
                is_dif=current_fragment.dif,
                length=current_fragment.num_bases(),
                coords=current_fragment.coords,
                blast_identity=br.identity if br else None,
                blast_alignment_length=br.alignment_length if br else None,
                blast_bit_score=br.bit_score if br else None,
                blast_e_value=br.e_value if br else None,
                blast_mismatches=br.mismatches if br else None,
                blast_subject=str(br.subject) if br else None,
            )
            fragment_results.append(fr)

        is_novel = not tg.has_previously_seen

        # Calculate confidence score
        confidence = calculate_confidence(
            fragment_results,
            tg.quality,
            is_frag_valid,
            is_novel=is_novel,
        )

        # Build analysis result (without QC flags initially)
        analysis_result = AnalysisResult(
            genome_file=input_file,
            genome_length=genome_length,
            is_circular=self.is_circular,
            num_operons=len(boundries),
            gs_type=gs_type,
            quality=tg.quality,
            is_novel=is_novel,
            fragment_pattern=str(gat_profile),
            orientation_binary=gat_profile.orientation_binary(),
            confidence_score=confidence,
            fragments=fragment_results,
            operons=operon_results,
            validation_passed=is_frag_valid,
            operon_direction_string=operon_directions_str,
        )

        # Generate QC flags
        analysis_result.qc_flags = generate_qc_flags(
            analysis_result,
            expected_fragment_count=p.num_fragments,
        )

        # Assess novelty if the profile is novel
        if is_novel:
            from dataclasses import asdict as _asdict
            known_profiles = [g.fragments for g in p.gats]
            blast_identities = [
                fr.blast_identity for fr in fragment_results
                if fr.blast_identity is not None
            ]
            novelty = assess_novelty(
                query_fragments=gat_profile.fragments,
                known_profiles=known_profiles,
                confidence_score=confidence,
                blast_identities=blast_identities,
            )
            analysis_result.novelty_assessment = _asdict(novelty)

        return type_output_string, analysis_result

    def _print_summary_table(self):
        """
        Print a summary table of all analysis results to stdout.

        Formatted as a fixed-width table suitable for terminal display.
        Only printed when multiple input files are processed.
        """
        if not self.analysis_results:
            return

        header = "{:<25s} {:<12s} {:>7s} {:>10s} {:>7s} {:>10s} {:>5s}  {}".format(
            'Assembly', 'GS Type', 'Quality', 'Confidence',
            'Operons', 'Fragments', 'Novel', 'Flags',
        )
        print("\n" + header)
        print("-" * len(header))

        for r in self.analysis_results:
            basename = os.path.basename(r.genome_file)
            if len(basename) > 24:
                basename = basename[:21] + "..."
            total_frags = len(r.fragments)
            matched_frags = sum(1 for f in r.fragments if str(f.number) != '?')
            frag_str = "{}/{}".format(matched_frags, total_frags)
            novel_str = "Yes" if r.is_novel else "No"
            flag_codes = ",".join(sorted({f.code for f in r.qc_flags}))

            print("{:<25s} {:<12s} {:>7s} {:>10.1f} {:>7d} {:>10s} {:>5s}  {}".format(
                basename, r.gs_type, r.quality, r.confidence_score,
                r.num_operons, frag_str, novel_str, flag_codes,
            ))

    def _generate_batch_outputs(self, batch_stats_dict):
        """Generate batch visualization files to the output directory.

        Creates type distribution, confidence heatmap, synteny, and
        per-assembly fragment quality SVGs, plus a batch_stats.json file.

        Args:
            batch_stats_dict: Pre-computed batch statistics dictionary,
                or None if not yet computed.
        """
        output_dir = self.output_dir
        os.makedirs(output_dir, exist_ok=True)

        result_dicts = [r.to_dict() for r in self.analysis_results]

        # Batch stats JSON
        if batch_stats_dict is not None:
            stats_path = os.path.join(output_dir, 'batch_stats.json')
            with open(stats_path, 'w') as fh:
                fh.write(json.dumps(batch_stats_dict, indent=2))
            logger.info("Batch stats written to:\t%s", stats_path)

        # Build assembly data for SVG generators
        assemblies = []
        for rd in result_dicts:
            frags = []
            for f in rd.get('fragments', []):
                fnum = f.get('number', 0)
                try:
                    fnum = int(fnum)
                except (ValueError, TypeError):
                    fnum = 0
                frags.append({
                    'number': fnum,
                    'reversed': f.get('reversed', False),
                    'length': f.get('length', 1),
                    'blast_identity': f.get('blast_identity'),
                    'blast_alignment_length': f.get('blast_alignment_length'),
                    'is_dnaA': f.get('is_dnaA', False),
                    'is_dif': f.get('is_dif', False),
                })
            assemblies.append({
                'name': os.path.basename(rd.get('genome_file', '')),
                'gs_type': rd.get('gs_type', ''),
                'quality': rd.get('quality', 'RED'),
                'fragments': frags,
            })

        # Type distribution SVG
        stats = BatchStats(result_dicts)
        type_dist = stats.type_distribution()
        svg_path = os.path.join(output_dir, 'type_distribution.svg')
        with open(svg_path, 'w') as fh:
            fh.write(generate_type_distribution_svg(type_dist))
        logger.info("Type distribution SVG written to:\t%s", svg_path)

        # Confidence heatmap SVG
        svg_path = os.path.join(output_dir, 'confidence_heatmap.svg')
        with open(svg_path, 'w') as fh:
            fh.write(generate_confidence_heatmap_svg(assemblies))
        logger.info("Confidence heatmap SVG written to:\t%s", svg_path)

        # Synteny SVG
        svg_path = os.path.join(output_dir, 'synteny.svg')
        with open(svg_path, 'w') as fh:
            fh.write(generate_synteny_svg(assemblies))
        logger.info("Synteny SVG written to:\t%s", svg_path)

        # Per-assembly fragment quality SVGs
        for asm in assemblies:
            name = asm.get('name', 'unknown')
            safe_name = name.replace('/', '_').replace(' ', '_')
            svg_path = os.path.join(output_dir, f'fragment_quality_{safe_name}.svg')
            with open(svg_path, 'w') as fh:
                fh.write(generate_fragment_quality_svg(
                    asm['fragments'], genome_name=name,
                ))
            logger.info("Fragment quality SVG written to:\t%s", svg_path)

    def output_operon_direction(self, input_file, operon_directions):
        """
        Write operon direction information to output file.

        Args:
            input_file (str): Name of input genome file
            operon_directions (str): String describing operon orientations
        """
        with open(self.output_operon_directions_file, "a+") as output_fh:
            output_fh.write(input_file + "\t" + operon_directions + "\n")

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.cleanup()
        return False

    def cleanup(self):
        """Clean up temporary files and directories."""
        for f in getattr(self, 'dirs_to_cleanup', []):
            if os.path.exists(f):
                shutil.rmtree(f, ignore_errors=True)
        self.dirs_to_cleanup = []

    def __del__(self):
        """Safety net cleanup -- prefer using as context manager."""
        self.cleanup()
