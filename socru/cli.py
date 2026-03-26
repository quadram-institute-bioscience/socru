"""
Command-line interface entry points for socru.

Each function in this module serves as a console_scripts entry point,
wrapping the logic previously found in the scripts/ directory.
"""

import logging
import sys
import argparse
import pkg_resources


def _get_version():
    """Get the installed package version, falling back to 'x.y.z'."""
    try:
        return pkg_resources.get_distribution("socru").version
    except pkg_resources.DistributionNotFound:
        return "x.y.z"


def socru_main():
    """Entry point for the main socru command."""
    from socru.Socru import Socru

    version = _get_version()

    parser = argparse.ArgumentParser(
        description=(
            'Please cite our paper, "Socru: Typing of genome level order and orientation in bacteria", '
            "Andrew J Page, Gemma Langridge, "
            "bioRxiv 543702; (2019) doi: https://doi.org/10.1101/543702"
        ),
        usage="socru [options] species assembly.fasta",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("species", help="Species name, use socru_species to see all available", type=str)
    parser.add_argument("input_files", help="Input FASTA files (optionally gzipped)", nargs="+", type=str)
    parser.add_argument("--db_dir", "-d", help="Base directory for species databases, defaults to bundled", type=str)
    parser.add_argument("--threads", "-t", help="No. of threads to use", type=int, default=1)
    parser.add_argument("--output_file", "-o", help="Output filename, defaults to STDOUT", type=str)
    parser.add_argument(
        "--output_plot_file",
        "-p",
        help="Filename of plot of genome structure",
        type=str,
        default="genome_structure.pdf",
    )
    parser.add_argument(
        "--novel_profiles", "-n", help="Filename for novel profiles", type=str, default="profile.txt.novel"
    )
    parser.add_argument(
        "--new_fragments", "-f", help="Filename for novel fragments", type=str, default="novel_fragments.fa"
    )
    parser.add_argument("--top_blast_hits", "-b", help="Filename for top blast hits", type=str)
    parser.add_argument(
        "--output_operon_directions_file",
        "-r",
        help="Filename of directions of operons",
        type=str,
        default="operon_directions.txt",
    )
    parser.add_argument(
        "--max_bases_from_ends",
        "-m",
        help="Only look at this number of bases from start and end of fragment",
        type=int,
    )
    parser.add_argument(
        "--not_circular", "-c", action="store_true", help="Assume chromosome is not circularised", default=False
    )
    parser.add_argument("--min_bit_score", help="Minimum bit score", type=int, default=100)
    parser.add_argument("--min_alignment_length", help="Minimum alignment length", type=int, default=100)
    parser.add_argument("--debug", action="store_true", help="Turn on debugging", default=False)
    parser.add_argument("--verbose", "-v", action="store_true", help="Turn on verbose output", default=False)
    parser.add_argument("--version", action="version", version=str(version))

    options = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if options.verbose else logging.WARNING,
        format='%(levelname)s: %(message)s',
    )

    if options.debug:
        import cProfile
        import io
        import pstats

        pr = cProfile.Profile()
        pr.enable()
        g = Socru(options)
        g.run()
        pr.disable()
        s = io.StringIO()
        sortby = "cumulative"
        ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
        ps.print_stats()
        print(s.getvalue())
    else:
        g = Socru(options)
        g.run()


def socru_create_main():
    """Entry point for socru_create command."""
    from socru.SocruCreate import SocruCreate

    version = _get_version()

    parser = argparse.ArgumentParser(
        description=(
            "Create genome arrangement type scheme. "
            'Please cite our paper, "Socru: Typing of genome level order and orientation in bacteria", '
            "Andrew J Page, Gemma Langridge, "
            "bioRxiv 543702; (2019) doi: https://doi.org/10.1101/543702"
        ),
        usage="socru_create [options] output_directory assembly.fasta",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("output_directory", help="Output directory", type=str)
    parser.add_argument("input_file", help="Input FASTA file (optionally gzipped)", type=str)
    parser.add_argument(
        "--max_bases_from_ends",
        "-m",
        help="Only look at this number of bases from start and end of fragment",
        type=int,
    )
    parser.add_argument("--threads", "-t", help="No. of threads to use", type=int, default=1)
    parser.add_argument(
        "--fragment_order",
        "-f",
        help="Order of fragments, you may need to change this, example 1-2-3-4-5-6-7",
        type=str,
    )
    parser.add_argument("--dnaa_fasta", "-d", help="Location of dnaA FASTA file, defaults to bundled", type=str)
    parser.add_argument("--dif_fasta", "-e", help="Location of dif FASTA file, defaults to bundled", type=str)
    parser.add_argument("--debug", action="store_true", help="Turn on debugging", default=False)
    parser.add_argument("--verbose", "-v", action="store_true", help="Turn on verbose output", default=False)
    parser.add_argument("--version", action="version", version=str(version))

    options = parser.parse_args()
    g = SocruCreate(options)
    g.run()


def socru_rebuild_profile_main():
    """Entry point for socru_rebuild_profile command."""
    from socru.SocruRebuildProfile import SocruRebuildProfile

    version = _get_version()

    parser = argparse.ArgumentParser(
        description=(
            "Admin utility which will take in an existing profile and rebuild it. "
            'Please cite our paper, "Socru: Typing of genome level order and orientation in bacteria", '
            "Andrew J Page, Gemma Langridge, "
            "bioRxiv 543702; (2019) doi: https://doi.org/10.1101/543702"
        ),
        usage="socru_rebuild_profile [options]",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("profile_filename", help="profile.txt from database", type=str)
    parser.add_argument("--output_file", "-o", help="Output filename", type=str, default="updated_profile.txt")
    parser.add_argument("--prefix", "-p", help="Prefix", type=str, default="GS")
    parser.add_argument("--debug", action="store_true", help="Turn on debugging", default=False)
    parser.add_argument("--verbose", "-v", action="store_true", help="Turn on verbose output", default=False)
    parser.add_argument("--version", action="version", version=str(version))

    options = parser.parse_args()
    g = SocruRebuildProfile(options)
    g.run()


def socru_shrink_database_main():
    """Entry point for socru_shrink_database command."""
    from socru.SocruShrinkDatabase import SocruShrinkDatabase

    version = _get_version()

    parser = argparse.ArgumentParser(
        description="Admin utility take a database and blast results and create a new shrunk database",
        usage="socru_shrink_database [options]",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "blast_results",
        help="Blast results file from running socru -b xxx against multiple assemblies",
        type=str,
    )
    parser.add_argument("input_database", help="Directory containing database to shrink", type=str)
    parser.add_argument(
        "output_database", help="Output directory for new database, it must not already exist", type=str
    )
    parser.add_argument("--min_fragment_size", help="Minimum fragment size in bases", type=int, default=100000)
    parser.add_argument("--debug", action="store_true", help="Turn on debugging", default=False)
    parser.add_argument("--verbose", "-v", action="store_true", help="Turn on verbose output", default=False)
    parser.add_argument("--version", action="version", version=str(version))

    options = parser.parse_args()
    g = SocruShrinkDatabase(options)
    print(g.run())


def socru_species_main():
    """Entry point for socru_species command."""
    from socru.Schemas import Schemas

    version = _get_version()

    parser = argparse.ArgumentParser(
        description="List all available species",
        usage="socru_species [options]",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--extended",
        "-e",
        action="store_true",
        help="Extended information about the species databases",
        default=False,
    )
    parser.add_argument("--debug", action="store_true", help="Turn on debugging", default=False)
    parser.add_argument("--verbose", "-v", action="store_true", help="Turn on verbose output", default=False)
    parser.add_argument("--version", action="version", version=str(version))

    options = parser.parse_args()

    if options.extended:
        Schemas(options.verbose).print_extended()
    else:
        Schemas(options.verbose).print_all()


def socru_update_profile_main():
    """Entry point for socru_update_profile command."""
    from socru.SocruUpdateProfile import SocruUpdateProfile

    version = _get_version()

    parser = argparse.ArgumentParser(
        description="Admin utility to take the novel GS results and update the profile for the database",
        usage="socru_update_profile [options]",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("socru_output_filename", help="Socru output file", type=str)
    parser.add_argument("profile_filename", help="profile.txt from database", type=str)
    parser.add_argument("--output_file", "-o", help="Output filename", type=str, default="updated_profile.txt")
    parser.add_argument("--debug", action="store_true", help="Turn on debugging", default=False)
    parser.add_argument("--verbose", "-v", action="store_true", help="Turn on verbose output", default=False)
    parser.add_argument("--version", action="version", version=str(version))

    options = parser.parse_args()
    g = SocruUpdateProfile(options)
    g.run()
