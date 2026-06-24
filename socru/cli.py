"""Command-line interface entry points for socru."""
from __future__ import annotations

import argparse
import logging
import sys
from importlib.metadata import PackageNotFoundError
from importlib.metadata import version as get_version


def _version():
    try:
        return get_version("socru")
    except PackageNotFoundError:
        return "unknown"


def socru_main():
    """Entry point for the main socru command."""
    from socru.Socru import Socru

    parser = argparse.ArgumentParser(
        description='Socru: typing of genome-level order and orientation around ribosomal operons in bacteria',
        usage='socru [options] species assembly.fasta',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument('species', help='Species name (use socru_species to list)', type=str)
    parser.add_argument('input_files', help='Input FASTA files (optionally gzipped)', nargs='+', type=str)
    parser.add_argument('--db_dir', '-d', help='Base directory for species databases', type=str)
    parser.add_argument('--threads', '-t', help='Number of threads', type=int, default=1)
    parser.add_argument('--output_file', '-o', help='Output filename (default: STDOUT)', type=str)
    parser.add_argument('--output_plot_file', '-p', help='PDF plot filename', type=str, default='genome_structure.pdf')
    parser.add_argument('--output_svg', '-s', help='SVG circular genome diagram filename', type=str)
    parser.add_argument('--output_json', '-j', help='JSON structured results filename', type=str)
    parser.add_argument('--output_html', help='Self-contained HTML report filename', type=str)
    parser.add_argument('--output_dir', help='Directory for batch visualization outputs', type=str)
    parser.add_argument('--novel_profiles', '-n', help='Novel profiles filename', type=str, default='profile.txt.novel')
    parser.add_argument('--new_fragments', '-f', help='Novel fragments filename', type=str, default='novel_fragments.fa')
    parser.add_argument('--top_blast_hits', '-b', help='Top BLAST hits filename', type=str)
    parser.add_argument('--output_operon_directions_file', '-r', help='Operon directions filename', type=str, default='operon_directions.txt')
    parser.add_argument('--max_bases_from_ends', '-m', help='Max bases from fragment ends', type=int)
    parser.add_argument('--not_circular', '-c', action='store_true', help='Assume chromosome is not circular', default=False)
    parser.add_argument('--min_bit_score', help='Minimum BLAST bit score', type=int, default=100)
    parser.add_argument('--min_alignment_length', help='Minimum BLAST alignment length', type=int, default=100)
    parser.add_argument('--debug', action='store_true', help='Enable profiling', default=False)
    parser.add_argument('--verbose', '-v', action='store_true', help='Verbose output', default=False)
    parser.add_argument('--version', action='version', version='%(prog)s ' + _version())

    options = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if options.verbose else logging.WARNING,
        format='%(levelname)s: %(message)s',
    )

    try:
        if options.debug:
            import cProfile
            import io
            import pstats
            pr = cProfile.Profile()
            pr.enable()
            with Socru(options) as g:
                g.run()
            pr.disable()
            s = io.StringIO()
            ps = pstats.Stats(pr, stream=s).sort_stats('cumulative')
            ps.print_stats()
            print(s.getvalue())
        else:
            with Socru(options) as g:
                g.run()
    except (FileNotFoundError, ValueError) as e:
        print(str(e), file=sys.stderr)
        sys.exit(1)


def socru_create_main():
    """Entry point for socru_create command."""
    from socru.SocruCreate import SocruCreate

    parser = argparse.ArgumentParser(
        description='Create a new Socru species database from a reference genome',
        usage='socru_create [options] output_directory assembly.fasta',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument('output_directory', help='Output directory', type=str)
    parser.add_argument('input_file', help='Input FASTA file (optionally gzipped)', type=str)
    parser.add_argument('--max_bases_from_ends', '-m', help='Max bases from fragment ends', type=int)
    parser.add_argument('--threads', '-t', help='Number of threads', type=int, default=1)
    parser.add_argument('--fragment_order', '-f', help='Custom fragment order (e.g. 1-2-3-4-5-6-7)', type=str)
    parser.add_argument('--dnaa_fasta', '-d', help='Custom dnaA FASTA file', type=str)
    parser.add_argument('--dif_fasta', '-e', help='Custom dif FASTA file', type=str)
    parser.add_argument('--debug', action='store_true', help='Enable debugging', default=False)
    parser.add_argument('--verbose', '-v', action='store_true', help='Verbose output', default=False)
    parser.add_argument('--version', action='version', version='%(prog)s ' + _version())

    options = parser.parse_args()

    try:
        with SocruCreate(options) as g:
            g.run()
    except (FileNotFoundError, FileExistsError, ValueError) as e:
        print(str(e), file=sys.stderr)
        sys.exit(1)


def socru_rebuild_profile_main():
    """Entry point for socru_rebuild_profile command."""
    from socru.SocruRebuildProfile import SocruRebuildProfile

    parser = argparse.ArgumentParser(
        description='Rebuild and renumber a Socru profile database',
        usage='socru_rebuild_profile [options] profile.txt',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument('profile_filename', help='profile.txt from database', type=str)
    parser.add_argument('--output_file', '-o', help='Output filename', type=str, default='updated_profile.txt')
    parser.add_argument('--prefix', '-p', help='Type prefix', type=str, default='GS')
    parser.add_argument('--debug', action='store_true', help='Enable debugging', default=False)
    parser.add_argument('--verbose', '-v', action='store_true', help='Verbose output', default=False)
    parser.add_argument('--version', action='version', version='%(prog)s ' + _version())

    options = parser.parse_args()

    with SocruRebuildProfile(options) as g:
        g.run()


def socru_shrink_database_main():
    """Entry point for socru_shrink_database command."""
    from socru.SocruShrinkDatabase import SocruShrinkDatabase

    parser = argparse.ArgumentParser(
        description='Shrink a Socru database using BLAST coverage analysis',
        usage='socru_shrink_database [options] blast_results input_db output_db',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument('blast_results', help='BLAST results file from socru -b', type=str)
    parser.add_argument('input_database', help='Input database directory', type=str)
    parser.add_argument('output_database', help='Output database directory (must not exist)', type=str)
    parser.add_argument('--min_fragment_size', help='Minimum fragment size in bases', type=int, default=100000)
    parser.add_argument('--debug', action='store_true', help='Enable debugging', default=False)
    parser.add_argument('--verbose', '-v', action='store_true', help='Verbose output', default=False)
    parser.add_argument('--version', action='version', version='%(prog)s ' + _version())

    options = parser.parse_args()
    g = SocruShrinkDatabase(options)
    print(g.run())


def socru_species_main():
    """Entry point for socru_species command."""
    from socru.DatabaseManager import DatabaseManager
    from socru.Schemas import Schemas

    parser = argparse.ArgumentParser(
        description='List all available Socru species databases',
        usage='socru_species [options]',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument('--extended', '-e', action='store_true', help='Extended database info', default=False)
    parser.add_argument('--detailed', '-d', action='store_true', help='Detailed info (fragments, types, path)', default=False)
    parser.add_argument('--debug', action='store_true', help='Enable debugging', default=False)
    parser.add_argument('--verbose', '-v', action='store_true', help='Verbose output', default=False)
    parser.add_argument('--version', action='version', version='%(prog)s ' + _version())

    options = parser.parse_args()

    if options.detailed:
        dm = DatabaseManager()
        species_list = dm.list_species()
        print("\t".join(['Species', 'Fragments', 'Known Types', 'Bundled', 'Path']))
        for species in species_list:
            info = dm.database_info(species)
            if info is not None:
                print("\t".join([
                    info['species'],
                    str(info['fragment_count']),
                    str(info.get('known_types', 'N/A')),
                    'yes' if info.get('is_bundled', False) else 'no',
                    info['path'],
                ]))
    elif options.extended:
        Schemas(options.verbose).print_extended()
    else:
        Schemas(options.verbose).print_all()


def socru_update_profile_main():
    """Entry point for socru_update_profile command."""
    from socru.SocruUpdateProfile import SocruUpdateProfile

    parser = argparse.ArgumentParser(
        description='Update a Socru profile database with novel GS types',
        usage='socru_update_profile [options] socru_output profile.txt',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument('socru_output_filename', help='Socru output file with novel profiles', type=str)
    parser.add_argument('profile_filename', help='profile.txt from database', type=str)
    parser.add_argument('--output_file', '-o', help='Output filename', type=str, default='updated_profile.txt')
    parser.add_argument('--debug', action='store_true', help='Enable debugging', default=False)
    parser.add_argument('--verbose', '-v', action='store_true', help='Verbose output', default=False)
    parser.add_argument('--version', action='version', version='%(prog)s ' + _version())

    options = parser.parse_args()
    g = SocruUpdateProfile(options)
    g.run()
