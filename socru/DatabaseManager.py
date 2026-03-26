"""
Database management for socru species databases.

This module provides a unified interface for discovering, installing, and
querying species databases. It supports both bundled databases (shipped with
the socru package) and user-installed databases stored in a configurable
data directory.

Search order for database resolution:
    1. User data directory (SOCRU_DATA_DIR or ~/.socru/data/)
    2. Bundled package data (socru/data/)

Classes:
    DatabaseManager: Manage socru species databases -- bundled and user-installed.
"""

import glob
import logging
import os
import shutil
from typing import Dict, List, Optional

import importlib.resources

logger = logging.getLogger(__name__)

# Default locations for database storage
DEFAULT_DATA_DIR = os.path.expanduser('~/.socru/data')


class DatabaseManager:
    """Manage socru species databases -- bundled and user-installed.

    Provides methods to locate, list, install, and inspect species databases.
    Databases are searched first in the user data directory, then in the
    bundled package data directory.

    Attributes:
        data_dir: Path to user database directory.
    """

    def __init__(self, data_dir: Optional[str] = None) -> None:
        """Initialize DatabaseManager.

        Args:
            data_dir: Custom database directory. Falls back to:
                1. SOCRU_DATA_DIR env var
                2. ~/.socru/data/
                3. Bundled package data
        """
        self.data_dir: str = data_dir or os.environ.get('SOCRU_DATA_DIR') or DEFAULT_DATA_DIR

    def _bundled_data_dir(self) -> Optional[str]:
        """Return path to bundled package data directory, or None."""
        try:
            bundled = str(importlib.resources.files('socru') / 'data')
            if os.path.isdir(bundled):
                return bundled
        except Exception:
            pass
        return None

    def get_database_dir(self, species: str) -> Optional[str]:
        """Find the database directory for a species.

        Search order:
            1. User data dir (self.data_dir / species)
            2. Bundled package data

        Args:
            species: Species name (directory name).

        Returns:
            Path to the database directory, or None if not found.
        """
        # Check user data dir first
        user_path = os.path.join(self.data_dir, species)
        if os.path.isdir(user_path):
            return user_path

        # Fall back to bundled
        bundled_dir = self._bundled_data_dir()
        if bundled_dir is not None:
            bundled_path = os.path.join(bundled_dir, species)
            if os.path.isdir(bundled_path):
                return bundled_path

        return None

    def list_species(self, include_bundled: bool = True) -> List[str]:
        """List all available species databases.

        Args:
            include_bundled: Whether to include bundled databases in the list.

        Returns:
            Sorted list of species names.
        """
        species: set = set()

        # User databases
        if os.path.isdir(self.data_dir):
            for name in os.listdir(self.data_dir):
                path = os.path.join(self.data_dir, name)
                if os.path.isdir(path) and not name.startswith('.'):
                    species.add(name)

        # Bundled databases
        if include_bundled:
            bundled_dir = self._bundled_data_dir()
            if bundled_dir is not None:
                for name in os.listdir(bundled_dir):
                    path = os.path.join(bundled_dir, name)
                    if os.path.isdir(path):
                        species.add(name)

        return sorted(species)

    def install_database(self, source_dir: str, species_name: Optional[str] = None) -> str:
        """Install a database from a directory.

        Copies all files from the source directory into the user data
        directory under the given species name.

        Args:
            source_dir: Path to the source database directory.
            species_name: Name for the installed database. Defaults to
                the basename of source_dir.

        Returns:
            Path to the installed database directory.
        """
        if species_name is None:
            species_name = os.path.basename(os.path.normpath(source_dir))

        dest = os.path.join(self.data_dir, species_name)
        os.makedirs(dest, exist_ok=True)

        for item in os.listdir(source_dir):
            src = os.path.join(source_dir, item)
            dst = os.path.join(dest, item)
            if os.path.isfile(src):
                shutil.copy2(src, dst)

        logger.info("Installed database for %s to %s", species_name, dest)
        return dest

    def database_info(self, species: str) -> Optional[Dict]:
        """Return info dict about a species database.

        Args:
            species: Species name.

        Returns:
            Dictionary with database metadata, or None if the species
            database is not found. Keys include:
                - species: species name
                - path: absolute path to database directory
                - fragment_count: number of fragment FASTA files
                - has_profile: whether profile.txt exists
                - is_bundled: whether the database is from the package
                - known_types: number of known types (if profile exists)
        """
        db_dir = self.get_database_dir(species)
        if db_dir is None:
            return None

        fa_files = glob.glob(os.path.join(db_dir, '*.fa*'))
        profile_path = os.path.join(db_dir, 'profile.txt')

        info: Dict = {
            'species': species,
            'path': db_dir,
            'fragment_count': len([f for f in fa_files if not f.endswith('.yml')]),
            'has_profile': os.path.exists(profile_path),
            'is_bundled': 'site-packages' in db_dir or os.sep + 'socru' + os.sep + 'data' + os.sep in db_dir,
        }

        if info['has_profile']:
            with open(profile_path) as f:
                lines = [line for line in f if line.strip() and not line.startswith('#')]
                info['known_types'] = len(lines)

        return info
