"""
External tool availability checks.

This module verifies that required external tools (barrnap, BLAST+) are
available on the system PATH before the pipeline attempts to use them.
Provides clear, actionable error messages when tools are missing.

Functions:
    check_tool: Verify a single tool is on PATH
    check_all_tools: Verify all required external tools
"""

import logging
import shutil

logger = logging.getLogger(__name__)


class MissingToolError(RuntimeError):
    """Raised when a required external tool is not found on PATH."""
    pass


def check_tool(name: str, hint: str = "") -> None:
    """Check if an external tool is on PATH. Raise MissingToolError if not.

    Args:
        name: Name of the executable to look up.
        hint: Optional installation hint appended to the error message.

    Raises:
        MissingToolError: If the tool cannot be found via ``shutil.which``.
    """
    if shutil.which(name) is None:
        msg = f"Required tool '{name}' not found on PATH."
        if hint:
            msg += f" {hint}"
        raise MissingToolError(msg)


def check_all_tools() -> None:
    """Check that barrnap and BLAST+ are available.

    Raises:
        MissingToolError: If any required tool is missing.
    """
    check_tool('barrnap', 'Install via: conda install -c bioconda barrnap')
    check_tool('blastn', 'Install via: conda install -c bioconda blast')
    check_tool('makeblastdb', 'Install via: conda install -c bioconda blast')
