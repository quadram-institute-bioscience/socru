"""Tests for ToolCheck module -- external tool availability checks."""

import unittest
from unittest.mock import patch

from socru.ToolCheck import MissingToolError, check_all_tools, check_tool


class TestCheckTool(unittest.TestCase):
    """Verify check_tool behaviour for present and missing executables."""

    def test_check_tool_finds_python(self):
        """python3 is always available in the test environment."""
        # Should not raise
        check_tool('python3')

    def test_check_tool_raises_for_missing(self):
        """A nonexistent tool name must raise MissingToolError."""
        with self.assertRaises(MissingToolError):
            check_tool('absolutely_nonexistent_tool_xyz_999')

    def test_missing_tool_error_message(self):
        """The error message must contain the tool name and optional hint."""
        tool_name = 'fake_tool_abc'
        hint = 'Install via: apt install fake_tool_abc'
        with self.assertRaises(MissingToolError) as ctx:
            check_tool(tool_name, hint=hint)
        self.assertIn(tool_name, str(ctx.exception))
        self.assertIn(hint, str(ctx.exception))

    def test_missing_tool_error_no_hint(self):
        """Without a hint the message should still be informative."""
        with self.assertRaises(MissingToolError) as ctx:
            check_tool('no_such_binary')
        self.assertIn("no_such_binary", str(ctx.exception))
        self.assertIn("not found on PATH", str(ctx.exception))

    def test_check_all_tools_when_available(self):
        """check_all_tools should pass when barrnap and BLAST+ are installed."""
        # This test will only pass in environments where the tools are present.
        # It is safe to run in CI with bioconda packages installed.
        try:
            check_all_tools()
        except MissingToolError:
            self.skipTest("barrnap/blast not installed in this environment")

    @patch('socru.ToolCheck.shutil.which', return_value=None)
    def test_check_all_tools_missing_barrnap(self, mock_which):
        """check_all_tools raises MissingToolError when barrnap is missing."""
        with self.assertRaises(MissingToolError) as ctx:
            check_all_tools()
        self.assertIn('barrnap', str(ctx.exception))
