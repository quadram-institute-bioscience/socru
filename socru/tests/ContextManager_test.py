"""Tests for context manager implementations across socru modules."""

import unittest
import tempfile
import os
import shutil

from socru.Barrnap import Barrnap
from socru.Blast import Blast
from socru.SocruCreate import SocruCreate
from socru.SocruRebuildProfile import SocruRebuildProfile


class TestBarrnapContextManager(unittest.TestCase):
    """Tests for Barrnap context manager cleanup."""

    def test_cleanup_removes_temp_files(self):
        """Barrnap.cleanup() should remove all tracked temporary files."""
        b = Barrnap("dummy_input.fa", 1, False)
        # Create actual temp files and register them
        fd1, path1 = tempfile.mkstemp()
        fd2, path2 = tempfile.mkstemp()
        os.close(fd1)
        os.close(fd2)
        b.files_to_cleanup = [path1, path2]

        b.cleanup()

        self.assertFalse(os.path.exists(path1))
        self.assertFalse(os.path.exists(path2))
        self.assertEqual(b.files_to_cleanup, [])

    def test_context_manager_calls_cleanup(self):
        """Using Barrnap as context manager should clean up on exit."""
        fd, path = tempfile.mkstemp()
        os.close(fd)

        with Barrnap("dummy_input.fa", 1, False) as b:
            b.files_to_cleanup.append(path)
            self.assertTrue(os.path.exists(path))

        self.assertFalse(os.path.exists(path))

    def test_context_manager_cleanup_on_exception(self):
        """Barrnap should clean up even if an exception occurs inside the with block."""
        fd, path = tempfile.mkstemp()
        os.close(fd)

        with self.assertRaises(ValueError):
            with Barrnap("dummy_input.fa", 1, False) as b:
                b.files_to_cleanup.append(path)
                raise ValueError("test error")

        self.assertFalse(os.path.exists(path))

    def test_context_manager_returns_self(self):
        """__enter__ should return the Barrnap instance."""
        b = Barrnap("dummy_input.fa", 1, False)
        with b as ctx:
            self.assertIs(ctx, b)

    def test_cleanup_tolerates_missing_files(self):
        """cleanup() should not raise if files are already gone."""
        b = Barrnap("dummy_input.fa", 1, False)
        b.files_to_cleanup = ["/tmp/nonexistent_file_12345"]
        # Should not raise
        b.cleanup()
        self.assertEqual(b.files_to_cleanup, [])

    def test_cleanup_is_idempotent(self):
        """Calling cleanup() multiple times should be safe."""
        fd, path = tempfile.mkstemp()
        os.close(fd)
        b = Barrnap("dummy_input.fa", 1, False)
        b.files_to_cleanup = [path]

        b.cleanup()
        b.cleanup()  # second call should not raise
        self.assertFalse(os.path.exists(path))


class TestBlastContextManager(unittest.TestCase):
    """Tests for Blast context manager cleanup."""

    def test_cleanup_removes_temp_files(self):
        """Blast.cleanup() should remove all tracked temporary files."""
        b = Blast("dummy_db", 1, False)
        fd1, path1 = tempfile.mkstemp()
        fd2, path2 = tempfile.mkstemp()
        os.close(fd1)
        os.close(fd2)
        b.files_to_cleanup = [path1, path2]

        b.cleanup()

        self.assertFalse(os.path.exists(path1))
        self.assertFalse(os.path.exists(path2))
        self.assertEqual(b.files_to_cleanup, [])

    def test_context_manager_calls_cleanup(self):
        """Using Blast as context manager should clean up on exit."""
        fd, path = tempfile.mkstemp()
        os.close(fd)

        with Blast("dummy_db", 1, False) as b:
            b.files_to_cleanup.append(path)
            self.assertTrue(os.path.exists(path))

        self.assertFalse(os.path.exists(path))

    def test_context_manager_cleanup_on_exception(self):
        """Blast should clean up even if an exception occurs."""
        fd, path = tempfile.mkstemp()
        os.close(fd)

        with self.assertRaises(RuntimeError):
            with Blast("dummy_db", 1, False) as b:
                b.files_to_cleanup.append(path)
                raise RuntimeError("test error")

        self.assertFalse(os.path.exists(path))

    def test_context_manager_returns_self(self):
        """__enter__ should return the Blast instance."""
        b = Blast("dummy_db", 1, False)
        with b as ctx:
            self.assertIs(ctx, b)


class TestSocruRebuildProfileContextManager(unittest.TestCase):
    """Tests for SocruRebuildProfile context manager cleanup."""

    def _make_options(self):
        """Create a minimal options object for SocruRebuildProfile."""

        class Options:
            profile_filename = "dummy_profile.txt"
            output_file = "dummy_output.txt"
            prefix = "GS"
            verbose = False

        return Options()

    def test_cleanup_removes_temp_files(self):
        """SocruRebuildProfile.cleanup() should remove tracked temp files."""
        opts = self._make_options()
        s = SocruRebuildProfile(opts)
        fd, path = tempfile.mkstemp()
        os.close(fd)
        s.files_to_cleanup = [path]

        s.cleanup()

        self.assertFalse(os.path.exists(path))
        self.assertEqual(s.files_to_cleanup, [])

    def test_context_manager_cleanup_on_exception(self):
        """SocruRebuildProfile should clean up even on exception."""
        opts = self._make_options()
        fd, path = tempfile.mkstemp()
        os.close(fd)

        with self.assertRaises(ValueError):
            with SocruRebuildProfile(opts) as s:
                s.files_to_cleanup.append(path)
                raise ValueError("test error")

        self.assertFalse(os.path.exists(path))

    def test_context_manager_returns_self(self):
        """__enter__ should return the SocruRebuildProfile instance."""
        opts = self._make_options()
        s = SocruRebuildProfile(opts)
        with s as ctx:
            self.assertIs(ctx, s)


class TestExitReturnValue(unittest.TestCase):
    """Verify that __exit__ returns False so exceptions propagate."""

    def test_barrnap_exit_returns_false(self):
        b = Barrnap("dummy.fa", 1, False)
        result = b.__exit__(None, None, None)
        self.assertFalse(result)

    def test_blast_exit_returns_false(self):
        b = Blast("dummy_db", 1, False)
        result = b.__exit__(None, None, None)
        self.assertFalse(result)

    def test_rebuild_profile_exit_returns_false(self):
        class Options:
            profile_filename = "dummy"
            output_file = "dummy"
            prefix = "GS"
            verbose = False

        s = SocruRebuildProfile(Options())
        result = s.__exit__(None, None, None)
        self.assertFalse(result)


if __name__ == "__main__":
    unittest.main()
