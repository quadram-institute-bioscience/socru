"""Pytest configuration: auto-mark integration tests that need barrnap/BLAST."""
from __future__ import annotations

import pytest

# Test modules that require barrnap and/or BLAST+ on PATH
_INTEGRATION_MODULES = {
    "EndToEnd_test",
    "EndToEndCreate_test",
    "Socru_test",
    "SocruCreate_test",
    "SocruRebuildProfile_test",
    "SocruUpdateProfile_test",
}


def pytest_collection_modifyitems(config, items):
    """Auto-apply the 'integration' marker to tests in integration modules."""
    for item in items:
        module_name = item.module.__name__.rsplit(".", 1)[-1]
        if module_name in _INTEGRATION_MODULES:
            item.add_marker(pytest.mark.integration)
