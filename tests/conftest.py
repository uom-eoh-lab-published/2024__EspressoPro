# -*- coding: utf-8 -*-
"""
pytest configuration for EspressoPro tests.
"""

import pytest
import sys
from pathlib import Path

# Add the src directory to the Python path
src_path = Path(__file__).parent.parent / "src"
sys.path.insert(0, str(src_path))

# Check for optional dependencies
try:
    import numpy
    import pandas
    import anndata
    import scanpy
    FULL_DEPS_AVAILABLE = True
except ImportError:
    FULL_DEPS_AVAILABLE = False

try:
    import missionbio
    MISSIONBIO_AVAILABLE = True
except ImportError:
    MISSIONBIO_AVAILABLE = False


@pytest.fixture(scope="session")
def test_data_available():
    """Check if test data is available."""
    project_root = Path(__file__).parent.parent
    data_folder = project_root / "Data" / "Pre_trained_models" / "TotalSeqD_Heme_Oncology_CAT399906"
    return data_folder.exists()


@pytest.fixture
def skip_if_no_deps():
    """Skip test if dependencies are not available."""
    if not FULL_DEPS_AVAILABLE:
        pytest.skip("Required dependencies (numpy, pandas, anndata, scanpy) not available")


@pytest.fixture
def skip_if_no_missionbio():
    """Skip test if MissionBio is not available."""
    if not MISSIONBIO_AVAILABLE:
        pytest.skip("MissionBio not available")


@pytest.fixture
def skip_if_no_test_data(test_data_available):
    """Skip test if test data is not available."""
    if not test_data_available:
        pytest.skip("Test data not available")


# Pytest markers
def pytest_configure(config):
    """Configure pytest markers."""
    config.addinivalue_line(
        "markers", "slow: marks tests as slow (deselect with '-m \"not slow\"')"
    )
    config.addinivalue_line(
        "markers", "integration: marks tests as integration tests"
    )
    config.addinivalue_line(
        "markers", "requires_data: marks tests that require actual model data"
    )
    config.addinivalue_line(
        "markers", "requires_missionbio: marks tests that require MissionBio"
    )
