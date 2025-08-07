# EspressoPro Tests

This directory contains all tests for the EspressoPro package.

## Test Files

### Core Tests
- **`test_automatic_paths.py`** - Tests for automatic path detection functionality
- **`test_annotation_auto_paths.py`** - Tests for annotation functions with auto paths
- **`test_simple.py`** - Simple path detection test
- **`test_espressopro.py`** - Main package functionality tests

### Legacy Tests  
- **`test_paths.py`** - Legacy path tests
- **`example_automatic_paths.py`** - Example demonstrating automatic path usage

### Test Infrastructure
- **`__init__.py`** - Test package initialization
- **`conftest.py`** - pytest configuration and fixtures
- **`run_tests.py`** - Simple test runner without external dependencies

## Running Tests

### Option 1: Simple Test Runner (No Dependencies)
```bash
# Run all core tests
python tests/run_tests.py

# Run individual tests
python tests/test_simple.py
python tests/test_automatic_paths.py
```

### Option 2: Using pytest (If Available)
```bash
# Install pytest first
pip install pytest

# Run all tests
pytest tests/ -v

# Run specific test files
pytest tests/test_automatic_paths.py -v

# Run tests with markers
pytest tests/ -m "not slow" -v
```

### Option 3: Use the Test Script
```bash
# Make executable and run
chmod +x run_tests.sh
./run_tests.sh
```

## Test Categories

### Unit Tests
- Path detection functionality
- Core module imports
- Data structure validation

### Integration Tests  
- Full annotation pipeline
- Automatic path detection with real data
- MissionBio integration (if dependencies available)

### Test Markers
- `@pytest.mark.slow` - Slow tests (can be skipped)
- `@pytest.mark.integration` - Integration tests
- `@pytest.mark.requires_data` - Tests that need actual model data
- `@pytest.mark.requires_missionbio` - Tests that need MissionBio

## Test Dependencies

### Required (Always)
- `pathlib` (built-in)
- `sys` (built-in)

### Optional (For Full Testing)
- `pytest` - Test framework
- `numpy` - Numerical computing
- `pandas` - Data manipulation  
- `anndata` - Single-cell data
- `scanpy` - Single-cell analysis
- `missionbio` - MissionBio platform integration

### Development Setup
```bash
# Install all test dependencies
pip install pytest numpy pandas anndata scanpy

# Or install development dependencies
pip install -e ".[dev]"
```

## Expected Test Results

✅ **All tests should pass if:**
- The `Data/` folder exists with proper structure
- All required model files are present
- Python path configuration is correct

⚠️ **Some tests may be skipped if:**
- Optional dependencies are not installed
- MissionBio package is not available
- Test data is not available

❌ **Tests will fail if:**
- Data structure is incorrect
- Model files are missing
- Import paths are broken

## Test Data Structure

Tests expect this data structure:
```
Data/
└── Pre_trained_models/
    └── TotalSeqD_Heme_Oncology_CAT399906/
        ├── Hao/
        ├── Zhang/
        ├── Triana/
        └── Luecken/
```

## Continuous Integration

For CI/CD pipelines, use:
```bash
# Basic tests (no optional dependencies)
python tests/run_tests.py

# Full test suite (with dependencies)
pip install -e ".[dev]"
pytest tests/ -v --tb=short
```
