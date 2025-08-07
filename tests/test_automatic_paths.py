#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test automatic path detection functionality.
"""

import pytest
import sys
from pathlib import Path

# Add the src directory to the path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

try:
    from src.core import (
        get_package_data_path,
        get_default_models_path, 
        get_default_data_path,
        _get_data_path_relative
    )
except ImportError:
    # Fallback for different import scenarios
    import core
    get_package_data_path = core.get_package_data_path
    get_default_models_path = core.get_default_models_path
    get_default_data_path = core.get_default_data_path
    _get_data_path_relative = core._get_data_path_relative


class TestAutomaticPathDetection:
    """Test automatic path detection functionality."""
    
    def test_relative_path_detection(self):
        """Test that relative path detection finds the Data folder."""
        data_path = _get_data_path_relative()
        assert data_path is not None, "Should find relative data path"
        assert data_path.exists(), f"Data path should exist: {data_path}"
        assert data_path.name in ["Data", "data"], f"Should find Data or data folder, got: {data_path.name}"
    
    def test_package_data_path(self):
        """Test package data path detection."""
        try:
            data_path = get_package_data_path()
            assert data_path is not None, "Package data path should not be None"
            assert data_path.exists(), f"Package data path should exist: {data_path}"
        except FileNotFoundError:
            pytest.skip("Package data path not found (expected in development)")
    
    def test_default_models_path(self):
        """Test default models path detection."""
        try:
            models_path = get_default_models_path()
            assert models_path is not None, "Models path should not be None"
            assert models_path.exists(), f"Models path should exist: {models_path}"
            
            # Check that it contains the expected structure
            expected_atlases = ["Hao", "Zhang", "Triana", "Luecken"]
            existing_atlases = [d.name for d in models_path.iterdir() if d.is_dir() and not d.name.startswith('.')]
            
            # At least one atlas should exist
            assert len(existing_atlases) > 0, f"Should have at least one atlas directory in {models_path}"
            
            # Check if any expected atlases exist
            found_atlases = [atlas for atlas in expected_atlases if atlas in existing_atlases]
            assert len(found_atlases) > 0, f"Should find at least one expected atlas. Expected: {expected_atlases}, Found: {existing_atlases}"
            
        except FileNotFoundError as e:
            pytest.fail(f"Default models path should be found: {e}")
    
    def test_default_data_path(self):
        """Test default data path detection."""
        try:
            data_path = get_default_data_path()
            assert data_path is not None, "Data path should not be None"
            assert data_path.exists(), f"Data path should exist: {data_path}"
            
        except FileNotFoundError as e:
            pytest.fail(f"Default data path should be found: {e}")
    
    def test_data_structure_exists(self):
        """Test that the expected data structure exists."""
        project_root = Path(__file__).parent.parent
        
        # Check Data folder structure
        data_folder = project_root / "Data"
        assert data_folder.exists(), f"Data folder should exist at {data_folder}"
        
        models_folder = data_folder / "Pre_trained_models" / "TotalSeqD_Heme_Oncology_CAT399906"
        assert models_folder.exists(), f"Models folder should exist at {models_folder}"
        
        # Check for atlas directories
        expected_atlases = ["Hao", "Zhang", "Triana", "Luecken"]
        existing_atlases = [d.name for d in models_folder.iterdir() if d.is_dir() and not d.name.startswith('.')]
        
        found_atlases = [atlas for atlas in expected_atlases if atlas in existing_atlases]
        assert len(found_atlases) > 0, f"Should find at least one atlas directory. Expected: {expected_atlases}, Found: {existing_atlases}"


def test_path_detection_integration():
    """Integration test for path detection."""
    try:
        # This should work without any explicit paths
        models_path = get_default_models_path()
        data_path = get_default_data_path()
        
        print(f"✅ Models path detected: {models_path}")
        print(f"✅ Data path detected: {data_path}")
        
        assert models_path.exists()
        assert data_path.exists()
        
    except Exception as e:
        pytest.fail(f"Path detection integration failed: {e}")


if __name__ == "__main__":
    # Run tests when script is executed directly
    print("🧪 Testing EspressoPro Automatic Path Detection")
    print("=" * 50)
    
    test_instance = TestAutomaticPathDetection()
    
    tests = [
        ("Relative path detection", test_instance.test_relative_path_detection),
        ("Package data path", test_instance.test_package_data_path),
        ("Default models path", test_instance.test_default_models_path),
        ("Default data path", test_instance.test_default_data_path),
        ("Data structure exists", test_instance.test_data_structure_exists),
        ("Integration test", test_path_detection_integration),
    ]
    
    passed = 0
    failed = 0
    
    for test_name, test_func in tests:
        try:
            print(f"\n🔍 Testing: {test_name}")
            test_func()
            print(f"✅ PASSED: {test_name}")
            passed += 1
        except Exception as e:
            print(f"❌ FAILED: {test_name} - {e}")
            failed += 1
    
    print(f"\n📊 Results: {passed} passed, {failed} failed")
    
    if failed == 0:
        print("🎉 All path detection tests passed!")
    else:
        print("⚠️ Some tests failed. Check the data structure.")
