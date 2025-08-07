#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test automatic path detection functionality without pytest.
"""

import sys
from pathlib import Path

# Add the src directory to the path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

try:
    import core
    get_package_data_path = core.get_package_data_path
    get_default_models_path = core.get_default_models_path
    get_default_data_path = core.get_default_data_path
    _get_data_path_relative = core._get_data_path_relative
except ImportError as e:
    print(f"❌ Import error: {e}")
    sys.exit(1)


def test_relative_path_detection():
    """Test that relative path detection finds the Data folder."""
    print("🔍 Testing relative path detection...")
    
    data_path = _get_data_path_relative()
    assert data_path is not None, "Should find relative data path"
    assert data_path.exists(), f"Data path should exist: {data_path}"
    assert data_path.name in ["Data", "data"], f"Should find Data or data folder, got: {data_path.name}"
    
    print(f"✅ Relative path detection works: {data_path}")


def test_package_data_path():
    """Test package data path detection."""
    print("🔍 Testing package data path...")
    
    try:
        data_path = get_package_data_path()
        assert data_path is not None, "Package data path should not be None"
        assert data_path.exists(), f"Package data path should exist: {data_path}"
        print(f"✅ Package data path works: {data_path}")
    except FileNotFoundError:
        print("⚠️ Package data path not found (expected in development)")


def test_default_models_path():
    """Test default models path detection."""
    print("🔍 Testing default models path...")
    
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
        
        print(f"✅ Models path works: {models_path}")
        print(f"   Found atlases: {found_atlases}")
        
    except Exception as e:
        raise AssertionError(f"Default models path should be found: {e}")


def test_default_data_path():
    """Test default data path detection."""
    print("🔍 Testing default data path...")
    
    try:
        data_path = get_default_data_path()
        assert data_path is not None, "Data path should not be None"
        assert data_path.exists(), f"Data path should exist: {data_path}"
        
        print(f"✅ Data path works: {data_path}")
        
    except Exception as e:
        raise AssertionError(f"Default data path should be found: {e}")


def test_data_structure_exists():
    """Test that the expected data structure exists."""
    print("🔍 Testing data structure...")
    
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
    
    print(f"✅ Data structure exists: {models_folder}")
    print(f"   Found atlases: {found_atlases}")


def test_path_detection_integration():
    """Integration test for path detection."""
    print("🔍 Testing path detection integration...")
    
    try:
        # This should work without any explicit paths
        models_path = get_default_models_path()
        data_path = get_default_data_path()
        
        print(f"✅ Models path detected: {models_path}")
        print(f"✅ Data path detected: {data_path}")
        
        assert models_path.exists()
        assert data_path.exists()
        
        print("✅ Integration test passed!")
        
    except Exception as e:
        raise AssertionError(f"Path detection integration failed: {e}")


def main():
    """Run all tests."""
    print("🧪 Testing EspressoPro Automatic Path Detection")
    print("=" * 50)
    
    tests = [
        ("Relative path detection", test_relative_path_detection),
        ("Package data path", test_package_data_path),
        ("Default models path", test_default_models_path),
        ("Default data path", test_default_data_path),
        ("Data structure exists", test_data_structure_exists),
        ("Integration test", test_path_detection_integration),
    ]
    
    passed = 0
    failed = 0
    
    for test_name, test_func in tests:
        try:
            print(f"\n🔄 Running: {test_name}")
            test_func()
            passed += 1
        except Exception as e:
            print(f"❌ FAILED: {test_name} - {e}")
            failed += 1
    
    print(f"\n📊 Results: {passed} passed, {failed} failed")
    
    if failed == 0:
        print("🎉 All path detection tests passed!")
        return True
    else:
        print("⚠️ Some tests failed. Check the data structure.")
        return False


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
