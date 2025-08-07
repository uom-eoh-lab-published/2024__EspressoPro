#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Simple test runner for EspressoPro without external dependencies.
"""

import sys
from pathlib import Path

# Add the src directory to the Python path
src_path = Path(__file__).parent.parent / "src"
sys.path.insert(0, str(src_path))


def run_all_tests():
    """Run all tests in the tests directory."""
    print("🧪 EspressoPro Test Suite")
    print("=" * 50)
    
    test_files = [
        "test_automatic_paths.py",
        "test_annotation_auto_paths.py",
        "test_simple.py"
    ]
    
    tests_dir = Path(__file__).parent
    passed = 0
    failed = 0
    
    for test_file in test_files:
        test_path = tests_dir / test_file
        if not test_path.exists():
            print(f"⚠️  Test file not found: {test_file}")
            continue
            
        print(f"\n📋 Running: {test_file}")
        print("-" * 30)
        
        try:
            # Import and run the test
            exec(open(test_path).read())
            print(f"✅ PASSED: {test_file}")
            passed += 1
        except Exception as e:
            print(f"❌ FAILED: {test_file} - {e}")
            failed += 1
    
    print(f"\n📊 Test Summary: {passed} passed, {failed} failed")
    
    if failed == 0:
        print("🎉 All tests passed!")
        return True
    else:
        print("⚠️ Some tests failed.")
        return False


def test_import_core():
    """Test that core modules can be imported."""
    print("\n🔍 Testing core module imports...")
    
    try:
        sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
        import core
        
        get_package_data_path = core.get_package_data_path
        get_default_models_path = core.get_default_models_path
        get_default_data_path = core.get_default_data_path
        load_models = core.load_models
        
        print("✅ Core functions imported successfully")
        return True
    except ImportError as e:
        print(f"❌ Core import failed: {e}")
        return False


def test_path_detection():
    """Test path detection functionality."""
    print("\n🔍 Testing path detection...")
    
    try:
        sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
        import core
        
        get_default_models_path = core.get_default_models_path
        get_default_data_path = core.get_default_data_path
        
        models_path = get_default_models_path()
        data_path = get_default_data_path()
        
        print(f"🤖 Models path: {models_path}")
        print(f"📂 Data path: {data_path}")
        
        assert models_path.exists(), f"Models path should exist: {models_path}"
        assert data_path.exists(), f"Data path should exist: {data_path}"
        
        print("✅ Path detection successful!")
        return True
        
    except Exception as e:
        print(f"❌ Path detection failed: {e}")
        return False


def test_data_structure():
    """Test that the expected data structure exists."""
    print("\n🔍 Testing data structure...")
    
    project_root = Path(__file__).parent.parent
    
    # Check Data folder
    data_folder = project_root / "Data"
    if not data_folder.exists():
        print(f"❌ Data folder not found: {data_folder}")
        return False
    
    # Check models folder
    models_folder = data_folder / "Pre_trained_models" / "TotalSeqD_Heme_Oncology_CAT399906"
    if not models_folder.exists():
        print(f"❌ Models folder not found: {models_folder}")
        return False
    
    # Check for atlas directories
    expected_atlases = ["Hao", "Zhang", "Triana", "Luecken"]
    existing_atlases = [d.name for d in models_folder.iterdir() if d.is_dir() and not d.name.startswith('.')]
    
    found_atlases = [atlas for atlas in expected_atlases if atlas in existing_atlases]
    
    print(f"📂 Data folder: {data_folder} ✅")
    print(f"🤖 Models folder: {models_folder} ✅")
    print(f"🗂️  Found atlases: {existing_atlases}")
    print(f"✅ Expected atlases found: {found_atlases}")
    
    if len(found_atlases) == 0:
        print(f"❌ No expected atlases found. Expected: {expected_atlases}")
        return False
    
    print("✅ Data structure looks good!")
    return True


def main():
    """Main test runner."""
    all_passed = True
    
    # Run individual tests
    tests = [
        ("Core imports", test_import_core),
        ("Data structure", test_data_structure),
        ("Path detection", test_path_detection),
    ]
    
    for test_name, test_func in tests:
        try:
            if not test_func():
                all_passed = False
        except Exception as e:
            print(f"❌ {test_name} failed with exception: {e}")
            all_passed = False
    
    print("\n" + "=" * 50)
    
    if all_passed:
        print("🎉 All core tests passed! EspressoPro looks ready to go!")
    else:
        print("⚠️ Some tests failed. Please check the setup.")
        return 1
    
    return 0


if __name__ == "__main__":
    exit(main())
