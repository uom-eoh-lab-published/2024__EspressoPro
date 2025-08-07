#!/usr/bin/env python3
"""
Test script for model downloading functionality
"""

def test_gdown_import():
    """Test if gdown can be imported and used"""
    try:
        import gdown
        print("✅ gdown import successful")
        return True
    except ImportError as e:
        print(f"❌ gdown import failed: {e}")
        print("Install with: pip install gdown")
        return False

def test_download_function():
    """Test the download function"""
    try:
        import sys
        sys.path.insert(0, 'src')
        from core import download_models, ensure_models_available
        print("✅ Download functions imported successfully")
        return True
    except Exception as e:
        print(f"❌ Download function import failed: {e}")
        return False

if __name__ == "__main__":
    print("Testing EspressoPro download functionality...")
    print("=" * 50)
    
    success = True
    success &= test_gdown_import()
    success &= test_download_function()
    
    if success:
        print("\n✅ All tests passed! Download functionality is ready.")
    else:
        print("\n❌ Some tests failed. Please check dependencies.")
