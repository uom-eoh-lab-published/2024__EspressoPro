#!/usr/bin/env python3
"""
Test script to check automatic path detection logic.
"""

import sys
import os
from pathlib import Path

# Add src to path
sys.path.insert(0, 'src')

def test_path_detection():
    """Test the path detection functions."""
    try:
        # Import the core module
        from src.core import get_default_models_path, get_default_data_path, _get_data_path_relative
        
        print("🔍 Testing automatic path detection...")
        
        # Test relative path detection (for development)
        print("\n1. Testing relative path detection:")
        rel_path = _get_data_path_relative()
        if rel_path:
            print(f"   ✅ Found relative data path: {rel_path}")
            print(f"   📁 Path exists: {rel_path.exists()}")
        else:
            print("   ❌ No relative data path found")
        
        # Test models path
        print("\n2. Testing models path detection:")
        try:
            models_path = get_default_models_path()
            print(f"   ✅ Default models path: {models_path}")
            print(f"   📁 Path exists: {models_path.exists()}")
            
            # List contents if it exists
            if models_path.exists():
                contents = list(models_path.iterdir())
                print(f"   📂 Contains {len(contents)} items:")
                for item in contents[:5]:  # Show first 5
                    print(f"      - {item.name}")
                if len(contents) > 5:
                    print(f"      ... and {len(contents) - 5} more")
        except Exception as e:
            print(f"   ❌ Models path error: {e}")
        
        # Test data path
        print("\n3. Testing data path detection:")
        try:
            data_path = get_default_data_path()
            print(f"   ✅ Default data path: {data_path}")
            print(f"   📁 Path exists: {data_path.exists()}")
        except Exception as e:
            print(f"   ❌ Data path error: {e}")
        
        return True
        
    except ImportError as e:
        print(f"❌ Import error: {e}")
        return False
    except Exception as e:
        print(f"❌ Unexpected error: {e}")
        return False


if __name__ == "__main__":
    print("🧪 EspressoPro Path Detection Test")
    print("=" * 50)
    
    # Show current working directory
    print(f"📍 Current directory: {os.getcwd()}")
    
    # Check if Data folder exists
    data_folder = Path("Data")
    print(f"📂 Data folder exists: {data_folder.exists()}")
    if data_folder.exists():
        print(f"   Contents: {list(data_folder.iterdir())}")
    
    # Run tests
    success = test_path_detection()
    
    if success:
        print("\n✅ Path detection test completed!")
    else:
        print("\n❌ Path detection test failed!")
    
    print("\n💡 Summary:")
    print("The automatic path detection will:")
    print("1. First try importlib.resources (when package is installed)")
    print("2. Then try pkg_resources (alternative)")
    print("3. Finally try relative paths (for development)")
    print("4. Fall back to the Data/ folder in your current workspace")
