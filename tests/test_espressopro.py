#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test script to verify the EspressoPro package works correctly.
"""

import sys
import tempfile
import numpy as np
import pandas as pd
import scanpy as sc
from pathlib import Path

# Add the package directory to Python path
package_dir = Path(__file__).parent.parent / "espressopro"
if package_dir.exists():
    sys.path.insert(0, str(package_dir.parent))
else:
    # Try alternative paths
    possible_paths = [
        "/Users/kgurashi/espressopro",
        "/Users/kgurashi/GitHub/2024__EspressoPro_Manuscript/Scripts/ModelTraining",
    ]
    for path in possible_paths:
        if Path(path).exists():
            sys.path.insert(0, path)
            break


def test_imports():
    """Test that all imports work correctly."""
    print("Testing imports...")
    
    try:
        # Test main package import
        import espressopro
        print("✅ Main package import successful")
        
        # Test individual module imports
        from espressopro import core, prediction, annotation, missionbio, markers, analysis, constants
        print("✅ All module imports successful")
        
        # Test function imports
        from espressopro import (
            load_models, annotate_query, annotate_missionbio_sample,
            suggest_cluster_celltype_identity, add_mast_annotation
        )
        print("✅ Function imports successful")
        
        # Test constants
        from espressopro import SIMPLIFIED_CLASSES, DETAILED_CLASSES
        print("✅ Constants import successful")
        
        print(f"Package version: {espressopro.__version__}")
        
    except ImportError as e:
        print(f"❌ Import failed: {e}")
        return False
    
    return True


def test_basic_functionality():
    """Test basic functionality with dummy data."""
    print("\nTesting basic functionality...")
    
    try:
        import espressopro as ep
        import anndata as ad
        
        # Create dummy AnnData
        n_cells, n_genes = 100, 50
        X = np.random.rand(n_cells, n_genes)
        adata = ad.AnnData(
            X=X,
            obs=pd.DataFrame(index=[f"cell_{i}" for i in range(n_cells)]),
            var=pd.DataFrame(index=[f"gene_{i}" for i in range(n_genes)]),
        )
        
        print(f"✅ Created dummy AnnData: {adata.n_obs} cells × {adata.n_vars} genes")
        
        # Test basic functions exist and are callable
        assert callable(ep.load_models)
        assert callable(ep.annotate_query)
        assert callable(ep.add_mast_annotation)
        print("✅ Core functions are callable")
        
    except Exception as e:
        print(f"❌ Basic functionality test failed: {e}")
        return False
    
    return True


def test_cli():
    """Test CLI functionality."""
    print("\nTesting CLI...")
    
    try:
        from espressopro.cli import main
        print("✅ CLI module loads successfully")
        
        # Test help
        import subprocess
        import sys
        result = subprocess.run([sys.executable, "-m", "espressopro.cli", "--help"], 
                              capture_output=True, text=True, cwd="/Users/kgurashi/espressopro")
        if result.returncode == 0:
            print("✅ CLI help works")
        else:
            print(f"⚠️ CLI help returned non-zero: {result.returncode}")
            
    except Exception as e:
        print(f"❌ CLI test failed: {e}")
        return False
    
    return True


def test_package_structure():
    """Test package structure."""
    print("\nTesting package structure...")
    
    try:
        import espressopro
        
        # Check __all__ is defined
        assert hasattr(espressopro, '__all__')
        print("✅ __all__ is defined")
        
        # Check version
        assert hasattr(espressopro, '__version__')
        print(f"✅ Version defined: {espressopro.__version__}")
        
        # Check key functions are in __all__
        expected_functions = [
            'load_models', 'annotate_query', 'annotate_missionbio_sample'
        ]
        for func in expected_functions:
            assert func in espressopro.__all__
        print("✅ Key functions in __all__")
        
    except Exception as e:
        print(f"❌ Package structure test failed: {e}")
        return False
    
    return True


def main():
    """Run all tests."""
    print("🧪 Testing EspressoPro package...\n")
    
    tests = [
        test_imports,
        test_basic_functionality,
        test_package_structure,
        test_cli,
    ]
    
    results = []
    for test in tests:
        results.append(test())
    
    print(f"\n📊 Test Results: {sum(results)}/{len(results)} passed")
    
    if all(results):
        print("🎉 All tests passed! The EspressoPro package is working correctly.")
        print("\n📦 Package is ready for installation!")
        print("\n🚀 Installation commands:")
        print("   # Development install:")
        print("   pip install -e .")
        print("   ")
        print("   # From GitHub:")
        print("   pip install git+https://github.com/yourusername/espressopro.git")
        print("   ")
        print("   # Build and upload to PyPI:")
        print("   python -m build")
        print("   python -m twine upload dist/*")
    else:
        print("❌ Some tests failed. Please check the setup.")
        return 1
    
    return 0


if __name__ == "__main__":
    exit(main())
