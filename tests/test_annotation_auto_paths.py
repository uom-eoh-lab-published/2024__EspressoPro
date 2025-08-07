#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test core annotation functionality with automatic path detection.
"""

import pytest
import sys
import numpy as np
import pandas as pd
from pathlib import Path

# Add the src directory to the path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

try:
    import anndata as ad
    import scanpy as sc
    ANNDATA_AVAILABLE = True
except ImportError:
    ANNDATA_AVAILABLE = False
    pytest.skip("AnnData not available", allow_module_level=True)


class TestAnnotationWithAutoPaths:
    """Test annotation functions with automatic path detection."""
    
    @pytest.fixture
    def dummy_adata(self):
        """Create dummy AnnData for testing."""
        if not ANNDATA_AVAILABLE:
            pytest.skip("AnnData not available")
        
        # Create dummy protein expression data
        n_cells, n_proteins = 100, 30
        
        protein_names = [
            'CD3', 'CD4', 'CD8', 'CD19', 'CD20', 'CD56', 'CD16', 'CD14', 'CD11c',
            'CD123', 'CD38', 'CD27', 'CD45RA', 'CD45RO', 'CD25', 'CD127', 'FOXP3',
            'CD57', 'CD279', 'CD154', 'CD69', 'CD107a', 'CD62L', 'CCR7', 'CXCR5',
            'CD161', 'CD103', 'CD49d', 'CD29', 'CD44'
        ][:n_proteins]
        
        # Generate realistic expression data
        np.random.seed(42)
        X = np.random.negative_binomial(5, 0.3, size=(n_cells, n_proteins)).astype(float)
        
        # Create AnnData object
        adata = ad.AnnData(
            X=X,
            obs=pd.DataFrame(index=[f"cell_{i}" for i in range(n_cells)]),
            var=pd.DataFrame(index=protein_names)
        )
        
        return adata
    
    def test_import_core_functions(self):
        """Test that core functions can be imported."""
        try:
            from src.core import get_default_models_path, get_default_data_path
            from src.annotation import annotate_anndata
            
            assert callable(get_default_models_path)
            assert callable(get_default_data_path)
            assert callable(annotate_anndata)
            
        except ImportError as e:
            pytest.fail(f"Failed to import core functions: {e}")
    
    def test_annotation_with_auto_paths(self, dummy_adata):
        """Test annotation with automatic path detection."""
        if not ANNDATA_AVAILABLE:
            pytest.skip("AnnData not available")
            
        try:
            from src.annotation import annotate_anndata
            
            # This should work with automatic path detection
            annotated_adata = annotate_anndata(dummy_adata)
            
            # Check that annotation columns were added
            assert hasattr(annotated_adata, 'obs')
            
            # Check for some expected annotation columns
            expected_cols = [
                'CommonBroad.Celltype',
                'CommonSimplified.Celltype', 
                'CommonDetailed.Celltype'
            ]
            
            found_cols = [col for col in expected_cols if col in annotated_adata.obs.columns]
            
            # At least some annotation should have been added
            assert len(found_cols) > 0, f"Should find at least one annotation column. Available columns: {list(annotated_adata.obs.columns)}"
            
        except FileNotFoundError:
            pytest.skip("Models/data not found (expected in some environments)")
        except Exception as e:
            pytest.fail(f"Annotation with auto paths failed: {e}")
    
    def test_custom_paths_still_work(self, dummy_adata):
        """Test that custom paths still work (backward compatibility)."""
        if not ANNDATA_AVAILABLE:
            pytest.skip("AnnData not available")
            
        try:
            from src.annotation import annotate_anndata
            from src.core import get_default_models_path, get_default_data_path
            
            # Get the default paths
            models_path = get_default_models_path()
            data_path = get_default_data_path()
            
            # Use them explicitly (should still work)
            annotated_adata = annotate_anndata(
                dummy_adata,
                models_path=str(models_path),
                data_path=str(data_path)
            )
            
            assert hasattr(annotated_adata, 'obs')
            
        except FileNotFoundError:
            pytest.skip("Models/data not found (expected in some environments)")
        except Exception as e:
            pytest.fail(f"Annotation with custom paths failed: {e}")


def create_simple_test_data():
    """Create simple test data without external dependencies."""
    if not ANNDATA_AVAILABLE:
        return None
        
    n_cells, n_genes = 50, 20
    X = np.random.rand(n_cells, n_genes)
    
    adata = ad.AnnData(
        X=X,
        obs=pd.DataFrame(index=[f"cell_{i}" for i in range(n_cells)]),
        var=pd.DataFrame(index=[f"gene_{i}" for i in range(n_genes)])
    )
    
    return adata


def test_annotation_pipeline_integration():
    """Integration test for the full annotation pipeline."""
    if not ANNDATA_AVAILABLE:
        pytest.skip("AnnData not available")
    
    print("\n🧬 Testing EspressoPro annotation pipeline with auto paths...")
    
    try:
        # Create test data
        adata = create_simple_test_data()
        if adata is None:
            pytest.skip("Could not create test data")
        
        print(f"📊 Created test data: {adata.n_obs} cells × {adata.n_vars} features")
        
        # Test path detection
        from src.core import get_default_models_path, get_default_data_path
        
        models_path = get_default_models_path()
        data_path = get_default_data_path()
        
        print(f"🤖 Models path: {models_path}")
        print(f"📂 Data path: {data_path}")
        
        assert models_path.exists(), f"Models path should exist: {models_path}"
        assert data_path.exists(), f"Data path should exist: {data_path}"
        
        print("✅ Path detection successful!")
        
    except Exception as e:
        print(f"❌ Integration test failed: {e}")
        pytest.fail(f"Integration test failed: {e}")


if __name__ == "__main__":
    # Run tests when script is executed directly
    print("🧪 Testing EspressoPro Annotation with Auto Paths")
    print("=" * 55)
    
    # Run the integration test
    test_annotation_pipeline_integration()
    
    print("\n🎉 Annotation tests completed!")
