#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Example demonstrating automatic path detection in EspressoPro.

This example shows how to use EspressoPro without specifying 
models_path and data_path parameters.
"""

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

# Import EspressoPro
import sys
sys.path.append('src')
import espressopro as ep

def create_dummy_data():
    """Create dummy single-cell protein data for testing."""
    print("Creating dummy single-cell protein data...")
    
    # Create dummy protein expression data
    n_cells, n_proteins = 1000, 50
    
    # Generate some realistic protein names
    protein_names = [
        'CD3', 'CD4', 'CD8', 'CD19', 'CD20', 'CD56', 'CD16', 'CD14', 'CD11c',
        'CD123', 'CD38', 'CD27', 'CD45RA', 'CD45RO', 'CD25', 'CD127', 'FOXP3',
        'CD57', 'CD279', 'CD154', 'CD69', 'CD107a', 'CD62L', 'CCR7', 'CXCR5',
        'CD161', 'CD103', 'CD49d', 'CD29', 'CD44', 'CD95', 'CD134', 'CD137',
        'CD152', 'CD183', 'CD185', 'CD194', 'CD196', 'CD197', 'CD275', 'CD278',
        'CD366', 'CD223', 'CD39', 'CD73', 'CD96', 'CD226', 'CD244', 'CD160'
    ][:n_proteins]
    
    # Generate expression data with some structure
    np.random.seed(42)
    X = np.random.negative_binomial(5, 0.3, size=(n_cells, n_proteins)).astype(float)
    
    # Add some realistic cell type structure
    # Simulate T cells (high CD3, CD4 or CD8)
    t_cell_mask = np.random.choice([True, False], n_cells, p=[0.4, 0.6])
    if 'CD3' in protein_names:
        cd3_idx = protein_names.index('CD3')
        X[t_cell_mask, cd3_idx] *= 3
    
    # Simulate B cells (high CD19, CD20)
    b_cell_mask = np.random.choice([True, False], n_cells, p=[0.2, 0.8]) & ~t_cell_mask
    if 'CD19' in protein_names:
        cd19_idx = protein_names.index('CD19')
        X[b_cell_mask, cd19_idx] *= 4
    
    # Create AnnData object
    adata = ad.AnnData(
        X=X,
        obs=pd.DataFrame(index=[f"cell_{i}" for i in range(n_cells)]),
        var=pd.DataFrame(index=protein_names)
    )
    
    print(f"Created AnnData: {adata.n_obs} cells × {adata.n_vars} proteins")
    return adata


def example_automatic_annotation():
    """Example of automatic annotation with default paths."""
    print("\n" + "="*60)
    print("EXAMPLE 1: Automatic annotation with default paths")
    print("="*60)
    
    # Create dummy data
    adata = create_dummy_data()
    
    try:
        # This is the new way - no need to specify paths!
        print("\nAnnotating with automatic path detection...")
        annotated_adata = ep.annotate_anndata(adata)
        
        print(f"\n✅ Success! Annotated {annotated_adata.n_obs} cells")
        
        # Show results
        if "CommonDetailed.Celltype" in annotated_adata.obs:
            print("\nDetailed cell types found:")
            counts = annotated_adata.obs["CommonDetailed.Celltype"].value_counts()
            for celltype, count in counts.head(10).items():
                print(f"  {celltype}: {count}")
        
        return annotated_adata
        
    except FileNotFoundError as e:
        print(f"❌ Could not find package data: {e}")
        print("This is expected if the models/data aren't properly packaged yet.")
        return None
    except Exception as e:
        print(f"❌ Error during annotation: {e}")
        return None


def example_custom_paths():
    """Example showing backward compatibility with custom paths."""
    print("\n" + "="*60)
    print("EXAMPLE 2: Custom paths (backward compatibility)")
    print("="*60)
    
    # Create dummy data
    adata = create_dummy_data()
    
    try:
        # This still works - specify custom paths
        models_path = "Data/Pre_trained_models/TotalSeqD_Heme_Oncology_CAT399906"
        data_path = "Data/Pre_trained_models/TotalSeqD_Heme_Oncology_CAT399906"
        
        print(f"\nAnnotating with custom paths:")
        print(f"  models_path: {models_path}")
        print(f"  data_path: {data_path}")
        
        annotated_adata = ep.annotate_anndata(
            adata, 
            models_path=models_path,
            data_path=data_path
        )
        
        print(f"\n✅ Success! Annotated {annotated_adata.n_obs} cells")
        return annotated_adata
        
    except FileNotFoundError as e:
        print(f"❌ Could not find data at specified paths: {e}")
        return None
    except Exception as e:
        print(f"❌ Error during annotation: {e}")
        return None


def example_missionbio_automatic():
    """Example of MissionBio annotation with automatic paths."""
    print("\n" + "="*60)
    print("EXAMPLE 3: MissionBio annotation with automatic paths")
    print("="*60)
    
    print("Note: This example requires a MissionBio Sample object.")
    print("The function signature is now:")
    print("  annotate_missionbio_sample(sample)  # No paths needed!")
    print("  # OR with custom paths:")
    print("  annotate_missionbio_sample(sample, models_path='...', data_path='...')")


def example_path_detection():
    """Example showing path detection functions."""
    print("\n" + "="*60)
    print("EXAMPLE 4: Path detection functions")
    print("="*60)
    
    try:
        print("Testing automatic path detection...")
        
        # Test package data path
        try:
            package_data = ep.get_package_data_path()
            print(f"📂 Package data path: {package_data}")
        except Exception as e:
            print(f"❌ Package data path: {e}")
        
        # Test models path
        try:
            models_path = ep.get_default_models_path()
            print(f"🤖 Default models path: {models_path}")
        except Exception as e:
            print(f"❌ Default models path: {e}")
        
        # Test data path  
        try:
            data_path = ep.get_default_data_path()
            print(f"📊 Default data path: {data_path}")
        except Exception as e:
            print(f"❌ Default data path: {e}")
            
    except Exception as e:
        print(f"❌ Error testing paths: {e}")


def main():
    """Run all examples."""
    print("🧬 EspressoPro Automatic Path Detection Examples")
    print("This demonstrates how to use EspressoPro without specifying paths manually.")
    
    # Test path detection first
    example_path_detection()
    
    # Try automatic annotation
    example_automatic_annotation()
    
    # Show custom paths for comparison
    example_custom_paths()
    
    # Show MissionBio example
    example_missionbio_automatic()
    
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    print("✨ EspressoPro now supports automatic path detection!")
    print()
    print("NEW USAGE (recommended):")
    print("  import espressopro as ep")
    print("  annotated_adata = ep.annotate_anndata(adata)")
    print("  # That's it! No paths needed.")
    print()
    print("OLD USAGE (still supported):")
    print("  annotated_adata = ep.annotate_anndata(")
    print("      adata, ")
    print("      models_path='/path/to/models',")
    print("      data_path='/path/to/data'")
    print("  )")
    print()
    print("🚀 This makes EspressoPro much easier to use!")


if __name__ == "__main__":
    main()
