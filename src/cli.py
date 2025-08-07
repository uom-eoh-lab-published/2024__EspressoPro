#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Command-line interface for EspressoPro package.
"""

import argparse
import scanpy as sc

from . import annotate_query


def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        description="EspressoPro: Automated cell type annotation for single-cell protein data",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic annotation (automatic paths)
  espressopro --query data.h5ad --out annotated.h5ad
  
  # Custom paths
  espressopro --query data.h5ad --models ./models --data ./data --out annotated.h5ad
        """
    )
    
    parser.add_argument("--query", required=True, 
                       help="Path to query AnnData file (.h5ad)")
    parser.add_argument("--models", 
                       help="Path to pre-trained models directory (optional, uses package default)")
    parser.add_argument("--data", 
                       help="Path to shared features data directory (optional, uses package default)")
    parser.add_argument("--out", required=True, 
                       help="Path for output annotated file (.h5ad)")
    
    args = parser.parse_args()

    print("🧬 EspressoPro: Starting cell type annotation...")
    print(f"   Query: {args.query}")
    print(f"   Models: {args.models or 'auto-detect'}")
    print(f"   Data: {args.data or 'auto-detect'}")
    print(f"   Output: {args.out}")

    # Load and annotate
    try:
        adata = sc.read_h5ad(args.query)
        print(f"   Loaded {adata.n_obs} cells × {adata.n_vars} features")
        
        adata = annotate_query(adata, args.models, args.data)
        adata.write_h5ad(args.out)
        
        # Print summary
        print("\n📊 Annotation Summary:")
        if "CommonDetailed.Celltype" in adata.obs:
            detailed_counts = adata.obs["CommonDetailed.Celltype"].value_counts()
            print("   Detailed cell types:")
            for celltype, count in detailed_counts.items():
                print(f"     {celltype}: {count}")
        
        print(f"\n✅ Annotation finished → {args.out}")
        
    except Exception as e:
        print(f"❌ Error during annotation: {e}")
        return 1
    
    return 0


if __name__ == "__main__":
    exit(main())
