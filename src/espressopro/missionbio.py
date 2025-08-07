# -*- coding: utf-8 -*-
"""
MissionBio-specific integration functions.
"""

from __future__ import annotations
from pathlib import Path
from typing import Tuple, Optional, List, Union, Dict

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import anndata as ad
from anndata import AnnData
from scipy.sparse import csr_matrix, csc_matrix, isspmatrix, issparse
from sklearn import preprocessing
from warnings import warn

from .annotation import annotate_anndata
from .markers import add_mast_annotation, add_signature_annotation


def clr(adata: AnnData, inplace: bool = True, axis: int = 0) -> Union[None, AnnData]:
    """
    Apply the centered log ratio (CLR) transformation
    to normalize counts in adata.X.

    Args:
        adata: AnnData object with protein expression counts.
        inplace: Whether to update adata.X inplace.
        axis: Axis across which CLR is performed.
    """

    if axis not in [0, 1]:
        raise ValueError("Invalid value for `axis` provided. Admissible options are `0` and `1`.")

    if not inplace:
        adata = adata.copy()

    if issparse(adata.X) and axis == 0 and not isinstance(adata.X, csc_matrix):
        warn("adata.X is sparse but not in CSC format. Converting to CSC.")
        x = csc_matrix(adata.X)
    elif issparse(adata.X) and axis == 1 and not isinstance(adata.X, csr_matrix):
        warn("adata.X is sparse but not in CSR format. Converting to CSR.")
        x = csr_matrix(adata.X)
    else:
        x = adata.X

    if issparse(x):
        x.data /= np.repeat(
            np.exp(np.log1p(x).sum(axis=axis).A / x.shape[axis]), x.getnnz(axis=axis)
        )
        np.log1p(x.data, out=x.data)
    else:
        np.log1p(
            x / np.exp(np.log1p(x).sum(axis=axis, keepdims=True) / x.shape[axis]),
            out=x,
        )

    adata.X = x

    return None if inplace else adata


def protein_normalization(x):
    """
    CLR normalize the input data using improved CLR transformation, followed by standard scaling.
    
    This is a fallback implementation when SCUtils is not available.

    Parameters
    ----------
    x : array-like
        The input data to be normalized. It can be a numpy array, 
        a pandas DataFrame, or a sparse matrix.

    Returns
    -------
    numpy.ndarray
        The normalized and scaled data.
    """
    if isinstance(x, (pd.DataFrame, np.ndarray)) or isspmatrix(x):
        # Create temporary AnnData for CLR transformation
        temp_adata = ad.AnnData(X=x)
        
        # Apply improved CLR transformation
        clr(temp_adata, inplace=True, axis=1)  # CLR across features (axis=1)
        
        # Extract the CLR-transformed data
        normalised_counts = temp_adata.X
        
        # Convert to dense array if sparse
        if issparse(normalised_counts):
            normalised_counts = normalised_counts.toarray()
        
        # Apply StandardScaler for feature scaling
        scaler = preprocessing.StandardScaler()
        data_scaled_standard = scaler.fit_transform(normalised_counts)
        
        return data_scaled_standard
    else:
        raise ValueError("Input x must be a numpy array, a pandas DataFrame, or a sparse matrix")


def annotate_missionbio_sample(
    sample,
    *,
    protein_norm_fn=None,
    add_mast: bool = True,
    mast_thresh: Optional[float] = None,
    custom_signatures: Optional[List[Dict]] = None,
    umap_n_neighbors: int = 30,
    umap_min_dist: float = 0.3,
    label_field: str = "CommonDetailed.Celltype.Refined",
) -> Tuple[pd.DataFrame, AnnData]:
    """
    Annotate a MissionBio Sample object using protein read counts.
    
    Models are automatically downloaded on first use (~400MB).
    
    Parameters
    ----------
    sample : missionbio.mosaic.Sample
        MissionBio sample object with protein assay
    protein_norm_fn : callable, optional
        Function for protein normalization. If None, tries SCUtils.Protein_normalization,
        falls back to built-in improved CLR + StandardScaler normalization if SCUtils is not available
    add_mast : bool, default True
        Whether to add mast cell detection
    mast_thresh : float, optional
        Mast cell detection threshold. If None, auto-determined
    custom_signatures : List[Dict], optional
        List of custom signature dictionaries. Each dict should contain:
        - 'positive_markers': List[str] - positive marker genes/proteins
        - 'negative_markers': List[str] - negative marker genes/proteins  
        - 'cell_type_label': str - cell type label to assign
        - 'thresh': float, optional - manual threshold
        - 'show_plots': bool, optional - whether to show plots (default True)
    umap_n_neighbors : int, default 30
        Number of neighbors for UMAP
    umap_min_dist : float, default 0.3
        Minimum distance for UMAP
    label_field : str, default "CommonDetailed.Celltype.Refined"
        Field name for final cell type labels
        
    Returns
    -------
    Tuple[pd.DataFrame, AnnData]
        Annotated DataFrame with barcodes and cell types at all levels, and full AnnData object
        
    Examples
    --------
    >>> import missionbio.mosaic as ms
    >>> sample = ms.load("sample.h5")
    >>> 
    >>> # Basic annotation (uses default package models automatically)
    >>> annotations_df, adata = annotate_missionbio_sample(sample)
    >>> 
    >>> # Annotation with mast cell detection
    >>> annotations_df, adata = annotate_missionbio_sample(sample, add_mast=True)
    >>> 
    >>> # Annotation without mast cell detection
    >>> annotations_df, adata = annotate_missionbio_sample(
    ...     sample, 
    ...     add_mast=False
    ... )
    >>> 
    >>> # Annotation with custom signatures
    >>> custom_sigs = [
    ...     {
    ...         'positive_markers': ['CD56', 'CD16'],
    ...         'negative_markers': ['CD3', 'CD19'],
    ...         'cell_type_label': 'NK',
    ...         'show_plots': False
    ...     },
    ...     {
    ...         'positive_markers': ['FOXP3', 'CD25'],
    ...         'negative_markers': ['CD8'],
    ...         'cell_type_label': 'Treg'
    ...     }
    ... ]
    >>> annotations_df, adata = annotate_missionbio_sample(
    ...     sample, 
    ...     custom_signatures=custom_sigs
    ... )
    """
    # Extract data
    if not hasattr(sample, 'protein') or sample.protein is None:
        raise ValueError("Sample must have a protein assay")
    
    # Ensure models are available (download if needed)
    print("[annotate_missionbio_sample] Checking model availability...")
    try:
        from .core import ensure_models_available
        ensure_models_available()
    except Exception as e:
        print(f"[annotate_missionbio_sample] Warning: Could not ensure models available: {e}")
    
    print("[annotate_missionbio_sample] Extracting protein read counts...")
    
    # Extract read counts following the user's pattern
    df_raw = sample.protein.get_attribute('read_counts', constraint='row+col')
    cols = list(df_raw.columns.values)
    
    # Add labels from existing clustering (if available)
    try:
        existing_labels = sample.protein.get_labels()
        df_raw.loc[:, 'existing_label'] = existing_labels
        print(f"[annotate_missionbio_sample] Found existing labels: {set(existing_labels)}")
    except Exception as e:
        print(f"[annotate_missionbio_sample] No existing labels found: {e}")
        df_raw.loc[:, 'existing_label'] = 'Unknown'
    
    # Filter numeric columns
    numeric_cols = df_raw.select_dtypes(include=["number"]).columns
    non_numeric = [c for c in df_raw.columns if c not in numeric_cols and c != 'existing_label']
    
    if non_numeric:
        print(f"[annotate_missionbio_sample] Dropping non-numeric columns: {non_numeric}")
    
    df_counts = df_raw[numeric_cols].copy()
    
    if df_counts.empty:
        raise ValueError("No numeric protein columns found after filtering")
    
    # Clean data
    df_counts = df_counts.apply(pd.to_numeric, errors="coerce")
    df_counts = df_counts.replace([np.inf, -np.inf], np.nan)
    df_counts = df_counts.fillna(0.0)
    
    print(f"[annotate_missionbio_sample] Using {df_counts.shape[1]} protein markers for {df_counts.shape[0]} cells")
    
    # Create AnnData
    adata = ad.AnnData(
        X=df_counts.values.astype(float),
        obs=pd.DataFrame({
            'barcode': df_counts.index,
            'existing_label': df_raw.loc[df_counts.index, 'existing_label']
        }, index=df_counts.index),
        var=pd.DataFrame(index=df_counts.columns),
    )
    
    # Normalization
    print("[annotate_missionbio_sample] Applying protein normalization...")
    
    if protein_norm_fn is None:
        try:
            import SCUtils
            protein_norm_fn = SCUtils.Protein_normalization
            print("[annotate_missionbio_sample] Using SCUtils.Protein_normalization")
        except ImportError:
            print("[annotate_missionbio_sample] SCUtils not found, using built-in improved CLR + StandardScaler normalization")
            protein_norm_fn = protein_normalization
    
    try:
        adata.X = protein_norm_fn(adata.X)
        print("[annotate_missionbio_sample] Normalization completed successfully")
    except Exception as e:
        print(f"[annotate_missionbio_sample] WARNING: Normalization failed ({e}), using raw data")
        # Keep original data if normalization fails
    
    # Dimensionality reduction
    print("[annotate_missionbio_sample] Computing PCA and neighbors...")
    
    sc.pp.pca(adata, svd_solver="arpack", zero_center=True)
    sc.pp.neighbors(
        adata,
        n_neighbors=min(100, adata.n_obs - 1),
        n_pcs=min(10, adata.obsm["X_pca"].shape[1]),
        metric="cosine",
        use_rep="X_pca",
    )
    
    # ML annotation
    print("[annotate_missionbio_sample] Running ML cell type annotation...")
    
    adata = annotate_anndata(adata)
    
    # Mast cell detection
    if add_mast:
        print("[annotate_missionbio_sample] Adding mast cell detection...")
        adata = add_mast_annotation(adata, thresh=mast_thresh, field_out=label_field)
    else:
        print("[annotate_missionbio_sample] Skipping mast cell detection")
    
    # Custom signature annotations
    if custom_signatures:
        print(f"[annotate_missionbio_sample] Adding {len(custom_signatures)} custom signature(s)...")
        for i, sig_config in enumerate(custom_signatures):
            try:
                # Extract required parameters
                positive_markers = sig_config['positive_markers']
                negative_markers = sig_config['negative_markers']
                cell_type_label = sig_config['cell_type_label']
                
                # Extract optional parameters
                thresh = sig_config.get('thresh', None)
                show_plots = sig_config.get('show_plots', True)
                
                print(f"  Adding {cell_type_label} signature...")
                adata = add_signature_annotation(
                    adata=adata,
                    positive_markers=positive_markers,
                    negative_markers=negative_markers,
                    cell_type_label=cell_type_label,
                    thresh=thresh,
                    field_out=label_field,
                    show_plots=show_plots
                )
            except KeyError as e:
                print(f"  ERROR: Custom signature {i+1} missing required key: {e}")
            except Exception as e:
                print(f"  ERROR: Failed to add custom signature {i+1}: {e}")
    else:
        print("[annotate_missionbio_sample] No custom signatures provided")
    
    # UMAP for visualization
    print("[annotate_missionbio_sample] Computing UMAP...")
    
    if "X_umap" not in adata.obsm or adata.obsm["X_umap"].shape[1] != 2:
        sc.pp.neighbors(adata, n_neighbors=umap_n_neighbors)
        sc.tl.umap(adata, min_dist=umap_min_dist)
    
    # Prepare output
    broad_col = "CommonBroad.Celltype"
    simplified_col = "CommonSimplified.Celltype" 
    detailed_col = label_field if label_field in adata.obs else "CommonDetailed.Celltype"
    
    # Create comprehensive annotations DataFrame
    annotations_df = pd.DataFrame({
        'barcode': adata.obs_names,
        'existing_label': adata.obs['existing_label'].astype(str).values,
        
        # All annotation levels
        'celltype_broad': adata.obs.get(broad_col, 'Unknown').astype(str).values,
        'celltype_simplified': adata.obs.get(simplified_col, 'Unknown').astype(str).values,
        'celltype_detailed': adata.obs.get(detailed_col, 'Unknown').astype(str).values,
        
        # Primary celltype (most detailed available)
        'celltype': adata.obs.get(detailed_col, 
                                 adata.obs.get(simplified_col, 
                                             adata.obs.get(broad_col, 'Unknown'))).astype(str).values,
        
        # Confidence scores
        'confidence_broad': adata.obs.get('CommonBroad.Score', 0),
        'confidence_simplified': adata.obs.get('CommonSimplified.Score', 0),
        'confidence_detailed': adata.obs.get('CommonDetailed.Score', 0),
        
        # Low confidence flags
        'low_conf_broad': adata.obs.get('CommonBroad.LowConf', False),
        'low_conf_simplified': adata.obs.get('CommonSimplified.LowConf', False),
        'low_conf_detailed': adata.obs.get('CommonDetailed.LowConf', False),
    })
    
    # Print summary
    print(f"\n[annotate_missionbio_sample] Annotation Summary:")
    print(f"Total cells annotated: {len(annotations_df)}")
    
    print("\nBroad cell type distribution:")
    for celltype, count in annotations_df['celltype_broad'].value_counts().items():
        print(f"  {celltype}: {count}")
    
    print("\nSimplified cell type distribution:")
    for celltype, count in annotations_df['celltype_simplified'].value_counts().items():
        print(f"  {celltype}: {count}")
        
    print("\nDetailed cell type distribution:")
    for celltype, count in annotations_df['celltype_detailed'].value_counts().items():
        print(f"  {celltype}: {count}")
    
    # Optional UMAP plot
    try:
        sc.pl.umap(
            adata,
            color=[detailed_col, simplified_col, broad_col, 'existing_label'],
            legend_loc="on data",
            legend_fontsize=6,
            wspace=0.4,
            show=True,
        )
    except Exception as e:
        print(f"[annotate_missionbio_sample] Could not create UMAP plot: {e}")
    
    return annotations_df, adata


def suggest_cluster_celltype_identity(sample, celltype_detailed, min_frequency_threshold=0.3):
    """
    Suggest cell type identity for each cluster based on frequency of detailed cell types.
    
    Parameters
    ----------
    sample : missionbio.mosaic.Sample
        A MissionBio sample object containing protein assay
    celltype_detailed : array-like
        Array of detailed cell type annotations (e.g., from annotations_df['celltype_detailed'])
    min_frequency_threshold : float, default 0.3
        Minimum frequency threshold for a cell type to be considered dominant
    
    Returns
    -------
    dict
        Dictionary with cluster labels as keys and suggested cell types as values
    pd.DataFrame
        Detailed frequency table for each cluster
    """
    # Get cluster labels from Sample.protein
    cluster_labels = sample.protein.get_labels()
    
    # Create DataFrame with cluster and cell type information
    df = pd.DataFrame({
        'cluster': cluster_labels,
        'celltype_detailed': celltype_detailed
    })
    
    # Calculate frequency of each cell type within each cluster
    cluster_celltype_counts = df.groupby(['cluster', 'celltype_detailed']).size().reset_index(name='count')
    cluster_totals = df.groupby('cluster').size().reset_index(name='total')
    
    # Merge to get frequencies
    frequency_df = cluster_celltype_counts.merge(cluster_totals, on='cluster')
    frequency_df['frequency'] = frequency_df['count'] / frequency_df['total']
    
    # Create a pivot table for better visualization
    pivot_df = frequency_df.pivot(index='cluster', columns='celltype_detailed', values='frequency').fillna(0)
    
    # Suggest cell type identity for each cluster
    cluster_suggestions = {}
    
    for cluster in pivot_df.index:
        cluster_row = pivot_df.loc[cluster]
        max_frequency = cluster_row.max()
        dominant_celltype = cluster_row.idxmax()
        
        if max_frequency >= min_frequency_threshold:
            cluster_suggestions[cluster] = {
                'suggested_celltype': dominant_celltype,
                'frequency': max_frequency,
                'confidence': 'High' if max_frequency >= 0.7 else 'Medium',
                'total_cells': cluster_totals[cluster_totals['cluster'] == cluster]['total'].iloc[0]
            }
        else:
            # If no single cell type dominates, suggest "Mixed" or the top cell type with low confidence
            cluster_suggestions[cluster] = {
                'suggested_celltype': f"Mixed ({dominant_celltype})",
                'frequency': max_frequency,
                'confidence': 'Low',
                'total_cells': cluster_totals[cluster_totals['cluster'] == cluster]['total'].iloc[0]
            }
    
    return cluster_suggestions, pivot_df


def print_cluster_suggestions(cluster_suggestions):
    """
    Print cluster suggestions in a readable format.
    
    Parameters
    ----------
    cluster_suggestions : dict
        Dictionary from suggest_cluster_celltype_identity()
    """
    print("Cluster Cell Type Suggestions:")
    print("=" * 60)
    
    for cluster, info in cluster_suggestions.items():
        print(f"Cluster {cluster}:")
        print(f"  Suggested: {info['suggested_celltype']}")
        print(f"  Frequency: {info['frequency']:.2%}")
        print(f"  Confidence: {info['confidence']}")
        print(f"  Total cells: {info['total_cells']}")
        print()


def visualize_cluster_celltype_frequencies(pivot_df, figsize=(12, 8)):
    """
    Create a heatmap visualization of cell type frequencies per cluster.
    
    Parameters
    ----------
    pivot_df : pd.DataFrame
        Pivot table from suggest_cluster_celltype_identity()
    figsize : tuple, default (12, 8)
        Figure size for the plot
    """
    plt.figure(figsize=figsize)
    sns.heatmap(pivot_df, 
                annot=True, 
                fmt='.2f', 
                cmap='Blues',
                cbar_kws={'label': 'Frequency'},
                linewidths=0.5)
    
    plt.title('Cell Type Frequency per Cluster', fontsize=14, fontweight='bold')
    plt.xlabel('Cell Type (Detailed)', fontsize=12)
    plt.ylabel('Cluster', fontsize=12)
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.show()
