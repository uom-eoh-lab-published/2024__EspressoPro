# -*- coding: utf-8 -*-
"""
Marker-based cell type detection and annotation.
"""

from __future__ import annotations
from typing import List, Optional, Union

import numpy as np
import matplotlib.pyplot as plt
import scanpy as sc
from anndata import AnnData

from .constants import MAST_POS, MAST_NEG


def _mean_marker_score(a: AnnData, genes: List[str], key: str) -> np.ndarray:
    """
    Calculate mean marker score for given genes.
    
    Parameters
    ----------
    a : AnnData
        AnnData object
    genes : List[str]
        List of gene/marker names
    key : str
        Key to store result in .obs
        
    Returns
    -------
    np.ndarray
        Mean marker scores
    """
    keep = [g for g in genes if g in a.var_names]
    if not keep:
        a.obs[key] = 0.0
        return a.obs[key].values
    X = a[:, keep].X
    X = X.toarray() if hasattr(X, "toarray") else X
    a.obs[key] = X.mean(1)
    return a.obs[key].values


def _otsu_1d(x: np.ndarray) -> float:
    """
    1D Otsu thresholding for automatic threshold selection.
    
    Parameters
    ----------
    x : np.ndarray
        Input values to threshold
        
    Returns
    -------
    float
        Optimal threshold value
    """
    h, bins = np.histogram(x, 256)
    mids = (bins[:-1] + bins[1:]) / 2
    w1 = np.cumsum(h)
    w2 = np.cumsum(h[::-1])[::-1]
    m1 = np.cumsum(h * mids) / (w1 + 1e-9)
    m2 = (np.cumsum((h * mids)[::-1]) / (w2 + 1e-9))[::-1]
    var = w1[:-1] * w2[1:] * (m1[:-1] - m2[1:]) ** 2
    return mids[np.argmax(var)]


def add_signature_annotation(
    adata: AnnData,
    positive_markers: List[str],
    negative_markers: List[str],
    cell_type_label: str,
    *,
    thresh: Optional[float] = None,
    q: float = 0.99,
    k_neighbors: int = 15,
    field_out: str = "CommonDetailed.Celltype.Refined",
    signature_key: Optional[str] = None,
    show_plots: bool = True,
) -> AnnData:
    """
    Append a cell type category to field_out based on marker signature.
    
    Parameters
    ----------
    adata : AnnData
        AnnData object with expression data
    positive_markers : List[str]
        List of positive marker genes/proteins
    negative_markers : List[str]
        List of negative marker genes/proteins
    cell_type_label : str
        Cell type label to assign to high-signature cells
    thresh : float, optional
        Manual threshold for signature. If None, auto-determined
    q : float, default 0.99
        Quantile threshold for auto-threshold if Otsu fails
    k_neighbors : int, default 15
        Number of neighbors for UMAP if needed
    field_out : str, default "CommonDetailed.Celltype.Refined"
        Field to store refined cell type annotations
    signature_key : str, optional
        Key to store signature scores. If None, uses f"{cell_type_label}_signature"
    show_plots : bool, default True
        Whether to show histogram and UMAP plots
        
    Returns
    -------
    AnnData
        Updated AnnData with cell type annotations
    """
    # Set default signature key
    if signature_key is None:
        signature_key = f"{cell_type_label}_signature"
    
    # Calculate cell type signature
    pos_key = f"{cell_type_label}_pos"
    neg_key = f"{cell_type_label}_neg"
    
    pos = _mean_marker_score(adata, positive_markers, pos_key)
    neg = _mean_marker_score(adata, negative_markers, neg_key)
    sig = pos - neg
    adata.obs[signature_key] = sig

    # Determine threshold
    if thresh is None:
        thresh = np.quantile(sig, q) if np.unique(sig).size > 50 else _otsu_1d(sig)
    cell_mask = sig > thresh
    print(f"[{cell_type_label}] auto‑threshold {thresh:.3f} → {cell_mask.sum()} cells")

    # Update categorical field
    if field_out not in adata.obs:
        adata.obs[field_out] = adata.obs.get("CommonDetailed.Celltype", "Unknown")
    adata.obs[field_out] = adata.obs[field_out].astype("category")
    if cell_type_label not in adata.obs[field_out].cat.categories:
        adata.obs[field_out] = adata.obs[field_out].cat.add_categories(cell_type_label)
    adata.obs.loc[cell_mask, field_out] = cell_type_label

    # Remove obsolete color palette
    key = f"{field_out}_colors"
    if key in adata.uns:
        del adata.uns[key]

    if show_plots:
        # Optional signature histogram
        plt.figure(figsize=(6, 3))
        plt.hist(sig, bins=60, color="steelblue", ec="black")
        plt.axvline(thresh, color="red", ls="--")
        plt.title(f"{cell_type_label} signature")
        plt.tight_layout()
        plt.show()

        # Ensure UMAP exists for plotting
        if "X_umap" not in adata.obsm:
            sc.pp.neighbors(adata, n_neighbors=k_neighbors)
            sc.tl.umap(adata)

        # Plot results
        sc.pl.umap(
            adata,
            color=[field_out, signature_key],
            cmap="coolwarm",
            wspace=0.35,
            na_color="lightgrey",
            show=True,
        )
    
    return adata


def add_mast_annotation(
    adata: AnnData,
    *,
    thresh: Optional[float] = None,
    q: float = 0.99,
    k_neighbors: int = 15,
    field_out: str = "CommonDetailed.Celltype.Refined",
    show_plots: bool = True,
) -> AnnData:
    """
    Append a "Mast" category to field_out based on marker signature.
    
    Convenience wrapper around add_signature_annotation for mast cells.
    
    Parameters
    ----------
    adata : AnnData
        AnnData object with protein expression data
    thresh : float, optional
        Manual threshold for mast cell signature. If None, auto-determined
    q : float, default 0.99
        Quantile threshold for auto-threshold if Otsu fails
    k_neighbors : int, default 15
        Number of neighbors for UMAP if needed
    field_out : str, default "CommonDetailed.Celltype.Refined"
        Field to store refined cell type annotations
    show_plots : bool, default True
        Whether to show histogram and UMAP plots
        
    Returns
    -------
    AnnData
        Updated AnnData with mast cell annotations
    """
    return add_signature_annotation(
        adata=adata,
        positive_markers=MAST_POS,
        negative_markers=MAST_NEG,
        cell_type_label="Mast",
        thresh=thresh,
        q=q,
        k_neighbors=k_neighbors,
        field_out=field_out,
        show_plots=show_plots,
    )
