# -*- coding: utf-8 -*-
"""
Prediction and scoring utilities for ML models.
"""

from __future__ import annotations
from pathlib import Path
from typing import Dict, List, Mapping, Sequence, Tuple, Union, Optional

import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.base import ClassifierMixin
from sklearn.neighbors import KDTree
from anndata import AnnData

from .constants import SIMPLIFIED_CLASSES, _DETAILED_LABELS


def stack_prediction(model: ClassifierMixin, X: pd.DataFrame) -> np.ndarray:
    """
    Generate prediction scores from a stacked model.
    
    Parameters
    ----------
    model : ClassifierMixin
        Trained sklearn classifier
    X : pd.DataFrame
        Feature matrix
        
    Returns
    -------
    np.ndarray
        Prediction scores (binary: positive class probability; multiclass: max probability)
    """
    proba = model.predict_proba(X)
    return proba[:, 1] if proba.shape[1] == 2 else proba.max(axis=1)


def generate_predictions(
    models: Mapping,
    data_path: str,
    query_df: pd.DataFrame,
    adata: AnnData,
) -> AnnData:
    """
    Generate predictions for all loaded models and write scores to .obs.
    
    Parameters
    ----------
    models : Mapping
        Nested model dictionary from load_models()
    data_path : str
        Path to shared features directory
    query_df : pd.DataFrame
        Query feature matrix
    adata : AnnData
        AnnData object to store predictions in
        
    Returns
    -------
    AnnData
        Updated AnnData with prediction scores in .obs
    """
    feature_cache: Dict[str, List[str]] = {}
    
    for atlas, ref in models.items():
        feats_fp = Path(data_path) / atlas / "Shared_Features.csv"
        if atlas not in feature_cache:
            feature_cache[atlas] = pd.read_csv(feats_fp).squeeze("columns").tolist()
        feats = feature_cache[atlas]

        # Handle missing features by adding zeros
        missing = set(feats) - set(query_df.columns)
        if missing:
            query_df = query_df.assign(**{g: 0.0 for g in missing})
        X_sub = query_df[feats]

        # Generate predictions for each cell type and depth
        for depth, cell_dict in ref.items():
            for cell_label, mdl_dict in cell_dict.items():
                if "Stacked" not in mdl_dict:
                    continue
                score = stack_prediction(mdl_dict["Stacked"], X_sub)
                adata.obs[f"{atlas}.{depth}.{cell_label}.predscore"] = score
                
    return adata


def _local_dispersion(vals: np.ndarray, coords: np.ndarray,
                      k: int = 15, p: int = 2, eps: float = 1e-9) -> np.ndarray:
    """
    Calculate distance-weighted absolute deviation from neighbourhood mean.
    
    Parameters
    ----------
    vals : np.ndarray
        Values to analyze for local dispersion
    coords : np.ndarray
        Coordinate matrix (e.g., UMAP coordinates)
    k : int, default 15
        Number of neighbors to consider
    p : int, default 2
        Minkowski metric parameter
    eps : float, default 1e-9
        Small epsilon to avoid division by zero
        
    Returns
    -------
    np.ndarray
        Local dispersion scores
    """
    tree = KDTree(coords, metric="minkowski", p=p)
    dist, idx = tree.query(coords, k=k + 1)            # +1 include self
    weights = 1 / (dist + eps)
    neighb_means = (weights * vals[idx]).sum(1) / weights.sum(1)
    deviation_abs = np.abs(vals[idx] - neighb_means[:, None])
    return (weights * deviation_abs).sum(1) / weights.sum(1)


def _best_localised_score(
    adata: AnnData,
    depth: str,
    label: str,
    atlases: Sequence[str],
    q: float = 0.90,
    k: int = 15,
    p: int = 2,
) -> Tuple[Optional[str], Optional[np.ndarray]]:
    """
    Find the best localised score track for a given cell type across atlases.
    
    Parameters
    ----------
    adata : AnnData
        AnnData object with prediction scores and UMAP coordinates
    depth : str
        Annotation depth ("Broad", "Simplified", "Detailed")
    label : str
        Cell type label to analyze
    atlases : Sequence[str]
        Atlas names to compare
    q : float, default 0.90
        Quantile threshold for high-scoring cells
    k : int, default 15
        Number of neighbors for locality analysis
    p : int, default 2
        Minkowski metric parameter
        
    Returns
    -------
    Tuple[Optional[str], Optional[np.ndarray]]
        Best atlas name and normalized score vector, or (None, None) if none found
    """
    coords = adata.obsm["X_umap"]
    best_med = np.inf
    best_vec = best_atl = None

    for atl in atlases:
        col = f"{atl}.{depth}.{label}.predscore"
        if col not in adata.obs:
            continue

        x = adata.obs[col].to_numpy()
        if x.ptp() == 0:        # constant → useless
            continue
        x = (x - x.min()) / x.ptp()   # min‑max 0‑1

        hi_mask = x >= np.quantile(x, q)
        if hi_mask.sum() < k + 1:
            continue

        disp = _local_dispersion(x, coords, k=k, p=p)
        med = np.median(disp[hi_mask])

        if med < best_med:
            best_med, best_vec, best_atl = med, x, atl

    return best_atl, best_vec


def add_best_localised_tracks(
    adata: AnnData,
    *,
    atlases: Sequence[str],
    q: float = 0.90,
    k: int = 15,
    p: int = 2,
) -> None:
    """
    Create BestBroad.*, BestSimplified.*, BestDetailed.* score columns.
    
    Parameters
    ----------
    adata : AnnData
        AnnData object to add tracks to
    atlases : Sequence[str]
        Atlas names to compare
    q : float, default 0.90
        Quantile threshold for high-scoring cells
    k : int, default 15
        Number of neighbors for locality analysis 
    p : int, default 2
        Minkowski metric parameter
    """
    # Ensure a 2‑D UMAP exists
    if "X_umap" not in adata.obsm or adata.obsm["X_umap"].shape[1] != 2:
        sc.pp.neighbors(adata, n_neighbors=30, use_rep="X_pca")
        sc.tl.umap(adata, n_components=2, min_dist=0.3)

    depth_to_labels = {
        "Broad":      ["Immature", "Mature"],
        "Simplified": list(SIMPLIFIED_CLASSES.keys()),
        "Detailed":   _DETAILED_LABELS,
    }

    for depth, labels in depth_to_labels.items():
        for lbl in labels:
            atl, vec = _best_localised_score(
                adata, depth, lbl, atlases, q=q, k=k, p=p,
            )
            if vec is None:
                continue
            adata.obs[f"Best{depth}.{lbl}.predscore"] = vec
            adata.obs[f"Best{depth}.{lbl}.atlas"] = atl
