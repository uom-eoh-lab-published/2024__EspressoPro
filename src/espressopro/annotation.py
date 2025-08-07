# -*- coding: utf-8 -*-
"""
Cell type annotation workflows and voting systems.
"""

from __future__ import annotations
from pathlib import Path
from typing import Dict, List, Tuple, Union, Optional

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from anndata import AnnData
from scipy.sparse import csr_matrix, csc_matrix, isspmatrix, issparse
from sklearn import preprocessing
from warnings import warn

from .core import load_models
from .prediction import generate_predictions, add_best_localised_tracks
from .constants import SIMPLIFIED_CLASSES, DETAILED_CLASSES, SIMPLIFIED_PARENT_MAP, DETAILED_PARENT_MAP


def Normalise_protein_data(adata: AnnData, inplace: bool = True, axis: int = 0) -> Union[None, AnnData]:
    """
    Apply the centered log ratio (CLR) transformation
    to normalize counts in adata.X.
    From Muon package (credit: https://github.com/muonlab/muon)

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


def Scale_protein_data(x):
    '''
    CLR normalize the input data using improved CLR transformation, followed by standard scaling.

    Parameters:
    x (array-like): The input data to be normalized. It can be a numpy array, a pandas DataFrame, or a sparse matrix.

    Returns:
    numpy.ndarray: The normalized and scaled data.

    '''
    if isinstance(x, (pd.DataFrame, np.ndarray)) or isspmatrix(x):
        # Create temporary AnnData for CLR transformation
        temp_adata = ad.AnnData(X=x)
        
        # Apply improved CLR transformation
        Normalise_protein_data(temp_adata, inplace=True, axis=1)  # CLR across features (axis=1)
        
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


def _safe_mean(arrs: List[np.ndarray]) -> np.ndarray:
    """Safely compute mean of arrays."""
    return np.mean(np.stack(arrs, axis=0), axis=0)


def voting_annotator(
    adata: AnnData,
    level_name: str,
    class_to_sources: Dict[str, List[str]],
    parent_field: Optional[str] = None,
    parent_to_subset: Optional[Dict[str, List[str]]] = None,
    conf_threshold: float = 0.75,
) -> None:
    """
    Ensemble voting for cell type annotation.
    
    Parameters
    ----------
    adata : AnnData
        AnnData object with prediction scores
    level_name : str
        Name for this annotation level (e.g., "CommonBroad")
    class_to_sources : Dict[str, List[str]]
        Mapping from cell types to their source score columns
    parent_field : Optional[str]
        Parent annotation field for hierarchical constraints
    parent_to_subset : Optional[Dict[str, List[str]]]
        Mapping from parent labels to allowed child score columns
    conf_threshold : float, default 0.75
        Confidence threshold for low-confidence flagging
    """
    # Aggregate scores for each class
    for out_cls, cols in class_to_sources.items():
        present = [c for c in cols if c in adata.obs.columns]
        if present:
            adata.obs[f"{level_name}.{out_cls}.predscore"] = _safe_mean(
                [adata.obs[c].values for c in present]
            )

    score_cols = [c for c in adata.obs if c.startswith(f"{level_name}.") and c.endswith(".predscore")]
    if not score_cols:
        return
        
    score_mat = adata.obs[score_cols].copy()

    # Apply hierarchical constraints
    if parent_field and parent_to_subset:
        for parent_lbl, keep_cols in parent_to_subset.items():
            mask = adata.obs[parent_field] == parent_lbl
            drop = [c for c in score_cols if c not in keep_cols]
            adata.obs.loc[mask, drop] = 0.0
            score_mat.loc[mask, drop] = 0.0

    # Softmax normalization
    z = np.exp(score_mat.values)
    score_mat.iloc[:, :] = z / z.sum(axis=1, keepdims=True)
    adata.obs[score_cols] = score_mat

    # Winner determination
    winner_idx = score_mat.values.argmax(axis=1)
    winner_cols = np.array(score_cols)[winner_idx]
    winner_scores = score_mat.values.max(axis=1)

    adata.obs[f"{level_name}.Celltype"] = pd.Categorical([c.split(".")[1] for c in winner_cols])
    adata.obs[f"{level_name}.Score"] = winner_scores
    adata.obs[f"{level_name}.LowConf"] = winner_scores < conf_threshold


def Broad_Annotation(adata: AnnData, conf_threshold: float = 0.75) -> AnnData:
    """
    Perform broad-level cell type annotation (Mature vs Immature).
    
    Parameters
    ----------
    adata : AnnData
        AnnData object with BestBroad prediction scores
    conf_threshold : float, default 0.75
        Confidence threshold
        
    Returns
    -------
    AnnData
        Updated AnnData with CommonBroad annotations
    """
    REF = {
        "Mature":   ["BestBroad.Mature.predscore"],
        "Immature": ["BestBroad.Immature.predscore"],
    }
    voting_annotator(adata, "CommonBroad", REF, conf_threshold=conf_threshold)
    
    # 3‑way split based on max
    m = adata.obs["CommonBroad.Mature.predscore"].to_numpy()
    i = adata.obs["CommonBroad.Immature.predscore"].to_numpy()
    triple = np.where(np.abs(m - i) < 1e-4, "Immature/Mature",
                      np.where(i > m, "Immature", "Mature"))
    adata.obs["CommonBroad.Triple"] = pd.Categorical(triple)
    return adata


def Simplified_Annotation(adata: AnnData, conf_threshold: float = 0.75) -> AnnData:
    """
    Perform simplified-level cell type annotation.
    
    Parameters
    ----------
    adata : AnnData
        AnnData object with BestSimplified prediction scores
    conf_threshold : float, default 0.75
        Confidence threshold
        
    Returns
    -------
    AnnData
        Updated AnnData with CommonSimplified annotations
    """
    if "CommonBroad.Celltype" not in adata.obs:
        Broad_Annotation(adata)
    voting_annotator(
        adata, "CommonSimplified", SIMPLIFIED_CLASSES,
        parent_field="CommonBroad.Celltype",
        parent_to_subset=SIMPLIFIED_PARENT_MAP,
        conf_threshold=conf_threshold,
    )
    return adata


def Detailed_Annotation(adata: AnnData, conf_threshold: float = 0.6) -> AnnData:
    """
    Perform detailed-level cell type annotation.
    
    Parameters
    ----------
    adata : AnnData
        AnnData object with BestDetailed prediction scores
    conf_threshold : float, default 0.6
        Confidence threshold
        
    Returns
    -------
    AnnData
        Updated AnnData with CommonDetailed annotations
    """
    if "CommonSimplified.Celltype" not in adata.obs:
        Simplified_Annotation(adata)
    voting_annotator(
        adata, "CommonDetailed", DETAILED_CLASSES,
        parent_field="CommonSimplified.Celltype",
        parent_to_subset=DETAILED_PARENT_MAP,
        conf_threshold=conf_threshold,
    )
    return adata


def annotate_anndata(
    query_adata: AnnData, 
    models_path: Optional[Union[str, Path]] = None, 
    data_path: Optional[Union[str, Path]] = None
) -> AnnData:
    """
    Main annotation pipeline driver.
    
    Models are automatically downloaded on first use (~400MB).
    
    Parameters
    ----------
    query_adata : AnnData
        Query dataset to annotate
    models_path : Optional[Union[str, Path]], default None
        Path to pre-trained models. If None, uses package default.
    data_path : Optional[Union[str, Path]], default None
        Path to shared features data. If None, uses package default.
        
    Returns
    -------
    AnnData
        Fully annotated AnnData object
    """
    from .core import get_default_models_path, get_default_data_path, ensure_models_available
    
    # Ensure models are available (download if needed)
    print("[annotate_anndata] Checking model availability...")
    try:
        ensure_models_available()
    except Exception as e:
        print(f"[annotate_anndata] Warning: Could not ensure models available: {e}")
    
    # Use default paths if not provided
    if models_path is None:
        models_path = str(get_default_models_path())
        print(f"[annotate_anndata] Using default models path: {models_path}")
    
    if data_path is None:
        data_path = str(get_default_data_path())
        print(f"[annotate_anndata] Using default data path: {data_path}")

    models = load_models(models_path)

    X = query_adata.X.A if hasattr(query_adata.X, "A") else query_adata.X
    query_df = pd.DataFrame(X, index=query_adata.obs_names, columns=query_adata.var_names)

    generate_predictions(models, data_path, query_df, query_adata)

    # Pick best‑localised tracks
    add_best_localised_tracks(
        query_adata,
        atlases=("Hao", "Zhang", "Triana", "Luecken"),
        q=0.90, k=15, p=2,
    )

    # Hierarchical annotation
    Broad_Annotation(query_adata)
    Simplified_Annotation(query_adata)
    Detailed_Annotation(query_adata)
    return query_adata


def annotate_counts_matrix(
    csv: Union[str, Path],
    models_path: Optional[Union[str, Path]] = None,
    data_path: Optional[Union[str, Path]] = None,
    *,
    index_col: Union[int, str, None] = 0,
) -> None:
    """
    Annotate a CSV count matrix by adding cell type annotation columns directly to the file.
    
    Models are automatically downloaded on first use (~400MB).
    
    Parameters
    ----------
    csv : Union[str, Path]
        Input CSV file path with count data (will be modified in place)
    models_path : Optional[Union[str, Path]], default None
        Path to pre-trained models. If None, uses package default.
    data_path : Optional[Union[str, Path]], default None
        Path to shared features data. If None, uses package default.
    index_col : Union[int, str, None], default 0
        Column to use as index (barcodes/cell IDs)
        
    Returns
    -------
    None
        Modifies the input CSV file by adding annotation columns
        
    Examples
    --------
    >>> # Adds annotation columns to "protein_counts.csv" (automatic paths)
    >>> annotate_counts_matrix("protein_counts.csv")
    >>> 
    >>> # Adds annotation columns with custom paths
    >>> annotate_counts_matrix(
    ...     "protein_counts.csv", 
    ...     models_path="/path/to/models",
    ...     data_path="/path/to/data"
    ... )
    """
    from .core import get_default_models_path, get_default_data_path, ensure_models_available
    
    # Ensure models are available (download if needed)
    print("[annotate_counts_matrix] Checking model availability...")
    try:
        ensure_models_available()
    except Exception as e:
        print(f"[annotate_counts_matrix] Warning: Could not ensure models available: {e}")
    
    # Use default paths if not provided
    if models_path is None:
        models_path = str(get_default_models_path())
        print(f"[annotate_counts_matrix] Using default models path: {models_path}")
    
    if data_path is None:
        data_path = str(get_default_data_path())
        print(f"[annotate_counts_matrix] Using default data path: {data_path}")

    # Load and clean data
    df_raw = pd.read_csv(csv, index_col=index_col)
    num_cols = df_raw.select_dtypes(include=["number"]).columns
    dropped = [c for c in df_raw.columns if c not in num_cols]
    if dropped:
        print(f"[annotate_counts_matrix] Dropped non‑numeric columns: {dropped}")
    df = df_raw[num_cols].copy()

    if df.empty:
        raise ValueError("No numeric count columns found after filtering.")

    # Clean data
    df = df.apply(pd.to_numeric, errors="coerce")
    df = df.replace([np.inf, -np.inf], np.nan)
    df = df.fillna(0.0)

    print(f"[annotate_counts_matrix] Using {df.shape[1]} features for {df.shape[0]} cells")

    # Create AnnData
    adata = ad.AnnData(
        X=df.values.astype(float),
        obs=pd.DataFrame(index=df.index),
        var=pd.DataFrame(index=df.columns),
    )

    # Normalization - CLR transformation for protein data
    adata.X = Scale_protein_data(adata.X)

    # Basic preprocessing
    sc.pp.pca(adata, svd_solver="arpack", zero_center=True)
    sc.pp.neighbors(
        adata,
        n_neighbors=min(100, adata.n_obs - 1),
        n_pcs=min(10, adata.obsm["X_pca"].shape[1]),
        metric="cosine",
        use_rep="X_pca",
    )

    # Run annotation pipeline
    adata = annotate_anndata(adata, str(models_path), str(data_path))

    # Add annotation columns to original data
    df_annotated = df_raw.copy()
    
    # Add key annotation columns
    annotation_cols = [
        'CommonBroad.Celltype', 'CommonBroad.Score', 'CommonBroad.LowConf',
        'CommonSimplified.Celltype', 'CommonSimplified.Score', 'CommonSimplified.LowConf',
        'CommonDetailed.Celltype', 'CommonDetailed.Score', 'CommonDetailed.LowConf'
    ]
    
    for col in annotation_cols:
        if col in adata.obs.columns:
            df_annotated[col] = adata.obs[col].values
    
    # Save back to the original CSV file
    df_annotated.to_csv(csv, index=True)
    print(f"[annotate_counts_matrix] Added annotation columns to {csv}")


