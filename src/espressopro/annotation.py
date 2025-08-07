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


def Normalise_protein_data(data, inplace: bool = True, axis: int = 0):
    """
    Apply the centered log ratio (CLR) transformation to normalize counts.
    From Muon package (credit: https://github.com/muonlab/muon)

    Args:
        data: Can be either:
            - AnnData object with protein expression counts in .X
            - Raw matrix/DataFrame (numpy array, pandas DataFrame, or sparse matrix)
        inplace: Whether to update data inplace (only applies to AnnData objects)
        axis: Axis across which CLR is performed (0=samples, 1=features)
        
    Returns:
        - If data is AnnData: None (if inplace=True) or AnnData copy (if inplace=False)
        - If data is matrix/DataFrame: CLR-transformed array
    """
    from anndata import AnnData

    if axis not in [0, 1]:
        raise ValueError("Invalid value for `axis` provided. Admissible options are `0` and `1`.")

    # Handle AnnData input (original behavior)
    if isinstance(data, AnnData):
        adata = data
        if not inplace:
            adata = adata.copy()
        x = adata.X
    # Handle raw matrix/DataFrame input (new functionality)
    else:
        if isinstance(data, (pd.DataFrame, np.ndarray)) or isspmatrix(data):
            x = data
        else:
            raise ValueError("Input data must be an AnnData object, numpy array, pandas DataFrame, or sparse matrix")

    # Ensure proper sparse matrix format for the given axis
    if issparse(x) and axis == 0 and not isinstance(x, csc_matrix):
        warn("Input is sparse but not in CSC format. Converting to CSC.")
        x = csc_matrix(x)
    elif issparse(x) and axis == 1 and not isinstance(x, csr_matrix):
        warn("Input is sparse but not in CSR format. Converting to CSR.")
        x = csr_matrix(x)

    # Apply CLR transformation
    if issparse(x):
        x.data /= np.repeat(
            np.exp(np.log1p(x).sum(axis=axis).A / x.shape[axis]), x.getnnz(axis=axis)
        )
        np.log1p(x.data, out=x.data)
    else:
        # Handle pandas DataFrame separately due to keepdims parameter issue
        if isinstance(x, pd.DataFrame):
            # Convert to numpy for CLR calculation, ensuring float dtype
            x_values = x.values.astype(np.float64)
            x_values = np.log1p(
                x_values / np.exp(np.log1p(x_values).sum(axis=axis, keepdims=True) / x_values.shape[axis])
            )
            x = pd.DataFrame(x_values, index=x.index, columns=x.columns)
        else:
            # Standard numpy array handling - ensure float type
            if x.dtype.kind in ['i', 'u']:  # integer types
                x = x.astype(np.float64)
            np.log1p(
                x / np.exp(np.log1p(x).sum(axis=axis, keepdims=True) / x.shape[axis]),
                out=x,
            )

    # Return based on input type
    if isinstance(data, AnnData):
        data.X = x
        return None if inplace else data
    else:
        return x


def Scale_protein_data(data, inplace: bool = True):
    """
    Apply StandardScaler to normalized protein data.

    Parameters
    ----------
    data : AnnData or array-like
        Can be either:
            - AnnData object with normalized protein expression in .X
            - Raw matrix/DataFrame (numpy array, pandas DataFrame, or sparse matrix)
    inplace : bool, default True
        Whether to update data inplace (only applies to AnnData objects)

    Returns
    -------
        - If data is AnnData: None (if inplace=True) or AnnData copy (if inplace=False)
        - If data is matrix/DataFrame: scaled numpy array
    """
    from anndata import AnnData

    # Handle AnnData input
    if isinstance(data, AnnData):
        adata = data
        if not inplace:
            adata = adata.copy()
        x = adata.X
        
        # Convert to numpy array for StandardScaler
        if isinstance(x, pd.DataFrame):
            data_to_scale = x.values
        elif issparse(x):
            data_to_scale = x.toarray()
        else:
            data_to_scale = x
            
        # Apply StandardScaler
        scaler = preprocessing.StandardScaler()
        scaled_data = scaler.fit_transform(data_to_scale)
        
        # Update AnnData
        adata.X = scaled_data
        return None if inplace else adata
        
    # Handle raw matrix/DataFrame input
    else:
        if isinstance(data, (pd.DataFrame, np.ndarray)) or isspmatrix(data):
            # Convert to numpy array for StandardScaler (handles DataFrame column name issues)
            if isinstance(data, pd.DataFrame):
                data_to_scale = data.values
                # Apply StandardScaler for feature scaling
                scaler = preprocessing.StandardScaler()
                data_scaled_standard = scaler.fit_transform(data_to_scale)
                # Return as DataFrame with original index and columns
                return pd.DataFrame(data_scaled_standard, index=data.index, columns=data.columns)
            elif issparse(data):
                data_to_scale = data.toarray()
            else:
                data_to_scale = data
            
            # Apply StandardScaler for feature scaling (non-DataFrame inputs)
            scaler = preprocessing.StandardScaler()
            data_scaled_standard = scaler.fit_transform(data_to_scale)
            
            return data_scaled_standard
        else:
            raise ValueError("Input data must be an AnnData object, numpy array, pandas DataFrame, or sparse matrix")


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


