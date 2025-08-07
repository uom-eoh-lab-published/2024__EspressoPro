# -*- coding: utf-8 -*-
"""
EspressoPro Package - Modular cell type annotation pipeline

This package provides end-to-end helpers for:
1. Loading pre-trained models
2. Scoring query AnnData objects  
3. Picking the most localised score tracks per label
4. Hierarchical ensemble voting (Broad → Simplified → Detailed)
5. MissionBio integration
6. Analysis and visualization tools
"""

from __future__ import annotations
from typing import Optional

# Core functionality
from .core import (
    load_models, 
    get_package_data_path, 
    get_default_models_path, 
    get_default_data_path,
    download_models,
    ensure_models_available
)
from .prediction import (
    stack_prediction, 
    generate_predictions, 
    add_best_localised_tracks
)
from .annotation import (
    voting_annotator,
    Broad_Annotation,
    Simplified_Annotation, 
    Detailed_Annotation,
    annotate_anndata,
    annotate_counts_matrix,
    Normalise_protein_data,
    Scale_protein_data
)
from .missionbio import (
    annotate_missionbio_sample,
    suggest_cluster_celltype_identity,
    print_cluster_suggestions,
    visualize_cluster_celltype_frequencies,
    protein_normalization
)
from .markers import add_mast_annotation, add_signature_annotation

# Constants
from .constants import (
    SIMPLIFIED_CLASSES,
    DETAILED_CLASSES,
    SIMPLIFIED_PARENT_MAP,
    DETAILED_PARENT_MAP
)

__version__ = "1.0.1"

__all__ = [
    # Core
    "load_models",
    "get_package_data_path",
    "get_default_models_path", 
    "get_default_data_path",
    "download_models",
    "ensure_models_available",
    
    # Prediction
    "stack_prediction",
    "generate_predictions", 
    "add_best_localised_tracks",
    
    # Annotation
    "voting_annotator",
    "Broad_Annotation",
    "Simplified_Annotation",
    "Detailed_Annotation", 
    "annotate_anndata",
    "annotate_counts_matrix",
    "Normalise_protein_data",
    "Scale_protein_data",
    
    # MissionBio
    "annotate_missionbio_sample",
    "suggest_cluster_celltype_identity",
    "print_cluster_suggestions", 
    "visualize_cluster_celltype_frequencies",
    "protein_normalization",
    
    # Markers
    "add_mast_annotation",
    "add_signature_annotation",
    
    # Constants
    "SIMPLIFIED_CLASSES",
    "DETAILED_CLASSES", 
    "SIMPLIFIED_PARENT_MAP",
    "DETAILED_PARENT_MAP",
]
