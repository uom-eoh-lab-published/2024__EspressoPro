# -*- coding: utf-8 -*-
"""
Core utilities for model loading and package data management.
"""

from __future__ import annotations
from pathlib import Path
from collections import defaultdict
from typing import Dict, Mapping, Sequence, Union

import joblib


def load_models(
    models_path: Union[str, Path],
    model_names: Sequence[str] = ("Hao", "Zhang", "Triana", "Luecken"),
    annotation_depth: Sequence[str] = ("Broad", "Simplified", "Detailed"),
) -> Dict[str, Mapping]:
    """
    Load pre-trained ML models from hierarchical directory structure.
    
    Parameters
    ----------
    models_path : Union[str, Path]
        Root path containing model directories
    model_names : Sequence[str]
        Atlas names to load (e.g., "Hao", "Zhang", "Triana", "Luecken")
    annotation_depth : Sequence[str] 
        Annotation levels to load ("Broad", "Simplified", "Detailed")
        
    Returns
    -------
    Dict[str, Mapping]
        Nested dictionary: {atlas: {depth: {cell_type: {model_type: model}}}}
    """
    models: Dict[str, dict] = defaultdict(lambda: defaultdict(dict))
    root = Path(models_path)
    
    for atlas in model_names:
        atlas_dir = root / atlas / "Models"
        if not atlas_dir.exists():
            print(f"[load_models] ⚠ atlas folder missing: {atlas_dir}")
            continue
            
        for cell_dir in atlas_dir.iterdir():
            if not cell_dir.is_dir():
                continue
                
            try:
                depth, cell = cell_dir.name.split("_", 1)
            except ValueError:
                print(f"[load_models] 🔴 bad dir name: {cell_dir}")
                continue
                
            if depth not in annotation_depth:
                continue
                
            for jf in cell_dir.glob("*.joblib"):
                mtype = jf.stem.split("_")[-1]
                if mtype == "Stacked":
                    models[atlas][depth].setdefault(cell, {})[mtype] = joblib.load(jf)
                    
    return models


def get_package_data_path():
    """
    Get the path to package data directory with fallback strategies.
    
    Returns
    -------
    Path
        Path to package data directory
        
    Raises
    ------
    FileNotFoundError
        If no valid data path is found
    """
    # Try different strategies for finding package data
    strategies = [
        # 1. Try importlib.resources (Python 3.9+)
        lambda: _get_data_path_importlib_resources(),
        
        # 2. Try pkg_resources 
        lambda: _get_data_path_pkg_resources(),
        
        # 3. Try relative to this file (development)
        lambda: _get_data_path_relative(),
    ]
    
    for strategy in strategies:
        try:
            path = strategy()
            if path and path.exists():
                return path
        except Exception:
            continue
            
    raise FileNotFoundError("Could not locate package data directory")


def get_default_models_path():
    """
    Get the default path to pre-trained models.
    
    Returns
    -------
    Path
        Path to pre-trained models directory
        
    Raises
    ------
    FileNotFoundError
        If models directory is not found
    """
    try:
        data_root = get_package_data_path()
        models_path = data_root / "Pre_trained_models" / "TotalSeqD_Heme_Oncology_CAT399906"
        if models_path.exists():
            return models_path
        else:
            # Try alternative structure for development
            package_root = Path(__file__).parent.parent
            alt_models_path = package_root / "Data" / "Pre_trained_models" / "TotalSeqD_Heme_Oncology_CAT399906"
            if alt_models_path.exists():
                return alt_models_path
            raise FileNotFoundError(f"Models directory not found at {models_path} or {alt_models_path}")
    except Exception as e:
        raise FileNotFoundError(f"Could not locate default models path: {e}")


def get_default_data_path():
    """
    Get the default path to shared features data.
    
    Returns
    -------
    Path
        Path to shared features data directory
        
    Raises
    ------
    FileNotFoundError
        If data directory is not found
    """
    try:
        models_root = get_default_models_path()
        return models_root  # The shared features are in the same directory as models
    except Exception as e:
        raise FileNotFoundError(f"Could not locate default data path: {e}")


def _get_data_path_importlib_resources():
    """Try importlib.resources approach."""
    try:
        import importlib.resources as resources
        return resources.files('espressopro') / 'data'
    except (ImportError, AttributeError):
        return None


def _get_data_path_pkg_resources():
    """Try pkg_resources approach."""
    try:
        import pkg_resources
        return Path(pkg_resources.resource_filename('espressopro', 'data'))
    except ImportError:
        return None


def _get_data_path_relative():
    """Try relative to this file (development mode)."""
    current_file = Path(__file__).resolve()
    possible_paths = [
        current_file.parent / 'data',
        current_file.parent.parent / 'data', 
        current_file.parent.parent.parent / 'data',
        current_file.parent.parent / 'Data',  # For development structure
        current_file.parent.parent.parent / 'Data',  # For development structure
    ]
    
    for path in possible_paths:
        if path.exists():
            return path
    return None
