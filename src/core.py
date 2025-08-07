# -*- coding: utf-8 -*-
"""
Core utilities for model loading and package data management.
"""

from __future__ import annotations
from pathlib import Path
from collections import defaultdict
from typing import Dict, Mapping, Sequence, Union
import tarfile
import tempfile

import joblib


def download_models():
    """
    Download and extract pre-trained models from Google Drive.
    
    This function downloads the compressed model data and extracts it to the 
    package's Data directory. It's called automatically when models are needed
    but not found locally.
    """
    try:
        import gdown
    except ImportError:
        raise ImportError(
            "gdown is required for automatic model downloading. "
            "Please install it with: pip install gdown"
        )
    
    # Google Drive file ID
    file_id = "1NOsu3hQAgeBQ4CWkqDG_epkWD-2kGqTe"
    
    # Get package root directory (where Data should be extracted)
    package_root = Path(__file__).parent.parent
    data_dir = package_root / "Data"
    
    # Check if models already exist
    models_path = data_dir / "Pre_trained_models" / "TotalSeqD_Heme_Oncology_CAT399906"
    if models_path.exists() and any(models_path.rglob("*_Stacked.joblib")):
        print("[download_models] Models already exist, skipping download")
        return data_dir
    
    print("[download_models] Downloading pre-trained models from Google Drive...")
    print("[download_models] This may take a few minutes (~400MB)...")
    
    try:
        # Create temporary directory for download
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_file = Path(temp_dir) / "models.tar.xz"
            
            # Download using gdown
            print("[download_models] Starting download...")
            gdown.download(id=file_id, output=str(temp_file), quiet=False)
            
            # Verify file was downloaded
            if not temp_file.exists():
                raise RuntimeError("Download failed - file not found after download")
            
            # Ensure Data directory exists
            data_dir.mkdir(exist_ok=True)
            
            # Extract the archive
            print("[download_models] Extracting models...")
            with tarfile.open(temp_file, 'r:xz') as tar:
                # Extract all contents to the package root
                tar.extractall(package_root)
            
            print(f"[download_models] ✅ Models successfully downloaded and extracted to {data_dir}")
            
            # Verify extraction was successful
            if models_path.exists() and any(models_path.rglob("*_Stacked.joblib")):
                print(f"[download_models] ✅ Model verification successful")
                return data_dir
            else:
                raise RuntimeError("Model extraction failed - models not found after extraction")
                
    except Exception as e:
        print(f"[download_models] ❌ Failed to download models: {e}")
        print("[download_models] Please check your internet connection and try again")
        print("[download_models] You may also need to install gdown: pip install gdown")
        raise RuntimeError(f"Model download failed: {e}")


def ensure_models_available():
    """
    Ensure that pre-trained models are available, downloading if necessary.
    
    Returns
    -------
    Path
        Path to the Data directory containing models
    """
    package_root = Path(__file__).parent.parent
    data_dir = package_root / "Data"
    models_path = data_dir / "Pre_trained_models" / "TotalSeqD_Heme_Oncology_CAT399906"
    
    # Check if models exist
    if models_path.exists() and any(models_path.rglob("*_Stacked.joblib")):
        return data_dir
    
    # Models don't exist, download them
    return download_models()


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
        # First try to find existing models
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
            
            # Models don't exist, download them
            print("[get_default_models_path] Models not found locally, downloading...")
            data_dir = ensure_models_available()
            downloaded_models_path = data_dir / "Pre_trained_models" / "TotalSeqD_Heme_Oncology_CAT399906"
            if downloaded_models_path.exists():
                return downloaded_models_path
            
            raise FileNotFoundError(f"Models directory not found even after download attempt")
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


def get_package_data_path():
    """
    Get path to package data directory using multiple fallback strategies.
    
    Returns
    -------
    Path
        Path to package data directory
        
    Raises
    ------
    FileNotFoundError
        If no data directory can be located
    """
    # Try multiple approaches to find the data
    strategies = [
        _get_data_path_importlib_resources,
        _get_data_path_pkg_resources, 
        _get_data_path_relative
    ]
    
    for strategy in strategies:
        try:
            path = strategy()
            if path and path.exists():
                return path
        except Exception:
            continue
    
    # If all strategies fail, ensure models are available via download
    try:
        package_root = Path(__file__).parent.parent
        data_dir = package_root / "Data"
        if not data_dir.exists():
            print("[get_package_data_path] Data directory not found, downloading models...")
            data_dir = ensure_models_available()
        return data_dir
    except Exception:
        pass
    
    raise FileNotFoundError("Could not locate package data directory")


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
