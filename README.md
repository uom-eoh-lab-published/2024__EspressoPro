# EspressoPro

[![PyPI version](https://badge.fury.io/py/espressopro.svg)](https://badge.fury.io/py/espressopro)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/release/python-380/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**EspressoPro** is a modular Python package for automated cell type annotation of single-cell protein data. It provides end-to-end machine learning pipelines with hierarchical ensemble voting and specialized integration for MissionBio platform data.

## ✨ Key Features

🚀 **Zero-Configuration**: Automatic path detection - no manual path setup required!  
🧬 **Multi-level Annotation**: Hierarchical cell type annotation (Broad → Simplified → Detailed)  
🤖 **Machine Learning**: Ensemble voting with multiple atlas references (Hao, Zhang, Triana, Luecken)  
📊 **Localization-aware**: Best atlas selection based on spatial coherence  
🔬 **MissionBio Integration**: Native support for MissionBio Mosaic platform  
🎯 **Specialized Detection**: Mast cell detection using marker signatures  
📈 **Analysis Tools**: Cluster-celltype relationship analysis and visualization  
🧪 **Comprehensive Testing**: Professional test suite with multiple runners  
⚡ **Easy to Use**: Simple one-line annotations with sensible defaults  

## Installation

### From GitHub (recommended for latest features)

```bash
pip install git+https://github.com/EspressoKris/EspressoPro.git
```

### From PyPI (when released)

```bash
pip install EspressoPro
```

### Development Installation

```bash
git clone https://github.com/EspressoKris/EspressoPro.git
cd EspressoPro
pip install -e ".[dev]"
```

## 🚀 Quick Start

### ⚡ Super Simple Usage (Recommended)

```python
import espressopro as ep
import scanpy as sc

# Load your single-cell data
adata = sc.read_h5ad("your_data.h5ad")

# That's it! One line annotation with automatic path detection
annotated_adata = ep.annotate_anndata(adata)

# Access results at multiple levels
print("Broad cell types:")
print(annotated_adata.obs['CommonBroad.Celltype'].value_counts())

print("\nDetailed cell types:")
print(annotated_adata.obs['CommonDetailed.Celltype'].value_counts())
```

### 🎛️ Advanced Usage (Custom Paths)

```python
# You can still specify custom paths if needed
annotated_adata = ep.annotate_anndata(
    adata, 
    models_path="/path/to/models",
    data_path="/path/to/shared_features"
)
```

### 🔬 MissionBio Integration

```python
import missionbio.mosaic as ms
import espressopro as ep

# Load MissionBio sample
sample = ms.load("sample.h5")

# One-line annotation with automatic paths!
annotations_df, adata = ep.annotate_missionbio_sample(sample)

# View comprehensive results
print(annotations_df[['celltype_broad', 'celltype_simplified', 'celltype_detailed']].head())

# Optional: Add custom signatures
custom_signatures = [
    {
        'positive_markers': ['CD56', 'CD16'],
        'negative_markers': ['CD3', 'CD19'],
        'cell_type_label': 'NK',
        'show_plots': False
    }
]

annotations_df, adata = ep.annotate_missionbio_sample(
    sample, 
    custom_signatures=custom_signatures
)
```

### 📊 CSV Processing

```python
import espressopro as ep

# Direct CSV annotation (adds columns to original file)
ep.annotate_counts_matrix("protein_counts.csv")

# The CSV now contains new annotation columns:
# - CommonBroad.Celltype, CommonBroad.Score, CommonBroad.LowConf
# - CommonSimplified.Celltype, CommonSimplified.Score, CommonSimplified.LowConf  
# - CommonDetailed.Celltype, CommonDetailed.Score, CommonDetailed.LowConf
```

### 🎯 Specialized Cell Detection

```python
# Add mast cell detection to any AnnData
adata_with_mast = ep.add_mast_annotation(adata)

# Add custom marker signatures
adata_custom = ep.add_signature_annotation(
    adata,
    positive_markers=['FOXP3', 'CD25'],
    negative_markers=['CD8'],
    cell_type_label='Treg'
)
```

## 🔧 Core Functions

### 🎯 Main Annotation Functions
- **`annotate_anndata(adata)`**: Main annotation pipeline with automatic path detection
- **`annotate_missionbio_sample(sample)`**: MissionBio-specific workflow with automatic paths  
- **`annotate_counts_matrix("file.csv")`**: Direct CSV processing with automatic paths

### 🔍 Specialized Detection
- **`add_mast_annotation(adata)`**: Mast cell detection using marker signatures
- **`add_signature_annotation(adata, ...)`**: Custom cell type detection
- **`suggest_cluster_celltype_identity(sample, ...)`**: Cluster analysis and suggestions

### ⚙️ Utility Functions
- **`load_models(models_path)`**: Model loading utilities
- **`get_default_models_path()`**: Get package default models path  
- **`get_default_data_path()`**: Get package default data path

### 📊 Analysis Functions
- **`print_cluster_suggestions(suggestions)`**: Display cluster analysis results
- **`visualize_cluster_celltype_frequencies(pivot_df)`**: Heatmap visualization

## 📦 Package Structure

```
EspressoPro/
├── src/
│   ├── __init__.py          # Main package exports
│   ├── core.py              # Model loading & automatic path detection
│   ├── prediction.py        # Prediction & scoring algorithms
│   ├── annotation.py        # Cell type annotation workflows  
│   ├── missionbio.py        # MissionBio platform integration
│   ├── markers.py           # Marker-based cell detection
│   ├── constants.py         # Configuration & cell type mappings
│   └── cli.py               # Command-line interface
├── tests/                   # Comprehensive test suite
│   ├── run_tests.py         # Simple test runner (no dependencies)
│   ├── test_*.py            # Individual test modules
│   └── README.md            # Test documentation
├── Data/                    # Pre-trained models & shared features
│   └── Pre_trained_models/
│       └── TotalSeqD_Heme_Oncology_CAT399906/
│           ├── Hao/         # Hao et al. atlas models
│           ├── Zhang/       # Zhang et al. atlas models
│           ├── Triana/      # Triana et al. atlas models
│           └── Luecken/     # Luecken et al. atlas models
└── examples/                # Usage examples and tutorials
```

## 📋 Requirements

### Core Dependencies
- **Python**: ≥3.8 (tested on 3.8-3.11)
- **numpy**: ≥1.20.0 (tested with 1.26.4)
- **pandas**: ≥1.3.0 (tested with 1.5.3)
- **scikit-learn**: ≥1.0.0 (tested with 1.5.2)
- **scipy**: ≥1.7.0
- **joblib**: ≥1.0.0 (tested with 1.5.1)

### Single-Cell Analysis
- **scanpy**: ≥1.8.0 (tested with 1.11.3)
- **anndata**: ≥0.8.0 (tested with 0.11.4)
- **umap-learn**: ≥0.5.0 (tested with 0.5.6)
- **leidenalg**: ≥0.8.0

### Visualization & Analysis
- **matplotlib**: ≥3.3.0 (tested with 3.10.3)
- **seaborn**: ≥0.11.0 (tested with 0.13.2)

### Machine Learning
- **xgboost**: ≥1.5.0

### Platform Integration
- **missionbio.mosaic-base**: ≥3.12.0 (tested with 3.12.2)

### Optional Dependencies

Install with specific features:

```bash
# Development tools (testing, linting, docs)
pip install espressopro[dev]

# Documentation generation
pip install espressopro[docs]

# All optional dependencies
pip install espressopro[all]
```

## 💻 Command Line Interface

EspressoPro includes a convenient CLI for batch processing:

```bash
# Simple annotation with automatic paths
espressopro --query data.h5ad --out annotated.h5ad

# With custom paths (if needed)
espressopro --query data.h5ad --models ./models --data ./data --out annotated.h5ad

# Get help
espressopro --help
```

## 🧪 Testing

EspressoPro includes a comprehensive test suite:

### Quick Test (No Dependencies)
```bash
# Run core functionality tests
python tests/run_tests.py

# Run simple path detection test
python tests/test_simple.py
```

### Full Test Suite (With pytest)
```bash
# Install pytest
pip install pytest

# Run all tests
pytest tests/ -v

# Run specific test categories  
pytest tests/ -m "not slow" -v
```

### Test Categories
- **Unit tests**: Core functionality, path detection
- **Integration tests**: Full annotation pipelines
- **Platform tests**: MissionBio integration
- **Performance tests**: Large dataset handling

## 🎯 Automatic Path Detection

EspressoPro automatically finds your models and data - no configuration needed!

### How It Works
1. **Package Resources**: Checks installed package data
2. **Relative Paths**: Finds development data structure
3. **Fallback**: Uses provided custom paths

### Manual Path Check
```python
import espressopro as ep

# Check detected paths
print("Models:", ep.get_default_models_path())
print("Data:", ep.get_default_data_path())
```

## 📖 Documentation

- **Full Documentation**: [espressopro.readthedocs.io](https://espressopro.readthedocs.io)
- **API Reference**: Complete function documentation with examples
- **Tutorials**: Step-by-step guides for common workflows
- **Test Documentation**: [`tests/README.md`](tests/README.md)

## 🔬 Scientific Background

EspressoPro implements hierarchical cell type annotation using multiple reference atlases:

### Reference Atlases
- **Hao et al.**: CITE-seq reference with comprehensive immune cell types
- **Zhang et al.**: Multi-modal single-cell reference  
- **Triana et al.**: Specialized hematopoietic reference
- **Luecken et al.**: Benchmarked annotation reference

### Annotation Levels
1. **Broad**: Immature vs Mature cell distinction
2. **Simplified**: 14 major cell type categories
3. **Detailed**: 35+ specific cell subtypes

### Machine Learning Pipeline
1. **Model Loading**: Pre-trained classifiers per atlas
2. **Feature Matching**: Automatic protein marker alignment
3. **Ensemble Voting**: Best-localized track selection
4. **Hierarchical Constraints**: Parent-child annotation consistency
5. **Confidence Scoring**: Quality assessment and low-confidence flagging

## 📄 Citation

If you use EspressoPro in your research, please cite:

```bibtex
@software{espressopro2024,
  title={EspressoPro: Automated Cell Type Annotation for Single-Cell Protein Data},
  author={Gurashi, Kristian and [Additional Authors]},
  year={2024},
  url={https://github.com/uom-eoh-lab-published/2024__EspressoPro},
  version={1.0.0}
}
```

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

- **Reference Atlas Authors**: Hao, Zhang, Triana, Luecken et al.
- **MissionBio**: Platform integration and support
- **Single-Cell Community**: scanpy, anndata, and related tools