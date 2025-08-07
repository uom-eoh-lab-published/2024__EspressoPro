# EspressoPro

[![PyPI version](https://badge.fury.io/py/espressopro.svg)](https://badge.fury.io/py/espressopro)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/release/python-380/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**EspressoPro** is a modular Python package for automated cell type annotation of single-cell protein data. It provides end-to-end machine learning pipelines with hierarchical ensemble voting and specialized integration for MissionBio platform data.

## Features

🧬 **Multi-level Annotation**: Hierarchical cell type annotation (Broad → Simplified → Detailed)  
🤖 **Machine Learning**: Ensemble voting with multiple atlas references  
📊 **Localization-aware**: Best atlas selection based on spatial coherence  
🔬 **MissionBio Integration**: Native support for MissionBio Mosaic platform  
🎯 **Specialized Detection**: Mast cell detection using marker signatures  
📈 **Analysis Tools**: Cluster-celltype relationship analysis and visualization  

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

## Quick Start

### Basic Usage

```python
import EspressoPro as ep
import scanpy as sc

# Load your single-cell data
adata = sc.read_h5ad("your_data.h5ad")

# Annotate cell types
annotated_adata = ep.annotate_query(
    adata, 
    models_path="/path/to/models",
    data_path="/path/to/shared_features"
)

# Access results
print(annotated_adata.obs['CommonDetailed.Celltype'].value_counts())
```

### MissionBio Integration

```python
import missionbio.mosaic as ms
import espressopro as ep

# Load MissionBio sample
sample = ms.load("sample.h5")

# Annotate with EspressoPro
annotations_df, adata = ep.annotate_missionbio_sample(
    sample,
    models_path="/path/to/models", 
    data_path="/path/to/data"
)

# View results
print(annotations_df.head())
```

### From CSV Count Matrix

```python
import espressopro as ep

# Annotate directly from CSV
csv_out, umap_png = ep.annotate_counts_matrix(
    csv_in="protein_counts.csv",
    models_path="/path/to/models",
    data_path="/path/to/data"
)
```

## Core Functions

- **`annotate_query()`**: Main annotation pipeline
- **`annotate_missionbio_sample()`**: MissionBio-specific workflow  
- **`annotate_counts_matrix()`**: Direct CSV processing
- **`load_models()`**: Model loading utilities
- **`add_mast_annotation()`**: Mast cell detection
- **`suggest_cluster_celltype_identity()`**: Cluster analysis

## Package Structure

```
espressopro/
├── core.py              # Model loading & utilities
├── prediction.py        # Prediction & scoring  
├── annotation.py        # Cell type annotation workflows
├── missionbio.py        # MissionBio integration
├── markers.py           # Marker-based detection
├── analysis.py          # Analysis & visualization
├── constants.py         # Configuration & constants
└── cli.py               # Command-line interface
```

## Requirements

- Python = 3.10.18
- numpy = 1.26.4
- pandas = 1.5.3
- scikit-learn = 1.5.2
- scanpy = 1.11.3
- anndata = 0.11.4
- joblib = 1.5.1
- matplotlib = 3.10.3
- seaborn = 0.13.2
- missionbio.mosaic-base = 3.12.2
- umap-learn = 0.5.6

### Optional Dependencies

- **MissionBio support**: `pip install espressopro[missionbio]`
- **Development tools**: `pip install espressopro[dev]`
- **Documentation**: `pip install espressopro[docs]`

## Command Line Interface

```bash
# Annotate from command line
espressopro --query data.h5ad --models ./models --data ./data --out annotated.h5ad
```

## Documentation

Full documentation is available at [espressopro.readthedocs.io](https://espressopro.readthedocs.io)

## Citation

If you use EspressoPro in your research, please cite:

```bibtex
@article{espressopro2024,
  title={EspressoPro: Automated Cell Type Annotation for Single-Cell Protein Data},
  author={Your Name},
  journal={Your Journal},
  year={2024}
}
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## Support

- 📧 Email: your.email@example.com
- 🐛 Issues: [GitHub Issues](https://github.com/yourusername/espressopro/issues)
- 💬 Discussions: [GitHub Discussions](https://github.com/yourusername/espressopro/discussions)
