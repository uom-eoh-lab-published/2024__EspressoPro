# EspressoPro

[![PyPI version](https://badge.fury.io/py/espressopro.svg)](https://badge.fury.io/py/espressopro)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

**EspressoPro** is a modular toolkit for automated cell-type annotation of single-cell protein data. It ships pre-trained models, atlas consensus blending, and hierarchy-aware voting—plus optional MissionBio integration.

## Features

- **Hierarchical labels:** Broad → Simplified → Detailed  
- **Pre-trained models:** Stacked ensembles for multiple atlases (Hao, Zhang, Triana, Luecken)  
- **Atlas consensus:** Weighted blending by atlas consensus predictions, excluding erroneous predictors
- **Signatures:** Built-ins (e.g., Mast) and custom marker signatures  
- **Utilities:** Mixedness scoring, label refinement, quick visual checks  
- **One-liner UX:** Sensible defaults; models auto-download on first run

## Install

**Latest (GitHub):**
```bash
pip install git+https://github.com/uom-eoh-lab-published/2024__EspressoPro.git
```

**Dev mode:**
```bash
git clone https://github.com/uom-eoh-lab-published/2024__EspressoPro.git
cd EspressoPro
pip install -e ".[dev]"
```

**MissionBio (optional, via conda):**
```bash
conda create -n mosaic -c missionbio -c conda-forge   python=3.10 missionbio.mosaic-base=3.12.2 python-kaleido -y
conda activate mosaic
pip install git+https://github.com/uom-eoh-lab-published/2024__EspressoPro.git
```

## Quick start

### MissionBio Sample
```python
import missionbio.mosaic as ms
import espressopro as ep

# Load and preprocess data
sample = ms.load("sample.h5")  
ep.Normalise_protein_data(sample)  
ep.Scale_protein_data(sample)  

# Dimensionality reduction and clustering
sample.protein.run_pca(
    attribute='Normalized_reads', 
    components=8, show_plot=False, random_state=42, svd_solver='randomized')  
sample.protein.run_umap(
    attribute='pca', 
    random_state=42, n_neighbors=50, min_dist=0.1, spread=8, n_components=2)  
sample.protein.cluster(
    attribute='umap', 
    method='graph-community', k=5, random_state=42)  

# Cell type annotation
sample = ep.generate_predictions(obj=sample)  
sample = ep.annotate_data(obj=sample)  

# Quality control and refinement
sample = ep.mark_small_clusters(sample, "Simplified.Celltype", min_cells=3)
sample = ep.mark_small_clusters(sample, "Detailed.Celltype", min_cells=3)
sample = ep.mark_mixed_clusters(sample, "Simplified.Celltype")
sample = ep.mark_mixed_clusters(sample, "Detailed.Celltype")
sample = ep.refine_labels_by_knn_consensus(sample, label_col='Detailed.Celltype')

# Optional: add marker-based calls
sample = ep.add_mast_annotation(sample)

# Optional: expand cell type labels to clusters
sample = ep.suggest_cluster_celltype_identity(
    sample=sample,
    annotation="Detailed.Celltype_refined_consensus", rewrite=True)
```

**Full tutorial:** [MissionBio_Tapestri.ipynb](https://github.com/uom-eoh-lab-published/2024__EspressoPro/blob/main/tutorials/MissionBio_Tapestri.ipynb)

## Requirements (core)

Python ≥3.8; key libs: `numpy`, `pandas`, `scipy`, `scikit-learn`, `scanpy`, `anndata`, `gdown`.  
MissionBio features require `missionbio.mosaic` (conda).

## Documentation

- Docs: <https://espressopro.readthedocs.io>  
- Repo: <https://github.com/uom-eoh-lab-published/2024__EspressoPro>

## Citation
```bibtex
@software{espressopro,
  title   = {EspressoPro: Modular cell type annotation pipeline for PB/BM-MNCs single-cell protein data},
  author  = {Gurashi, Kristian},
  year    = {2025},
  url     = {https://github.com/uom-eoh-lab-published/2024__EspressoPro},
  version = {1.0.0}
}
```

## License

CC © Kristian Gurashi. See [LICENSE](LICENSE).
