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
pip install git+https://github.com/EspressoKris/2024__EspressoPro.git
```

**Dev mode:**
```bash
git clone https://github.com/EspresssoKris/2024__EspressoPro.git
cd EspressoPro
pip install -e ".[dev]"
```

**MissionBio (optional, via conda):**
```bash
conda create -n mosaic -c missionbio -c conda-forge   python=3.10 missionbio.mosaic-base=3.12.2 python-kaleido -y
conda activate mosaic
pip install git+https://github.com/EspressoKris/2024__EspressoPro.git
```

## Quick start

### AnnData
```python
import scanpy as sc
import espressopro as ep

adata = sc.read_h5ad("your_data.h5ad")
ep.annotate_data(adata)  # downloads models on first use

print(adata.obs["Broad.Celltype"].value_counts().head())
print(adata.obs["Detailed.Celltype"].value_counts().head())
```

### MissionBio Sample
```python
import missionbio.mosaic as ms
import espressopro as ep

sample = ms.load("sample.h5")
ep.annotate_data(sample)

# Optional: add marker-based calls
ep.add_mast_annotation(sample)
```

## Main API (high level)

- `
- `annotate_data(obj)` — end-to-end annotation (AnnData or MissionBio Sample)  
- `refine_labels_by_centroid_knn(obj, ...)` — label outlier cleanup  
- `add_mast_annotation(obj)` / `add_signature_annotation(obj, ...)` — marker signatures  
- `score_mixed_clusters(obj, clusters, labels)` — mixedness metrics

Results appear in:
- `Broad.Celltype`, `Simplified.Celltype`, `Detailed.Celltype`  
- Plus per-class `Averaged.*.<Label>.predscore` tracks

## Requirements (core)

Python ≥3.8; key libs: `numpy`, `pandas`, `scipy`, `scikit-learn`, `scanpy`, `anndata`, `gdown`.  
MissionBio features require `missionbio.mosaic` (conda).

## Documentation

- Docs: <https://espressopro.readthedocs.io>  
- Repo: <https://github.com/EspressoKris/2024__EspressoPro>

## Citation
```bibtex
@software{espressopro,
  title   = {EspressoPro: Modular cell type annotation pipeline for PB/BM-MNCs single-cell protein data},
  author  = {Gurashi, Kristian},
  year    = {2025},
  url     = {https://github.com/EspressoKris/2024__EspressoPro},
  version = {1.0.0}
}
```

## License

CC © Kristian Gurashi. See [LICENSE](LICENSE).
