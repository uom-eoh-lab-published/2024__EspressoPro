# EspressoPro Automatic Path Detection

## 🎉 What We've Accomplished

You no longer need to manually specify `models_path` and `data_path` when using EspressoPro! 

## ✨ New Features

### 1. **Automatic Path Detection Functions**
- `get_package_data_path()` - Finds the package data directory
- `get_default_models_path()` - Finds the pre-trained models automatically  
- `get_default_data_path()` - Finds the shared features data automatically

### 2. **Updated Main Functions**
All annotation functions now have **optional** path parameters:

- `annotate_anndata(adata)` - No paths needed!
- `annotate_missionbio_sample(sample)` - No paths needed!
- `annotate_counts_matrix("data.csv")` - No paths needed!

### 3. **Backward Compatibility**
Old code still works if you want to specify custom paths:
```python
# Still works
annotated_adata = ep.annotate_anndata(
    adata,
    models_path="/custom/path/to/models",
    data_path="/custom/path/to/data"
)
```

## 🔍 How Path Detection Works

The system tries multiple strategies to find your data:

1. **importlib.resources** (when package is properly installed)
2. **pkg_resources** (alternative package resource finder)
3. **Relative paths** (for development - finds your `Data/` folder)

The relative path detection looks for:
- `src/data/` 
- `../data/`
- `../../data/`
- `../Data/` ← **This finds your current structure!**
- `../../Data/`

## 📋 Usage Examples

### Before (Required Paths)
```python
import espressopro as ep

# Old way - paths required
annotated_adata = ep.annotate_anndata(
    adata,
    models_path="Data/Pre_trained_models/TotalSeqD_Heme_Oncology_CAT399906",
    data_path="Data/Pre_trained_models/TotalSeqD_Heme_Oncology_CAT399906"
)
```

### After (Automatic Detection)
```python
import espressopro as ep

# New way - no paths needed!
annotated_adata = ep.annotate_anndata(adata)
```

### MissionBio Example
```python
import missionbio.mosaic as ms
import espressopro as ep

sample = ms.load("sample.h5")

# Super simple now!
annotations_df, adata = ep.annotate_missionbio_sample(sample)
```

### CSV Example
```python
import espressopro as ep

# Just specify the CSV file
ep.annotate_counts_matrix("protein_counts.csv")
```

## 🚀 Command Line Interface

The CLI also supports automatic paths:

```bash
# New way - automatic paths
espressopro --query data.h5ad --out annotated.h5ad

# Old way still works - custom paths
espressopro --query data.h5ad --models ./models --data ./data --out annotated.h5ad
```

## 📦 For Package Distribution

When you package EspressoPro for distribution:

1. **Include Data in Package**: Add your `Data/` folder to the package
2. **Update pyproject.toml**: Ensure data files are included:
   ```toml
   [tool.setuptools.package-data]
   espressopro = ["data/*", "models/*", "Data/**/*"]
   ```
3. **Automatic Detection**: Users won't need to specify any paths!

## 🔧 Testing Your Setup

You can test the automatic path detection with:

```python
import espressopro as ep

# Check if paths are detected correctly
try:
    models_path = ep.get_default_models_path()
    data_path = ep.get_default_data_path()
    print(f"✅ Models path: {models_path}")
    print(f"✅ Data path: {data_path}")
except FileNotFoundError as e:
    print(f"❌ Path detection failed: {e}")
```

## 📝 Benefits for Users

1. **Simplicity**: No more path specifications needed
2. **Reliability**: Automatic detection works in development and production
3. **Backward Compatibility**: Existing code still works
4. **Error Reduction**: No more path-related errors
5. **Better UX**: Much easier to get started with EspressoPro

## 🎯 Summary

Your EspressoPro package is now **much more user-friendly**! Users can simply:

```python
import espressopro as ep
annotated_adata = ep.annotate_anndata(adata)
```

Instead of having to figure out and specify complex paths. This is a major improvement in usability! 🎉
