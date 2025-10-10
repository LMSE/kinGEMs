# Virtual Environment Setup Guide

## ✅ Virtual Environment Successfully Created!

Your virtual environment has been set up with all required packages for the kinGEMs project.

### Location
- **Path**: `/project/def-mahadeva/ranaab/kinGEMs_v2/venv`

### Python Version
- Python 3.11.4

## Activation

To activate the virtual environment, run:

```bash
source venv/bin/activate
```

To deactivate:

```bash
deactivate
```

## Installed Packages

All required packages have been installed, including:

### Core Dependencies
- ✅ kinGEMs (0.0.1) - installed in editable mode
- ✅ loguru - logging
- ✅ python-dotenv - environment variables
- ✅ tqdm - progress bars

### Scientific Computing
- ✅ numpy (2.2.2)
- ✅ pandas (2.2.3)
- ✅ scipy (1.15.2)
- ✅ matplotlib (3.10.1)
- ✅ seaborn (0.13.2)

### Bioinformatics
- ✅ bioservices (1.12.1)
- ✅ cobra (0.29.1)
- ✅ biopython (1.84)
- ✅ pubchempy (1.0.5)

### Optimization
- ✅ pyomo (6.9.4)

### Machine Learning (for validation_utils.py)
- ✅ scikit-learn (1.7.1)
- ✅ lightgbm (4.5.0)
- ✅ shap (0.48.0)
- ✅ statsmodels (0.14.5)

### Clustering & Parallel Computing
- ✅ fastcluster (1.3.0)
- ✅ dask (2025.9.1)
- ✅ distributed (2025.9.1)

### Development Tools
- ✅ ipykernel (6.30.1) - Jupyter notebook support
- ✅ requests (2.32.5)

## Testing the Installation

To verify all packages are working correctly:

```bash
source venv/bin/activate
python -c "import kinGEMs; import cobra; import numpy; import pandas; import pyomo; import scipy; import sklearn; import lightgbm; import shap; print('All packages imported successfully')"
```

## Running Scripts

With the virtual environment activated, you can now run any of the scripts:

```bash
source venv/bin/activate
python scripts/run_382_genome_cpd03198_GEM.py
```

Or run Jupyter notebooks:

```bash
source venv/bin/activate
jupyter notebook
```

## Requirements File

The updated `requirements.txt` now includes all necessary packages with proper organization:
- Core dependencies
- Scientific computing packages
- Bioinformatics tools
- Machine learning libraries
- Parallel computing support
- Development tools

All packages were analyzed from the codebase to ensure nothing is missing.

## Notes

- Installation was done with `--no-cache-dir` to save disk space on your HPC system
- The package is installed in editable mode (`-e .`), so changes to the code are immediately reflected
- All dependencies from `validation_utils.py`, `dataset.py`, `optimize.py`, `tuning.py`, and `fva.py` are included
