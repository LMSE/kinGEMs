# kinGEMs

<a target="_blank" href="https://cookiecutter-data-science.drivendata.org/">
    <img src="https://img.shields.io/badge/CCDS-Project%20template-328F97?logo=cookiecutter" />
</a>

A pipeline for automatic reconstruction of enzyme-constrained genome-scale models using CPI-Pred for kinetic parameter annotation.

## Project Organization

```
├── LICENSE            <- Open-source license if one is chosen
├── Makefile           <- Makefile with convenience commands like `make data` or `make train`
├── README.md          <- The top-level README for developers using this project.
├── data
│   ├── external       <- Data from third party sources.
│   ├── interim        <- Intermediate data that has been transformed.
│   ├── processed      <- The final, canonical data sets for modeling.
│   └── raw            <- The original, immutable data dump.
│
├── docs               <- A default mkdocs project; see www.mkdocs.org for details
│
├── models             <- Trained and serialized models, model predictions, or model summaries
│
├── notebooks          <- Jupyter notebooks. Naming convention is a number (for ordering),
│                         the creator's initials, and a short `-` delimited description, e.g.
│                         `1.0-jqp-initial-data-exploration`.
│
├── pyproject.toml     <- Project configuration file with package metadata for
│                         kinGEMs and configuration for tools like black
│
├── references         <- Data dictionaries, manuals, and all other explanatory materials.
│
├── reports            <- Generated analysis as HTML, PDF, LaTeX, etc.
│   └── figures        <- Generated graphics and figures to be used in reporting
│
├── requirements.txt   <- The requirements file for reproducing the analysis environment, e.g.
│                         generated with `pip freeze > requirements.txt`
│
├── setup.cfg          <- Configuration file for flake8
│
└── kinGEMs   <- Source code for use in this project.
    │
    ├── __init__.py             <- Makes kinGEMs a Python module
    │
    ├── config.py               <- Store useful variables and configuration
    │
    ├── dataset.py              <- Scripts to download or generate data
    │
    ├── features.py             <- Code to create features for modeling
    │
    ├── modeling
    │   ├── __init__.py
    │   ├── predict.py          <- Code to run model inference with trained models
    │   └── train.py            <- Code to train models
    │
    └── plots.py                <- Code to create visualizations
```

## Installation

### Prerequisites

- Python 3.9
- [Anaconda](https://www.anaconda.com/products/distribution) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html)
- Git (for cloning the repository)

### Quick Start

We recommend using [mamba](https://github.com/mamba-org/mamba) since the linear programming libraries
require non-python dependencies (i.g. [glpk](https://www.gnu.org/software/glpk/) & glpsol.exe)

```bash
   git clone https://github.com/FILL_ME/kinGEMs.git
   cd kinGEMs
   mamba env create -f enviroment.yml
```

1. **Clone the repository:**

   ```bash
   git clone https://github.com/your-username/kinGEMs.git
   cd kinGEMs
   ```

2. **Create and activate a conda environment:**

   ```bash
   conda create -n kingems python=3.9 -y
   conda activate kingems
   ```

3. **Install dependencies:**

   ```bash
   pip install -r requirements.txt
   ```

4. **Install System-Level Dependencies:**

   **GLPK Solver**:
   The GLPK (GNU Linear Programming Kit) solver is required for optimization operations. There are two options for installation:

   **_Option A: Using Conda (Recommended)_**

   1. Run:

   ```bash
   conda install -c conda-forge glpk
   ```

   2. (On Windows) Ensure the following directory is added to your system PATH environment variable:

   ```bash
   C:\Users\<your-username>\anaconda3\envs\kingems\Library\bin
   ```

   This step is required so that Pyomo can locate the `glpsol.exe` executable.

   **_Option B: Manual Installation (Windows)_**

   1. Get the precompiled binary for Windows from [GLPK for Windows](https://sourceforge.net/projects/winglpk/).
   2. Extract the downloaded files to a directory on your system (e.g., `C:\glpk`).
   3. Add the directory containing `glpsol.exe` to your system's PATH environment variable.

   **_Verify GLPK Installation_**

   After installing GLPK, verify that it is accessible by running:

   ```bash
   glpsol --version
   ```
