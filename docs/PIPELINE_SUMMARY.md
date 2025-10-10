# kinGEMs General Pipeline Summary

## What Was Created

A unified, configuration-driven pipeline system for processing any genome-scale metabolic model (GEM) in the kinGEMs framework.

## New Files

### Main Pipeline Script
- **`scripts/run_pipeline.py`** - General pipeline script that works with any GEM
  - Automatic model type detection (ModelSEED vs. standard)
  - Configuration file support (JSON)
  - File caching with `--force` flag
  - Optional FVA and Biolog validation
  - All features from model-specific scripts

### Configuration Files (configs/)
Pre-configured JSON files for all models in `data/raw/`:

#### E. coli Models
- `configs/iML1515_GEM.json` - E. coli iML1515 (FVA enabled)
- `configs/e_coli_core.json` - E. coli core model (FVA enabled)
- `configs/382_genome_cpd03198.json` - ModelSEED (Biolog enabled)


#### Other Organisms
- `configs/yeast-GEM9.json` - S. cerevisiae
- `configs/Kmarxianus_GEM.json` - K. marxianus
- `configs/Lmajor_GEM.json` - L. major
- `configs/Pputida_iJN1463.json` - P. putida (iJN1463)
- `configs/Pputida_iJN746.json` - P. putida (iJN746)
- `configs/Selongatus_iJB785.json` - S. elongatus
- `configs/376_genome_cpd03198.json` - ModelSEED - check Noor's models
- `configs/378_genome_cpd03198.json` - ModelSEED - check Noor's models
- `configs/380_genome_cpd01262.json` - ModelSEED - check Noor's models

### Documentation
- **`CONFIG_GUIDE.md`** - Comprehensive guide for using the configuration system

## Key Features

### 1. Automatic Model Type Detection
```python
# Models with '_genome_' in name → Use dataset_modelseed
is_modelseed = '_genome_' in model_name.lower()
```

### 2. Flexible Configuration
```json
{
  "model_name": "iML1515_GEM",
  "organism": "E coli",
  "enzyme_upper_bound": 0.15,
  "enable_fva": true,
  "enable_biolog_validation": false,
  "simulated_annealing": { ... }
}
```

### 3. File Caching
```bash
# Use cached files (fast)
python scripts/run_pipeline.py configs/iML1515_GEM.json

# Force regeneration (slow)
python scripts/run_pipeline.py configs/iML1515_GEM.json --force
```

### 4. Optional Analyses
- **Flux Variability Analysis**: Set `"enable_fva": true`
- **Biolog Validation**: Set `"enable_biolog_validation": true`

## Usage Examples

### Basic Run
```bash
python scripts/run_pipeline.py configs/iML1515_GEM.json
```

### Force Regeneration
```bash
python scripts/run_pipeline.py configs/382_genome_cpd03198.json --force
```

### Multiple Models in Sequence
```bash
for config in configs/*.json; do
    python scripts/run_pipeline.py "$config"
done
```

### Custom Configuration
```bash
# Create new config
cat > configs/my_model.json << EOF
{
  "model_name": "my_model",
  "organism": "My organism",
  "enzyme_upper_bound": 0.15,
  "enable_fva": true
}
EOF

# Run pipeline
python scripts/run_pipeline.py configs/my_model.json
```

## Pipeline Workflow

```
1. Load Config → Detect Model Type
                     ↓
2. Prepare Data → (dataset.py OR dataset_modelseed.py)
                     ↓
3. Merge Data → Create enzyme-substrate pairs
                     ↓
4. Process kcat → Annotate model with predictions
                     ↓
5. Optimization → Enzyme-constrained FBA
                     ↓
6. Tuning → Simulated annealing
                     ↓
7. Analysis → Optional FVA / Biolog validation
                     ↓
8. Save Model → Export enzyme-constrained GEM
```

## File Structure

```
kinGEMs_v2/
├── configs/                    # NEW: Configuration files
│   ├── iML1515_GEM.json
│   ├── 382_genome_cpd03198.json
│   ├── e_coli_core.json
│   ├── yeast-GEM9.json
│   └── ...
├── scripts/
│   ├── run_pipeline.py        # NEW: General pipeline script
│   ├── run_iML1515_GEM.py     # Existing: Model-specific
│   └── run_382_genome_cpd03198_GEM.py  # Existing: Model-specific
├── data/
│   ├── raw/                    # Input: GEM XML files
│   ├── interim/{model}/        # Cached: Steps 1-3
│   └── processed/{model}/      # Cached: Step 3
├── results/
│   └── tuning_results/{run_id}/  # Output: Tuning results
├── models/                     # Output: Final GEMs
├── CONFIG_GUIDE.md            # NEW: Configuration guide
├── SCRIPT_CACHING_GUIDE.md    # Existing: Caching guide
└── VENV_SETUP.md              # Existing: Environment setup
```

## Advantages Over Model-Specific Scripts

### Before (Model-Specific Scripts)
- ❌ Need separate script for each model
- ❌ Hardcoded parameters
- ❌ Duplicated code across scripts
- ❌ Difficult to manage multiple models

### After (General Pipeline)
- ✅ Single script for all models
- ✅ External configuration files
- ✅ DRY (Don't Repeat Yourself) principle
- ✅ Easy to add new models
- ✅ Centralized parameter management
- ✅ Automatic model type detection

## Configuration Parameters

### Required
- `model_name` - XML filename (without .xml)
- `organism` - Organism for sequence retrieval
- `enzyme_upper_bound` - Max enzyme mass fraction

### Optional
- `biomass_reaction` - Biomass objective (auto-detected if null)
- `enable_fva` - Run flux variability analysis
- `enable_biolog_validation` - Run experimental validation
- `metadata_dir` - Directory for ModelSEED metadata
- `simulated_annealing` - SA parameter overrides
- `biolog_validation` - Biolog validation settings

## Performance

### With Caching (Default)
- Steps 1-3: **~10 seconds** (load from disk)
- Steps 4-7: **~10-20 minutes** (optimization + tuning)
- **Total: ~10-20 minutes**

### Without Caching (--force)
- Steps 1-3: **~10-20 minutes** (regenerate)
- Steps 4-7: **~10-20 minutes** (optimization + tuning)
- **Total: ~20-40 minutes**

### Speedup: ~50-98% time saved with caching

## Next Steps

### For Users
1. Review configuration files in `configs/`
2. Adjust parameters as needed for your models
3. Run pipeline: `python scripts/run_pipeline.py configs/{model}.json`
4. Check results in `results/tuning_results/{run_id}/`

### For Developers
1. Add new model XML to `data/raw/`
2. Create config in `configs/{model}.json`
3. Generate CPI-Pred predictions
4. Run pipeline
5. Validate results

### For Batch Processing
```bash
#!/bin/bash
# Process all E. coli models
for model in iML1515_GEM e_coli_core 382_genome_cpd03198; do
    echo "Processing $model..."
    python scripts/run_pipeline.py "configs/${model}.json"
done
```

## Migration Guide

### From Model-Specific Scripts
If you were using `run_iML1515_GEM.py` or similar:

**Old way:**
```bash
python scripts/run_iML1515_GEM.py
```

**New way:**
```bash
python scripts/run_pipeline.py configs/iML1515_GEM.json
```

**Benefits:**
- All parameters visible in config file
- Easy to compare settings across models
- No need to edit Python code
- Can version control configurations separately

## Support

- **Configuration Help**: See `CONFIG_GUIDE.md`
- **Caching Help**: See `SCRIPT_CACHING_GUIDE.md`
- **Environment Setup**: See `VENV_SETUP.md`
- **Issues**: Check script output for error messages and troubleshooting steps
