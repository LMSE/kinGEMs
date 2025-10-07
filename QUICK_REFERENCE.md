# kinGEMs Pipeline Quick Reference

## Single Model Runs

```bash
# Run a single model (use cached files)
python scripts/run_pipeline.py configs/iML1515_GEM.json

# Run with forced regeneration
python scripts/run_pipeline.py configs/382_genome_cpd03198.json --force
```

## Batch Processing

```bash
# Run all E. coli models
bash scripts/batch_run.sh ecoli

# Run all ModelSEED models
bash scripts/batch_run.sh genome

# Run all standard models
bash scripts/batch_run.sh standard

# Run all models
bash scripts/batch_run.sh all

# Run specific models
bash scripts/batch_run.sh iML1515_GEM yeast-GEM9

# Run with forced regeneration
bash scripts/batch_run.sh ecoli --force
```

## Available Configurations

### E. coli Models
| Config File | Model Type | Features |
|-------------|------------|----------|
| `iML1515_GEM.json` | Standard | FVA enabled |
| `e_coli_core.json` | Standard | FVA enabled |
| `382_genome_cpd03198.json` | ModelSEED | Biolog validation |
| `376_genome_cpd03198.json` | ModelSEED | Basic |
| `378_genome_cpd03198.json` | ModelSEED | Basic |
| `380_genome_cpd01262.json` | ModelSEED | Basic |

### Other Organisms
| Config File | Organism | Features |
|-------------|----------|----------|
| `yeast-GEM9.json` | S. cerevisiae | FVA enabled |
| `Kmarxianus_GEM.json` | K. marxianus | FVA enabled |
| `Lmajor_GEM.json` | L. major | FVA enabled |
| `Pputida_iJN1463.json` | P. putida | FVA enabled |
| `Pputida_iJN746.json` | P. putida | FVA enabled |
| `Selongatus_iJB785.json` | S. elongatus | FVA enabled |

## Model Type Detection

**Automatic detection based on filename:**
- Models with `_genome_` → Use `dataset_modelseed` functions
- Other models → Use standard `dataset` functions

## Configuration Template

```json
{
  "model_name": "my_model",
  "organism": "My organism",
  "biomass_reaction": null,
  "enzyme_upper_bound": 0.15,
  "enable_fva": false,
  "enable_biolog_validation": false,
  "simulated_annealing": {
    "temperature": 1.0,
    "cooling_rate": 0.95,
    "min_temperature": 0.01,
    "max_iterations": 100,
    "max_unchanged_iterations": 5,
    "change_threshold": 0.009,
    "biomass_goal": 0.5
  }
}
```

## Output Locations

### Results
```
results/tuning_results/{run_id}/
├── df_new.csv                 # Enzyme masses
├── kcat_dict.csv              # Tuned kcats
├── final_model_info.csv       # Merged results
├── *_fva_results.csv          # FVA (if enabled)
├── *_fva_*_plot.png           # FVA plots
└── biolog_comparison.*        # Biolog (if enabled)
```

### Models
```
models/
└── {model_name}_{date}_{random}.xml
```

### Cached Files
```
data/interim/{model_name}/
├── {model}_substrates.csv
├── {model}_sequences.csv
└── {model}_merged_data.csv

data/processed/{model_name}/
└── {model}_processed_data.csv
```

## Pipeline Steps

1. **Prepare Data** - Extract substrates & sequences
2. **Merge Data** - Combine enzyme-substrate pairs
3. **Process kcat** - Annotate with predictions
4. **Optimize** - Enzyme-constrained FBA
5. **Tune** - Simulated annealing
6. **Analyze** - Optional FVA/Biolog
7. **Save** - Export final GEM

## Performance

| Mode | Steps 1-3 | Steps 4-7 | Total |
|------|-----------|-----------|-------|
| Cached (default) | ~10 sec | ~10-20 min | ~10-20 min |
| Force (--force) | ~10-20 min | ~10-20 min | ~20-40 min |

**Speedup: 50-98% with caching**

## Common Tasks

### Add New Model
1. Place XML in `data/raw/my_model.xml`
2. Create `configs/my_model.json`
3. Run: `python scripts/run_pipeline.py configs/my_model.json`

### Modify Parameters
1. Edit config file: `configs/my_model.json`
2. Run: `python scripts/run_pipeline.py configs/my_model.json`

### Compare Results
```bash
# Run same model with different parameters
python scripts/run_pipeline.py configs/iML1515_GEM.json
# Edit config...
python scripts/run_pipeline.py configs/iML1515_GEM.json --force
```

### Regenerate After Data Update
```bash
# Force regeneration when inputs change
python scripts/run_pipeline.py configs/my_model.json --force
```

## Troubleshooting

| Issue | Solution |
|-------|----------|
| Config not found | Check filename matches model name |
| Model XML not found | Ensure file in `data/raw/` |
| No CPI-Pred predictions | Generate predictions first |
| Biomass detection fails | Set `biomass_reaction` in config |
| Memory error | Reduce `max_iterations` in SA config |

## Documentation

- **`CONFIG_GUIDE.md`** - Full configuration guide
- **`SCRIPT_CACHING_GUIDE.md`** - Caching system details
- **`PIPELINE_SUMMARY.md`** - Complete overview
- **`VENV_SETUP.md`** - Environment setup

## Command Cheat Sheet

```bash
# Single runs
python scripts/run_pipeline.py configs/{model}.json
python scripts/run_pipeline.py configs/{model}.json --force

# Batch runs
bash scripts/batch_run.sh all
bash scripts/batch_run.sh ecoli
bash scripts/batch_run.sh genome
bash scripts/batch_run.sh standard
bash scripts/batch_run.sh {model1} {model2} ...
bash scripts/batch_run.sh all --force

# Check results
ls -lh results/tuning_results/
ls -lh models/

# Clean cache for specific model
rm -rf data/interim/{model}/ data/processed/{model}/
```
