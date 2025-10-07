# Results Directory

This directory contains output files from pipeline runs. **These files are not tracked by git** due to their large size.

## Directory Structure

```
results/
├── tuning_results/           # Simulated annealing results (NOT in git)
│   └── {model}_{date}_{id}/
│       ├── df_new.csv        # Enzyme mass contributions
│       ├── df_FBA.csv         # Flux balance analysis results
│       ├── kcat_dict.csv      # Tuned kcat values
│       ├── final_model_info.csv  # Merged results
│       ├── iterations.csv     # Annealing iterations
│       ├── annealing_progress.png  # Progress visualization
│       └── *_fva_*.csv/png    # FVA results (if enabled)
│
└── validation/               # Validation results (NOT in git)
    └── ...
```

## Important Notes

### Files Are Excluded from Git

Result files can be **very large** (100+ MB) and are automatically excluded from version control via `.gitignore`:

```
results/tuning_results/
results/validation/
```

### Disk Space

- Each pipeline run can generate 300-500 MB of results
- Monitor your disk quota on HPC systems
- Clean up old runs periodically

### Cleaning Up Old Results

To remove old tuning results:

```bash
# List results by size
du -sh results/tuning_results/*/ | sort -h

# Remove specific run
rm -rf results/tuning_results/{model}_{date}_{id}/

# Remove all tuning results (careful!)
rm -rf results/tuning_results/*/
```

### Backing Up Results

For important results, copy them to a permanent storage location:

```bash
# On Compute Canada
cp -r results/tuning_results/{run_id}/ ~/projects/def-mahadeva/results/

# Or archive and compress
tar -czf results_{run_id}.tar.gz results/tuning_results/{run_id}/
```

## Result File Descriptions

### df_new.csv
- Enzyme mass contributions for each reaction-gene pair
- Contains kcat values and enzyme masses
- Used for final model annotation
- **Size: 100-400 MB**

### final_model_info.csv
- Merged dataframe with tuned kcat values
- Includes kcat_tuned column from simulated annealing
- **Size: 100-400 MB**

### kcat_dict.csv
- Mapping of reaction_gene pairs to tuned kcat values
- Compact format (< 1 MB)
- Can be used to recreate tuned models

### df_FBA.csv
- Flux balance analysis results
- Reaction fluxes and metabolite concentrations
- **Size: 1-10 MB**

### iterations.csv
- Simulated annealing iteration history
- Biomass values at each iteration
- Useful for convergence analysis
- **Size: < 1 MB**

### annealing_progress.png
- Visualization of annealing convergence
- Shows biomass improvement over iterations
- **Size: 100-500 KB**

## Sharing Results

To share results with collaborators:

1. **Archive the run directory:**
   ```bash
   tar -czf run_{date}.tar.gz results/tuning_results/{run_id}/
   ```

2. **Upload to shared storage:**
   ```bash
   # Compute Canada project space
   cp run_{date}.tar.gz ~/projects/def-mahadeva/shared/
   ```

3. **Or use external services:**
   - Globus for large file transfers
   - Google Drive / Dropbox for smaller archives

## Troubleshooting

### "Disk quota exceeded"

**Solution:** Clean up old results or move to project space

```bash
# Check disk usage
du -sh results/tuning_results/

# Move to project space
mv results/tuning_results/* ~/projects/def-mahadeva/kinGEMs_results/
```

### "Git won't push - files too large"

**Solution:** Results should already be in `.gitignore`. If you see this error:

```bash
# Check if files are tracked
git ls-files results/

# If files are tracked, remove them
git rm --cached -r results/tuning_results/
git commit -m "Remove large result files from git"
```

## See Also

- `CONFIG_GUIDE.md` - Pipeline configuration
- `PIPELINE_SUMMARY.md` - Complete pipeline overview
- `SCRIPT_CACHING_GUIDE.md` - Intermediate file caching
