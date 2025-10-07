
# kinGEMs Pipeline Script for E. coli iML1515
#
# Usage:
#   python scripts/run_iML1515_GEM.py           # Use cached intermediate files if available
#   python scripts/run_iML1515_GEM.py --force   # Force regeneration of all intermediate files
#   python scripts/run_iML1515_GEM.py -f        # Same as --force
#
# The script will automatically skip Steps 1-3 if intermediate files exist:
#   - Step 1: substrates.csv, sequences.csv
#   - Step 2: merged_data.csv
#   - Step 3: processed_data.csv
#
# Use --force to regenerate these files (useful if input data has changed)
#
from datetime import datetime
import logging
import os
import random
import sys
import warnings

import cobra
from cobra.flux_analysis import flux_variability_analysis as cobra_fva
from cobra.io import write_sbml_model
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Add parent directory to Python path before kinGEMs imports
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from kinGEMs.dataset import (
    annotate_model_with_kcat_and_gpr,
    assign_kcats_to_model,
    load_model,
    merge_substrate_sequences,
    prepare_model_data,
    process_kcat_predictions,
)
from kinGEMs.modeling.fva import (
    flux_variability_analysis,
    plot_cumulative_fvi_distribution,
    plot_flux_variability,
)
from kinGEMs.modeling.optimize import run_optimization_with_dataframe
from kinGEMs.modeling.tuning import simulated_annealing

# Suppress warnings and configure logging
warnings.filterwarnings('ignore')
logging.getLogger('distributed').setLevel(logging.ERROR)
try:
    import gurobipy
    gurobipy.setParam('OutputFlag', 0)
except ImportError:
    pass

# Check for --force flag to regenerate all intermediate files
FORCE_REGENERATE = '--force' in sys.argv or '-f' in sys.argv

# === Configuration ===
organism_strain_GEMname = "ecoli_iML1515"
organism = "E coli"
run_id = f"{organism_strain_GEMname}_{datetime.today().strftime('%Y%m%d')}_{random.randint(1000, 9999)}"

project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
data_dir = os.path.join(project_root, "data")
raw_data_dir = os.path.join(data_dir, "raw")
interim_data_dir = os.path.join(data_dir, "interim")
interim_data_dir_spec = os.path.join(interim_data_dir, f"{organism_strain_GEMname}")
processed_data_dir = os.path.join(data_dir, "processed")
processed_data_dir_spec = os.path.join(processed_data_dir, f"{organism_strain_GEMname}")
CPIPred_data_dir = os.path.join(interim_data_dir, "CPI-Pred predictions")
results_dir = os.path.join(project_root, "results")
tuning_results_dir = os.path.join(results_dir, "tuning_results", run_id)
os.makedirs(tuning_results_dir, exist_ok=True)

# Input files
model_path = os.path.join(raw_data_dir, "iML1515_GEM.xml")
predictions_csv_path = os.path.join(CPIPred_data_dir, f"X06A_kinGEMs_{organism_strain_GEMname}_predictions.csv")

# Output files
substrates_output = os.path.join(interim_data_dir_spec, f"{organism_strain_GEMname}_substrates.csv")
sequences_output = os.path.join(interim_data_dir_spec, f"{organism_strain_GEMname}_sequences.csv")
merged_data_output = os.path.join(interim_data_dir_spec, f"{organism_strain_GEMname}_merged_data.csv")
processed_data_output = os.path.join(processed_data_dir_spec, f"{organism_strain_GEMname}_processed_data.csv")

# Simulation parameters
biomass_reaction = 'BIOMASS_Ec_iML1515_core_75p37M'
enzyme_upper_bound = 0.15

print(f"=== kinGEMs Pipeline for {organism_strain_GEMname} ===")
print(f"Run ID: {run_id}")
print(f"Results will be saved to: {tuning_results_dir}")
if FORCE_REGENERATE:
    print("⚠️  Force regenerate mode: will regenerate all intermediate files")

# === Step 1: Preparing model data ===
print("\n=== Step 1: Preparing model data ===")
if not FORCE_REGENERATE and os.path.exists(substrates_output) and os.path.exists(sequences_output):
    print("  Found existing substrate and sequence files. Loading from disk...")
    substrate_df = pd.read_csv(substrates_output)
    sequences_df = pd.read_csv(sequences_output)
    from kinGEMs.dataset import load_model
    model = load_model(model_path)
    print(f"  Loaded {len(substrate_df)} substrates and {len(sequences_df)} sequences")
else:
    print("  Generating new substrate and sequence data...")
    model, substrate_df, sequences_df = prepare_model_data(
        model_path=model_path,
        substrates_output=substrates_output,
        sequences_output=sequences_output,
        organism=organism
    )
    print("  Generated and saved substrate and sequence data")

print(f"  Genes in model: {len(model.genes)}")
print(f"  Reactions in model: {len(model.reactions)}")

# === Step 2: Merging substrate and sequence data ===
print("\n=== Step 2: Merging substrate and sequence data ===")
if not FORCE_REGENERATE and os.path.exists(merged_data_output):
    print("  Found existing merged data file. Loading from disk...")
    merged_data = pd.read_csv(merged_data_output)
    print(f"  Loaded {len(merged_data)} merged rows")
else:
    print("  Generating new merged data...")
    merged_data = merge_substrate_sequences(
        substrate_df=substrate_df,
        sequences_df=sequences_df,
        model=model,
        output_path=merged_data_output
    )
    print(f"  Generated and saved {len(merged_data)} merged rows")

# === Step 3: Processing CPI-Pred kcat values ===
print("\n=== Step 3: Processing CPI-Pred kcat values & annotating model ===")
if not FORCE_REGENERATE and os.path.exists(processed_data_output):
    print("  Found existing processed data file. Loading from disk...")
    processed_data = pd.read_csv(processed_data_output)
    print(f"  Loaded {len(processed_data)} processed rows")
else:
    print("  Generating new processed kcat data...")
    processed_data = process_kcat_predictions(
        merged_df=merged_data,
        predictions_csv_path=predictions_csv_path,
        output_path=processed_data_output
    )
    print(f"  Generated and saved {len(processed_data)} processed rows")

# Annotate model with kcat values
model = annotate_model_with_kcat_and_gpr(
    model=model,
    df=processed_data
)

# Check kcat coverage
rxn_with_kcat = sum(1 for rxn in model.reactions 
                    if hasattr(rxn, 'annotation') and 'kcat' in rxn.annotation 
                    and rxn.annotation['kcat'] not in [None, '', 0, '0'])
print(f"  Reactions with kcat: {rxn_with_kcat}/{len(model.reactions)}")

# === Step 4: Optimization (FBA Sanity Check + kinGEMs) ===
print("\n=== Step 4: Running optimization with kcat constraints ===")
(solution_value, df_FBA, gene_sequences_dict, _) = run_optimization_with_dataframe(
    model=model,
    processed_df=processed_data,
    objective_reaction=biomass_reaction,
    enzyme_upper_bound=enzyme_upper_bound,
    enzyme_ratio=True,
    maximization=True,
    multi_enzyme_off=False,
    isoenzymes_off=False,
    promiscuous_off=False,
    complexes_off=False,
    output_dir=None,
    save_results=False,
    print_reaction_conditions=True,
    verbose=False
)

print(f"  Initial biomass value: {solution_value:.4f}")

# === Step 5: Simulated Annealing ===
print("\n=== Step 5: Running simulated annealing ===")
temperature = 1.0
cooling_rate = 0.975
min_temperature = 0.001
max_iterations = 100
max_unchanged_iterations = 5
change_threshold = 0.009
biomass_goal = 0.5

kcat_dict, top_targets, df_new, iterations, biomasses, df_FBA = simulated_annealing(
    model=model,
    processed_data=processed_data,
    biomass_reaction=biomass_reaction,
    objective_value=biomass_goal,
    gene_sequences_dict=gene_sequences_dict,
    output_dir=tuning_results_dir,
    enzyme_fraction=enzyme_upper_bound,
    temperature=temperature,
    cooling_rate=cooling_rate,
    min_temperature=min_temperature,
    max_iterations=max_iterations,
    max_unchanged_iterations=max_unchanged_iterations,
    change_threshold=change_threshold
)

print(f"  Final biomass: {biomasses[-1]:.4f}")
print(f"  Improvement: {(biomasses[-1] - biomasses[0]) / biomasses[0] * 100:.1f}%")
print("\n  Top 10 enzymes by mass contribution:")
print(top_targets[['Reactions', 'Single_gene', 'enzyme_mass']].head(10))

# Save tuning results
df_new_path = os.path.join(tuning_results_dir, "df_new.csv")
df_new.to_csv(df_new_path, index=False)
kcat_dict_path = os.path.join(tuning_results_dir, "kcat_dict.csv")

# Merge kcat_dict into df_new for final model info
kcat_dict_df = pd.read_csv(kcat_dict_path)
if 'reaction_gene' not in kcat_dict_df.columns:
    kcat_dict_df.columns = ['reaction_gene', 'kcat_value']

df_new['reaction_gene'] = df_new['Reactions'].astype(str) + '_' + df_new['Single_gene'].astype(str)
df_new = df_new.merge(kcat_dict_df, on='reaction_gene', how='left')
df_new.rename(columns={'kcat_value': 'kcat_tuned'}, inplace=True)

final_info_path = os.path.join(tuning_results_dir, "final_model_info.csv")
df_new.to_csv(final_info_path, index=False)
print(f"\n  Saved merged DataFrame with kcat_tuned to: {final_info_path}")

# === Step 6: Flux Variability Analysis ===
print("\n=== Step 6: Running Flux Variability Analysis ===")
fva_results_path = os.path.join(tuning_results_dir, f"{organism_strain_GEMname}_fva_results.csv")
fva_plot_path = os.path.join(tuning_results_dir, f"{organism_strain_GEMname}_fva_flux_range_plot.png")
fva_cumulative_path = os.path.join(tuning_results_dir, f"{organism_strain_GEMname}_fva_cumulative_plot.png")

# Run kinGEMs FVA
fva_results, _, _ = flux_variability_analysis(
    model=model,
    processed_df=df_new,
    biomass_reaction=biomass_reaction,
    output_file=fva_results_path,
    enzyme_upper_bound=enzyme_upper_bound
)

print(f"  kinGEMs FVA completed: {len(fva_results)} reactions analyzed")

# Run COBRApy FVA for comparison
print("  Running COBRApy FVA for comparison...")
cobra_fva_results = cobra_fva(model, fraction_of_optimum=0.9)
cobra_fva_df = pd.DataFrame({
    "Reactions": cobra_fva_results.index,
    "Min Solutions": cobra_fva_results['minimum'],
    "Max Solutions": cobra_fva_results['maximum'],
    "Solution Biomass": [model.slim_optimize()] * len(cobra_fva_results)
})

print(f"  COBRApy FVA completed: {len(cobra_fva_df)} reactions analyzed")

# Plot standard FVA range
print("  Generating FVA plots...")
fig = plot_flux_variability(fva_results, output_file=fva_plot_path)
print(f"  Saved FVA flux range plot to: {fva_plot_path}")

# Plot cumulative distribution comparing kinGEMs and COBRApy
plot_cumulative_fvi_distribution(
    dfs=[fva_results, cobra_fva_df],
    labels=["kinGEMs FVA", "COBRApy FVA"],
    output_file=fva_cumulative_path
)
print(f"  Saved FVA cumulative plot to: {fva_cumulative_path}")

# === Step 7: Save Final Model ===
print("\n=== Step 7: Saving final model ===")
model_output_dir = os.path.join(project_root, "models")
os.makedirs(model_output_dir, exist_ok=True)
model_output_path = os.path.join(model_output_dir, f"{run_id}.xml")

# Assign kcats to model
model_with_kcats = assign_kcats_to_model(model, df_new)

# Clean up reaction annotations before SBML export
def clean_annotations(model):
    """Convert float values in annotations to strings for SBML compatibility."""
    for rxn in model.reactions:
        ann = rxn.annotation
        if not isinstance(ann, dict):
            rxn.annotation = {}
        else:
            new_ann = {}
            for k, v in ann.items():
                if isinstance(v, float):
                    new_ann[k] = str(v)
                elif isinstance(v, (list, tuple)):
                    new_ann[k] = [str(item) if isinstance(item, float) else item for item in v]
                elif isinstance(v, (str, dict)):
                    new_ann[k] = v
            rxn.annotation = new_ann
    return model

model_with_kcats = clean_annotations(model_with_kcats)

# Save the final model
write_sbml_model(model_with_kcats, model_output_path)
print(f"  Final GEM saved to: {model_output_path}")

# === Summary ===
print("\n" + "="*60)
print("=== Pipeline Complete ===")
print("="*60)
print(f"Run ID: {run_id}")
print(f"Initial biomass: {biomasses[0]:.4f}")
print(f"Final biomass: {biomasses[-1]:.4f}")
print(f"Improvement: {(biomasses[-1] - biomasses[0]) / biomasses[0] * 100:.1f}%")
print(f"Iterations: {iterations}")
print(f"\nResults directory: {tuning_results_dir}")
print(f"Final model: {model_output_path}")
print("="*60)
