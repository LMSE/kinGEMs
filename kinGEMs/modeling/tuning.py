"""
Parameter tuning module for kinGEMs.

This module provides simulated annealing functionality to tune kcat parameters
and optimize the model's performance.
"""
import copy  # noqa: F401
import math  # noqa: F401
import os
import random  # noqa: F401
import warnings

from Bio.Data.IUPACData import protein_letters
from Bio.SeqUtils import molecular_weight
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import matplotlib.pyplot as plt
import pandas as pd

from ..config import ensure_dir_exists
from ..dataset import annotate_model_with_kcat_and_gpr
from ..plots import plot_annealing_progress
from .optimize import run_optimization_with_dataframe

warnings.filterwarnings('ignore')
import logging
logging.getLogger('distributed').setLevel(logging.ERROR)
try:
    import gurobipy
    gurobipy.setParam('OutputFlag', 0)
except ImportError:
    pass




def simulated_annealing(
    model,
    processed_data,
    biomass_reaction,
    objective_value,
    gene_sequences_dict,
    output_dir=None,
    enzyme_fraction=0.15,
    temperature=1.0,
    cooling_rate=0.98,
    min_temperature=0.01,
    max_iterations=250,
    max_unchanged_iterations=3,
    change_threshold=0.001,
    verbose=False,
):
    """
    Use simulated annealing to tune kcat values for improved biomass production.
    """

    def acceptance_probability(old_cost, new_cost, temperature):
        if new_cost > old_cost:
            return 1.0
        return math.exp((new_cost - old_cost) / temperature)

    def get_neighbor(kcat_value, std):
        k_val_hr = kcat_value * 3600
        std_hr = std * 3600
        new_kcat = k_val_hr * (1.5 + random.uniform(-1, 3.5))
        if std_hr == 0:
            std_hr = k_val_hr * 0.1
        ub = k_val_hr + std_hr
        lb = k_val_hr - std_hr
        # Also limit to max reported kcat (per hour)
        max_hr = 4.6e9  # carbonic anhydrase per hour
        return min(max(new_kcat, lb), min(ub, max_hr))

    def update_kcat(df, reaction_id, gene_id, new_kcat_value):
        updated_df = df.copy()
        cond = (
            (updated_df['Reactions'] == reaction_id) &
            (updated_df['Single_gene'] == gene_id)
        )
        # convert back to per-second
        updated_df.loc[cond, 'kcat_mean'] = new_kcat_value / 3600
        return updated_df

    def calculate_molecular_weight(seq):
        return ProteinAnalysis(seq).molecular_weight()

    # Precompute MWs
    # mw_dict = {
    #     gene: calculate_molecular_weight(seq)
    #     for gene, seq in gene_sequences_dict.items() if seq
    # }
    
    def safe_mw(seq: str) -> float:
        # keep only standard amino acids
        cleaned = "".join([aa for aa in seq if aa in protein_letters])
        if not cleaned:
            cleaned = "A"  # fallback to alanine
        try:
            return molecular_weight(cleaned, seq_type="protein")
        except Exception:
            return 1e5  # large default if something still weird

    mw_dict = {
        gene: safe_mw(seq)
        for gene, seq in gene_sequences_dict.items() if seq
    }

    # INITIAL FBA
    biomass, df_FBA, _, _ = run_optimization_with_dataframe(
        model=model,
        processed_df=processed_data,
        objective_reaction=biomass_reaction,
        enzyme_upper_bound=enzyme_fraction,
        output_dir=output_dir,
        save_results=False,
        verbose=False
    )

    # pick top 65
    enzyme_df = df_FBA[df_FBA['Variable']=='enzyme'].copy()
    enzyme_df['MW'] = enzyme_df['Index'].map(mw_dict).fillna(0)
    enzyme_df['enzyme_mass'] = enzyme_df['Value'] * enzyme_df['MW'] * 1e-3
    top65 = enzyme_df.nlargest(65, 'enzyme_mass')
    top_targets = (
        top65[['Index','enzyme_mass']]
        .rename(columns={'Index':'Single_gene'})
        .merge(processed_data, on='Single_gene')
        [['Reactions','Single_gene','enzyme_mass','kcat_mean','kcat_std']]
        .reset_index(drop=True)
    )

    largest_rxn_id  = top_targets['Reactions'].tolist()
    largest_gene_id = top_targets['Single_gene'].tolist()
    current_solution = top_targets['kcat_mean'].tolist()
    stds             = top_targets['kcat_std'].fillna(0.1).tolist()

    df_new = processed_data.copy()
    current_biomass = biomass
    best_solution   = current_solution[:]
    best_biomass    = current_biomass
    best_df         = df_new.copy()

    iteration = 1
    no_change_counter = 0
    iterations = [0]
    biomasses  = [biomass]

    # ANNEALING
    while (temperature > min_temperature
           and iteration < max_iterations
           and current_biomass < objective_value):
        if verbose:
            print(f"\n--- Iteration {iteration} ---")
            print(f"Current biomass = {current_biomass:.6e}")

        # PROPOSE & print old vs new kcats
        updated_df = df_new.copy()
        for i, (rxn, gene) in enumerate(zip(largest_rxn_id, largest_gene_id)):
            old_k = current_solution[i]
            new_k_hr = get_neighbor(old_k, stds[i])
            new_k = new_k_hr / 3600.0
            if verbose:
                print(f"  {rxn}_{gene}: old kcat = {old_k:.3e}  →  new kcat = {new_k:.3e}")
            updated_df = update_kcat(updated_df, rxn, gene, new_k_hr)

        # ANNOTATE & EVALUATE
        temp_model = annotate_model_with_kcat_and_gpr(model=model, df=updated_df)

        new_biomass, temp_df_FBA, _, _ = run_optimization_with_dataframe(
            model=temp_model,
            processed_df=updated_df,
            objective_reaction=biomass_reaction,
            enzyme_upper_bound=enzyme_fraction,
            output_dir=None,
            save_results=False,
            verbose=False
        )
        if verbose:
            print(f"Proposed biomass = {new_biomass:.6e}")

        # compute change *before* accepting new solution
        change = abs(new_biomass - current_biomass) / max(current_biomass, 1e-6)

        # ACCEPT or REJECT
        prob = acceptance_probability(current_biomass, new_biomass, temperature)
        if prob > random.random():
            if verbose:
                print(f"Iteration {iteration}: ACCEPTED (Δ = {new_biomass-current_biomass:.2e})")
            model = temp_model
            df_FBA = temp_df_FBA
            df_new = updated_df.copy()
            current_biomass = new_biomass
            current_solution = [
                df_new.loc[
                    (df_new['Reactions']==rxn)&(df_new['Single_gene']==gene),
                    'kcat_mean'
                ].iat[0]
                for rxn,gene in zip(largest_rxn_id, largest_gene_id)
            ]
            if new_biomass > best_biomass:
                best_biomass  = new_biomass
                best_solution = current_solution[:]
                best_df       = df_new.copy()
        else:
            if verbose:    
                print(f"Iteration {iteration}: REJECTED (Δ = {new_biomass-current_biomass:.2e})")

        iterations.append(iteration)
        biomasses.append(new_biomass)

        # STAGNATION check uses `change` from above
        if change < change_threshold:
            no_change_counter += 1
            if no_change_counter >= max_unchanged_iterations:
                if verbose:
                    print(f"No significant change for {max_unchanged_iterations} iterations; stopping early.")
                break
        else:
            no_change_counter = 0

        temperature *= cooling_rate
        iteration += 1

    # FINALIZE: build kcat_dict from best_solution
    kcat_dict = {
        f"{rxn}_{gene}": k
        for (rxn, gene), k in zip(zip(largest_rxn_id, largest_gene_id), best_solution)
    }

    if output_dir:
        save_annealing_results(
            output_dir,
            kcat_dict,
            top_targets,
            best_df,
            iterations,
            biomasses,
            df_FBA
        )

    return kcat_dict, top_targets, best_df, iterations, biomasses, df_FBA


def save_annealing_results(output_dir, kcat_dict, df_enzyme_sorted, df_new, iterations, biomasses, df_FBA, prefix=""):
    """
    Save the results of the simulated annealing process.
    
    Parameters
    ----------
    output_dir : str
        Directory to save output files
    kcat_dict : dict
        Dictionary of optimized kcat values
    df_enzyme_sorted : pandas.DataFrame
        DataFrame with sorted enzyme data
    df_new : pandas.DataFrame
        DataFrame with updated kcat values
    iterations : list
        List of iteration numbers
    biomasses : list
        List of biomass values at each iteration
    df_FBA : pandas.DataFrame
        DataFrame with FBA results
    prefix : str, optional
        Prefix for output filenames
    """
    # Ensure directory exists
    ensure_dir_exists(output_dir)
    
    # Save kcat dictionary
    kcat_dict_df = pd.DataFrame(list(kcat_dict.items()), columns=['Key', 'Value'])
    kcat_dict_df.to_csv(os.path.join(output_dir, f"{prefix}kcat_dict.csv"), index=False)
    
    # Save sorted enzyme data
    df_enzyme_sorted.to_csv(os.path.join(output_dir, f"{prefix}df_enzyme_sorted.csv"), index=False)
    
    # Save updated data
    df_new.to_csv(os.path.join(output_dir, f"{prefix}df_new.csv"), index=False)
    
    # Save FBA results
    df_FBA.to_csv(os.path.join(output_dir, f"{prefix}df_FBA.csv"), index=False)
    
    # Save iterations data
    df_iterations = pd.DataFrame({"Iteration": iterations, "Biomass": biomasses})
    df_iterations.to_csv(os.path.join(output_dir, f"{prefix}iterations.csv"), index=False)
    
    # Create and save plot
    plot_annealing_progress(iterations, biomasses, 
                           output_path=os.path.join(output_dir, f"{prefix}annealing_progress.png"))
