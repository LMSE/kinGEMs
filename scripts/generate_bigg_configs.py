#!/usr/bin/env python3
"""
Generate configuration files for all BiGG models.

This script creates config JSON files for each model in the BiGG_models folder,
using the default biomass from COBRApy's slim_optimize as the simulated annealing goal.
"""

import json
import os
import sys
from pathlib import Path
import requests
import time

import cobra

# Add parent directory to path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


def get_bigg_organism(model_bigg_id):
    """Get organism name from BiGG database API."""
    try:
        url = f"http://bigg.ucsd.edu/api/v2/models/{model_bigg_id}"
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        data = response.json()
        return data.get('organism', None)
    except Exception as e:
        print(f"  Warning: Could not fetch organism from BiGG API: {e}")
        return None


def get_default_biomass(model_path, solver='glpk'):
    """Get the default biomass value using slim_optimize."""
    try:
        model = cobra.io.read_sbml_model(model_path)
        model.solver = solver
        biomass = model.slim_optimize()

        # Get biomass reaction ID
        obj_rxns = [rxn.id for rxn in model.reactions if rxn.objective_coefficient != 0]
        biomass_rxn = obj_rxns[0] if obj_rxns else None

        return biomass, biomass_rxn, len(model.genes), len(model.reactions)
    except Exception as e:
        print(f"  Error loading {os.path.basename(model_path)}: {e}")
        return None, None, None, None


def create_config(model_name, biomass_value, biomass_rxn, n_genes, n_reactions, bigg_organism=None):
    """Create a config dictionary following iML1515_GEM.json format."""

    # Map BiGG organism names to our standardized names in TAXONOMY_IDS
    organism_name_map = {
        # Bacteria - Proteobacteria - Enterobacteriaceae
        'Escherichia coli': 'E coli',
        'Shigella': 'E coli',
        'Salmonella': 'Salmonella',
        'Klebsiella': 'Klebsiella',
        'Yersinia': 'Yersinia',
        
        # Bacteria - Proteobacteria - Other
        'Pseudomonas': 'Pseudomonas',
        'Acinetobacter': 'Acinetobacter',
        'Geobacter': 'Geobacter',
        'Burkholderia': 'Burkholderia',
        
        # Bacteria - Firmicutes
        'Bacillus subtilis': 'B subtilis',
        'Staphylococcus aureus': 'S aureus',
        'Lactococcus lactis': 'L lactis',
        'Clostridium': 'Clostridium',
        'Clostridioides': 'Clostridioides',
        'Thermoanaerobacter': 'Thermoanaerobacter',
        
        # Bacteria - Actinobacteria
        'Mycobacterium tuberculosis': 'M tuberculosis',
        'Streptomyces': 'Streptomyces',
        
        # Bacteria - Cyanobacteria
        'Synechococcus': 'Synechococcus',
        'Synechocystis': 'Synechocystis',
        
        # Bacteria - Other
        'Helicobacter pylori': 'H pylori',
        'Thermotoga': 'Thermotoga',
        'Bacteroides': 'Bacteroides',
        
        # Archaea
        'Methanosarcina': 'Methanosarcina',
        'Methanobacterium': 'Methanobacterium',
        
        # Eukaryotes - Fungi
        'Saccharomyces cerevisiae': 'Yeast',
        
        # Eukaryotes - Mammals
        'Homo sapiens': 'Human',
        'Mus musculus': 'Mouse',
        'Cricetulus griseus': 'Chinese Hamster Ovary',
        
        # Eukaryotes - Parasites
        'Trypanosoma cruzi': 'T cruzi',
        'Leishmania': 'Leishmania',
        'Plasmodium falciparum': 'P falciparum',
        'Plasmodium berghei': 'P berghei',
        'Plasmodium chabaudi': 'P chabaudi',
        'Plasmodium knowlesi': 'P knowlesi',
        'Plasmodium vivax': 'P vivax',
        'Plasmodium cynomolgi': 'P cynomolgi',
        
        # Eukaryotes - Algae/Diatoms
        'Phaeodactylum': 'Phaeodactylum',
        'Chlamydomonas': 'Chlamydomonas',
    }

    # Try to map BiGG organism to our standardized name
    organism = 'E coli'  # Default fallback
    if bigg_organism:
        # Try exact match first
        if bigg_organism in organism_name_map:
            organism = organism_name_map[bigg_organism]
        else:
            # Try partial match (e.g., "Escherichia coli str. K-12" contains "Escherichia coli")
            for bigg_name, our_name in organism_name_map.items():
                if bigg_name.lower() in bigg_organism.lower():
                    organism = our_name
                    break

    # Set n_top_enzymes based on model size
    if n_reactions < 1000:
        n_top_enzymes = 50
    elif n_reactions < 2000:
        n_top_enzymes = 100
    elif n_reactions < 3000:
        n_top_enzymes = 200
    else:
        n_top_enzymes = 500

    config = {
        "model_name": model_name,
        "organism": organism,
        "biomass_reaction": biomass_rxn,
        "enzyme_upper_bound": 0.25,
        "solver": "glpk",
        "enable_fva": True,
        "enable_biolog_validation": False,
        "enable_maintenance_sweep": False,
        "edit_ngam": False,
        "ngam_rxn_id": "ATPM",
        "results_subdir": "BiGG_models",
        "_comment_results_subdir": "Subdirectory under tuning_results/ for organizing results",
        "_comment_edit_ngam": "Set to true to set NGAM (ATPM) lower bound to 0 during optimization",
        "_comment_ngam_rxn_id": "Reaction ID for NGAM/ATPM (default: 'ATPM')",
        "_comment_enable_maintenance_sweep": "Set to true to run NGAM/GAM parameter sweep",

        "fva": {
            "parallel": True,
            "workers": 4,
            "chunk_size": 50,
            "method": "dask",
            "opt_ratio": 0.9,
            "_comment_opt_ratio": "Fraction of optimal biomass to constrain during FVA (0.9 = 90% of optimal)",
            "_comment_method": "Options: 'dask' (distributed, scalable) or 'multiprocessing' (simpler, single-machine)",
            "_comment_workers": "Set to null for auto-detection (uses all CPU cores)",
            "_comment_chunk_size": "Number of reactions per task. null = auto-calculate. Larger = less overhead but less load balancing"
        },

        "simulated_annealing": {
            "temperature": 100,
            "cooling_rate": 0.95,
            "min_temperature": 0.001,
            "max_iterations": 1000,
            "max_unchanged_iterations": 50,
            "change_threshold": 0.001,
            "biomass_goal": round(biomass_value, 4) if biomass_value else 0.5,
            "n_top_enzymes": n_top_enzymes,
            "verbose": False,
            "_comment_aggressive": "AGGRESSIVE SETTINGS - larger perturbations and longer search",
            "_comment_temperature": "High starting temp (100) for large initial exploration",
            "_comment_cooling_rate": "Moderate cooling (0.95) for thorough exploration",
            "_comment_change_threshold": "Standard threshold (0.1%) for meaningful improvements",
            "_comment_biomass_goal": f"Target goal based on COBRApy slim_optimize: {biomass_value:.4f}" if biomass_value else "Target biomass goal",
            "_comment_n_top_enzymes": f"Tune top {n_top_enzymes} enzymes (scaled by model size: {n_reactions} reactions)"
        },

        "maintenance_sweep": {
            "ngam_range": [0, 1, 2, 3, 4, 5, 6, 7],
            "gam_range": [40, 45, 50, 55, 60, 70, 75.55, 80, 100],
            "verbose": False,
            "_comment_ngam_range": "NGAM values to test (ATPM lower bound in mmol/gDW/h). Higher values (10-100) to reduce biomass toward goal.",
            "_comment_gam_range": "GAM values to test (ATP in biomass, mmol/gDW). Higher values (75.55-300) to reduce biomass."
        }
    }

    return config


def main():
    """Generate config files for all BiGG models."""

    project_root = Path(__file__).parent.parent
    bigg_models_dir = project_root / "data" / "raw" / "BiGG_models"
    configs_dir = project_root / "configs" / "BiGG_models"

    # Create configs directory
    configs_dir.mkdir(parents=True, exist_ok=True)

    # Get all XML files
    model_files = sorted(bigg_models_dir.glob("*.xml"))

    print(f"Found {len(model_files)} models in {bigg_models_dir}")
    print(f"Generating config files in {configs_dir}\n")

    successful = 0
    failed = 0

    for model_file in model_files:
        model_name = model_file.stem  # Filename without .xml extension
        print(f"Processing {model_name}...")

        # Get organism from BiGG API
        bigg_organism = get_bigg_organism(model_name)
        if bigg_organism:
            print(f"  BiGG organism: {bigg_organism}")
        time.sleep(0.1)  # Be nice to the API

        # Get default biomass
        biomass, biomass_rxn, n_genes, n_reactions = get_default_biomass(str(model_file))

        if biomass is None:
            print(f"  ✗ Failed to load model")
            failed += 1
            continue

        print(f"  Genes: {n_genes}, Reactions: {n_reactions}")
        print(f"  Biomass reaction: {biomass_rxn}")
        print(f"  Default biomass: {biomass:.4f}")

        # Create config
        config = create_config(model_name, biomass, biomass_rxn, n_genes, n_reactions, bigg_organism)

        # Save config file
        config_path = configs_dir / f"{model_name}.json"
        with open(config_path, 'w') as f:
            json.dump(config, f, indent=2)

        print(f"  ✓ Saved config to {config_path.name}\n")
        successful += 1

    print("="*70)
    print(f"Config generation complete!")
    print(f"  Successful: {successful}")
    print(f"  Failed: {failed}")
    print(f"  Total: {len(model_files)}")
    print(f"\nConfig files saved to: {configs_dir}")
    print("="*70)


if __name__ == '__main__':
    main()
