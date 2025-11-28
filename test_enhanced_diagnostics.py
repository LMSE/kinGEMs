#!/usr/bin/env python3
"""
Test script to run optimization with enhanced diagnostics
"""

import pandas as pd
import cobra as cb
from kinGEMs.modeling.optimize import run_optimization_with_dataframe

def main():
    print("=== Testing Enhanced Diagnostics with iML1515_GEM ===\n")
    
    # Load model
    model_path = "/project/6000502/ranaab/kinGEMs_v2/models/ecoli_iML1515_20250826_4941.xml"
    print(f"Loading model: {model_path}")
    model = cb.io.read_sbml_model(model_path)
    print(f"Model loaded: {len(model.genes)} genes, {len(model.reactions)} reactions")
    
    # Load processed data
    data_path = "/project/6000502/ranaab/kinGEMs_v2/data/processed/iML1515_GEM/iML1515_GEM_processed_data.csv"
    print(f"Loading processed data: {data_path}")
    df = pd.read_csv(data_path)
    print(f"Data loaded: {df.shape[0]} rows")
    print(f"Direction counts: {dict(df['Direction'].value_counts())}")
    
    # Run optimization with enhanced diagnostics
    print("\n=== Running optimization with bidirectional constraints and verbose output ===")
    
    # Test with a smaller enzyme upper bound to ensure constraints are active
    sol_val, df_FBA, gene_seq, output_path = run_optimization_with_dataframe(
        model=model,
        processed_df=df,
        objective_reaction='BIOMASS_Ec_iML1515_core_75p37M',
        enzyme_upper_bound=0.05,  # Smaller bound to ensure constraints are active
        enzyme_ratio=True,
        maximization=True,
        multi_enzyme_off=False,
        isoenzymes_off=False,
        promiscuous_off=False,
        complexes_off=False,
        output_dir=None,
        save_results=False,
        verbose=True,  # Enable verbose output to see enhanced diagnostics
        solver_name='glpk',
        bidirectional_constraints=True  # Enable bidirectional constraints
    )
    
    print("\n=== Results ===")
    print(f"Objective value: {sol_val:.6f}")
    print(f"Result shape: {df_FBA.shape}")
    
    # Show some enzyme allocations
    enzyme_data = df_FBA[df_FBA['Variable'] == 'enzyme']
    nonzero_enzymes = enzyme_data[enzyme_data['Value'] > 1e-6]
    print(f"Non-zero enzyme allocations: {len(nonzero_enzymes)}")
    
    if len(nonzero_enzymes) > 0:
        print("Top 10 enzyme allocations:")
        print(nonzero_enzymes.nlargest(10, 'Value')[['Index', 'Value']])
    
    # Show bidirectional flux components if available
    if 'flux_fwd' in df_FBA['Variable'].values:
        flux_fwd = df_FBA[df_FBA['Variable'] == 'flux_fwd']
        nonzero_fwd = flux_fwd[flux_fwd['Value'] > 1e-6]
        print(f"\nNon-zero forward fluxes: {len(nonzero_fwd)}")
        
        flux_rev = df_FBA[df_FBA['Variable'] == 'flux_rev']
        nonzero_rev = flux_rev[flux_rev['Value'] > 1e-6]
        print(f"Non-zero reverse fluxes: {len(nonzero_rev)}")

if __name__ == "__main__":
    main()
