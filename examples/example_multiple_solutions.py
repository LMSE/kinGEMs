"""
Example: Generate Multiple Alternative Solutions using Solution Pools

This script demonstrates how to use the solution pool feature to find
multiple near-optimal solutions for enzyme-constrained FBA.

This is particularly useful for:
1. Exploring solution space diversity
2. Finding alternative metabolic strategies
3. Understanding flux flexibility
4. Identifying robust vs flexible reactions

Only works with Gurobi or CPLEX solvers.
"""

import os
import sys
import pandas as pd

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from kinGEMs.modeling.optimize import run_optimization_with_dataframe

def compare_solutions(solutions):
    """
    Compare multiple solutions to identify differences.
    
    Parameters
    ----------
    solutions : list of (sol_val, df_FBA) tuples
        List of solution results
    
    Returns
    -------
    pandas.DataFrame
        DataFrame showing flux variability across solutions
    """
    if len(solutions) < 2:
        print("Need at least 2 solutions to compare")
        return None
    
    print(f"\n{'='*60}")
    print(f"COMPARING {len(solutions)} SOLUTIONS")
    print(f"{'='*60}\n")
    
    # Extract objective values
    obj_values = [sol_val for sol_val, _ in solutions]
    print(f"Objective values: {[f'{v:.6f}' for v in obj_values]}")
    print(f"Range: {min(obj_values):.6f} to {max(obj_values):.6f}")
    print(f"Spread: {(max(obj_values) - min(obj_values)):.6f} ({(max(obj_values) - min(obj_values))/max(obj_values)*100:.2f}%)\n")
    
    # Extract flux data for comparison
    flux_dfs = []
    for idx, (_, df_fba) in enumerate(solutions):
        flux_data = df_fba[df_fba['Variable'] == 'flux'][['Index', 'Value']].copy()
        flux_data = flux_data.rename(columns={'Value': f'Sol_{idx+1}'})
        flux_dfs.append(flux_data)
    
    # Merge all flux data
    comparison = flux_dfs[0]
    for df in flux_dfs[1:]:
        comparison = comparison.merge(df, on='Index', how='outer')
    
    # Fill NaN with 0
    comparison = comparison.fillna(0)
    
    # Calculate statistics
    flux_cols = [col for col in comparison.columns if col.startswith('Sol_')]
    comparison['Mean'] = comparison[flux_cols].mean(axis=1)
    comparison['Std'] = comparison[flux_cols].std(axis=1)
    comparison['Range'] = comparison[flux_cols].max(axis=1) - comparison[flux_cols].min(axis=1)
    
    # Identify reactions with high variability
    # (excluding reactions that are zero in all solutions)
    non_zero = comparison[comparison[flux_cols].abs().sum(axis=1) > 1e-6]
    variable_rxns = non_zero[non_zero['Std'] > 0.01].sort_values('Std', ascending=False)
    
    print(f"Total reactions: {len(comparison)}")
    print(f"Active reactions (non-zero): {len(non_zero)}")
    print(f"Variable reactions (std > 0.01): {len(variable_rxns)}\n")
    
    if len(variable_rxns) > 0:
        print("Top 10 most variable reactions:")
        print(variable_rxns[['Index', 'Mean', 'Std', 'Range']].head(10).to_string(index=False))
    
    return comparison


def main():
    """
    Example: Find multiple solutions for E. coli core model
    """
    
    # Configuration
    model_path = "path/to/your/model.json"  # Update this
    processed_data_path = "path/to/your/processed_data.csv"  # Update this
    
    print("="*60)
    print("SOLUTION POOL EXAMPLE")
    print("="*60)
    print("\nThis example demonstrates finding multiple alternative")
    print("optimal/near-optimal solutions using Gurobi's solution pool.\n")
    
    # Load data
    print("Loading model and data...")
    processed_df = pd.read_csv(processed_data_path)
    
    # Run optimization with solution pool
    num_solutions = 10  # Request 10 alternative solutions
    solution_pool_gap = 0.02  # Allow solutions within 2% of optimal
    
    print(f"\nRunning optimization:")
    print(f"  - Requesting {num_solutions} solutions")
    print(f"  - Allowing gap of {solution_pool_gap*100:.1f}% from optimal")
    print(f"  - Using Gurobi solver\n")
    
    try:
        result = run_optimization_with_dataframe(
            model=model_path,
            processed_df=processed_df,
            objective_reaction='BIOMASS_Ecoli_core_w_GAM',
            enzyme_upper_bound=0.125,
            maximization=True,
            solver_name='gurobi',
            num_solutions=num_solutions,
            solution_pool_gap=solution_pool_gap,
            output_dir='results/solution_pool',
            save_results=True,
            verbose=True
        )
        
        # Check if multiple solutions were returned
        if isinstance(result[0], list):
            solutions, gene_sequences_dict, output_paths = result
            print(f"\n✓ Successfully found {len(solutions)} solutions!")
            
            # Compare solutions
            comparison_df = compare_solutions(solutions)
            
            # Save comparison
            if comparison_df is not None:
                comparison_path = 'results/solution_pool/flux_comparison.csv'
                comparison_df.to_csv(comparison_path, index=False)
                print(f"\n✓ Flux comparison saved to: {comparison_path}")
            
        else:
            print("\n⚠️  Only single solution returned (solver may not support solution pools)")
            solution_value, df_FBA, gene_sequences_dict, output_path = result
            print(f"Objective value: {solution_value:.6f}")
    
    except Exception as e:
        print(f"\n✗ Error: {e}")
        import traceback
        traceback.print_exc()


if __name__ == '__main__':
    main()
