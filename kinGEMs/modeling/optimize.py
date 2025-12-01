"""
Optimization module for kinGEMs.

This module provides functions for enzyme-constrained flux balance analysis,
including core optimization functionality from the original KG03b module.
"""

from collections import Counter  # noqa: F401
from itertools import product
import logging
import logging as pyomo_logging
import math
import os

import psutil  # For memory monitoring
import re
import time  # For timing measurements
import warnings

from Bio.SeqUtils import molecular_weight
import cobra as cb
from cobra.util.array import create_stoichiometric_matrix
import numpy as np  # noqa: F401

import pandas as pd
import pyomo.environ as pyo
from pyomo.environ import *  # noqa: F403
from pyomo.environ import value  # Explicit import to fix lint warning
from pyomo.opt import SolverFactory

from ..config import ensure_dir_exists  # noqa: F401

warnings.filterwarnings('ignore')

# Suppress Pyomo warnings about setting variables slightly outside bounds due to numerical precision
pyomo_logging.getLogger('pyomo.core').setLevel(logging.ERROR)


logging.getLogger('distributed').setLevel(logging.ERROR)
try:
    import gurobipy
    gurobipy.setParam('OutputFlag', 0)
except ImportError:
    pass

def _tokenize_gpr(rule):
    """Split a GPR rule into tokens: parentheses, 'and', 'or', gene IDs."""
    # Use case-insensitive matching for 'and'/'or', but preserve gene ID case
    token_spec = r'\(|\)|(?i:and)|(?i:or)|[^()\s]+'
    tokens = re.findall(token_spec, rule)
    # Normalize 'and'/'or' to lowercase for parser, keep gene IDs as-is
    return [t.lower() if t.lower() in ('and', 'or') else t for t in tokens]


def _get_memory_usage():
    """Get current memory usage in MB."""
    try:
        process = psutil.Process(os.getpid())
        return process.memory_info().rss / 1024 / 1024  # Convert bytes to MB
    except Exception:
        return -1  # Return -1 if memory monitoring fails


def _print_step_timing(step_name, start_time, verbose=True):
    """Print timing and memory info for a step."""
    if not verbose:
        return

    elapsed = time.time() - start_time
    memory_mb = _get_memory_usage()

    if memory_mb > 0:
        print(f"  ✓ {step_name}: {elapsed:.2f}s (Memory: {memory_mb:.1f} MB)")
    else:
        print(f"  ✓ {step_name}: {elapsed:.2f}s")


def _parse_gpr_to_dnf(tokens):
    """
    Parse tokens into Disjunctive Normal Form (DNF):
    returns a list of clauses, each clause is a list of gene IDs.
    e.g. "g1 and (g2 or g3)" -> [['g1','g2'], ['g1','g3']]
    """
    def parse_expression(idx=0):
        clauses, idx = parse_term(idx)
        while idx < len(tokens) and tokens[idx] == 'or':
            right, idx = parse_term(idx+1)
            clauses += right
        return clauses, idx

    def parse_term(idx):
        clauses, idx = parse_factor(idx)
        while idx < len(tokens) and tokens[idx] == 'and':
            right, idx = parse_factor(idx+1)
            clauses = [c1 + c2 for c1, c2 in product(clauses, right)]
        return clauses, idx

    def parse_factor(idx):
        tok = tokens[idx]
        if tok == '(':
            clauses, idx = parse_expression(idx+1)
            if tokens[idx] != ')':
                raise ValueError(f"Mismatched parentheses in GPR: {tokens}")
            return clauses, idx+1
        else:
            return [[tok]], idx+1

    dnf, _ = parse_expression(0)
    # dedupe and sort
    unique = []
    for clause in dnf:
        cl = sorted(set(clause))
        if cl not in unique:
            unique.append(cl)
    return unique


def run_optimization(
    model,
    kcat_dict,
    objective_reaction,
    gene_sequences_dict=None,
    enzyme_upper_bound=0.125,
    enzyme_ratio=True,
    maximization=True,
    multi_enzyme_off=False,
    isoenzymes_off=False,
    promiscuous_off=False,
    complexes_off=False,
    solver_name='gurobi',
    tee=False,
    verbose=False,
    medium=None,
    medium_upper_bound=False,
    bidirectional_constraints=True,
):
    """
    Enzyme-constrained FBA via Pyomo, handling:
      - single enzymes (max kcat)
      - enzyme complexes (avg kcat)
      - isoenzymes (OR-GPR)
      - promiscuous enzymes
      - bidirectional constraints (substrate-specific kcats for reversible reactions)
    Returns: sol_val, df_FBA, gene_sequences_dict, model
    """

    # Initialize timing
    overall_start = time.time()
    step_start = time.time()

    if verbose:
        initial_memory = _get_memory_usage()
        print(f"\n{'='*60}")
        print("ENZYME-CONSTRAINED OPTIMIZATION - PERFORMANCE MONITORING")
        print(f"{'='*60}")
        if initial_memory > 0:
            print(f"Initial memory usage: {initial_memory:.1f} MB")
        print(f"Solver: {solver_name}")
        print(f"Enzyme upper bound: {enzyme_upper_bound}")
        print(f"Bidirectional constraints: {bidirectional_constraints}")
        print()

    # 1) Load COBRA model
    step_start = time.time()
    if isinstance(model, str):
        mod = (
            cb.io.read_sbml_model(model)
            if model.endswith(('.xml', '.sbml'))
            else cb.io.load_json_model(model)
        )
    else:
        mod = model
    _print_step_timing("Step 1: Load COBRA model", step_start, verbose)

    # 1a) Apply medium conditions if provided
    if medium is not None:
        step_start = time.time()
        for rxn_id, flux_value in medium.items():
            try:
                rxn = mod.reactions.get_by_id(rxn_id)
                rxn.lower_bound = flux_value
                if medium_upper_bound:
                    rxn.upper_bound = flux_value
                if verbose:
                    if medium_upper_bound:
                        print(f"  Fixed {rxn_id}: lower={flux_value}, upper={flux_value}")
                    else:
                        print(f"  Set {rxn_id}: lower={flux_value}, upper={rxn.upper_bound}")
            except KeyError:
                print(f"  Warning: Reaction provided '{rxn_id}' was not found in model")
        _print_step_timing("Step 1a: Apply medium conditions", step_start, verbose)

    # 2) Initial flux guess
    step_start = time.time()
    mod.objective = objective_reaction
    if verbose:
        print("Model objective:", mod.objective)
    cobra_sol = mod.optimize()
    flux0 = cobra_sol.fluxes.to_dict()
    _print_step_timing("Step 2: Get initial flux guess (COBRApy optimization)", step_start, verbose)

    # 3) Load & normalize kcat_dict
    step_start = time.time()
    if isinstance(kcat_dict, str):
        df = pd.read_csv(kcat_dict)
        tmp = {}
        for r, g, k in zip(df.reaction, df.gene, df.kcat):
            tmp.setdefault((r, g), []).append(k)
        kcat_dict = tmp
    # convert all values to hr^-1 lists
    # IMPORTANT: We assume all incoming kcats are in s⁻¹ and need conversion to hr⁻¹
    # The old heuristic (if v < 1000) was unreliable during simulated annealing
    for key, vals in list(kcat_dict.items()):
        kcat_dict[key] = [v * 3600 for v in vals]
    _print_step_timing(f"Step 3: Process kcat dictionary ({len(kcat_dict)} entries)", step_start, verbose)

    # 4) Build stoichiometry, bounds, objective
    step_start = time.time()
    S = create_stoichiometric_matrix(mod)
    mets = [m.id for m in mod.metabolites]
    rxns = [r.id for r in mod.reactions]
    genes = [g.id for g in mod.genes]
    lb = {r.id: r.lower_bound for r in mod.reactions}
    ub = {r.id: r.upper_bound for r in mod.reactions}

    # DIAGNOSTIC: Print the bounds that will be used in optimization
    # if medium is not None:
    #     print("\n=== DIAGNOSTIC: Checking captured bounds ===")
    #     for rxn_id in medium.keys():
    #         if rxn_id in lb:
    #             print(f"  {rxn_id} lower bound in lb dict: {lb[rxn_id]:.4f}")
    #         if rxn_id in ub:
    #             print(f"  {rxn_id} upper bound in ub dict: {ub[rxn_id]:.4f}")
    #         else:
    #             print(f"  {rxn_id} NOT FOUND in lb or ub dict!")
    #     print("=== END DIAGNOSTIC ===\n")
    obj_coef = {r.id: (1.0 if r.id == objective_reaction else 0.0) for r in mod.reactions} #obj_coef = {r.id: r.objective_coefficient for r in mod.reactions}
    met_index = {m: i for i, m in enumerate(mets)}
    rxn_index = {r: j for j, r in enumerate(rxns)}
    _print_step_timing(f"Step 4: Build model data structures ({len(mets)} mets, {len(rxns)} rxns, {len(genes)} genes)", step_start, verbose)

    # 5) Generate DNF clauses
    step_start = time.time()
    dnf_clauses = {}
    for r in mod.reactions:
        rule = (r.gene_reaction_rule or '').strip()
        if not rule:
            continue
        tokens = _tokenize_gpr(rule)
        dnf_clauses[r.id] = _parse_gpr_to_dnf(tokens)
    _print_step_timing(f"Step 5: Parse gene-protein-reaction rules ({len(dnf_clauses)} rules)", step_start, verbose)

    # 6) Build Pyomo model
    step_start = time.time()
    m = ConcreteModel()  # noqa: F405
    m.M = Set(initialize=mets)  # noqa: F405
    m.R = Set(initialize=rxns)  # noqa: F405
    m.G = Set(initialize=genes)  # noqa: F405
    _print_step_timing("Step 6a: Initialize Pyomo model sets", step_start, verbose)

    # For bidirectional constraints, identify reversible reactions
    if bidirectional_constraints:
        step_start = time.time()
        # Find reversible reactions with both forward and reverse kcat data
        reversible_reactions = set()
        irreversible_reactions = set(rxns)  # Start with all reactions

        for key in kcat_dict.keys():
            if len(key) == 3:  # (reaction, gene, direction) format
                rxn_id, gene_id, direction = key
                if direction in ['forward', 'reverse']:
                    reversible_reactions.add(rxn_id)

        # Filter to only those with both directions
        true_reversible = set()
        for rxn_id in reversible_reactions:
            has_forward = any(key[0] == rxn_id and key[2] == 'forward' for key in kcat_dict.keys() if len(key) == 3)
            has_reverse = any(key[0] == rxn_id and key[2] == 'reverse' for key in kcat_dict.keys() if len(key) == 3)
            if has_forward and has_reverse:
                true_reversible.add(rxn_id)
                irreversible_reactions.discard(rxn_id)

        m.R_rev = Set(initialize=list(true_reversible))  # noqa: F405
        m.R_irr = Set(initialize=list(irreversible_reactions))  # noqa: F405

        if verbose:
            print(f"Identified {len(true_reversible)} truly reversible reactions with bidirectional kcat data")
            print(f"Treating {len(irreversible_reactions)} reactions as irreversible")
        _print_step_timing(f"Step 6b: Analyze bidirectional reactions ({len(true_reversible)} reversible)", step_start, verbose)
    else:
        m.R_rev = Set(initialize=[])  # noqa: F405
        m.R_irr = Set(initialize=rxns)  # noqa: F405

    # OPTIMIZED: Only include (rxn, gene) pairs that have kcat data
    # This avoids creating millions of unnecessary constraint checks
    if bidirectional_constraints:
        # For bidirectional constraints, we have 3-tuples (rxn, gene, direction)
        m.K = Set(initialize=list(kcat_dict.keys()), dimen=3)  # noqa: F405
    else:
        # For standard constraints, we have 2-tuples (rxn, gene)
        m.K = Set(initialize=list(kcat_dict.keys()), dimen=2)  # noqa: F405

    # Variables
    step_start = time.time()
    if bidirectional_constraints:
        # Optimize bounds computation - calculate once, reuse
        rxn_ub_pos = {r: max(0, ub[r]) for r in m.R_rev}
        rxn_lb_neg = {r: max(0, -lb[r]) for r in m.R_rev}

        # Create separate forward and reverse variables for reversible reactions
        def v_fwd_bounds(mo, j):
            return (0, rxn_ub_pos.get(j, 0))

        def v_rev_bounds(mo, j):
            return (0, rxn_lb_neg.get(j, 0))

        def v_irr_bounds(mo, j):
            if j in m.R_irr:
                return (lb[j], ub[j])
            else:
                return (0, 0)

        # Forward and reverse flux variables for reversible reactions
        m.v_fwd = Var(m.R_rev, domain=NonNegativeReals, bounds=v_fwd_bounds,  # noqa: F405
                      initialize=lambda mo, j: max(0, flux0.get(j, 0.0)))
        m.v_rev = Var(m.R_rev, domain=NonNegativeReals, bounds=v_rev_bounds,  # noqa: F405
                      initialize=lambda mo, j: max(0, -flux0.get(j, 0.0)))

        # Standard flux variables for irreversible reactions
        m.v_irr = Var(m.R_irr, domain=Reals, bounds=v_irr_bounds,  # noqa: F405
                      initialize=lambda mo, j: flux0.get(j, 0.0))

        # OPTION C: Direction-specific enzyme allocation for reversible reactions
        # Each gene can allocate enzyme to forward OR reverse direction, but not both simultaneously
        m.E_fwd = Var(m.G, m.R_rev, domain=NonNegativeReals, initialize=0.001)  # noqa: F405
        m.E_rev = Var(m.G, m.R_rev, domain=NonNegativeReals, initialize=0.001)  # noqa: F405

        # Standard enzyme variables
        m.E = Var(m.G, domain=NonNegativeReals, initialize=0.001)  # noqa: F405

        # Constraint: direction-specific enzyme allocation cannot exceed total enzyme
        def enzyme_allocation_constraint(mo, g, j):
            if j in mo.R_rev:
                return mo.E_fwd[g, j] + mo.E_rev[g, j] <= mo.E[g]
            else:
                return Constraint.Feasible  # noqa: F405

        m.enzyme_allocation = Constraint(m.G, m.R, rule=enzyme_allocation_constraint)  # noqa: F405

    else:
        # Standard single flux variable for all reactions
        m.v = Var(  # noqa: F405
            m.R,
            domain=Reals,  # noqa: F405
            bounds=lambda mo, j: (lb[j], ub[j]),
            initialize=lambda mo, j: flux0.get(j, 0.0)
        )
        # Standard enzyme variables
        m.E = Var(m.G, domain=NonNegativeReals, initialize=0.001)  # noqa: F405

    _print_step_timing("Step 6c: Create variables (flux and enzyme)", step_start, verbose)

    # Mass balance
    step_start = time.time()
    # OPTIMIZED mass balance

    # Pre-compute which reactions involve each metabolite
    met_reactions = {}
    for met in mets:
        met_reactions[met] = []
        for rxn in rxns:
            if S[met_index[met], rxn_index[rxn]] != 0:
                met_reactions[met].append(rxn)

    # print(f"Sparse matrix analysis: avg {sum(len(v) for v in met_reactions.values()) / len(mets):.1f} reactions per metabolite")

    def mass_balance_sparse(mo, met):
        # Only iterate over reactions that involve this metabolite
        relevant_rxns = met_reactions[met]
        if not relevant_rxns:
            return Constraint.Feasible  # noqa: F405

        i = met_index[met]

        if bidirectional_constraints:
            # Sum contributions from reversible and irreversible reactions
            flux_terms = []
            for r in relevant_rxns:
                stoich = S[i, rxn_index[r]]
                if r in mo.R_rev:
                    # For reversible reactions: net flux = v_fwd - v_rev
                    flux_terms.append(stoich * (mo.v_fwd[r] - mo.v_rev[r]))
                elif r in mo.R_irr:
                    # For irreversible reactions: use standard flux variable
                    flux_terms.append(stoich * mo.v_irr[r])
            return sum(flux_terms) == 0
        else:
            # Standard mass balance
            return sum(S[i, rxn_index[r]] * mo.v[r] for r in relevant_rxns) == 0

    m.mass_balance = Constraint(m.M, rule=mass_balance_sparse)  # noqa: F405
    _print_step_timing(f"Step 6d: Create mass balance constraints ({len(mets)} constraints)", step_start, verbose)

    # Objective
    step_start = time.time()
    sense = maximize if maximization else minimize  # noqa: F405
    if bidirectional_constraints:
        # Objective using net flux for reversible reactions and standard flux for irreversible
        obj_terms = []
        for r in rxns:
            coeff = obj_coef[r]
            if coeff != 0:  # Only include reactions with non-zero objective coefficient
                if r in m.R_rev:
                    # Net flux for reversible reactions
                    obj_terms.append(coeff * (m.v_fwd[r] - m.v_rev[r]))
                elif r in m.R_irr:
                    # Standard flux for irreversible reactions
                    obj_terms.append(coeff * m.v_irr[r])
        m.obj = Objective(expr=sum(obj_terms) if obj_terms else 0, sense=sense)  # noqa: F405
    else:
        # Standard objective
        m.obj = Objective(expr=sum(obj_coef[r] * m.v[r] for r in m.R), sense=sense)  # noqa: F405
    _print_step_timing("Step 6e: Create objective function", step_start, verbose)

    # Initialize constraint counters
    step_start = time.time()
    constraints_added = 0
    constraints_skipped = 0
    skip_reasons = {'no_clauses': 0, 'gene_mismatch': 0, 'no_kcat': 0, 'multiple_clauses': 0}
    iso_added = 0
    iso_skipped = 0
    iso_skip_reasons = {'single_clause': 0, 'has_complex': 0, 'bidirectional_rev': 0, 'no_terms': 0}
    mixed_added = 0
    mixed_skipped = 0
    mixed_skip_reasons = {'single_clause': 0, 'no_complex': 0, 'bidirectional_rev': 0, 'no_capacities': 0}
    promis_added = 0
    promis_skipped = 0

    # Bidirectional constraint counters
    bidirectional_iso_added = 0
    bidirectional_iso_skipped = 0
    bidirectional_mixed_added = 0
    bidirectional_mixed_skipped = 0

    if verbose:
        print("\n  Starting enzyme constraint generation...")
        print(f"  Available constraint types: AND={not multi_enzyme_off}, ISO={not isoenzymes_off}, MIXED={not (isoenzymes_off or complexes_off)}, PROMISCUOUS={not promiscuous_off}")
        print(f"  kcat entries: {len(kcat_dict)}, GPR rules: {len(dnf_clauses)}")
    _print_step_timing("Step 7a: Initialize enzyme constraint generation", step_start, verbose)

    # 6a) AND‐GPR: single or complex (simple cases only)
    def and_rule(mo, *args):
        nonlocal constraints_added, constraints_skipped

        if bidirectional_constraints:
            # For bidirectional: args = (rxn_id, gene_id, direction)
            if len(args) != 3:
                constraints_skipped += 1
                return Constraint.Skip  # noqa: F405
            rxn_id, gene_id, direction = args
            kcat_key = (rxn_id, gene_id, direction)
        else:
            # For standard: args = (rxn_id, gene_id)
            if len(args) != 2:
                constraints_skipped += 1
                return Constraint.Skip  # noqa: F405
            rxn_id, gene_id = args
            kcat_key = (rxn_id, gene_id)

        clauses = dnf_clauses.get(rxn_id, [])
        if not clauses:
            constraints_skipped += 1
            skip_reasons['no_clauses'] += 1
            return Constraint.Skip  # noqa: F405

        # Get kcat values
        k_list = kcat_dict.get(kcat_key, [])
        if not k_list:
            constraints_skipped += 1
            skip_reasons['no_kcat'] += 1
            return Constraint.Skip  # noqa: F405

        # single enzyme: max kcat
        if len(clauses) == 1 and len(clauses[0]) == 1:
            g = clauses[0][0]
            if gene_id != g:
                constraints_skipped += 1
                skip_reasons['gene_mismatch'] += 1
                return Constraint.Skip  # noqa: F405

            k_val = max(k_list)
            constraints_added += 1

            # Choose the appropriate flux variable based on reaction type and direction
            if bidirectional_constraints and rxn_id in mo.R_rev:
                # For reversible reactions with bidirectional constraints
                # Use direction-specific enzyme allocation (Option C)
                if direction == 'forward':
                    return mo.v_fwd[rxn_id] <= k_val * mo.E_fwd[g, rxn_id]
                elif direction == 'reverse':
                    return mo.v_rev[rxn_id] <= k_val * mo.E_rev[g, rxn_id]
                else:
                    constraints_skipped += 1
                    return Constraint.Skip  # noqa: F405
            elif bidirectional_constraints and rxn_id in mo.R_irr:
                # For irreversible reactions in bidirectional mode
                return mo.v_irr[rxn_id] <= k_val * mo.E[g]
            else:
                # Standard mode
                return mo.v[rxn_id] <= k_val * mo.E[g]
        # enzyme complex: avg kcat
        if len(clauses) == 1 and len(clauses[0]) > 1:
            clause = clauses[0]
            if gene_id not in clause:
                constraints_skipped += 1
                skip_reasons['gene_mismatch'] += 1
                return Constraint.Skip  # noqa: F405
            # Skip if complexes are disabled
            if complexes_off:
                constraints_skipped += 1
                skip_reasons['no_kcat'] += 1  # Reuse this category
                return Constraint.Skip  # noqa: F405

            # For complexes in bidirectional mode, we need to handle each direction separately
            if bidirectional_constraints:
                # Get kcats for this specific direction and gene
                all_ks = k_list  # Already filtered for this direction and gene
            else:
                # Standard mode: collect all kcats for all genes in complex
                all_ks = []
                for g in clause:
                    all_ks.extend(kcat_dict.get((rxn_id, g), []))

            if not all_ks:
                constraints_skipped += 1
                skip_reasons['no_kcat'] += 1
                return Constraint.Skip  # noqa: F405

            k_val = sum(all_ks) / len(all_ks)
            constraints_added += 1

            # Choose the appropriate flux variable
            if bidirectional_constraints and rxn_id in mo.R_rev:
                # For enzyme complexes in reversible reactions
                # Use direction-specific enzyme allocation (Option C)
                if direction == 'forward':
                    return mo.v_fwd[rxn_id] <= k_val * mo.E_fwd[gene_id, rxn_id]
                elif direction == 'reverse':
                    return mo.v_rev[rxn_id] <= k_val * mo.E_rev[gene_id, rxn_id]
                else:
                    constraints_skipped += 1
                    return Constraint.Skip  # noqa: F405
            elif bidirectional_constraints and rxn_id in mo.R_irr:
                return mo.v_irr[rxn_id] <= k_val * mo.E[gene_id]
            else:
                return mo.v[rxn_id] <= k_val * mo.E[gene_id]
        # Multiple clauses - will be handled by ISO or MIXED constraints
        constraints_skipped += 1
        skip_reasons['multiple_clauses'] += 1
        return Constraint.Skip  # noqa: F405

    if not multi_enzyme_off:
        step_start = time.time()
        m.kcat_and = Constraint(m.K, rule=and_rule)  # noqa: F405
        _print_step_timing(f"Step 7b: AND constraints ({constraints_added} added, {constraints_skipped} skipped)", step_start, verbose)

    # 6b) OR‐GPR: pure isoenzymes (no AND within clauses)
    def iso_rule(mo, rxn_id):
        nonlocal iso_added, iso_skipped, iso_skip_reasons
        clauses = dnf_clauses.get(rxn_id, [])
        if len(clauses) <= 1:
            iso_skipped += 1
            iso_skip_reasons['single_clause'] += 1
            return Constraint.Skip  # noqa: F405

        # Check if this is pure isoenzymes (all clauses have single genes)
        # If any clause has multiple genes, it's a mixed case - skip for now
        has_complex_clause = any(len(clause) > 1 for clause in clauses)
        if has_complex_clause:
            iso_skipped += 1
            iso_skip_reasons['has_complex'] += 1
            return Constraint.Skip  # noqa: F405

        # For bidirectional constraints, we need to create separate constraints for forward and reverse
        if bidirectional_constraints and rxn_id in mo.R_rev:
            # Handle separately for forward and reverse directions
            iso_skipped += 1
            iso_skip_reasons['bidirectional_rev'] += 1
            return Constraint.Skip  # noqa: F405  # Will be handled by bidirectional iso rules

        # Pure isoenzymes: g1 or g2 or g3 (all single genes)
        terms = []
        for clause in clauses:
            g = clause[0]  # Single gene in this clause
            if bidirectional_constraints:
                # For irreversible reactions in bidirectional mode, use appropriate key format
                # We need to check if we have direction-specific data or standard data
                forward_key = (rxn_id, g, 'forward')
                standard_key = (rxn_id, g)

                if forward_key in kcat_dict:
                    kl = kcat_dict[forward_key]
                elif standard_key in kcat_dict:
                    kl = kcat_dict[standard_key]
                else:
                    kl = []
            else:
                kl = kcat_dict.get((rxn_id, g), [])

            if kl:
                kmin = min(kl)
                if bidirectional_constraints and rxn_id in mo.R_irr:
                    terms.append(kmin * mo.E[g])
                elif not bidirectional_constraints:
                    terms.append(kmin * mo.E[g])

        if not terms:
            iso_skipped += 1
            iso_skip_reasons['no_terms'] += 1
            return Constraint.Skip  # noqa: F405

        iso_added += 1
        if bidirectional_constraints and rxn_id in mo.R_irr:
            return mo.v_irr[rxn_id] <= sum(terms)
        else:
            return mo.v[rxn_id] <= sum(terms)

    if not isoenzymes_off:
        step_start = time.time()
        m.kcat_iso = Constraint(m.R, rule=iso_rule)  # noqa: F405
        _print_step_timing(f"Step 7c: ISO constraints ({iso_added} added, {iso_skipped} skipped)", step_start, verbose)

    # 6b2) Bidirectional isoenzyme constraints for reversible reactions
    if bidirectional_constraints and not isoenzymes_off:
        def iso_bidirectional_rule(mo, rxn_id, direction):
            nonlocal bidirectional_iso_added, bidirectional_iso_skipped
            clauses = dnf_clauses.get(rxn_id, [])
            if len(clauses) <= 1:
                bidirectional_iso_skipped += 1
                return Constraint.Skip  # noqa: F405

            # Check if this is pure isoenzymes
            has_complex_clause = any(len(clause) > 1 for clause in clauses)
            if has_complex_clause:
                bidirectional_iso_skipped += 1
                return Constraint.Skip  # noqa: F405

            # Pure isoenzymes for specific direction
            terms = []
            for clause in clauses:
                g = clause[0]  # Single gene
                kl = kcat_dict.get((rxn_id, g, direction), [])
                if kl:
                    kmin = min(kl)
                    # Use direction-specific enzyme allocation (Option C)
                    if direction == 'forward':
                        terms.append(kmin * mo.E_fwd[g, rxn_id])
                    else:  # reverse
                        terms.append(kmin * mo.E_rev[g, rxn_id])

            if not terms:
                bidirectional_iso_skipped += 1
                return Constraint.Skip  # noqa: F405

            bidirectional_iso_added += 1
            if direction == 'forward':
                return mo.v_fwd[rxn_id] <= sum(terms)
            else:  # reverse
                return mo.v_rev[rxn_id] <= sum(terms)

        # Create constraints for both directions of reversible reactions
        bidirectional_iso_keys = []
        for rxn_id in m.R_rev:
            # Check if we have isoenzyme data for both directions
            clauses = dnf_clauses.get(rxn_id, [])
            if len(clauses) > 1 and all(len(clause) == 1 for clause in clauses):
                has_forward = any((rxn_id, clause[0], 'forward') in kcat_dict for clause in clauses)
                has_reverse = any((rxn_id, clause[0], 'reverse') in kcat_dict for clause in clauses)
                if has_forward:
                    bidirectional_iso_keys.append((rxn_id, 'forward'))
                if has_reverse:
                    bidirectional_iso_keys.append((rxn_id, 'reverse'))

        if bidirectional_iso_keys:
            m.BI = Set(initialize=bidirectional_iso_keys, dimen=2)  # noqa: F405
            m.kcat_iso_bidirectional = Constraint(m.BI, rule=iso_bidirectional_rule)  # noqa: F405

    # 6c) MIXED OR+AND: Handle complex nested cases like (g1 and g2) or (g3 and g4) or g5
    def mixed_rule(mo, rxn_id):
        nonlocal mixed_added, mixed_skipped, mixed_skip_reasons
        clauses = dnf_clauses.get(rxn_id, [])
        if len(clauses) <= 1:
            mixed_skipped += 1
            mixed_skip_reasons['single_clause'] += 1
            return Constraint.Skip  # noqa: F405

        # Only handle if at least one clause has multiple genes (complex)
        has_complex_clause = any(len(clause) > 1 for clause in clauses)
        if not has_complex_clause:
            mixed_skipped += 1
            mixed_skip_reasons['no_complex'] += 1
            return Constraint.Skip  # noqa: F405

        # Choose appropriate flux variable based on bidirectional constraints
        if bidirectional_constraints and rxn_id in mo.R_rev:
            # For reversible reactions in bidirectional mode, we need separate constraints
            # Skip here and let bidirectional mixed constraints handle this
            mixed_skipped += 1
            mixed_skip_reasons['bidirectional_rev'] += 1
            return Constraint.Skip  # noqa: F405  # Let bidirectional constraints handle this
        elif bidirectional_constraints and rxn_id in mo.R_irr:
            flux_var = mo.v_irr[rxn_id]
        else:
            flux_var = mo.v[rxn_id]

        # NESTED APPROACH: Handle each clause based on its structure
        # 1. Single genes (isoenzymes): use their individual capacity
        # 2. Multi-gene clauses (complexes): handle as enzyme complexes
        # 3. Sum all alternative capacities (OR relationship between clauses)

        def get_clause_capacity(clause):
            """Calculate the capacity contribution of a single clause"""
            if len(clause) == 1:
                # Single gene - this is an isoenzyme alternative
                g = clause[0]
                if bidirectional_constraints:
                    # Check for direction-specific data first
                    forward_key = (rxn_id, g, 'forward')
                    reverse_key = (rxn_id, g, 'reverse')
                    standard_key = (rxn_id, g)

                    if forward_key in kcat_dict:
                        kl = kcat_dict[forward_key]
                    elif reverse_key in kcat_dict:
                        kl = kcat_dict[reverse_key]
                    elif standard_key in kcat_dict:
                        kl = kcat_dict[standard_key]
                    else:
                        kl = []
                else:
                    kl = kcat_dict.get((rxn_id, g), [])

                if kl:
                    # For isoenzymes, use minimum kcat (most conservative)
                    kcat_val = min(kl)
                    return kcat_val, mo.E[g]  # Return (kcat, enzyme_expr) tuple
                else:
                    return None  # No kcat data available

            else:
                # Multi-gene clause - this is an enzyme complex
                # All genes must be present for the complex to function
                complex_kcats = []

                for g in clause:
                    if bidirectional_constraints:
                        # Check for direction-specific data
                        forward_key = (rxn_id, g, 'forward')
                        reverse_key = (rxn_id, g, 'reverse')
                        standard_key = (rxn_id, g)

                        if forward_key in kcat_dict:
                            kl = kcat_dict[forward_key]
                        elif reverse_key in kcat_dict:
                            kl = kcat_dict[reverse_key]
                        elif standard_key in kcat_dict:
                            kl = kcat_dict[standard_key]
                        else:
                            kl = []
                    else:
                        kl = kcat_dict.get((rxn_id, g), [])

                    if kl:
                        complex_kcats.extend(kl)

                if not complex_kcats:
                    return None  # No kcat data for this complex

                # For enzyme complexes:
                # - Use average kcat across all subunits
                # - The limiting factor is the minimum enzyme amount among subunits
                # - But we need to handle this in a linearizable way

                avg_kcat = sum(complex_kcats) / len(complex_kcats)

                # Option C implementation: Direction-specific enzyme allocation
                # For now, use a simplified approach that approximates the complex constraint
                # The capacity is limited by the average enzyme availability
                # This is an approximation: in reality, we'd need min(E[g] for g in clause)
                avg_enzyme = sum(mo.E[g] for g in clause) / len(clause)

                return avg_kcat, avg_enzyme  # Return (kcat, enzyme_expr) tuple

        # Calculate total capacity as sum of all clause capacities (OR relationship)
        clause_capacities = []
        for clause in clauses:
            capacity_result = get_clause_capacity(clause)
            if capacity_result is not None:  # Only include clauses with valid kcat data
                kcat_val, enzyme_expr = capacity_result
                clause_capacities.append(kcat_val * enzyme_expr)

        if not clause_capacities:
            mixed_skipped += 1
            mixed_skip_reasons['no_capacities'] += 1
            return Constraint.Skip  # noqa: F405

        mixed_added += 1
        return flux_var <= sum(clause_capacities)

    # Only add mixed constraints if neither isoenzymes nor complexes are disabled
    if not isoenzymes_off and not complexes_off:
        step_start = time.time()
        m.kcat_mixed = Constraint(m.R, rule=mixed_rule)  # noqa: F405
        _print_step_timing(f"Step 7d: MIXED constraints ({mixed_added} added, {mixed_skipped} skipped)", step_start, verbose)

    # 6c2) Bidirectional MIXED constraints for reversible reactions
    if bidirectional_constraints and not isoenzymes_off and not complexes_off:
        def mixed_bidirectional_rule(mo, rxn_id, direction):
            nonlocal bidirectional_mixed_added, bidirectional_mixed_skipped
            clauses = dnf_clauses.get(rxn_id, [])
            if len(clauses) <= 1:
                bidirectional_mixed_skipped += 1
                return Constraint.Skip  # noqa: F405

            # Only handle if at least one clause has multiple genes (complex)
            has_complex_clause = any(len(clause) > 1 for clause in clauses)
            if not has_complex_clause:
                bidirectional_mixed_skipped += 1
                return Constraint.Skip  # noqa: F405

            def get_bidirectional_clause_capacity(clause):
                """Calculate capacity for a specific direction"""
                if len(clause) == 1:
                    # Single gene - isoenzyme alternative
                    g = clause[0]
                    kl = kcat_dict.get((rxn_id, g, direction), [])
                    if kl:
                        kcat_val = min(kl)
                        # Use direction-specific enzyme allocation (Option C)
                        if direction == 'forward':
                            return kcat_val, mo.E_fwd[g, rxn_id]
                        else:  # reverse
                            return kcat_val, mo.E_rev[g, rxn_id]
                    return None
                else:
                    # Multi-gene clause - enzyme complex
                    complex_kcats = []
                    for g in clause:
                        kl = kcat_dict.get((rxn_id, g, direction), [])
                        if kl:
                            complex_kcats.extend(kl)

                    if not complex_kcats:
                        return None

                    avg_kcat = sum(complex_kcats) / len(complex_kcats)

                    # For complexes, use average of direction-specific allocations
                    if direction == 'forward':
                        avg_enzyme = sum(mo.E_fwd[g, rxn_id] for g in clause) / len(clause)
                    else:  # reverse
                        avg_enzyme = sum(mo.E_rev[g, rxn_id] for g in clause) / len(clause)

                    return avg_kcat, avg_enzyme

            # Calculate capacity for this direction
            clause_capacities = []
            for clause in clauses:
                capacity_result = get_bidirectional_clause_capacity(clause)
                if capacity_result is not None:
                    kcat_val, enzyme_expr = capacity_result
                    clause_capacities.append(kcat_val * enzyme_expr)

            if not clause_capacities:
                bidirectional_mixed_skipped += 1
                return Constraint.Skip  # noqa: F405

            bidirectional_mixed_added += 1
            if direction == 'forward':
                return mo.v_fwd[rxn_id] <= sum(clause_capacities)
            else:  # reverse
                return mo.v_rev[rxn_id] <= sum(clause_capacities)

        # Create bidirectional mixed constraints
        bidirectional_mixed_keys = []
        for rxn_id in m.R_rev:
            clauses = dnf_clauses.get(rxn_id, [])
            if len(clauses) > 1 and any(len(clause) > 1 for clause in clauses):
                # Check if we have mixed constraint data for both directions
                has_data = {'forward': False, 'reverse': False}
                for clause in clauses:
                    for g in clause:
                        if (rxn_id, g, 'forward') in kcat_dict:
                            has_data['forward'] = True
                        if (rxn_id, g, 'reverse') in kcat_dict:
                            has_data['reverse'] = True

                if has_data['forward']:
                    bidirectional_mixed_keys.append((rxn_id, 'forward'))
                if has_data['reverse']:
                    bidirectional_mixed_keys.append((rxn_id, 'reverse'))

        if bidirectional_mixed_keys:
            m.BM = Set(initialize=bidirectional_mixed_keys, dimen=2)  # noqa: F405
            m.kcat_mixed_bidirectional = Constraint(m.BM, rule=mixed_bidirectional_rule)  # noqa: F405

    # 6d) Promiscuous enzymes
    def promis_rule(mo, g_id):
        nonlocal promis_added, promis_skipped
        usage = []
        for r_id, clauses in dnf_clauses.items():
            for clause in clauses:
                if g_id not in clause:
                    continue

                # Handle bidirectional vs standard constraints
                if bidirectional_constraints:
                    # For bidirectional constraints with Option C, use direction-specific allocation
                    # We need to account for usage in both forward and reverse directions

                    # Forward direction usage
                    if r_id in mo.R_rev:
                        # Reversible reaction: use v_fwd
                        for k in kcat_dict.get((r_id, g_id, 'forward'), []):
                            usage.append(mo.v_fwd[r_id] / k)
                    elif r_id in mo.R_irr:
                        # Irreversible reaction: use v_irr for forward direction
                        for k in kcat_dict.get((r_id, g_id, 'forward'), []):
                            usage.append(mo.v_irr[r_id] / k)
                        # Also check for standard kcat data (non-directional)
                        for k in kcat_dict.get((r_id, g_id), []):
                            usage.append(mo.v_irr[r_id] / k)

                    # Reverse direction usage (only for reversible reactions)
                    if r_id in mo.R_rev:
                        for k in kcat_dict.get((r_id, g_id, 'reverse'), []):
                            usage.append(mo.v_rev[r_id] / k)

                else:
                    # Standard constraints: use standard flux variable
                    for k in kcat_dict.get((r_id, g_id), []):
                        usage.append(mo.v[r_id] / k)

        if not usage:
            promis_skipped += 1
            return Constraint.Skip  # noqa: F405
        promis_added += 1
        return sum(usage) <= mo.E[g_id]

    if not promiscuous_off:
        step_start = time.time()
        m.promiscuous = Constraint(m.G, rule=promis_rule)  # noqa: F405
        _print_step_timing(f"Step 7e: Promiscuous constraints ({promis_added} added, {promis_skipped} skipped)", step_start, verbose)

    # Note: Bidirectional constraints are now handled in the main and_rule and iso_bidirectional_rule functions
    # This legacy bidirectional constraint code has been removed to avoid conflicts    # Print constraint summary (only if verbose)
    if verbose:
        print(f"\n{'='*60}")
        print("ENZYME CONSTRAINT SUMMARY")
        print(f"{'='*60}")
        if multi_enzyme_off:
            print("AND constraints (single/complex): DISABLED")
        else:
            print(f"AND constraints (single/complex):  {constraints_added:4d} added, {constraints_skipped:4d} skipped")
            if constraints_skipped > 0:
                print(f"  Skip reasons: {skip_reasons}")

        if isoenzymes_off:
            print("OR/ISO constraints (pure isozymes): DISABLED")
        else:
            print(f"OR/ISO constraints (pure isozymes): {iso_added:4d} added, {iso_skipped:4d} skipped")
            if iso_skipped > 0:
                print(f"  Skip reasons: {iso_skip_reasons}")

        if isoenzymes_off or complexes_off:
            print("MIXED OR+AND constraints: DISABLED")
        else:
            print(f"MIXED OR+AND constraints:          {mixed_added:4d} added, {mixed_skipped:4d} skipped")
            if mixed_skipped > 0:
                print(f"  Skip reasons: {mixed_skip_reasons}")

        if promiscuous_off:
            print("Promiscuous enzyme constraints: DISABLED")
        else:
            print(f"Promiscuous enzyme constraints:    {promis_added:4d} added, {promis_skipped:4d} skipped")

        if bidirectional_constraints:
            print("Bidirectional constraints: ENABLED (integrated with and_rule and iso constraints)")
            print(f"Bidirectional ISO constraints:     {bidirectional_iso_added:4d} added, {bidirectional_iso_skipped:4d} skipped")
            print(f"Bidirectional MIXED constraints:   {bidirectional_mixed_added:4d} added, {bidirectional_mixed_skipped:4d} skipped")

            # Calculate total bidirectional coverage
            total_reversible = len(m.R_rev) if hasattr(m, 'R_rev') else 0
            bidirectional_covered = bidirectional_iso_added + bidirectional_mixed_added
            print(f"Bidirectional coverage:            {bidirectional_covered}/{total_reversible*2} direction-constraints")
            print(f"  (Expected max: {total_reversible} reversible × 2 directions)")
        else:
            print("Bidirectional constraints: DISABLED")

        total_active = (
            (constraints_added if not multi_enzyme_off else 0) +
            (iso_added if not isoenzymes_off else 0) +
            (mixed_added if not (isoenzymes_off or complexes_off) else 0) +
            (promis_added if not promiscuous_off else 0) +
            (bidirectional_iso_added if bidirectional_constraints else 0) +
            (bidirectional_mixed_added if bidirectional_constraints else 0)
        )
        print(f"{'='*60}")
        print(f"Total active enzyme constraints:   {total_active:4d}")
        print(f"{'='*60}\n")

        # Enhanced diagnostic for understanding gaps
        if bidirectional_constraints:
            print("BIDIRECTIONAL CONSTRAINT ANALYSIS:")
            print(f"  Reversible reactions identified:   {len(m.R_rev)}")
            print(f"  Irreversible reactions:            {len(m.R_irr)}")

            # Analyze what types of GPR structures exist
            gpr_types = {'single': 0, 'pure_iso': 0, 'mixed': 0, 'none': 0}
            for rxn_id in m.R_rev:
                clauses = dnf_clauses.get(rxn_id, [])
                if not clauses:
                    gpr_types['none'] += 1
                elif len(clauses) == 1 and len(clauses[0]) == 1:
                    gpr_types['single'] += 1
                elif len(clauses) > 1 and all(len(c) == 1 for c in clauses):
                    gpr_types['pure_iso'] += 1
                elif any(len(c) > 1 for c in clauses):
                    gpr_types['mixed'] += 1

            print("  GPR structure breakdown for reversible reactions:")
            print(f"    Single gene:     {gpr_types['single']:4d} (handled by AND constraints)")
            print(f"    Pure isoenzymes: {gpr_types['pure_iso']:4d} (should be bidirectional ISO)")
            print(f"    Mixed complexes: {gpr_types['mixed']:4d} (should be bidirectional MIXED)")
            print(f"    No GPR:          {gpr_types['none']:4d} (no constraints possible)")

            expected_bidirectional = gpr_types['pure_iso'] * 2 + gpr_types['mixed'] * 2  # *2 for forward+reverse
            actual_bidirectional = bidirectional_iso_added + bidirectional_mixed_added
            print(f"  Expected bidirectional constraints: {expected_bidirectional}")
            print(f"  Actual bidirectional constraints:   {actual_bidirectional}")
            if expected_bidirectional > actual_bidirectional:
                print(f"  ⚠️  Gap: {expected_bidirectional - actual_bidirectional} missing bidirectional constraints!")
            print()

        # Diagnostic: Check if we actually have any constraints
        if total_active == 0:
            print("WARNING: NO ENZYME CONSTRAINTS WERE CREATED!")
            print("This means the optimization is unconstrained by enzymes.")
            print("Checking constraint generation issues:")
            print(f"  - kcat_dict entries: {len(kcat_dict)}")
            print(f"  - gene_sequences_dict entries: {len(gene_sequences_dict) if gene_sequences_dict else 0}")
            print(f"  - dnf_clauses entries: {len(dnf_clauses)}")
            if len(kcat_dict) > 0:
                print(f"  - Sample kcat_dict key: {list(kcat_dict.keys())[0]}")
            if len(dnf_clauses) > 0:
                print(f"  - Sample dnf_clauses: {list(dnf_clauses.items())[0]}")
        else:
            print(f"✓ Successfully created {total_active} enzyme constraints")

    # 7) Total enzyme pool / ratio
    step_start = time.time()
    if enzyme_ratio:
        if gene_sequences_dict is None:
            gene_sequences_dict = {}
        mw = {}
        num_default_mw = 0
        for g in genes:
            seq = gene_sequences_dict.get(g, '')
            try:
                mw_val = molecular_weight(seq, seq_type='protein')
                if not mw_val:
                    mw_val = 1e5
            except Exception as e:
                if verbose:
                    print(f"[MW ERROR] Gene: {g} | Sequence: '{seq}' | Error: {e}")
                mw_val = 1e5
            if mw_val == 1e5:
                num_default_mw += 1
            mw[g] = mw_val
        if verbose:
            print(f"[DIAG] {num_default_mw} out of {len(genes)} genes are using the default molecular weight (likely due to invalid/missing sequences).")
        m.E_ratio = Var(domain=NonNegativeReals, bounds=(0, enzyme_upper_bound))  # noqa: F405
        m.total_enzyme = Constraint(  # noqa: F405
            expr=sum(m.E[g] * mw[g] for g in m.G) * 1e-3 <= m.E_ratio
        )
    else:
        m.E_total = Var(domain=NonNegativeReals, bounds=(0, enzyme_upper_bound))  # noqa: F405
        m.total_enzyme = Constraint(expr=sum(m.E[g] for g in m.G) <= m.E_total)  # noqa: F405
    _print_step_timing(f"Step 8: Create enzyme pool constraints (ratio={enzyme_ratio})", step_start, verbose)

    # 8) Solve
    step_start = time.time()
    solver = SolverFactory(solver_name)

    # Debug: Check if solver is available
    if not solver.available():
        if verbose:
            print(f"⚠️  WARNING: Solver '{solver_name}' not available in Pyomo, falling back to GLPK")
        solver = SolverFactory('glpk')

    # Set solver tolerances to avoid numerical precision warnings
    if solver_name.lower() == 'glpk':
        # Use minimal GLPK options to avoid conflicts
        solver.options['tmlim'] = 600  # time limit in seconds (increased)
    elif solver_name.lower() == 'gurobi':
        solver.options['FeasibilityTol'] = 1e-9
        solver.options['OptimalityTol'] = 1e-9
        solver.options['MemLimit'] = 4.0  # 4GB memory limit
    # elif solver_name.lower() == 'cplex':
        # CPLEX options for better performance
        # solver.options['timelimit'] = 300  # time limit in seconds

    #solver.options['threads'] = 4
    # print("Solver:", solver)
    # print(f"Solver options: {dict(solver.options)}")
    _print_step_timing("Step 9a: Setup solver", step_start, verbose)

    # Actually solve the optimization problem
    step_start = time.time()
    solver.solve(m, tee=False, load_solutions=True)  # Always suppress solver output
    _print_step_timing("Step 9b: Solve optimization problem", step_start, verbose)
    # Print total enzyme usage after optimization if possible
    if enzyme_ratio and verbose:
        try:
            # Count non-zero enzyme values
            non_zero_enzymes = [g for g in m.G if value(m.E[g]) > 1e-6]
            total_enzyme = value(sum(m.E[g] * mw[g] for g in m.G) * 1e-3)
            e_ratio_value = value(m.E_ratio) if hasattr(m, 'E_ratio') else None
            print(f"[DIAG] Number of genes with enzyme allocation: {len(non_zero_enzymes)}/{len(m.G)}")
            print(f"[DIAG] Total enzyme usage (g/gDW): {total_enzyme:.6g} (upper bound: {enzyme_upper_bound})")
            if e_ratio_value is not None:
                print(f"[DIAG] Final E_ratio value: {e_ratio_value:.6g}")
        except Exception as e:
            print(f"[DIAG] Could not compute total enzyme usage: {e}")

    # 9) Post-process - handle numerical precision issues
    step_start = time.time()
    # Clamp very small values to zero to avoid floating-point errors
    tolerance = 1e-10

    if bidirectional_constraints:
        # Handle forward and reverse flux variables
        for r in m.R_rev:
            # Forward flux
            val_fwd = m.v_fwd[r].value if m.v_fwd[r].value is not None else 0.0
            if abs(val_fwd) < tolerance:
                val_fwd = 0.0
            m.v_fwd[r].value = max(val_fwd, 0.0)

            # Reverse flux
            val_rev = m.v_rev[r].value if m.v_rev[r].value is not None else 0.0
            if abs(val_rev) < tolerance:
                val_rev = 0.0
            m.v_rev[r].value = max(val_rev, 0.0)

        # Handle irreversible flux variables
        for r in m.R_irr:
            val = m.v_irr[r].value if m.v_irr[r].value is not None else lb[r]
            if abs(val) < tolerance:
                val = 0.0
            m.v_irr[r].value = max(min(val, ub[r]), lb[r])
    else:
        # Standard flux variables
        for r in m.R:
            val = m.v[r].value if m.v[r].value is not None else lb[r]
            if abs(val) < tolerance:
                val = 0.0
            m.v[r].value = max(min(val, ub[r]), lb[r])

    # Handle enzyme variables
    for g in m.G:
        val = m.E[g].value if m.E[g].value is not None else 0.0
        if abs(val) < tolerance:
            val = 0.0
        m.E[g].value = max(val, 0.0)
    _print_step_timing("Step 10a: Post-process numerical precision", step_start, verbose)

    # 10) Collect results
    step_start = time.time()
    sol_val = value(m.obj)  # noqa: F405

    # Collect flux results more efficiently
    records = []

    # Process results efficiently
    if bidirectional_constraints:
        # For reversible reactions, compute net flux and store both forward/reverse components
        for r in m.R_rev:
            v_fwd = m.v_fwd[r].value if m.v_fwd[r].value is not None else 0.0
            v_rev = m.v_rev[r].value if m.v_rev[r].value is not None else 0.0
            net_flux = v_fwd - v_rev
            records.extend([
                ('flux', r, net_flux),
                ('flux_fwd', r, v_fwd),
                ('flux_rev', r, v_rev)
            ])

        # For irreversible reactions, use standard flux
        for r in m.R_irr:
            flux_val = m.v_irr[r].value if m.v_irr[r].value is not None else 0.0
            records.append(('flux', r, flux_val))
    else:
        # Standard flux collection
        records.extend([('flux', r, m.v[r].value if m.v[r].value is not None else 0.0) for r in m.R])

    # Collect enzyme results
    records.extend([('enzyme', g, m.E[g].value if m.E[g].value is not None else 0.0) for g in m.G])

    df_FBA = pd.DataFrame(records, columns=['Variable','Index','Value'])
    _print_step_timing(f"Step 10b: Collect results ({len(records)} variables)", step_start, verbose)

    # Final timing summary
    if verbose:
        total_time = time.time() - overall_start
        final_memory = _get_memory_usage()
        print(f"\n{'='*60}")
        print("OPTIMIZATION COMPLETED")
        print(f"{'='*60}")
        print(f"Total optimization time: {total_time:.2f}s")
        if final_memory > 0:
            print(f"Final memory usage: {final_memory:.1f} MB")
            if initial_memory > 0:
                memory_increase = final_memory - initial_memory
                print(f"Memory increase: {memory_increase:.1f} MB")
        print(f"Objective value: {sol_val:.6f}")
        print(f"Solution variables: {len(records)}")
        print(f"{'='*60}\n")

    # DIAGNOSTIC: Print exchange reaction fluxes if medium was provided
    # if medium is not None:
    #     print("\n=== DIAGNOSTIC: Exchange reaction fluxes after optimization ===")
    #     for rxn_id in medium.keys():
    #         flux_val = m.v[rxn_id].value if rxn_id in m.R else None
    #         if flux_val is not None:
    #             print(f"  {rxn_id}: flux = {flux_val:.4f} (bound was {medium[rxn_id]:.4f})")
    #         else:
    #             print(f"  {rxn_id}: NOT FOUND in optimization results")
    #     print(f"Objective value: {sol_val:.6f}")
    #     print("=== END DIAGNOSTIC ===\n")

    return sol_val, df_FBA, gene_sequences_dict, m

def create_descriptive_filename(objective_reaction, enzyme_upper_bound, maximization,
                         multi_enzyme_off, isoenzymes_off, promiscuous_off, complexes_off,
                         output_dir=None, extension='.csv'):
    """
    Create a descriptive filename based on optimization parameters.

    Parameters
    ----------
    objective_reaction : str
        The reaction ID used as objective
    enzyme_upper_bound : float
        Upper bound for total enzyme concentration
    maximization : bool
        Whether maximization was used
    multi_enzyme_off : bool
        Whether multi-enzyme reactions were disabled
    isoenzymes_off : bool
        Whether isoenzyme handling was disabled
    promiscuous_off : bool
        Whether promiscuous enzyme handling was disabled
    complexes_off : bool
        Whether enzyme complex handling was disabled
    output_dir : str, optional
        Directory to save the file in
    extension : str, optional
        File extension to use

    Returns
    -------
    str
        The complete filepath
    """
    import os

    # Shorten the objective reaction name if it's too long
    if len(objective_reaction) > 20:
        obj_short = objective_reaction[:20]
    else:
        obj_short = objective_reaction

    # Create descriptive part for optimization direction
    opt_dir = "max" if maximization else "min"

    # Create descriptive parts for enzyme constraints
    multi = "noMulti" if multi_enzyme_off else "Multi"
    iso = "noIso" if isoenzymes_off else "Iso"
    promis = "noPromis" if promiscuous_off else "Promis"
    complex_str = "noComplex" if complexes_off else "Complex"

    # Create the filename
    filename = f"FBA_{obj_short}_{opt_dir}_E{enzyme_upper_bound}_{multi}_{iso}_{promis}_{complex_str}{extension}"

    # Clean up any characters that might cause issues in filenames
    filename = filename.replace(' ', '_').replace('/', '_').replace('\\', '_')

    # Join with output directory if provided
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
        return os.path.join(output_dir, filename)
    else:
        return filename


def run_optimization_with_dataframe(model, processed_df, objective_reaction,
                    enzyme_upper_bound=0.125, enzyme_ratio=True, maximization=True,
                    multi_enzyme_off=False, isoenzymes_off=False,
                    promiscuous_off=False, complexes_off=False,
                    output_dir=None, save_results=True, print_reaction_conditions=True, verbose=True,
                    solver_name='glpk', medium=None, medium_upper_bound=False,
                    bidirectional_constraints=True):
    """
    Run enzyme-constrained flux balance analysis using a processed dataframe.

    Parameters
    ----------
    model : cobra.Model or str
        COBRA model object or path to model file
    processed_df : pandas.DataFrame
        DataFrame containing Reactions, Single_gene, SEQ, SMILES, and kcat_mean columns.
        For bidirectional constraints, must also contain Direction column.
    objective_reaction : str
        Reaction ID to maximize/minimize
    enzyme_upper_bound : float, optional
        Upper bound for total enzyme concentration
    enzyme_ratio : bool, optional
        Whether to use enzyme ratio constraint
    maximization : bool, optional
        Whether to maximize (True) or minimize (False) the objective
    multi_enzyme_off : bool, optional
        Whether to disable multi-enzyme reactions
    isoenzymes_off : bool, optional
        Whether to disable isoenzyme handling
    promiscuous_off : bool, optional
        Whether to disable promiscuous enzyme handling
    complexes_off : bool, optional
        Whether to disable enzyme complex handling
    output_dir : str, optional
        Directory to save results in
    save_results : bool, optional
        Whether to automatically save results to a file
    medium : dict, optional
        Dictionary mapping exchange reaction IDs to their flux values.
        Example: {"EX_glc__D_e": -10, "EX_o2_e": -14.49}
    medium_upper_bound : bool, optional
        If True, set both lower and upper bounds equal.
        If False, only set lower bound.
    bidirectional_constraints : bool, optional
        Whether to use substrate-specific bidirectional constraints for reversible reactions.
        Requires 'Direction' column in processed_df.

    Returns
    -------
    tuple
        (solution_value, df_FBA, gene_sequences_dict, output_filepath)
    """


    # Check if required columns exist in processed_df
    required_cols = ['Reactions', 'Single_gene', 'SEQ', 'kcat_mean']
    missing_cols = [col for col in required_cols if col not in processed_df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns in processed_df: {missing_cols}")

    # Check if bidirectional constraints are requested but Direction column is missing
    if bidirectional_constraints and 'Direction' not in processed_df.columns:
        print("WARNING: bidirectional_constraints=True but 'Direction' column not found. Falling back to standard constraints.")
        bidirectional_constraints = False

    # Extract kcat dictionary and gene sequences from processed_df
    kcat_dict = {}
    gene_sequences_dict = {}

    # Determine which kcat column to use (prefer 'kcat' if available, fallback to 'kcat_mean')
    kcat_col = 'kcat' if 'kcat' in processed_df.columns else 'kcat_mean'

    # OPTIMIZED: Filter valid rows first, then use vectorized operations
    # This is MUCH faster than iterrows() for large DataFrames (398K rows)
    valid_rows = processed_df[processed_df[kcat_col].notna() & processed_df['SEQ'].notna()].copy()

    if bidirectional_constraints:
        # Build direction-aware kcat_dict using (Reactions, Single_gene, Direction) as keys
        grouped = valid_rows.groupby(['Reactions', 'Single_gene', 'Direction'])[kcat_col].max()
        kcat_dict = {key: [value] for key, value in grouped.items()}

        if verbose:
            forward_count = sum(1 for key in kcat_dict if key[2] == 'forward')
            reverse_count = sum(1 for key in kcat_dict if key[2] == 'reverse')
            print(f"Created bidirectional kcat dictionary with {len(kcat_dict)} entries:")
            print(f"  Forward direction: {forward_count}")
            print(f"  Reverse direction: {reverse_count}")
            if len(kcat_dict) > 0:
                print(f"  Sample kcat keys: {list(kcat_dict.keys())[:3]}")
    else:
        # Build standard kcat_dict using (Reactions, Single_gene) as keys
        # Group by (Reactions, Single_gene) and take the max kcat if duplicates exist
        grouped = valid_rows.groupby(['Reactions', 'Single_gene'])[kcat_col].max()
        kcat_dict = {key: [value] for key, value in grouped.items()}

        if verbose:
            print(f"Created standard kcat dictionary with {len(kcat_dict)} entries")
            if len(kcat_dict) > 0:
                print(f"  Sample kcat keys: {list(kcat_dict.keys())[:3]}")

    # Build gene_sequences_dict (take first sequence for each gene)
    gene_seq_grouped = valid_rows.groupby('Single_gene')['SEQ'].first()
    gene_sequences_dict = gene_seq_grouped.to_dict()

    if verbose:
        if len(kcat_dict) > 0:
            print(f"Created kcat dictionary with {len(kcat_dict)} entries")
            if not bidirectional_constraints and len(kcat_dict) > 0:
                print(f"  Sample kcat keys: {list(kcat_dict.keys())[:3]}")
        else:
            print("WARNING: No valid kcat data found")
        print(f"Created gene_sequences_dict with {len(gene_sequences_dict)} entries")

    solution_value, df_FBA, gene_sequences_dict, m = run_optimization(
        model=model,
        kcat_dict=kcat_dict,
        objective_reaction=objective_reaction,
        gene_sequences_dict=gene_sequences_dict,
        enzyme_upper_bound=enzyme_upper_bound,
        maximization=maximization,
        multi_enzyme_off=multi_enzyme_off,
        isoenzymes_off=isoenzymes_off,
        promiscuous_off=promiscuous_off,
        complexes_off=complexes_off,
        enzyme_ratio=True,
        tee=verbose,
        solver_name=solver_name,
        medium=medium,
        medium_upper_bound=medium_upper_bound,
        bidirectional_constraints=bidirectional_constraints,
        verbose=verbose  # Make sure verbose is passed through
    )

    # Create descriptive filename and save results if requested
    output_filepath = create_descriptive_filename(
        objective_reaction,
        enzyme_upper_bound,
        maximization,
        multi_enzyme_off,
        isoenzymes_off,
        promiscuous_off,
        complexes_off,
        output_dir=output_dir
    )

    if save_results and solution_value is not None:
        # Create directory if it doesn't exist
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)

        # Save the results
        df_FBA.to_csv(output_filepath, index=False)
        print(f"Results saved to: {output_filepath}")

    return solution_value, df_FBA, gene_sequences_dict, output_filepath

# Create a version of run_optimization with fixed unit conversions



def debug_enzyme_constraints_detailed(model, processed_data, objective_reaction):
    """
    Debug the specific enzyme constraint implementation
    """
    print("=== ENZYME CONSTRAINT IMPLEMENTATION DEBUG ===\n")

    # Create kcat dictionary and gene sequences
    kcat_dict = {}
    gene_sequences_dict = {}

    # Process data - exactly as in run_optimization_with_dataframe
    print("1. Creating kcat dictionary and gene sequences...")
    for _, row in processed_data.iterrows():
        if pd.notna(row['kcat_mean']) and pd.notna(row['SEQ']):
            reaction_id = row['Reactions']
            gene_id = row['Single_gene']

            # Store the kcat value as a list
            kcat_dict[(reaction_id, gene_id)] = [row['kcat_mean']]

            # Store gene sequence
            if gene_id not in gene_sequences_dict and pd.notna(row['SEQ']):
                gene_sequences_dict[gene_id] = row['SEQ']

    print(f"Created kcat_dict with {len(kcat_dict)} entries")
    print(f"Created gene_sequences_dict with {len(gene_sequences_dict)} entries")

    # Sample some entries
    print("\nSample kcat_dict entries:")
    for i, (key, kcat_value) in enumerate(list(kcat_dict.items())[:3]):
        print(f"  {key}: {kcat_value}")

    # Check for problematic values
    print("\n2. Checking for problematic values...")

    problematic_kcats = []
    problematic_sequences = []

    for key, kcat_list in kcat_dict.items():
        if len(kcat_list) == 0 or kcat_list[0] <= 0 or math.isnan(kcat_list[0]) or math.isinf(kcat_list[0]):
            problematic_kcats.append((key, kcat_list))

    for gene, seq in gene_sequences_dict.items():
        if not seq or seq == '' or seq == 'None':
            problematic_sequences.append((gene, seq))

    print(f"Problematic kcat values: {len(problematic_kcats)}")
    print(f"Problematic sequences: {len(problematic_sequences)}")

    # Show examples
    if problematic_kcats:
        print("Sample problematic kcats:")
        for item in problematic_kcats[:3]:
            print(f"  {item}")

    if problematic_sequences:
        print("Sample problematic sequences:")
        for item in problematic_sequences[:3]:
            print(f"  {item}")

    # Test molecular weight calculation
    print("\n3. Testing molecular weight calculation...")

    test_sequences = list(gene_sequences_dict.values())[:5]
    for seq in test_sequences:
        try:
            mw = molecular_weight(seq, seq_type='protein')
            print(f"Sequence length {len(seq)}: MW = {mw}")
        except Exception as e:
            print(f"Error calculating MW for sequence: {e}")

    # Check enzyme concentration calculation
    print("\n4. Checking enzyme concentration calculation...")

    # Use just one reaction-gene pair for testing
    test_key = list(kcat_dict.keys())[0]
    test_kcat = kcat_dict[test_key][0]
    test_gene = test_key[1]
    test_seq = gene_sequences_dict.get(test_gene, '')

    print(f"Test case: {test_key}")
    print(f"kcat: {test_kcat}")
    print(f"Gene: {test_gene}")
    print(f"Sequence length: {len(test_seq)}")

    # Calculate enzyme requirement for a typical flux
    typical_flux = 1.0  # mmol/gDCW/hr
    enzyme_mmol = typical_flux / test_kcat

    if test_seq:
        mw = molecular_weight(test_seq, seq_type='protein')
        enzyme_g = enzyme_mmol * mw / 1000 / 1000  # Convert to g/gDCW

        print(f"For flux {typical_flux}, requires {enzyme_mmol} mmol enzyme")
        print(f"MW: {mw}")
        print(f"Enzyme requirement: {enzyme_g} g/gDCW")

        if enzyme_g > 1.0:
            print("WARNING: This reaction requires more than cell mass!")

    # Test the enzyme constraint creation
    print("\n5. Testing enzyme constraint creation...")
    # Create minimal Pyomo model to test constraints
    from pyomo.environ import (
        ConcreteModel,
        Constraint,
        NonNegativeReals,
        Objective,
        Suffix,  # noqa: F401
        Var,
        maximize,
        value,
    )
    from pyomo.opt import SolverFactory

    test_model = ConcreteModel()
    # Add variables
    test_model.reaction = Var(within=NonNegativeReals, bounds=(0, 10))
    test_model.enzyme = Var(within=NonNegativeReals, bounds=(0, 1))

    # Add enzyme constraint
    def enzyme_constraint(m):
        return m.reaction <= test_kcat * m.enzyme

    test_model.enzyme_con = Constraint(rule=enzyme_constraint)

    # Add enzyme bound
    def enzyme_bound(m):
        return m.enzyme <= 0.1  # 10% of cell mass

    test_model.enzyme_bound_con = Constraint(rule=enzyme_bound)

    # Add objective
    test_model.obj = Objective(expr=test_model.reaction, sense=maximize)

    # Solve
    solver = SolverFactory('glpk')
    # solver.options['max_iter'] = 100

    try:
        results = solver.solve(test_model, tee=False)
        print(f"Test model status: {results.solver.status}")
        print(f"Test reaction flux: {value(test_model.reaction)}")
        print(f"Test enzyme amount: {value(test_model.enzyme)}")

        # Calculate theoretical maximum
        theoretical_max = test_kcat * 0.1
        print(f"Theoretical maximum flux: {theoretical_max}")

    except Exception as e:
        print(f"Error solving test model: {e}")

    # Test with enzyme ratio constraint
    print("\n6. Testing enzyme ratio constraint...")

    test_model_ratio = ConcreteModel()

    # Add variables
    test_model_ratio.reaction = Var(within=NonNegativeReals, bounds=(0, 10))
    test_model_ratio.enzyme = Var(within=NonNegativeReals)
    test_model_ratio.E_ratio = Var(within=NonNegativeReals, bounds=(0, 0.125))

    # Add enzyme constraint
    def enzyme_constraint_ratio(m):
        return m.reaction <= test_kcat * m.enzyme

    test_model_ratio.enzyme_con = Constraint(rule=enzyme_constraint_ratio)

    # Add enzyme ratio constraint
    if test_seq:
        mw = molecular_weight(test_seq, seq_type='protein')
    else:
        mw = 50000  # Default MW

    def enzyme_ratio_constraint(m):
        return m.enzyme * mw * 0.001 <= m.E_ratio

    test_model_ratio.ratio_con = Constraint(rule=enzyme_ratio_constraint)

    # Add objective
    test_model_ratio.obj = Objective(expr=test_model_ratio.reaction, sense=maximize)

    # Solve
    try:
        results_ratio = solver.solve(test_model_ratio, tee=False)
        print(f"Test ratio model status: {results_ratio.solver.status}")
        print(f"Test reaction flux: {value(test_model_ratio.reaction)}")
        print(f"Test enzyme amount: {value(test_model_ratio.enzyme)}")
        print(f"Test E_ratio: {value(test_model_ratio.E_ratio)}")

        # Calculate theoretical limits
        max_enzyme_mmol = 0.125 / (mw / 1000 / 1000)
        max_flux = test_kcat * max_enzyme_mmol
        print(f"Theoretical max enzyme: {max_enzyme_mmol} mmol")
        print(f"Theoretical max flux: {max_flux}")

    except Exception as e:
        print(f"Error solving test ratio model: {e}")

    # Check for unit conversion issues
    print("\n7. Checking unit conversions...")

    # Example calculation
    flux_mmol_per_hr = 1.0  # mmol/gDCW/hr
    kcat_per_s = test_kcat  # 1/s

    # Convert flux to 1/s
    flux_per_s = flux_mmol_per_hr / 3600  # Convert hr to s

    # Calculate required enzyme
    enzyme_required_mmol = flux_per_s / kcat_per_s

    print(f"Flux: {flux_mmol_per_hr} mmol/gDCW/hr = {flux_per_s} mmol/gDCW/s")
    print(f"kcat: {kcat_per_s} 1/s")
    print(f"Required enzyme: {enzyme_required_mmol} mmol/gDCW")

    if test_seq:
        mw = molecular_weight(test_seq, seq_type='protein')
        enzyme_required_g = enzyme_required_mmol * mw / 1000 / 1000
        print(f"Required enzyme: {enzyme_required_g} g/gDCW")

    return kcat_dict, gene_sequences_dict

def validate_enzyme_constraints(df_FBA,
                                kcat_dict_hr,
                                gene_sequences_dict,
                                reaction_gene_list,
                                gpr_dict,
                                enzyme_ratio,            # True or False
                                enzyme_upper_bound,      # the same upper bound you passed
                                enzyme_mw_dict,          # {gene: molecular_weight_in_dalton}
                                S_mat):
    """
    Given the DataFrame df_FBA (with columns ['Variable','Index','Value'])
    and the same dictionaries used inside run_optimization, check:
      1) For every (reaction, gene) that should be constrained, v_j ≤ kcat_j * e_i.
      2) Total‐enzyme constraint: sum(e_i * MW_i)*0.001 ≤ enzyme_upper_bound (if enzyme_ratio=True),
         or sum(e_i) ≤ enzyme_upper_bound (if enzyme_ratio=False).
      3) Steady‐state (S · v = 0) up to a small numerical tolerance.
    """

    # 1) Build solution dicts for v_j and e_i
    #    df_FBA rows look like: Variable='reaction', Index='R_EX_glc__D_e', Value=4.5
    flux_sol = {}   # reaction_id → flux_value
    enz_sol  = {}   # gene_id     → enzyme_amount

    for _, row in df_FBA.iterrows():
        varname = row['Variable']
        idx     = row['Index']
        val     = float(row['Value'])
        if varname == 'reaction':
            flux_sol[idx] = val
        elif varname == 'enzyme':
            enz_sol[idx] = val
    # 1a) Check that we really extracted them
    if not flux_sol:
        raise RuntimeError("Could not find any 'reaction' entries in df_FBA.")
    if not enz_sol:
        raise RuntimeError("Could not find any 'enzyme' entries in df_FBA.")

    # 2) Check each (reaction, gene) in reaction_gene_list
    violations = []
    for (j, i) in reaction_gene_list:
        # Was (j,i) supposed to be constrained?  It will be constrained if:
        #   • j in kcat_dict_hr (or in reaction.annotation['kcat']) AND
        #   • i in gene_sequences_dict
        #   • AND inside your rule_kcat pipeline, it did NOT fall into the “missing data → Feasible” case.
        #
        # For simplicity, assume: if j in kcat_dict_hr and i in gene_sequences_dict, then it should be constrained.
        if j not in kcat_dict_hr or i not in gene_sequences_dict or not gene_sequences_dict[i]:
            continue

        vj = flux_sol.get(j, 0.0)
        # there might be multiple possible kcats for this reaction—your code picks the first (or does a max).
        # Here we assume kcat_dict_hr[j] is a float or a one-element list.
        kc = kcat_dict_hr[j]
        if isinstance(kc, list):
            kc_val = kc[0]
        else:
            kc_val = kc

        ei = enz_sol.get(i, 0.0)
        if vj > kc_val * ei + 1e-6:    # small tolerance for numerical solver noise
            violations.append((j, i, vj, kc_val*ei))

    if violations:
        print("FOUND violations of (v_j ≤ kcat·e_i):")
        for j, i, vj, bound in violations[:5]:
            print(f"  • Reaction {j}, gene {i}: flux {vj:.4g} > kcat·e = {bound:.4g}")
        print(f"...plus {len(violations)-5} more." if len(violations)>5 else "")
    else:
        print("All (reaction, gene) v_j ≤ kcat·e_i constraints are satisfied (within 1e-6 tolerance).")

    # 3) Check total‐enzyme constraint
    if enzyme_ratio:
        # Reconstruct: total_weight_grams_per_gDCW = ∑ (e_i [mmol/gDCW] * MW_i [g/mol]) * 1e−3 (to get g/gDCW).
        total_weight = 0.0
        for gene, ei in enz_sol.items():
            mw = enzyme_mw_dict.get(gene, 1e5)  # fallback if missing
            total_weight += ei * mw * 1e-3
        if total_weight > enzyme_upper_bound + 1e-6:
            print(f"Total‐enzyme weight constraint violated: {total_weight:.6f} > {enzyme_upper_bound:.6f}")
        else:
            print(f"Total‐enzyme weight = {total_weight:.6f} ≤ {enzyme_upper_bound:.6f} (OK)")

    else:
        total_e = sum(enz_sol.values())
        if total_e > enzyme_upper_bound + 1e-6:
            print(f"Total‐enzyme mmol constraint violated: {total_e:.6f} > {enzyme_upper_bound:.6f}")
        else:
            print(f"Sum(e_i) = {total_e:.6f} ≤ {enzyme_upper_bound:.6f} (OK)")

    # 4) Check steady‐state: S·v ≈ 0 for every metabolite t
    ss_violations = []
    for (met_id, rxn_id), coeff in S_mat.items():
        # S_mat keys are tuples (met_id, rxn_id) → stoich.  We want ∑_j S[t,j]·v_j = 0
        # so accumulate into a per‐met residual.
        pass
    # Instead, build a metabolite → sum(S[t,j]*v_j)
    met_residual = {t: 0.0 for (t, _) in S_mat.keys()}
    for (t, j), coeff in S_mat.items():
        vj = flux_sol.get(j, 0.0)
        met_residual[t] += coeff * vj

    # Now check if any |residual| > tol
    for t, resid in met_residual.items():
        if abs(resid) > 1e-6:
            ss_violations.append((t, resid))

    if ss_violations:
        print("Steady‐state (S·v=0) violations (|residual| > 1e−6):")
        for t, resid in ss_violations[:5]:
            print(f"  • Metabolite {t}: ∑ S[{t},j]·v_j = {resid:.4g}")
        print(f"...plus {len(ss_violations)-5} more." if len(ss_violations)>5 else "")
    else:
        print("All metabolites satisfy steady‐state (S·v≈0).")