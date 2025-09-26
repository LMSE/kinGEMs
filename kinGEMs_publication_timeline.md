# kinGEMs Publication Timeline (Updated – Aug 2025)

Target journal: **Nature Biotechnology**  
Target submission: **By Sep 30, 2025**

---

## ✅ Major Updates
- Noor's Biolog phenotype data → **classification only** (growth vs. no growth).
- Add **flux validation method (Lya’s analysis)**.
- Debugging & validation will run **in parallel with writing**.
- Manuscript writing starts **next week (Aug 25, 2025)**.

---

## 📅 Timeline & To-Do List

### Aug 19–24 (This Week)
- [X] Debug simulation issues (using E. coli runs).
- [ ] Coordinate with Lya on flux validation outputs.
- [X] Finalize plan for Biolog classification metrics.

### Aug 25 – Sep 7 (Writing + Validation in Parallel)
- Writing:
  - [ ] Draft **Introduction** + **Methods** sections.
  - [X] Start figure assembly (schematics, pipeline, Biolog classification workflow).
- Validation:
  - [ ] Benchmark kinGEMs growth predictions vs Biolog (classification only).
  - [ ] ~~Benchmark kinGEMs vs black-box model (scDCA from Genentech)~~.
  - .
  - [ ] Begin scaling tests on 10–20 species AGORA subset.

### Sep 8 – 14
- Writing:
  - [ ] Draft **Results (Part 1)**: P. putida & S. elongatus benchmarks, Biolog classification results.
  - [ ] Add figures: growth classification ROC/PR curves, confusion matrices.
- Validation:
  - [ ] Benchmark kinGEMs growth predictions vs Biolog classification
  - [ ] Rerun E. coli genetic perturbation prediction on new pipeline
  - [ ] Collate FVA precision + knockout accuracy metrics.
  - [ ] Continue AGORA scaling.

### Sep 15 – 21
- Writing:
  - [ ] Draft **Results (Part 2)**: Flux validation (Lya), scaling to multi-species, comparisons with black-box.
  - [ ] Draft **Discussion** section.
- Validation:
  - [ ] Run P. putida and S. elongatus on kinGEMs
  - [ ] Run P. putida and S. elongatus genetic perturbation
  - [ ] Finalize all stats/plots for flux validation & scaling.
  - [ ] QA check all figures and tables.
  - [ ] Integrate Lya’s flux validation results

### Sep 22 – 26
- [ ] Internal review with co-authors.
- [ ] Revise manuscript text, polish writing, captions, figure consistency.
- [ ] Prepare Supplementary Info (pipelines, extra tables).

### Sep 27 – 30
- [ ] Final checks: references, cover letter, journal formatting.
- [ ] Package and **submit to PNAS** by Sep 30, 2025.


### Sep 26 meeting
- **New approach**: utilize FBA data as input
- If gene essentiality doesn't improve, discard
- Sim annealing + media changes -Lya
- Revising gene essentiality (validation nb) -Lya
- Use Noor's Biolog for growth (binary) -Rana

--> New deadline: All validation resutls should be compeltefd by **Oct 13th**



---

## 🔑 Parallel Work Principles
- **Debug/validate daily** (AM sessions).
- **Write daily** (PM sessions) starting Aug 25.
- Keep figures and results modular → drop-in replacements possible as debugging/validation improves.

