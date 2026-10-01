# BASiCS workflow

Run stages in order:

1. `01_prepare.R` freezes tasks, settings, checksum, and session information.
2. `02_fit_cell_type.R` fits one task locally or under SLURM and writes the chain
   before post-processing.
3. `03_diagnose.R` calculates epsilon ESS and Geweke statistics.
4. `04_detect_lvg.R` reuses the saved chain and diagnostics.
5. `05_assess_experiment_support.R` measures cross-experiment detection support.
6. `06_build_master.R` joins HVG/LVG calls and diagnostics and removes the
   immutable 199-gene stress list from analysis eligibility.
7. `scripts/analysis` creates mappings and enrichment results.

Each cell type is stored beneath `RUN/cell_types/NNN_NAME/`; `_COMPLETE` and
`_LVG_COMPLETE` are written last. Production settings are 20,000 iterations,
10,000 burn-in, thinning by 10, HVG upper percentile 0.90, LVG lower percentile
0.10, EFDR 0.10, and minimum epsilon ESS 100.

For SLURM, count rows in `cell_types.tsv`, submit `run_basics_array.sh` with
`--array=1-N`, run diagnostics once, then submit `run_lvg_array.sh`.
