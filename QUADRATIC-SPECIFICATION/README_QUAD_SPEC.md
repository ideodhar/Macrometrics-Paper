# MATLAB to R translation

This folder translates the uploaded MATLAB replication files for the fiscal/Ramey levels branch used by `FINAL_OUTPUT.m`.

## Files

- `FINAL_OUTPUT.R` — main runner corresponding to `FINAL_OUTPUT.m`.
- `CROSS_VALIDATION_BAYESIAN_RAMEY_LEVELS.R` — translation of the levels cross-validation / Bayesian IRF routine.
- `generate_unobserved_shock_ramey_levels.R` — translation of `GENERATE_UNOBSERVED_SHOCK_RAMEY_LEVELS.m` for the fiscal branch.
- `nwest_balanced_panel_bayesian.R` — translation of the Bayesian Newey-West regression helper.
- `polorth.R` — translation of the polynomial basis helper.
- `bayes_utils.R` — utility functions for inverse-Wishart and multivariate-normal draws.

## Run

Install dependencies once:

```r
install.packages(c("R.matlab", "MASS"))
```

Then run:

```r
setwd("/path/to/matlab_to_R_translation")
source("FINAL_OUTPUT.R")
run_final_output(data_file = "../FISCAL_DATA_SHOCK.mat", out_dir = "R_OUTPUT", rep = 2000, seed = 123)
```

Outputs are written to `R_OUTPUT/`:

- `FISCAL_NEWS_FIGURE_NONLINEAR_LEVELS11.pdf`
- `FISCAL_NEWS_FIGURE_LINEAR_LEVELS11.pdf`
- `FIS_MC_INPUT_LEVELS_1.RData`

## Scope note

The uploaded `FINAL_OUTPUT.m` calls the levels fiscal branch (`state=0`, `levels=1`, `aic=1`, `lag_ramey=1`). This translation implements that executable path. The unused branches in the MATLAB file that call non-uploaded functions such as `CROSS_VALIDATION_BAYESIAN_RAMEY`, `CROSS_VALIDATION_BAYESIAN_RAMEY_STATE`, and other monetary/oil/uncertainty data files are not implemented here.
