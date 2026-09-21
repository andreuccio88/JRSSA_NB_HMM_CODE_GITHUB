# JRSSA_NB_HMM_CODE_GITHUB

R and Stan code for the statistical analyses in the paper on generational mortality change, negative-binomial mortality modelling and latent cohort transitions.

This repository is an analysis-only release derived directly from the project used to obtain the reported empirical results. Manuscript sources, publication-formatting code and pre-computed result files are not included.

## Contents

- `data/` — HMD death-count and exposure input objects used in the analysis.
- `config/` — data-selection and model configuration used for the reported fits.
- `stan/` — the exact Stan programs used for the hierarchical change-point model and smooth-cohort comparator.
- `R/` — R functions used for data preparation, estimation and diagnostics.
- `00_preflight.R` to `06_compare_masked_models.R` — the analysis stages used to obtain the reported results.
- `RUN_ALL_ANALYSIS.R` — convenience runner executing those stages in the original analysis order.

## Reproduction

From the project root:

```r
source("install_dependencies.R")
source("RUN_ALL_ANALYSIS.R")
```

The stages are computationally intensive. See `RUN_ORDER.md` for the exact manual run order.

## Analysis design

The primary model uses negative-binomial death counts with exposure offsets, smooth country-sex age profiles, country-sex net drift, nonlinear period curvature, and a hierarchical single cohort transition marginalized exactly over candidate transition cohorts. The comparator replaces the discrete cohort representation with a regularized smooth-cohort component. The two specifications are evaluated with support-shift diagnostics and a targeted 1925--1935 masked reconstruction exercise.

## Data

The included RData objects are the HMD-derived inputs used for the reported analysis. Users should consult the Human Mortality Database for current data access, citation and redistribution terms: https://www.mortality.org/.
