# Exact analysis run order

Open `JRSSA_NB_HMM_CODE_GITHUB.Rproj` and work from the project root.

## One-time setup

```r
source("install_dependencies.R")
```

## Reproduce the statistical analyses

Run:

```r
source("RUN_ALL_ANALYSIS.R")
```

or execute the stages manually in this exact order:

```r
source("00_preflight.R")
source("01_run_v7_primary.R")
source("02_run_v7_support_shift.R")
source("03_run_v7_masked_target.R")
source("04_run_smooth_primary.R")
source("05_run_smooth_masked_target.R")
source("06_compare_masked_models.R")
```

The workflow requires R, CmdStan and a working C++ toolchain. The Stan fits are computationally expensive.
