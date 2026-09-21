# Reproduce the statistical analyses used in the paper.
# Run from the project root after install_dependencies.R.
# Each stage runs in a fresh R session and stops immediately on failure.

scripts <- c(
  "00_preflight.R",
  "01_run_v7_primary.R",
  "02_run_v7_support_shift.R",
  "03_run_v7_masked_target.R",
  "04_run_smooth_primary.R",
  "05_run_smooth_masked_target.R",
  "06_compare_masked_models.R"
)

rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
if (!file.exists(rscript)) stop("Rscript executable not found: ", rscript)

for (s in scripts) {
  cat("\n============================================================\n")
  cat("RUNNING: ", s, "\n", sep = "")
  cat("============================================================\n")
  status <- system2(rscript, s)
  if (!identical(status, 0L)) stop("Stage failed: ", s, " (exit status ", status, ")")
}

cat("\nSTATISTICAL ANALYSIS COMPLETE.\n")
