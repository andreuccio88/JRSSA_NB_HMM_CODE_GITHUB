# Complete analysis runner.
# Each stage is executed in a fresh R session and the workflow stops on failure.

scripts <- c(
  "scripts/01_preflight.R",
  "scripts/02_fit_changepoint.R",
  "scripts/03_support_sensitivity.R",
  "scripts/04_masked_changepoint.R",
  "scripts/05_fit_smooth.R",
  "scripts/06_masked_smooth.R",
  "scripts/07_compare_models.R",
  "scripts/99_session_info.R"
)

rscript <- file.path(
  R.home("bin"),
  if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript"
)
if (!file.exists(rscript)) stop("Rscript executable not found: ", rscript)

for (s in scripts) {
  cat("\n============================================================\n")
  cat("RUNNING: ", s, "\n", sep = "")
  cat("============================================================\n")
  status <- system2(rscript, s)
  if (!identical(status, 0L)) {
    stop("Stage failed: ", s, " (exit status ", status, ")")
  }
}

cat("\nANALYSIS COMPLETE.\n")
cat("Generated files are in: output/\n")
