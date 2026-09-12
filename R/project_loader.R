source_project_files <- function(project_root) {
  files <- c(
    "R/packages.R",
    "R/config.R",
    "R/country_selection.R",
    "R/data.R",
    "R/model_utils.R",
    "R/changepoint.R",
    "R/smooth_cohort.R"
  )
  missing <- files[!file.exists(file.path(project_root, files))]
  if (length(missing)) stop("Project incomplete. Missing: ", paste(missing, collapse = ", "), call. = FALSE)
  for (f in files) source(file.path(project_root, f), local = FALSE)
  invisible(files)
}

assert_project_api <- function() {
  required <- c(
    "country_screening_manifest", "prepare_mortality_data", "truncate_supported_cohort_window",
    "build_changepoint_data", "fit_changepoint_model", "summarise_changepoint_fit",
    "build_smooth_data", "fit_smooth_cohort_model", "summarise_smooth_diagnostics",
    "summarise_smooth_masked"
  )
  missing <- required[!vapply(required, exists, logical(1), mode = "function", inherits = TRUE)]
  if (length(missing)) stop("Project API incomplete after sourcing: ", paste(missing, collapse = ", "), call. = FALSE)
  invisible(TRUE)
}
