`%||%` <- function(x, y) if (is.null(x)) y else x

is_absolute_path <- function(x) {
  x <- path.expand(x)
  grepl("^/", x) || grepl("^[A-Za-z]:[/\\\\]", x)
}

resolve_project_path <- function(root, x) {
  if (is_absolute_path(x)) normalizePath(x, winslash = "/", mustWork = FALSE)
  else normalizePath(file.path(root, x), winslash = "/", mustWork = FALSE)
}

resolve_stan_threads <- function(cfg) {
  requested <- as.integer(cfg$sampling$threads_per_chain %||% 0L)
  if (requested > 0L) return(requested)
  physical <- suppressWarnings(parallel::detectCores(logical = FALSE))
  if (!is.finite(physical) || physical < 1L) physical <- suppressWarnings(parallel::detectCores(logical = TRUE))
  if (!is.finite(physical) || physical < 1L) physical <- 4L
  parallel_chains <- max(1L, as.integer(cfg$sampling$parallel_chains %||% cfg$sampling$chains))
  cap <- max(1L, as.integer(cfg$runtime$stan_threads_max_per_chain %||% 4L))
  max(1L, min(cap, floor(physical / parallel_chains)))
}

read_project_config <- function(config_path, project_root) {
  cfg <- yaml::read_yaml(config_path)
  cfg$project_root <- normalizePath(project_root, winslash = "/", mustWork = TRUE)
  cfg$sampling$threads_per_chain <- resolve_stan_threads(cfg)

  for (nm in c("dx_file", "ex_file", "series_metadata", "oecd_universe")) {
    cfg$data[[nm]] <- resolve_project_path(cfg$project_root, cfg$data[[nm]])
  }
  cfg$paths$analysis <- resolve_project_path(cfg$project_root, cfg$paths$analysis)
  cfg$output_dir <- cfg$paths$analysis
  cfg$stan_changepoint_file <- file.path(cfg$project_root, "stan", "changepoint_model.stan")
  cfg$stan_smooth_file <- file.path(cfg$project_root, "stan", "smooth_cohort_model.stan")
  validate_project_config(cfg)
  cfg
}

validate_project_config <- function(cfg) {
  if (as.integer(cfg$data$expected_selected_populations) != 21L)
    stop("The application is pre-specified for 21 populations.", call. = FALSE)
  if (as.integer(cfg$data$cohort_step) != 5L)
    stop("The analysis is configured for five-year cohort groups.", call. = FALSE)

  positive_scales <- c(
    cfg$model$baseline_country_level_prior_sd,
    cfg$model$baseline_country_shape_prior_sd,
    cfg$model$baseline_global_rw2_prior_sd,
    cfg$model$drift_global_prior_sd,
    cfg$model$drift_country_deviation_prior_sd,
    cfg$model$drift_country_sd_prior,
    cfg$model$period_global_prior_sd,
    cfg$model$period_country_deviation_prior_sd,
    cfg$model$regime_range_prior_sd,
    cfg$model$regime_global_shape_prior_sd,
    cfg$model$country_regime_level_prior_sd,
    cfg$model$cohort_global_prior_sd,
    cfg$model$cohort_country_deviation_prior_sd,
    cfg$model$cohort_smoothness_sigma,
    cfg$model$switch_center_prior_sd_intervals,
    cfg$model$switch_sd_prior_intervals,
    cfg$model$switch_sd_floor_intervals,
    cfg$model$numerical_floor_phi,
    cfg$model$numerical_floor_baseline_rw2_sd
  )
  if (any(!is.finite(positive_scales)) || any(positive_scales <= 0))
    stop("All prior and numerical scales must be finite and positive.", call. = FALSE)
  invisible(TRUE)
}

configure_runtime <- function(cfg) {
  if (isTRUE(cfg$runtime$prevent_nested_threading)) {
    Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
  }
  invisible(TRUE)
}

make_output_directories <- function(cfg) {
  dirs <- c(
    cfg$paths$analysis,
    file.path(cfg$paths$analysis, "changepoint"),
    file.path(cfg$paths$analysis, "smooth"),
    file.path(cfg$paths$analysis, "comparison")
  )
  for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  invisible(dirs)
}

compile_project_stan_model <- function(stan_file, cfg, force_recompile = FALSE) {
  check_cmdstan_installation()
  cmdstanr::cmdstan_model(
    stan_file,
    cpp_options = list(stan_threads = TRUE),
    stanc_options = list("O1"),
    force_recompile = force_recompile
  )
}

smooth_config <- function(cfg, output_subdir = "primary") {
  x <- cfg
  x$output_dir <- file.path(cfg$paths$analysis, "smooth", output_subdir)
  x$sampling$current <- cfg$sampling$smooth
  dir.create(file.path(x$output_dir, "fit", "csv"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(x$output_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
  x
}
