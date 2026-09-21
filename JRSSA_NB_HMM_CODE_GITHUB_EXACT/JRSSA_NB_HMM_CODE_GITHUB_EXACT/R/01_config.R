`%||%` <- function(x, y) if (is.null(x)) y else x

is_absolute_path <- function(x) {
  x <- path.expand(x)
  grepl('^/', x) || grepl('^[A-Za-z]:[/\\\\]', x)
}

resolve_project_path <- function(root, x) {
  if (is_absolute_path(x)) normalizePath(x, winslash='/', mustWork=FALSE)
  else normalizePath(file.path(root, x), winslash='/', mustWork=FALSE)
}

resolve_stan_threads <- function(cfg) {
  requested <- as.integer(cfg$sampling$threads_per_chain %||% 0L)
  if (requested > 0L) return(requested)
  phys <- suppressWarnings(parallel::detectCores(logical=FALSE))
  if (!is.finite(phys) || phys < 1L) phys <- suppressWarnings(parallel::detectCores(logical=TRUE))
  if (!is.finite(phys) || phys < 1L) phys <- 4L
  pc <- max(1L, as.integer(cfg$sampling$parallel_chains %||% cfg$sampling$chains))
  cap <- max(1L, as.integer(cfg$runtime$stan_threads_max_per_chain %||% 4L))
  max(1L, min(cap, floor(phys / pc)))
}

read_project_config <- function(config_path, project_root) {
  cfg <- yaml::read_yaml(config_path)
  cfg$project_root <- normalizePath(project_root, winslash='/', mustWork=TRUE)
  cfg$sampling$threads_per_chain <- resolve_stan_threads(cfg)
  for (nm in c('dx_file','ex_file','series_metadata','oecd_universe'))
    cfg$data[[nm]] <- resolve_project_path(cfg$project_root, cfg$data[[nm]])
  for (nm in c('real_full','real_share'))
    cfg$paths[[nm]] <- resolve_project_path(cfg$project_root, cfg$paths[[nm]])
  cfg$output_dir <- cfg$paths$real_full
  cfg$stan_changepoint_file <- file.path(cfg$project_root, 'stan', 'hierarchical_changepoint_v7_light.stan')
  cfg$stan_smooth_file <- file.path(cfg$project_root, 'stan', 'smooth_cohort_benchmark_v7.stan')
  validate_project_config(cfg)
  cfg
}

validate_project_config <- function(cfg) {
  if (as.integer(cfg$data$expected_selected_populations) != 21L)
    stop('The frozen application is pre-specified for 21 populations.', call.=FALSE)
  if (as.integer(cfg$data$cohort_step) != 5L)
    stop('v7-light is configured for 5-year cohort groups.', call.=FALSE)
  positives <- c(
    cfg$model$baseline_country_level_prior_sd,
    cfg$model$baseline_country_shape_prior_sd,
    cfg$model$baseline_global_rw2_prior_sd,
    cfg$model$drift_global_prior_sd,
    cfg$model$drift_country_deviation_prior_sd,
    cfg$model$period_global_prior_sd,
    cfg$model$period_country_deviation_prior_sd,
    cfg$model$regime_range_prior_sd,
    cfg$model$regime_global_shape_prior_sd,
    cfg$model$country_regime_level_prior_sd,
    cfg$model$switch_center_prior_sd_intervals,
    cfg$model$switch_sd_prior_intervals,
    cfg$model$switch_sd_floor_intervals,
    cfg$model$numerical_floor_phi,
    cfg$model$numerical_floor_baseline_rw2_sd
  )
  if (any(!is.finite(positives)) || any(positives <= 0))
    stop('All v7 prior/floor scales must be finite and positive.', call.=FALSE)
  invisible(TRUE)
}

compile_project_stan_model <- function(stan_file, cfg, force_recompile=FALSE, allow_no_range_checks=FALSE) {
  check_cmdstan_installation()
  cmdstanr::cmdstan_model(
    stan_file,
    cpp_options=list(stan_threads=TRUE),
    stanc_options=list('O1'),
    force_recompile=force_recompile
  )
}

configure_runtime <- function(cfg) {
  if (isTRUE(cfg$runtime$prevent_nested_threading))
    Sys.setenv(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
  options(future.globals.maxSize=max(
    getOption('future.globals.maxSize', 0),
    as.numeric(cfg$runtime$future_globals_max_gb %||% 8) * 1024^3
  ))
  invisible(TRUE)
}

make_output_directories <- function(cfg) {
  dirs <- c(
    cfg$paths$real_full,
    file.path(cfg$paths$real_full, 'V7_Light_CP'),
    file.path(cfg$paths$real_full, 'Smooth'),
    file.path(cfg$paths$real_full, 'comparison'),
    file.path(cfg$paths$real_full, 'publication'),
    cfg$paths$real_share
  )
  for (d in dirs) dir.create(d, recursive=TRUE, showWarnings=FALSE)
  invisible(dirs)
}

smooth_config <- function(cfg, output_subdir='primary_21countries') {
  x <- cfg
  x$model$K <- 2L
  x$model$B_cohort <- as.integer(cfg$analysis_spec$smooth_B_cohort %||% cfg$model$B_cohort)
  x$model$cohort_smoothness_sigma <- as.numeric(cfg$analysis_spec$cohort_smoothness_sigma)
  x$output_dir <- file.path(cfg$paths$real_full, 'Smooth', output_subdir)
  x$sampling$iter_warmup <- as.integer(cfg$sampling_by_model$Smooth$iter_warmup)
  x$sampling$iter_sampling <- as.integer(cfg$sampling_by_model$Smooth$iter_sampling)
  dir.create(file.path(x$output_dir, 'fit', 'csv'), recursive=TRUE, showWarnings=FALSE)
  dir.create(file.path(x$output_dir, 'tables'), recursive=TRUE, showWarnings=FALSE)
  x
}
