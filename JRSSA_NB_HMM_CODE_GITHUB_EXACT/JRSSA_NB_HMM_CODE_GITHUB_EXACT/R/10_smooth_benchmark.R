# ======================================================================
# v6 Smooth benchmark: identical APC observation layer to the HSMM,
# without the latent regime component.
# ======================================================================

build_smooth_benchmark_data <- function(pack, cfg) {
  dat <- pack$metadata$dat_long
  sd <- pack$stan_data
  stan_data <- list(
    I = as.integer(sd$I),
    reduce_grainsize = as.integer(cfg$runtime$stan_smooth_reduce_grainsize %||% 64L),
    A = as.integer(sd$A), B_age = as.integer(sd$B_age), B_age_dev = as.integer(sd$B_age_dev),
    X_age_baseline = sd$X_age_baseline, X_age_deviation = sd$X_age_deviation,
    P = as.integer(sd$P), B_period = as.integer(sd$B_period), X_period = sd$X_period, drift_index = sd$drift_index,
    C = as.integer(sd$C), B_cohort = as.integer(sd$B_cohort), X_cohort = sd$X_cohort, D2_cohort = sd$D2_cohort,
    N = as.integer(nrow(dat)), country_idx = as.integer(dat$country_id), age_idx = as.integer(dat$age_idx),
    period_idx = as.integer(dat$period_idx), cohort_idx = as.integer(dat$global_c), include_likelihood = as.integer(dat$include_likelihood),
    has_F = as.integer(dat$has_F), has_M = as.integer(dat$has_M), D_F = as.integer(dat$Deaths_F), D_M = as.integer(dat$Deaths_M),
    logE_F = as.numeric(dat$logE_F), logE_M = as.numeric(dat$logE_M),
    baseline_country_level_prior_sd = as.numeric(cfg$model$baseline_country_level_prior_sd),
    baseline_country_shape_prior_sd = as.numeric(cfg$model$baseline_country_shape_prior_sd),
    baseline_global_rw2_prior_sd = as.numeric(cfg$model$baseline_global_rw2_prior_sd),
    drift_global_prior_sd = as.numeric(cfg$model$drift_global_prior_sd), drift_country_sd_prior = as.numeric(cfg$model$drift_country_sd_prior),
    period_global_prior_sd = as.numeric(cfg$model$period_global_prior_sd),
    period_country_deviation_prior_sd = as.numeric(cfg$model$period_country_deviation_prior_sd),
    cohort_global_prior_sd = as.numeric(cfg$model$cohort_global_prior_sd),
    cohort_country_deviation_prior_sd = as.numeric(cfg$model$cohort_country_deviation_prior_sd),
    cohort_smoothness_sigma = as.numeric(cfg$model$cohort_smoothness_sigma),
    numerical_floor_phi = as.numeric(cfg$model$numerical_floor_phi),
    numerical_floor_baseline_rw2_sd = as.numeric(cfg$model$numerical_floor_baseline_rw2_sd)
  )
  list(stan_data=stan_data, metadata=pack$metadata, model_label='Smooth', model_K=0L)
}

make_smooth_initialization <- function(pack, cfg) {
  sd <- pack$stan_data; meta <- pack$metadata; I <- sd$I
  beta_F <- beta_M <- matrix(0,I,sd$B_age)
  for(i in seq_len(I)) {
    d_i <- dplyr::filter(meta$dat_long,country_id==i)
    beta_F[i,] <- fit_age_baseline(d_i,meta$ages,sd$X_age_baseline,'Deaths_F','Exp_F')
    beta_M[i,] <- fit_age_baseline(d_i,meta$ages,sd$X_age_baseline,'Deaths_M','Exp_M')
  }
  global_F <- colMeans(beta_F); global_M <- colMeans(beta_M)
  function(chain_id=1L) {
    set.seed(as.integer(cfg$sampling$seed+70000L+100L*chain_id))
    list(
      baseline_age_global_F=global_F+rnorm(sd$B_age,0,.01), baseline_age_global_M=global_M+rnorm(sd$B_age,0,.01),
      baseline_global_rw2_sd_F=.08, baseline_global_rw2_sd_M=.08,
      baseline_country_level_raw_F=rnorm(I,0,.03), baseline_country_level_raw_M=rnorm(I,0,.03),
      baseline_country_shape_raw_F=matrix(rnorm(I*sd$B_age_dev,0,.02),I,sd$B_age_dev), baseline_country_shape_raw_M=matrix(rnorm(I*sd$B_age_dev,0,.02),I,sd$B_age_dev),
      drift_global_F=-.08, drift_global_M=-.07, drift_country_sd_F=.025, drift_country_sd_M=.025,
      drift_country_raw_F=rnorm(I,0,.02), drift_country_raw_M=rnorm(I,0,.02),
      period_global_coef_F=rep(0,sd$B_period), period_global_coef_M=rep(0,sd$B_period),
      period_country_raw_F=matrix(0,I,sd$B_period), period_country_raw_M=matrix(0,I,sd$B_period),
      cohort_global_coef_F=rep(0,sd$B_cohort), cohort_global_coef_M=rep(0,sd$B_cohort),
      cohort_country_raw_F=matrix(0,I,sd$B_cohort), cohort_country_raw_M=matrix(0,I,sd$B_cohort),
      phi_F=8, phi_M=8
    )
  }
}

fit_smooth_cohort_benchmark <- function(pack,cfg) {
  model <- compile_project_stan_model(cfg$stan_smooth_file,cfg,force_recompile=FALSE,allow_no_range_checks=TRUE)
  final_csv <- file.path(cfg$output_dir,'fit','csv'); clean_directory(final_csv)
  tune <- cfg$sampling_by_model$Smooth
  message('v6 Smooth APC benchmark: ',cfg$sampling$chains,' chains; ',cfg$sampling$iter_warmup,' warmup + ',cfg$sampling$iter_sampling,' sampling.')
  fit <- model$sample(
    data=pack$stan_data, seed=as.integer(cfg$sampling$seed+70000L), chains=cfg$sampling$chains,
    parallel_chains=min(cfg$sampling$parallel_chains,cfg$sampling$chains), threads_per_chain=as.integer(cfg$sampling$threads_per_chain),
    iter_warmup=cfg$sampling$iter_warmup, iter_sampling=cfg$sampling$iter_sampling,
    adapt_delta=as.numeric(tune$adapt_delta), step_size=as.numeric(tune$initial_step_size), max_treedepth=as.integer(tune$max_treedepth),
    refresh=cfg$sampling$refresh, init=make_smooth_initialization(pack,cfg), output_basename='sm', save_metric=TRUE
  )
  if(any(fit$return_codes()!=0L)) stop('One or more v6 Smooth Stan chains failed.',call.=FALSE)
  if(isTRUE(cfg$sampling$save_cmdstan_csv)) { dir.create(final_csv,recursive=TRUE,showWarnings=FALSE); fit$save_output_files(dir=final_csv,basename='sm',timestamp=FALSE,random=FALSE) }
  fit$save_object(file.path(cfg$output_dir,'fit','smooth_cohort_fit_cmdstanr.rds'))
  fit
}

smooth_core_parameter_mask <- function(parameters, pack) {
  grepl("^cohort_global_coef_[FM]\\[|^cohort_country_coef_[FM]\\[|^drift_global_[FM]$|^drift_country_sd_[FM]$|^drift_country_[FM]\\[|^phi_[FM]$", parameters)
}

summarise_smooth_diagnostics <- function(fit, cfg, pack, prefix = "smooth") {
  q025 <- function(x) stats::quantile(x, 0.025)
  q975 <- function(x) stats::quantile(x, 0.975)

  diagnostics <- posterior::summarise_draws(
    fit$draws(), mean, median, sd, q025, q975,
    posterior::rhat, posterior::ess_bulk, posterior::ess_tail
  ) |>
    as.data.frame() |>
    normalise_diagnostic_column_names()

  names(diagnostics)[names(diagnostics) == "variable"] <- "parameter"
  core <- smooth_core_parameter_mask(diagnostics$parameter, pack)
  diagnostics$parameter_group <- ifelse(core, "supported_substantive", "nuisance_or_generated")
  diagnostics$core_supported_parameter <- core

  report_rhat <- as.numeric(cfg$quality_control$reporting_rhat_threshold)
  report_bulk <- as.numeric(cfg$quality_control$reporting_min_bulk_ess)
  report_tail <- as.numeric(cfg$quality_control$reporting_min_tail_ess)
  hard_rhat <- as.numeric(cfg$quality_control$hard_rhat_threshold)
  hard_bulk <- as.numeric(cfg$quality_control$hard_min_bulk_ess)
  hard_tail <- as.numeric(cfg$quality_control$hard_min_tail_ess)

  diagnostics$report_bad_rhat <- is.finite(diagnostics$rhat) & diagnostics$rhat > report_rhat
  diagnostics$report_low_bulk_ess <- is.finite(diagnostics$ess_bulk) & diagnostics$ess_bulk < report_bulk
  diagnostics$report_low_tail_ess <- is.finite(diagnostics$ess_tail) & diagnostics$ess_tail < report_tail
  diagnostics$hard_bad_rhat <- is.finite(diagnostics$rhat) & diagnostics$rhat > hard_rhat
  diagnostics$hard_low_bulk_ess <- is.finite(diagnostics$ess_bulk) & diagnostics$ess_bulk < hard_bulk
  diagnostics$hard_low_tail_ess <- is.finite(diagnostics$ess_tail) & diagnostics$ess_tail < hard_tail

  sampler <- posterior::as_draws_matrix(fit$sampler_diagnostics())
  divergent <- if ("divergent__" %in% colnames(sampler)) sum(sampler[, "divergent__"]) else NA_real_
  tune <- cfg$sampling_by_model$Smooth
  treedepth <- if ("treedepth__" %in% colnames(sampler)) {
    sum(sampler[, "treedepth__"] >= as.integer(tune$max_treedepth))
  } else NA_real_

  cmdstan_diag <- suppressWarnings(fit$diagnostic_summary()) |> as.data.frame()
  cmdstan_diag$chain <- seq_len(nrow(cmdstan_diag))

  bfmi_name <- intersect(c("ebfmi", "e_bfmi", "E-BFMI", "bfmi"), names(cmdstan_diag))
  low_bfmi <- if (length(bfmi_name)) {
    sum(cmdstan_diag[[bfmi_name[1]]] < cfg$quality_control$bfmi_threshold, na.rm = TRUE)
  } else NA_integer_

  hard_fail <- diagnostics$hard_bad_rhat |
    diagnostics$hard_low_bulk_ess |
    diagnostics$hard_low_tail_ess
  report_fail <- diagnostics$report_bad_rhat |
    diagnostics$report_low_bulk_ess |
    diagnostics$report_low_tail_ess

  sampler_pass <-
    (is.na(divergent) || divergent == 0L) &&
    (is.na(treedepth) || treedepth <= cfg$quality_control$max_treedepth_allowed) &&
    (is.na(low_bfmi) || low_bfmi == 0L)

  core_pass <- !any(hard_fail[core], na.rm = TRUE)
  eligible <- sampler_pass && core_pass
  reporting_pass <- sampler_pass && !any(report_fail, na.rm = TRUE)
  status <- if (!eligible) "FAIL" else if (reporting_pass) "PASS" else "PASS_WITH_WARNINGS"

  qc <- data.frame(
    diagnostic_status = status,
    scientifically_eligible = eligible,
    strict_reporting_targets_passed = reporting_pass,
    full_parameters_above_reporting_rhat = sum(diagnostics$report_bad_rhat, na.rm = TRUE),
    full_parameters_below_reporting_bulk_ess = sum(diagnostics$report_low_bulk_ess, na.rm = TRUE),
    full_parameters_below_reporting_tail_ess = sum(diagnostics$report_low_tail_ess, na.rm = TRUE),
    supported_parameters_above_hard_rhat = sum(diagnostics$hard_bad_rhat[core], na.rm = TRUE),
    supported_parameters_below_hard_bulk_ess = sum(diagnostics$hard_low_bulk_ess[core], na.rm = TRUE),
    supported_parameters_below_hard_tail_ess = sum(diagnostics$hard_low_tail_ess[core], na.rm = TRUE),
    divergent_transitions = divergent,
    max_treedepth_hits = treedepth,
    chains_below_bfmi_threshold = low_bfmi,
    passed = eligible
  )

  readr::write_csv(
    diagnostics,
    file.path(cfg$output_dir, "tables", paste0(prefix, "_mcmc_parameter_diagnostics.csv"))
  )
  readr::write_csv(
    cmdstan_diag,
    file.path(cfg$output_dir, "tables", paste0(prefix, "_sampler_diagnostics_by_chain.csv"))
  )
  readr::write_csv(
    qc,
    file.path(cfg$output_dir, "tables", paste0(prefix, "_mcmc_quality_control.csv"))
  )
  readr::write_csv(
    diagnostics[report_fail, , drop = FALSE],
    file.path(cfg$output_dir, "tables", paste0(prefix, "_parameters_warning_reporting_targets.csv"))
  )
  readr::write_csv(
    diagnostics[core & hard_fail, , drop = FALSE],
    file.path(cfg$output_dir, "tables", paste0(prefix, "_parameters_failing_hard_supported_qc.csv"))
  )

  list(
    parameter_summary = diagnostics,
    sampler_by_chain = cmdstan_diag,
    quality_control = qc,
    status = status,
    eligible = eligible
  )
}
