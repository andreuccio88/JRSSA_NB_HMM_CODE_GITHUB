# v7-light change-point helpers.
# Primary v7-light helper used by this standalone project.
# It reuses prepare_mortality_data(), orthogonal_spline_basis(), fit_age_baseline(),
# and the project configuration/package infrastructure already present there.

`%||%` <- function(x, y) if (is.null(x)) y else x

v7_light_defaults <- function(cfg) {
  cfg$model$B_age <- cfg$model$B_age %||% 8L
  cfg$model$B_age_deviation <- cfg$model$B_age_deviation %||% 5L
  cfg$model$B_period <- cfg$model$B_period %||% 6L
  cfg$model$B_regime <- cfg$model$B_regime %||% 6L

  # Keep v6 APC priors, but FIX nuisance cross-country scales rather than estimating them.
  cfg$model$drift_country_deviation_prior_sd <- 0.06
  cfg$model$country_regime_level_prior_sd <- 0.08

  # Arbitrary numerical origin, deliberately NOT the 1925-1935 empirical window.
  cfg$model$switch_reference_year <- 1900
  cfg$model$switch_center_prior_sd_intervals <- 20   # 100 years with 5y cohorts
  cfg$model$switch_sd_prior_intervals <- 8           # 40-year half-normal scale
  cfg$model$switch_sd_floor_intervals <- 0.25        # 1.25 years numerical floor

  cfg$sampling$v7_iter_warmup <- 600L
  cfg$sampling$v7_iter_sampling <- 600L
  cfg$sampling$v7_adapt_delta <- 0.93
  cfg$sampling$v7_max_treedepth <- 11L
  cfg
}

build_changepoint_stan_data <- function(dat, cfg, cohorts_global = NULL) {
  cfg <- v7_light_defaults(cfg)
  countries <- sort(unique(as.character(dat$country)))
  I <- length(countries)
  ages <- sort(unique(as.integer(dat$Age)))
  periods <- sort(unique(as.integer(dat$Period)))
  step <- as.integer(cfg$data$cohort_step)

  if (is.null(cohorts_global)) {
    cohorts_global <- seq(min(dat$Cohort), max(dat$Cohort), by = step)
  } else {
    cohorts_global <- sort(unique(as.integer(cohorts_global)))
  }
  if (any(diff(cohorts_global) != step)) {
    stop("cohorts_global must be a contiguous grid with the configured cohort step.", call. = FALSE)
  }
  if (!all(unique(dat$Cohort) %in% cohorts_global)) {
    stop("Some retained cohorts fall outside the frozen global switch grid.", call. = FALSE)
  }

  A <- length(ages)
  P <- length(periods)
  C_global <- length(cohorts_global)
  if (A < 4L || P < 5L || C_global < 5L) stop("Lexis grid too small for v7-light.")

  X_age_baseline <- as.matrix(splines::bs(
    ages, df = cfg$model$B_age, degree = 3, intercept = TRUE
  ))
  X_age_deviation <- orthogonal_spline_basis(
    ages, cfg$model$B_age_deviation, remove_linear = FALSE
  )
  X_age_regime <- orthogonal_spline_basis(
    ages, cfg$model$B_regime, remove_linear = FALSE
  )
  X_period <- orthogonal_spline_basis(
    periods, cfg$model$B_period, remove_linear = TRUE
  )
  drift_index <- as.numeric((periods - mean(periods)) / as.numeric(cfg$data$period_step))

  dat_indexed <- dat |>
    dplyr::mutate(
      country_id = match(country, countries),
      age_idx = match(Age, ages),
      period_idx = match(Period, periods),
      global_c = match(Cohort, cohorts_global)
    )
  if (anyNA(dat_indexed$country_id) || anyNA(dat_indexed$age_idx) ||
      anyNA(dat_indexed$period_idx) || anyNA(dat_indexed$global_c)) {
    stop("Internal v7-light indexing failure.", call. = FALSE)
  }

  mask <- mask_window_from_cfg(cfg)
  dat_indexed$include_likelihood <- 1L
  if (!is.null(mask)) {
    dat_indexed$include_likelihood[
      dat_indexed$Cohort >= mask[1] & dat_indexed$Cohort <= mask[2]
    ] <- 0L
  }

  active_start <- integer(I)
  active_end <- integer(I)
  local_grid <- vector("list", I)
  T_vec <- integer(I)
  for (i in seq_len(I)) {
    obs_c <- sort(unique(dat_indexed$global_c[dat_indexed$country_id == i]))
    active_start[i] <- min(obs_c)
    active_end[i] <- max(obs_c)
    if (any(diff(obs_c) != 1L)) {
      stop("v7-light requires contiguous supported cohort grids within country: ",
           countries[i], call. = FALSE)
    }
    local_grid[[i]] <- seq.int(active_start[i], active_end[i])
    T_vec[i] <- length(local_grid[[i]])
  }
  maxT <- max(T_vec)

  n_obs <- matrix(0L, I, maxT)
  start_idx <- matrix(1L, I, maxT)
  end_idx <- matrix(1L, I, maxT)
  pieces <- vector("list", I)
  cohort_map <- vector("list", I)
  cursor <- 0L

  for (i in seq_len(I)) {
    d_i <- dat_indexed |>
      dplyr::filter(country_id == i) |>
      dplyr::mutate(t = global_c - active_start[i] + 1L) |>
      dplyr::arrange(t, Period, Age)

    for (t in seq_len(T_vec[i])) {
      local_rows <- which(d_i$t == t)
      if (length(local_rows)) {
        n_obs[i, t] <- length(local_rows)
        start_idx[i, t] <- cursor + min(local_rows)
        end_idx[i, t] <- cursor + max(local_rows)
      }
    }
    pieces[[i]] <- d_i
    cohort_map[[i]] <- data.frame(
      country = countries[i], t = seq_len(T_vec[i]),
      global_c = local_grid[[i]], Cohort = cohorts_global[local_grid[[i]]]
    )
    cursor <- cursor + nrow(d_i)
  }

  dat_long <- dplyr::bind_rows(pieces) |>
    dplyr::mutate(
      has_F = as.integer(Exp_F > 0),
      has_M = as.integer(Exp_M > 0),
      logE_F = dplyr::if_else(Exp_F > 0, log(Exp_F), 0),
      logE_M = dplyr::if_else(Exp_M > 0, log(Exp_M), 0)
    )
  N <- nrow(dat_long)

  # g is the first cohort assigned to post-switch state 1.
  # C+1 includes the all-pre possibility (switch after the final cohort).
  switch_grid_year <- c(cohorts_global, max(cohorts_global) + step)
  switch_ref <- as.numeric(cfg$model$switch_reference_year)
  switch_grid_std <- (switch_grid_year - switch_ref) / step

  stan_data <- list(
    I = as.integer(I),
    reduce_grainsize = as.integer(cfg$runtime$stan_reduce_grainsize %||% 1L),
    A = as.integer(A),
    B_age = as.integer(ncol(X_age_baseline)),
    B_age_dev = as.integer(ncol(X_age_deviation)),
    B_regime = as.integer(ncol(X_age_regime)),
    X_age_baseline = unname(X_age_baseline),
    X_age_deviation = unname(X_age_deviation),
    X_age_regime = unname(X_age_regime),
    P = as.integer(P),
    B_period = as.integer(ncol(X_period)),
    X_period = unname(X_period),
    drift_index = as.numeric(drift_index),
    C = as.integer(C_global),
    J_switch = as.integer(length(switch_grid_year)),
    switch_grid_std = as.numeric(switch_grid_std),
    switch_reference_year = switch_ref,
    cohort_step_years = as.numeric(step),
    maxT = as.integer(maxT),
    T = as.integer(T_vec),
    active_start = as.integer(active_start),
    N = as.integer(N),
    age_idx = as.integer(dat_long$age_idx),
    period_idx = as.integer(dat_long$period_idx),
    include_likelihood = as.integer(dat_long$include_likelihood),
    has_F = as.integer(dat_long$has_F),
    has_M = as.integer(dat_long$has_M),
    D_F = as.integer(dat_long$Deaths_F),
    D_M = as.integer(dat_long$Deaths_M),
    logE_F = as.numeric(dat_long$logE_F),
    logE_M = as.numeric(dat_long$logE_M),
    n_obs = unname(n_obs),
    start_idx = unname(start_idx),
    end_idx = unname(end_idx),
    baseline_country_level_prior_sd = as.numeric(cfg$model$baseline_country_level_prior_sd),
    baseline_country_shape_prior_sd = as.numeric(cfg$model$baseline_country_shape_prior_sd),
    baseline_global_rw2_prior_sd = as.numeric(cfg$model$baseline_global_rw2_prior_sd),
    drift_global_prior_sd = as.numeric(cfg$model$drift_global_prior_sd),
    drift_country_deviation_prior_sd = as.numeric(cfg$model$drift_country_deviation_prior_sd),
    period_global_prior_sd = as.numeric(cfg$model$period_global_prior_sd),
    period_country_deviation_prior_sd = as.numeric(cfg$model$period_country_deviation_prior_sd),
    regime_range_prior_sd = as.numeric(cfg$model$regime_range_prior_sd),
    regime_global_shape_prior_sd = as.numeric(cfg$model$regime_global_shape_prior_sd),
    country_regime_level_prior_sd = as.numeric(cfg$model$country_regime_level_prior_sd),
    switch_center_prior_sd_intervals = as.numeric(cfg$model$switch_center_prior_sd_intervals),
    switch_sd_prior_intervals = as.numeric(cfg$model$switch_sd_prior_intervals),
    switch_sd_floor_intervals = as.numeric(cfg$model$switch_sd_floor_intervals),
    numerical_floor_phi = as.numeric(cfg$model$numerical_floor_phi),
    numerical_floor_baseline_rw2_sd = as.numeric(cfg$model$numerical_floor_baseline_rw2_sd)
  )

  metadata <- list(
    countries = countries,
    ages = ages,
    periods = periods,
    cohorts_global = cohorts_global,
    switch_grid_year = switch_grid_year,
    active_start = active_start,
    active_end = active_end,
    T = T_vec,
    cohort_map = dplyr::bind_rows(cohort_map),
    dat_long = dat_long,
    cohort_step = step,
    mask_window = mask
  )
  list(stan_data = stan_data, metadata = metadata)
}

make_changepoint_initialization <- function(pack, cfg) {
  cfg <- v7_light_defaults(cfg)
  sd <- pack$stan_data
  I <- sd$I
  meta <- pack$metadata
  betaF <- betaM <- matrix(0, I, sd$B_age)

  for (i in seq_len(I)) {
    di <- dplyr::filter(meta$dat_long, country_id == i)
    betaF[i, ] <- fit_age_baseline(
      di, meta$ages, sd$X_age_baseline, "Deaths_F", "Exp_F"
    )
    betaM[i, ] <- fit_age_baseline(
      di, meta$ages, sd$X_age_baseline, "Deaths_M", "Exp_M"
    )
  }
  globalF <- colMeans(betaF)
  globalM <- colMeans(betaM)

  function(chain_id = 1L) {
    set.seed(cfg$sampling$seed + 701L * chain_id)
    list(
      baseline_age_global_F = globalF + rnorm(sd$B_age, 0, .01),
      baseline_age_global_M = globalM + rnorm(sd$B_age, 0, .01),
      baseline_global_rw2_sd_F = .08,
      baseline_global_rw2_sd_M = .08,
      baseline_country_level_raw_F = rnorm(I, 0, .03),
      baseline_country_level_raw_M = rnorm(I, 0, .03),
      baseline_country_shape_raw_F = matrix(rnorm(I * sd$B_age_dev, 0, .02), I, sd$B_age_dev),
      baseline_country_shape_raw_M = matrix(rnorm(I * sd$B_age_dev, 0, .02), I, sd$B_age_dev),
      drift_global_F = -.08,
      drift_global_M = -.07,
      drift_country_raw_F = rnorm(I, 0, .02),
      drift_country_raw_M = rnorm(I, 0, .02),
      period_global_coef_F = rep(0, sd$B_period),
      period_global_coef_M = rep(0, sd$B_period),
      period_country_raw_F = matrix(0, I, sd$B_period),
      period_country_raw_M = matrix(0, I, sd$B_period),
      regime_range_F = .20,
      regime_range_M = .20,
      regime_shape_F = rep(0, sd$B_regime),
      regime_shape_M = rep(0, sd$B_regime),
      country_regime_level_raw_F = rep(0, I),
      country_regime_level_raw_M = rep(0, I),
      phi_F = 8,
      phi_M = 8,
      # Deliberately dispersed initial locations across chains; no chain is forced to 1930.
      switch_center_std = c(0, 4, 8, 12)[((chain_id - 1L) %% 4L) + 1L] + rnorm(1, 0, .25),
      switch_sd_raw = 2
    )
  }
}

fit_changepoint_v7 <- function(pack, cfg, stan_file, output_dir,
                               iter_warmup = NULL, iter_sampling = NULL,
                               adapt_delta = NULL, max_treedepth = NULL) {
  cfg <- v7_light_defaults(cfg)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  csv_dir <- file.path(output_dir, "csv")
  dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)

  model <- cmdstanr::cmdstan_model(
    stan_file,
    cpp_options = list(stan_threads = TRUE),
    stanc_options = list("O1")
  )
  init <- make_changepoint_initialization(pack, cfg)
  fit <- model$sample(
    data = pack$stan_data,
    seed = cfg$sampling$seed,
    chains = cfg$sampling$chains,
    parallel_chains = min(cfg$sampling$parallel_chains, cfg$sampling$chains),
    threads_per_chain = as.integer(cfg$sampling$threads_per_chain %||% 2L),
    iter_warmup = as.integer(iter_warmup %||% cfg$sampling$v7_iter_warmup),
    iter_sampling = as.integer(iter_sampling %||% cfg$sampling$v7_iter_sampling),
    adapt_delta = as.numeric(adapt_delta %||% cfg$sampling$v7_adapt_delta),
    max_treedepth = as.integer(max_treedepth %||% cfg$sampling$v7_max_treedepth),
    refresh = cfg$sampling$refresh,
    init = init,
    output_basename = "v7_light_cp"
  )
  if (any(fit$return_codes() != 0L)) stop("One or more v7-light Stan chains failed.", call. = FALSE)
  fit$save_output_files(dir = csv_dir, basename = "v7_light_cp", timestamp = FALSE, random = FALSE)
  fit$save_object(file.path(output_dir, "fit_cmdstanr.rds"))
  fit
}

summarise_changepoint_v7 <- function(fit, pack, output_dir) {
  dir.create(file.path(output_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
  countries <- pack$metadata$countries
  grid <- pack$metadata$switch_grid_year

  pars <- c(
    "switch_center_year", "switch_sd_years", "regime_range_F", "regime_range_M",
    "phi_F", "phi_M", "drift_global_F", "drift_global_M"
  )
  sm <- fit$summary(variables = pars)
  readr::write_csv(sm, file.path(output_dir, "tables", "key_parameter_summary.csv"))

  # Each Stan draw contains a full posterior switch distribution for every country.
  probs <- fit$draws("switch_probability", format = "draws_matrix")
  ans <- lapply(seq_along(countries), function(i) {
    cols <- grep(paste0("^switch_probability\\[", i, ","), colnames(probs))
    P <- as.matrix(probs[, cols, drop = FALSE])
    pbar <- colMeans(P)
    data.frame(
      country = countries[i],
      switch_year = grid,
      posterior_probability = pbar
    )
  }) |>
    dplyr::bind_rows()
  readr::write_csv(ans, file.path(output_dir, "tables", "country_switch_probabilities.csv"))

  sw <- ans |>
    dplyr::group_by(country) |>
    dplyr::summarise(
      posterior_mean_switch = sum(switch_year * posterior_probability),
      map_switch = switch_year[which.max(posterior_probability)],
      p_1925_1935 = sum(posterior_probability[switch_year >= 1925 & switch_year <= 1935]),
      .groups = "drop"
    )
  readr::write_csv(sw, file.path(output_dir, "tables", "country_switch_summary.csv"))

  # The masked kernel is integrated over posterior draws with log-mean-exp.
  k <- fit$draws("log_masked_predictive_kernel_country", format = "draws_matrix")
  log_mean_exp <- function(x) {
    m <- max(x)
    m + log(mean(exp(x - m)))
  }
  masked <- data.frame(
    country = countries,
    log_predictive_density = vapply(seq_along(countries), function(i) {
      log_mean_exp(k[, i])
    }, numeric(1))
  )
  readr::write_csv(masked, file.path(output_dir, "tables", "masked_predictive_density_by_country.csv"))

  diag <- posterior::summarise_draws(
    fit$draws(), posterior::rhat, posterior::ess_bulk, posterior::ess_tail
  ) |> as.data.frame()
  readr::write_csv(diag, file.path(output_dir, "tables", "mcmc_diagnostics_all_parameters.csv"))
  sampler <- tryCatch(as.data.frame(fit$diagnostic_summary()), error = function(e) data.frame())
  if (nrow(sampler)) readr::write_csv(sampler, file.path(output_dir, "tables", "sampler_diagnostics_by_chain.csv"))

  invisible(list(key = sm, switch_probabilities = ans, switch_summary = sw, masked = masked, diagnostics = diag))
}
