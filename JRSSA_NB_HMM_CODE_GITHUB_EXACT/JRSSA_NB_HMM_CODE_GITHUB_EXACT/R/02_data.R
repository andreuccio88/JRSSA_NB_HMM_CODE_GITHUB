read_hmd_inputs <- function(cfg) {
  src <- read_hmd_sources_separately(cfg)
  dplyr::inner_join(src$dx, src$ex, by = c("country", "Year", "Age"))
}

longest_contiguous_block <- function(cohorts, step = 5L) {
  cohorts <- sort(unique(as.integer(cohorts)))
  if (!length(cohorts)) return(integer())
  run_id <- cumsum(c(TRUE, diff(cohorts) != as.integer(step)))
  groups <- split(cohorts, run_id)
  lens <- lengths(groups)
  groups[[which.max(lens)]]
}

prepare_mortality_data <- function(cfg) {
  selected <- selected_country_codes(cfg, ensure_manifest = TRUE)
  raw <- read_hmd_inputs(cfg)
  missing <- setdiff(selected, intersect(unique(raw$country), selected))
  if (length(missing)) stop("Selected HMD codes missing after Dx/Ex join: ", paste(missing, collapse = ", "), call. = FALSE)

  dat <- raw |>
    dplyr::filter(
      country %in% selected,
      Year >= cfg$data$year_min, Year <= cfg$data$year_max,
      Age >= cfg$data$age_min, Age <= cfg$data$age_max
    ) |>
    dplyr::mutate(
      Age_lower = 5L * floor(Age / 5L),
      Period_lower = 5L * floor(Year / 5L),
      Age = as.integer(Age_lower + 2L),
      Period = as.integer(Period_lower + 2L),
      BirthCohort = as.integer(Period - Age),
      Cohort = as.integer(5L * floor(BirthCohort / 5L)),
      Deaths_F = dplyr::if_else(is.finite(Deaths_F), pmax(Deaths_F, 0), 0),
      Deaths_M = dplyr::if_else(is.finite(Deaths_M), pmax(Deaths_M, 0), 0),
      Exp_F = dplyr::if_else(is.finite(Exp_F), pmax(Exp_F, 0), 0),
      Exp_M = dplyr::if_else(is.finite(Exp_M), pmax(Exp_M, 0), 0)
    ) |>
    dplyr::group_by(country, Period, Age, Cohort) |>
    dplyr::summarise(
      Deaths_F = sum(Deaths_F, na.rm = TRUE), Deaths_M = sum(Deaths_M, na.rm = TRUE),
      Exp_F = sum(Exp_F, na.rm = TRUE), Exp_M = sum(Exp_M, na.rm = TRUE), .groups = "drop"
    ) |>
    dplyr::mutate(
      Age_group = sprintf("%d-%d", Age - 2L, Age + 2L),
      Period_group = sprintf("%d-%d", Period - 2L, Period + 2L),
      Cohort_group = sprintf("%d-%d", Cohort, Cohort + 4L),
      Deaths_F = dplyr::if_else(Exp_F > 0, as.integer(round(Deaths_F)), 0L),
      Deaths_M = dplyr::if_else(Exp_M > 0, as.integer(round(Deaths_M)), 0L),
      Exp_F = dplyr::if_else(Exp_F > 0, Exp_F, 0),
      Exp_M = dplyr::if_else(Exp_M > 0, Exp_M, 0)
    ) |>
    dplyr::filter(Exp_F > 0 | Exp_M > 0) |>
    dplyr::arrange(country, Cohort, Period, Age)

  if (!nrow(dat)) stop("No mortality cells remain after objective sample selection.", call. = FALSE)

  support <- dat |>
    dplyr::group_by(country, Cohort) |>
    dplyr::summarise(n_age_groups = dplyr::n_distinct(Age), exposure = sum(Exp_F + Exp_M), .groups = "drop") |>
    dplyr::filter(n_age_groups >= as.integer(cfg$data$min_age_groups_per_supported_cohort))

  kept <- lapply(split(support, support$country), function(z) {
    block <- longest_contiguous_block(z$Cohort, cfg$data$cohort_step)
    data.frame(country = unique(z$country), Cohort = block)
  }) |>
    dplyr::bind_rows()

  counts <- kept |> dplyr::count(country, name = "supported_cohort_groups")
  bad <- counts |> dplyr::filter(supported_cohort_groups < as.integer(cfg$data$min_supported_cohorts)) |> dplyr::pull(country)
  if (length(bad)) stop("Populations fail the contiguous supported-cohort safeguard: ", paste(bad, collapse = ", "), call. = FALSE)

  # v6: the support criterion filters the actual HSMM cells, not just country eligibility.
  dat <- dplyr::semi_join(dat, kept, by = c("country", "Cohort")) |>
    dplyr::arrange(country, Cohort, Period, Age)

  expected <- as.integer(cfg$data$expected_selected_populations %||% 21L)
  if (dplyr::n_distinct(dat$country) != expected) {
    stop("The v7-light application is locked to ", expected, " selected populations; found ", dplyr::n_distinct(dat$country), ".", call. = FALSE)
  }
  attr(dat, "supported_blocks") <- kept
  dat
}

truncate_supported_cohort_window <- function(dat, cfg, direction = c("drop_oldest", "drop_newest"), n_groups = NULL) {
  direction <- match.arg(direction)
  n_groups <- as.integer(n_groups %||% cfg$analysis_spec$support_shift_cohorts %||% 5L)
  step <- as.integer(cfg$data$cohort_step)
  out <- lapply(split(dat, dat$country), function(z) {
    cs <- sort(unique(z$Cohort))
    if (length(cs) <= n_groups + 3L) stop("Too few cohorts for support-shift test in ", unique(z$country), call. = FALSE)
    keep <- if (direction == "drop_oldest") cs[(n_groups + 1L):length(cs)] else cs[seq_len(length(cs) - n_groups)]
    z[z$Cohort %in% keep, , drop = FALSE]
  })
  dplyr::bind_rows(out) |> dplyr::arrange(country, Cohort, Period, Age)
}

select_pilot_countries <- function(dat, n = 5L) {
  tab <- dat |>
    dplyr::group_by(country) |>
    dplyr::summarise(exposure = sum(Exp_F + Exp_M), n_cohorts = dplyr::n_distinct(Cohort), .groups = "drop") |>
    dplyr::arrange(exposure, country)
  n <- min(as.integer(n), nrow(tab))
  idx <- unique(pmax(1L, pmin(nrow(tab), round(seq(1, nrow(tab), length.out = n)))))
  while (length(idx) < n) idx <- sort(unique(c(idx, setdiff(seq_len(nrow(tab)), idx)[1])))
  tab$country[idx[seq_len(n)]]
}

orthogonal_spline_basis <- function(x, df, remove_linear = FALSE) {
  x <- as.numeric(x)
  raw <- as.matrix(splines::bs(x, df = as.integer(df), degree = 3, intercept = TRUE))
  nuisance <- if (isTRUE(remove_linear)) cbind(1, as.numeric(scale(x, center = TRUE, scale = FALSE))) else matrix(1, nrow = length(x), ncol = 1)
  residual <- qr.resid(qr(nuisance), raw)
  q <- qr(residual)
  if (q$rank < 1L) stop("Spline basis has zero residual rank.", call. = FALSE)
  keep <- min(as.integer(df), q$rank)
  qr.Q(q, complete = FALSE)[, seq_len(keep), drop = FALSE]
}

second_difference_matrix <- function(n) {
  n <- as.integer(n)
  if (n < 3L) stop("Need at least 3 points for second differences.")
  D <- matrix(0, n - 2L, n)
  for (j in seq_len(n - 2L)) D[j, j:(j + 2L)] <- c(1, -2, 1)
  D
}

cohort_local_linear_leakage <- function(dat, cohorts_global, X_cohort) {
  ans <- lapply(split(dat, dat$country), function(z) {
    rows <- match(sort(unique(z$Cohort)), cohorts_global)
    X <- X_cohort[rows, , drop = FALSE]
    x <- sort(unique(z$Cohort))
    Q <- cbind(1, as.numeric(scale(x, center = TRUE, scale = FALSE)))
    fit <- qr.fitted(qr(Q), X)
    ratios <- sqrt(colSums(fit^2) / pmax(colSums(X^2), .Machine$double.eps))
    data.frame(country = unique(z$country), n_cohorts = length(rows), max_linear_leakage = max(ratios), rms_linear_leakage = sqrt(mean(ratios^2)))
  })
  dplyr::bind_rows(ans)
}

mask_window_from_cfg <- function(cfg) {
  w <- cfg$mask_window %||% NULL
  if (is.null(w) || !isTRUE(w$enabled)) return(NULL)
  c(as.integer(w$start), as.integer(w$end))
}

window_design_features <- function(dat, start, end) {
  z <- dat |> dplyr::filter(Cohort >= start, Cohort <= end)
  if (!nrow(z)) return(NULL)
  ex <- z$Exp_F + z$Exp_M
  wmean <- if (sum(ex) > 0) weighted.mean(z$Age, ex) else mean(z$Age)
  wsd <- if (sum(ex) > 0) sqrt(weighted.mean((z$Age - wmean)^2, ex)) else stats::sd(z$Age)
  data.frame(
    start = start, end = end, n_countries = dplyr::n_distinct(z$country), n_cells = nrow(z), total_exposure = sum(ex),
    mean_age = wmean, sd_age = wsd,
    share_age_70plus = if (sum(ex) > 0) sum(ex[z$Age >= 70]) / sum(ex) else mean(z$Age >= 70),
    share_age_80plus = if (sum(ex) > 0) sum(ex[z$Age >= 80]) / sum(ex) else mean(z$Age >= 80)
  )
}

choose_matched_placebo_window <- function(dat, cfg) {
  step <- as.integer(cfg$data$cohort_step)
  target_start <- as.integer(cfg$analysis_spec$target_mask_start)
  target_end <- as.integer(cfg$analysis_spec$target_mask_end)
  width_n <- as.integer((target_end - target_start) / step + 1L)
  buffer <- as.integer(cfg$analysis_spec$placebo_buffer_cohorts %||% 4L)
  target <- window_design_features(dat, target_start, target_end)
  starts <- seq(min(dat$Cohort), max(dat$Cohort) - step * (width_n - 1L), by = step)
  candidates <- lapply(starts, function(s) {
    e <- s + step * (width_n - 1L)
    far_enough <- (e <= target_start - buffer * step) || (s >= target_end + buffer * step)
    if (!far_enough) return(NULL)
    window_design_features(dat, s, e)
  }) |> dplyr::bind_rows()
  if (!nrow(candidates)) stop("No placebo window satisfies the pre-specified buffer.", call. = FALSE)
  feats <- c("n_countries", "n_cells", "total_exposure", "mean_age", "sd_age", "share_age_70plus", "share_age_80plus")
  log_feats <- c("n_cells", "total_exposure")
  transform_feature <- function(x, nm) if (nm %in% log_feats) log1p(x) else x
  allx <- dplyr::bind_rows(target, candidates)
  scales <- vapply(feats, function(nm) {
    z <- transform_feature(allx[[nm]], nm)
    stats::sd(z[is.finite(z)])
  }, numeric(1))
  scales[!is.finite(scales) | scales < 1e-8] <- 1
  dist <- rep(0, nrow(candidates))
  for (nm in feats) {
    a <- transform_feature(candidates[[nm]], nm)
    b <- transform_feature(target[[nm]], nm)
    dist <- dist + ((a - b) / scales[[nm]])^2
  }
  candidates$design_distance <- sqrt(dist)
  candidates$target_start <- target_start
  candidates$target_end <- target_end
  candidates$buffer_cohort_groups <- buffer
  candidates <- candidates |> dplyr::arrange(design_distance, start)
  list(selected = candidates[1, , drop = FALSE], candidates = candidates, target = target)
}

build_stan_data <- function(dat, cfg) {
  countries <- sort(unique(as.character(dat$country)))
  I <- length(countries)
  ages <- sort(unique(as.integer(dat$Age)))
  periods <- sort(unique(as.integer(dat$Period)))
  cohorts_global <- seq(min(dat$Cohort), max(dat$Cohort), by = cfg$data$cohort_step)
  A <- length(ages); P <- length(periods); C_global <- length(cohorts_global)
  if (A < 4L || P < 5L || C_global < 5L) stop("The retained Lexis grid is too small for the model.")

  X_age_baseline <- as.matrix(splines::bs(ages, df = cfg$model$B_age, degree = 3, intercept = TRUE))
  X_age_deviation <- orthogonal_spline_basis(ages, cfg$model$B_age_deviation, remove_linear = FALSE)
  X_age_regime <- orthogonal_spline_basis(ages, cfg$model$B_regime, remove_linear = FALSE)
  X_regime_deviation <- orthogonal_spline_basis(ages, cfg$model$B_regime_deviation, remove_linear = FALSE)
  X_period <- orthogonal_spline_basis(periods, cfg$model$B_period, remove_linear = TRUE)
  X_cohort <- orthogonal_spline_basis(cohorts_global, cfg$model$B_cohort, remove_linear = TRUE)
  D2_cohort <- second_difference_matrix(C_global)
  drift_index <- as.numeric((periods - mean(periods)) / as.numeric(cfg$data$period_step))

  dat_indexed <- dat |>
    dplyr::mutate(
      country_id = match(country, countries), age_idx = match(Age, ages),
      period_idx = match(Period, periods), global_c = match(Cohort, cohorts_global)
    )
  if (anyNA(dat_indexed$country_id) || anyNA(dat_indexed$age_idx) || anyNA(dat_indexed$period_idx) || anyNA(dat_indexed$global_c)) stop("Internal indexing failure.", call. = FALSE)

  mask <- mask_window_from_cfg(cfg)
  dat_indexed$include_likelihood <- 1L
  if (!is.null(mask)) dat_indexed$include_likelihood[dat_indexed$Cohort >= mask[1] & dat_indexed$Cohort <= mask[2]] <- 0L

  active_start <- integer(I); active_end <- integer(I); local_grid <- vector("list", I); T_vec <- integer(I)
  for (i in seq_len(I)) {
    obs_c <- sort(unique(dat_indexed$global_c[dat_indexed$country_id == i]))
    active_start[i] <- min(obs_c); active_end[i] <- max(obs_c)
    local_grid[[i]] <- seq.int(active_start[i], active_end[i]); T_vec[i] <- length(local_grid[[i]])
  }
  maxT <- max(T_vec)
  D_max <- max(maxT, as.integer(cfg$model$duration_support_min_intervals), as.integer(ceiling(cfg$model$duration_support_multiplier * maxT)))

  n_obs <- matrix(0L, I, maxT); start_idx <- matrix(1L, I, maxT); end_idx <- matrix(1L, I, maxT)
  pieces <- vector("list", I); rows_by_country_cohort <- vector("list", I); cohort_map <- vector("list", I); cursor <- 0L
  for (i in seq_len(I)) {
    d_i <- dat_indexed |> dplyr::filter(country_id == i) |> dplyr::mutate(t = global_c - active_start[i] + 1L) |> dplyr::arrange(t, Period, Age)
    for (t in seq_len(T_vec[i])) {
      local_rows <- which(d_i$t == t)
      if (length(local_rows)) {
        n_obs[i, t] <- length(local_rows); start_idx[i, t] <- cursor + min(local_rows); end_idx[i, t] <- cursor + max(local_rows)
      }
    }
    pieces[[i]] <- d_i
    cohort_map[[i]] <- data.frame(country = countries[i], t = seq_len(T_vec[i]), global_c = local_grid[[i]], Cohort = cohorts_global[local_grid[[i]]])
    cursor <- cursor + nrow(d_i)
  }
  dat_long <- dplyr::bind_rows(pieces) |>
    dplyr::mutate(
      has_F = as.integer(Exp_F > 0), has_M = as.integer(Exp_M > 0),
      logE_F = dplyr::if_else(Exp_F > 0, log(Exp_F), 0), logE_M = dplyr::if_else(Exp_M > 0, log(Exp_M), 0)
    )
  N <- nrow(dat_long)
  for (i in seq_len(I)) rows_by_country_cohort[[i]] <- lapply(seq_len(T_vec[i]), function(t) which(dat_long$country_id == i & dat_long$t == t))

  stan_data <- list(
    I = as.integer(I), K = as.integer(cfg$model$K), reduce_grainsize = as.integer(cfg$runtime$stan_reduce_grainsize %||% 1L),
    A = as.integer(A), B_age = as.integer(ncol(X_age_baseline)), B_age_dev = as.integer(ncol(X_age_deviation)),
    B_regime = as.integer(ncol(X_age_regime)), B_regime_dev = as.integer(ncol(X_regime_deviation)),
    X_age_baseline = unname(X_age_baseline), X_age_deviation = unname(X_age_deviation), X_age_regime = unname(X_age_regime), X_regime_deviation = unname(X_regime_deviation),
    P = as.integer(P), B_period = as.integer(ncol(X_period)), X_period = unname(X_period), drift_index = drift_index,
    C = as.integer(C_global), B_cohort = as.integer(ncol(X_cohort)), X_cohort = unname(X_cohort), D2_cohort = unname(D2_cohort),
    maxT = as.integer(maxT), D_max = as.integer(D_max), T = as.integer(T_vec), N = as.integer(N),
    age_idx = as.integer(dat_long$age_idx), period_idx = as.integer(dat_long$period_idx), cohort_idx = as.integer(dat_long$global_c), include_likelihood = as.integer(dat_long$include_likelihood),
    has_F = as.integer(dat_long$has_F), has_M = as.integer(dat_long$has_M), D_F = as.integer(dat_long$Deaths_F), D_M = as.integer(dat_long$Deaths_M),
    logE_F = as.numeric(dat_long$logE_F), logE_M = as.numeric(dat_long$logE_M), n_obs = unname(n_obs), start_idx = unname(start_idx), end_idx = unname(end_idx),
    baseline_country_level_prior_sd = as.numeric(cfg$model$baseline_country_level_prior_sd), baseline_country_shape_prior_sd = as.numeric(cfg$model$baseline_country_shape_prior_sd), baseline_global_rw2_prior_sd = as.numeric(cfg$model$baseline_global_rw2_prior_sd),
    drift_global_prior_sd = as.numeric(cfg$model$drift_global_prior_sd), drift_country_sd_prior = as.numeric(cfg$model$drift_country_sd_prior),
    period_global_prior_sd = as.numeric(cfg$model$period_global_prior_sd), period_country_deviation_prior_sd = as.numeric(cfg$model$period_country_deviation_prior_sd),
    cohort_global_prior_sd = as.numeric(cfg$model$cohort_global_prior_sd), cohort_country_deviation_prior_sd = as.numeric(cfg$model$cohort_country_deviation_prior_sd), cohort_smoothness_sigma = as.numeric(cfg$model$cohort_smoothness_sigma),
    regime_range_prior_sd = as.numeric(cfg$model$regime_range_prior_sd), regime_global_shape_prior_sd = as.numeric(cfg$model$regime_global_shape_prior_sd), country_regime_level_sd_prior = as.numeric(cfg$model$country_regime_level_sd_prior), country_regime_shape_sd_prior = as.numeric(cfg$model$country_regime_shape_sd_prior),
    initial_occupancy_global_prior_sd = as.numeric(cfg$model$initial_occupancy_global_prior_sd), initial_occupancy_country_sd_prior = as.numeric(cfg$model$initial_occupancy_country_sd_prior),
    numerical_floor_phi = as.numeric(cfg$model$numerical_floor_phi), numerical_floor_baseline_rw2_sd = as.numeric(cfg$model$numerical_floor_baseline_rw2_sd), numerical_floor_duration_country_sd = as.numeric(cfg$model$numerical_floor_duration_country_sd), numerical_floor_duration_dispersion = as.numeric(cfg$model$numerical_floor_duration_dispersion),
    duration_prior_mean = as.numeric(cfg$model$duration_prior_mean_intervals), duration_prior_log_sd = as.numeric(cfg$model$duration_prior_log_sd), duration_country_sd_prior = as.numeric(cfg$model$duration_country_sd_prior)
  )

  leakage <- cohort_local_linear_leakage(dat_long, cohorts_global, X_cohort)
  metadata <- list(
    countries = countries, ages = ages, periods = periods, cohorts_global = cohorts_global, C_global = C_global,
    active_start = active_start, active_end = active_end, T = T_vec, cohort_map = dplyr::bind_rows(cohort_map), dat_long = dat_long,
    rows_by_country_cohort = rows_by_country_cohort, cohort_step = cfg$data$cohort_step, D_max = D_max, duration_support_years = D_max * cfg$data$cohort_step,
    mask_window = mask, cohort_local_linear_leakage = leakage
  )
  list(stan_data = stan_data, metadata = metadata)
}

write_data_artifacts <- function(pack, cfg) {
  dir.create(file.path(cfg$output_dir, "data"), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(pack$metadata$dat_long, file.path(cfg$output_dir, "data", "analysis_cells.csv"))
  readr::write_csv(pack$metadata$cohort_map, file.path(cfg$output_dir, "data", "cohort_map.csv"))
  readr::write_csv(pack$metadata$cohort_local_linear_leakage, file.path(cfg$output_dir, "data", "cohort_local_linear_leakage.csv"))
  coverage <- pack$metadata$dat_long |>
    dplyr::group_by(country) |>
    dplyr::summarise(first_period = min(Period), last_period = max(Period), first_cohort = min(Cohort), last_cohort = max(Cohort), min_age = min(Age), max_age = max(Age), n_cohort_groups = dplyr::n_distinct(Cohort), n_cells = dplyr::n(), deaths_F = sum(Deaths_F), deaths_M = sum(Deaths_M), exposure_F = sum(Exp_F), exposure_M = sum(Exp_M), .groups = "drop")
  readr::write_csv(coverage, file.path(cfg$output_dir, "data", "country_coverage.csv"))
  support_meta <- data.frame(C_global = pack$metadata$C_global, max_national_T = max(pack$metadata$T), D_max_intervals = pack$metadata$D_max, cohort_step_years = pack$metadata$cohort_step, D_max_years = pack$metadata$duration_support_years, mask_start = ifelse(is.null(pack$metadata$mask_window), NA, pack$metadata$mask_window[1]), mask_end = ifelse(is.null(pack$metadata$mask_window), NA, pack$metadata$mask_window[2]))
  readr::write_csv(support_meta, file.path(cfg$output_dir, "data", "duration_support_metadata.csv"))
  saveRDS(pack, file.path(cfg$output_dir, "data", "stan_pack.rds"))
  invisible(coverage)
}
