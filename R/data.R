read_hmd_inputs <- function(cfg) {
  src <- read_hmd_sources_separately(cfg)
  dplyr::inner_join(src$dx, src$ex, by = c("country", "Year", "Age"))
}

longest_contiguous_block <- function(cohorts, step = 5L) {
  cohorts <- sort(unique(as.integer(cohorts)))
  if (!length(cohorts)) return(integer())
  run_id <- cumsum(c(TRUE, diff(cohorts) != as.integer(step)))
  groups <- split(cohorts, run_id)
  groups[[which.max(lengths(groups))]]
}

prepare_mortality_data <- function(cfg) {
  selected <- selected_country_codes(cfg, ensure_manifest = TRUE)
  raw <- read_hmd_inputs(cfg)
  missing <- setdiff(selected, intersect(unique(raw$country), selected))
  if (length(missing)) stop("Selected HMD codes missing after deaths/exposures join: ", paste(missing, collapse = ", "), call. = FALSE)

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

  if (!nrow(dat)) stop("No mortality cells remain after sample selection.", call. = FALSE)

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
  bad <- counts |>
    dplyr::filter(supported_cohort_groups < as.integer(cfg$data$min_supported_cohorts)) |>
    dplyr::pull(country)
  if (length(bad)) stop("Populations fail the contiguous supported-cohort safeguard: ", paste(bad, collapse = ", "), call. = FALSE)

  dat <- dplyr::semi_join(dat, kept, by = c("country", "Cohort")) |>
    dplyr::arrange(country, Cohort, Period, Age)

  expected <- as.integer(cfg$data$expected_selected_populations %||% 21L)
  if (dplyr::n_distinct(dat$country) != expected) {
    stop("The application is locked to ", expected, " selected populations; found ", dplyr::n_distinct(dat$country), ".", call. = FALSE)
  }
  attr(dat, "supported_blocks") <- kept
  dat
}

truncate_supported_cohort_window <- function(dat, cfg, direction = c("drop_oldest", "drop_newest"), n_groups = NULL) {
  direction <- match.arg(direction)
  n_groups <- as.integer(n_groups %||% cfg$analysis$support_shift_cohorts %||% 5L)
  out <- lapply(split(dat, dat$country), function(z) {
    cs <- sort(unique(z$Cohort))
    if (length(cs) <= n_groups + 3L) stop("Too few cohorts for support-shift diagnostic in ", unique(z$country), call. = FALSE)
    keep <- if (direction == "drop_oldest") cs[(n_groups + 1L):length(cs)] else cs[seq_len(length(cs) - n_groups)]
    z[z$Cohort %in% keep, , drop = FALSE]
  })
  dplyr::bind_rows(out) |> dplyr::arrange(country, Cohort, Period, Age)
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

mask_window_from_cfg <- function(cfg) {
  w <- cfg$mask_window %||% NULL
  if (is.null(w) || !isTRUE(w$enabled)) return(NULL)
  c(as.integer(w$start), as.integer(w$end))
}

# Direct builder for the smooth comparator. It reproduces the same design matrices,
# indexing, row ordering and priors as the final analysis, while avoiding
# construction of objects needed only by earlier latent-state prototypes.
build_smooth_data <- function(dat, cfg) {
  countries <- sort(unique(as.character(dat$country)))
  I <- length(countries)
  ages <- sort(unique(as.integer(dat$Age)))
  periods <- sort(unique(as.integer(dat$Period)))
  cohorts_global <- seq(min(dat$Cohort), max(dat$Cohort), by = cfg$data$cohort_step)
  A <- length(ages); P <- length(periods); C_global <- length(cohorts_global)
  if (A < 4L || P < 5L || C_global < 5L) stop("The retained Lexis grid is too small for the model.")

  X_age_baseline <- as.matrix(splines::bs(ages, df = cfg$model$B_age, degree = 3, intercept = TRUE))
  X_age_deviation <- orthogonal_spline_basis(ages, cfg$model$B_age_deviation, remove_linear = FALSE)
  X_period <- orthogonal_spline_basis(periods, cfg$model$B_period, remove_linear = TRUE)
  X_cohort <- orthogonal_spline_basis(cohorts_global, cfg$model$B_cohort, remove_linear = TRUE)
  D2_cohort <- second_difference_matrix(C_global)
  drift_index <- as.numeric((periods - mean(periods)) / as.numeric(cfg$data$period_step))

  dat_indexed <- dat |>
    dplyr::mutate(
      country_id = match(country, countries),
      age_idx = match(Age, ages),
      period_idx = match(Period, periods),
      global_c = match(Cohort, cohorts_global)
    )
  if (anyNA(dat_indexed$country_id) || anyNA(dat_indexed$age_idx) || anyNA(dat_indexed$period_idx) || anyNA(dat_indexed$global_c))
    stop("Internal smooth-model indexing failure.", call. = FALSE)

  mask <- mask_window_from_cfg(cfg)
  dat_indexed$include_likelihood <- 1L
  if (!is.null(mask)) {
    dat_indexed$include_likelihood[dat_indexed$Cohort >= mask[1] & dat_indexed$Cohort <= mask[2]] <- 0L
  }

  pieces <- vector("list", I)
  for (i in seq_len(I)) {
    d_i <- dat_indexed |>
      dplyr::filter(country_id == i) |>
      dplyr::mutate(t = global_c - min(global_c) + 1L) |>
      dplyr::arrange(t, Period, Age)
    pieces[[i]] <- d_i
  }
  dat_long <- dplyr::bind_rows(pieces) |>
    dplyr::mutate(
      has_F = as.integer(Exp_F > 0),
      has_M = as.integer(Exp_M > 0),
      logE_F = dplyr::if_else(Exp_F > 0, log(Exp_F), 0),
      logE_M = dplyr::if_else(Exp_M > 0, log(Exp_M), 0)
    )

  stan_data <- list(
    I = as.integer(I),
    reduce_grainsize = as.integer(cfg$runtime$stan_smooth_reduce_grainsize %||% 64L),
    A = as.integer(A),
    B_age = as.integer(ncol(X_age_baseline)),
    B_age_dev = as.integer(ncol(X_age_deviation)),
    X_age_baseline = unname(X_age_baseline),
    X_age_deviation = unname(X_age_deviation),
    P = as.integer(P),
    B_period = as.integer(ncol(X_period)),
    X_period = unname(X_period),
    drift_index = drift_index,
    C = as.integer(C_global),
    B_cohort = as.integer(ncol(X_cohort)),
    X_cohort = unname(X_cohort),
    D2_cohort = unname(D2_cohort),
    N = as.integer(nrow(dat_long)),
    country_idx = as.integer(dat_long$country_id),
    age_idx = as.integer(dat_long$age_idx),
    period_idx = as.integer(dat_long$period_idx),
    cohort_idx = as.integer(dat_long$global_c),
    include_likelihood = as.integer(dat_long$include_likelihood),
    has_F = as.integer(dat_long$has_F),
    has_M = as.integer(dat_long$has_M),
    D_F = as.integer(dat_long$Deaths_F),
    D_M = as.integer(dat_long$Deaths_M),
    logE_F = as.numeric(dat_long$logE_F),
    logE_M = as.numeric(dat_long$logE_M),
    baseline_country_level_prior_sd = as.numeric(cfg$model$baseline_country_level_prior_sd),
    baseline_country_shape_prior_sd = as.numeric(cfg$model$baseline_country_shape_prior_sd),
    baseline_global_rw2_prior_sd = as.numeric(cfg$model$baseline_global_rw2_prior_sd),
    drift_global_prior_sd = as.numeric(cfg$model$drift_global_prior_sd),
    drift_country_sd_prior = as.numeric(cfg$model$drift_country_sd_prior),
    period_global_prior_sd = as.numeric(cfg$model$period_global_prior_sd),
    period_country_deviation_prior_sd = as.numeric(cfg$model$period_country_deviation_prior_sd),
    cohort_global_prior_sd = as.numeric(cfg$model$cohort_global_prior_sd),
    cohort_country_deviation_prior_sd = as.numeric(cfg$model$cohort_country_deviation_prior_sd),
    cohort_smoothness_sigma = as.numeric(cfg$model$cohort_smoothness_sigma),
    numerical_floor_phi = as.numeric(cfg$model$numerical_floor_phi),
    numerical_floor_baseline_rw2_sd = as.numeric(cfg$model$numerical_floor_baseline_rw2_sd)
  )

  metadata <- list(
    countries = countries,
    ages = ages,
    periods = periods,
    cohorts_global = cohorts_global,
    dat_long = dat_long,
    cohort_step = cfg$data$cohort_step,
    mask_window = mask
  )
  list(stan_data = stan_data, metadata = metadata, model_label = "Smooth", model_K = 0L)
}
