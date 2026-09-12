fit_age_baseline <- function(dat_i, ages, X, death_col, exposure_col) {
  agg <- dat_i |>
    dplyr::group_by(Age) |>
    dplyr::summarise(
      deaths = sum(.data[[death_col]], na.rm = TRUE),
      exposure = sum(.data[[exposure_col]], na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::arrange(Age)
  y <- rep(0, length(ages))
  m <- match(ages, agg$Age)
  valid <- !is.na(m) & agg$exposure[m] > 0
  y[valid] <- log((agg$deaths[m[valid]] + .5) / (agg$exposure[m[valid]] + 1))
  fallback <- if (any(valid)) mean(y[valid]) else -6
  if (sum(valid) < ncol(X)) return(rep(fallback, ncol(X)))
  tryCatch(as.numeric(qr.solve(X[valid, , drop = FALSE], y[valid])), error = function(e) rep(fallback, ncol(X)))
}

clean_directory <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  x <- list.files(path, full.names = TRUE, all.files = TRUE, no.. = TRUE)
  if (length(x)) unlink(x, recursive = TRUE, force = TRUE)
  invisible(path)
}

normalise_diagnostic_column_names <- function(x) {
  nm <- names(x)
  for (target in c("rhat", "ess_bulk", "ess_tail")) {
    exact <- which(nm == target)
    if (length(exact)) {
      names(x)[exact[1]] <- target
      next
    }
    suffix <- which(grepl(paste0("(^|::)", target, "$"), nm))
    if (length(suffix)) names(x)[suffix[1]] <- target
  }
  x
}
