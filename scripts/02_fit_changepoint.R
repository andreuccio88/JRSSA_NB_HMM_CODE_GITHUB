rm(list = ls(all.names = TRUE)); gc()
ROOT <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
source(file.path(ROOT, "R", "project_loader.R"), local = FALSE)
source_project_files(ROOT); assert_project_api(); load_project_packages()
cfg <- read_project_config(file.path(ROOT, "config", "default.yml"), ROOT); configure_runtime(cfg)

dat <- prepare_mortality_data(cfg)
full_grid <- seq(min(dat$Cohort), max(dat$Cohort), by = cfg$data$cohort_step)
pack <- build_changepoint_data(dat, cfg, cohorts_global = full_grid)
out <- file.path(cfg$paths$analysis, "changepoint", "primary")
tune <- cfg$sampling$changepoint_primary
fit <- fit_changepoint_model(pack, cfg, cfg$stan_changepoint_file, out,
  iter_warmup = tune$iter_warmup, iter_sampling = tune$iter_sampling,
  adapt_delta = tune$adapt_delta, max_treedepth = tune$max_treedepth)
saveRDS(pack, file.path(out, "stan_pack.rds"))
summarise_changepoint_fit(fit, pack, out)
print(fit$diagnostic_summary())
print(fit$summary(variables = c("switch_center_year", "switch_sd_years", "regime_range_F", "regime_range_M", "drift_global_F", "drift_global_M", "phi_F", "phi_M")))
