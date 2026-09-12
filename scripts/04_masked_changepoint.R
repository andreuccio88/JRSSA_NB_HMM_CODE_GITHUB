rm(list = ls(all.names = TRUE)); gc()
ROOT <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
source(file.path(ROOT, "R", "project_loader.R"), local = FALSE)
source_project_files(ROOT); assert_project_api(); load_project_packages()
cfg <- read_project_config(file.path(ROOT, "config", "default.yml"), ROOT); configure_runtime(cfg)
full <- prepare_mortality_data(cfg)
full_grid <- seq(min(full$Cohort), max(full$Cohort), by = cfg$data$cohort_step)
cfg$mask_window <- list(enabled = TRUE, start = as.integer(cfg$analysis$target_mask_start), end = as.integer(cfg$analysis$target_mask_end))
pack <- build_changepoint_data(full, cfg, cohorts_global = full_grid)
out <- file.path(cfg$paths$analysis, "changepoint", "masked_1925_1935")
tune <- cfg$sampling$changepoint_sensitivity
fit <- fit_changepoint_model(pack, cfg, cfg$stan_changepoint_file, out,
  iter_warmup = tune$iter_warmup, iter_sampling = tune$iter_sampling,
  adapt_delta = tune$adapt_delta, max_treedepth = tune$max_treedepth)
saveRDS(pack, file.path(out, "stan_pack.rds"))
summarise_changepoint_fit(fit, pack, out)
