rm(list = ls(all.names = TRUE)); gc()
ROOT <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
source(file.path(ROOT, "R", "project_loader.R"), local = FALSE)
source_project_files(ROOT); assert_project_api(); load_project_packages()
cfg <- read_project_config(file.path(ROOT, "config", "default.yml"), ROOT); configure_runtime(cfg)
full <- prepare_mortality_data(cfg)
full_grid <- seq(min(full$Cohort), max(full$Cohort), by = cfg$data$cohort_step)
tune <- cfg$sampling$changepoint_sensitivity

for (direction in c("drop_oldest", "drop_newest")) {
  dat <- truncate_supported_cohort_window(full, cfg, direction = direction, n_groups = cfg$analysis$support_shift_cohorts)
  pack <- build_changepoint_data(dat, cfg, cohorts_global = full_grid)
  out <- file.path(cfg$paths$analysis, "changepoint", paste0("support_", direction))
  fit <- fit_changepoint_model(pack, cfg, cfg$stan_changepoint_file, out,
    iter_warmup = tune$iter_warmup, iter_sampling = tune$iter_sampling,
    adapt_delta = tune$adapt_delta, max_treedepth = tune$max_treedepth)
  saveRDS(pack, file.path(out, "stan_pack.rds"))
  summarise_changepoint_fit(fit, pack, out)
}
