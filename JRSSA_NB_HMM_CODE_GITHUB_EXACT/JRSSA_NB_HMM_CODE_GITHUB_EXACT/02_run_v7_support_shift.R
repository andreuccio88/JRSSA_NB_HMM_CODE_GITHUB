rm(list=ls(all.names=TRUE)); gc()
ROOT <- normalizePath(getwd(), winslash='/', mustWork=TRUE)
source(file.path(ROOT,'R','project_loader.R'), local=FALSE)
source_project_files(ROOT); assert_project_api(); load_project_packages()
cfg <- read_project_config(file.path(ROOT,'config','default.yml'), ROOT)
cfg <- v7_light_defaults(cfg); configure_runtime(cfg)
full <- prepare_mortality_data(cfg)
full_grid <- seq(min(full$Cohort), max(full$Cohort), by=cfg$data$cohort_step)

for (direction in c('drop_oldest','drop_newest')) {
  dat <- truncate_supported_cohort_window(full, cfg, direction=direction,
                                          n_groups=cfg$analysis_spec$support_shift_cohorts)
  pack <- build_changepoint_stan_data(dat, cfg, cohorts_global=full_grid)
  out <- file.path(cfg$paths$real_full,'V7_Light_CP',paste0('support_shift_',direction))
  fit <- fit_changepoint_v7(pack, cfg, cfg$stan_changepoint_file, out,
                            iter_warmup=500L, iter_sampling=500L)
  saveRDS(pack, file.path(out,'stan_pack.rds'))
  summarise_changepoint_v7(fit, pack, out)
}
