rm(list=ls(all.names=TRUE)); gc()
ROOT <- normalizePath(getwd(), winslash='/', mustWork=TRUE)
if (!file.exists(file.path(ROOT,'config','default.yml'))) stop('Open the JRSS v7-light project root first.')
source(file.path(ROOT,'R','project_loader.R'), local=FALSE)
source_project_files(ROOT); assert_project_api(); load_project_packages()
cfg <- read_project_config(file.path(ROOT,'config','default.yml'), ROOT)
configure_runtime(cfg); make_output_directories(cfg)
check_cmdstan_installation()

manifest <- country_screening_manifest(cfg, write=TRUE)
if (sum(manifest$included == 1L) != cfg$data$expected_selected_populations)
  stop('Population screen did not return the expected 21 populations.', call.=FALSE)

dat <- prepare_mortality_data(cfg)
full_grid <- seq(min(dat$Cohort), max(dat$Cohort), by=cfg$data$cohort_step)
cp_pack <- build_changepoint_stan_data(dat, cfg, cohorts_global=full_grid)
if (cp_pack$stan_data$I != 21L) stop('v7 pack does not contain 21 populations.')
if (any(cp_pack$stan_data$include_likelihood != 1L)) stop('Unmasked preflight unexpectedly contains masked cells.')

cat('Compiling v7-light change-point model...\n')
cmdstanr::cmdstan_model(cfg$stan_changepoint_file, cpp_options=list(stan_threads=TRUE), stanc_options=list('O1'))
cat('Compiling smooth-cohort comparator...\n')
cmdstanr::cmdstan_model(cfg$stan_smooth_file, cpp_options=list(stan_threads=TRUE), stanc_options=list('O1'))

summary_lines <- c(
  paste0('version=', cfg$project$version),
  paste0('n_populations=', cp_pack$stan_data$I),
  paste0('age_groups=', cp_pack$stan_data$A),
  paste0('period_groups=', cp_pack$stan_data$P),
  paste0('global_cohort_groups=', cp_pack$stan_data$C),
  paste0('switch_grid_points=', cp_pack$stan_data$J_switch),
  'stan_compile=PASS'
)
dir.create(file.path(cfg$paths$real_full,'preflight'), recursive=TRUE, showWarnings=FALSE)
writeLines(summary_lines, file.path(cfg$paths$real_full,'preflight','preflight_summary.txt'))
cat('\nPREFLIGHT PASS. Next run 01_run_v7_primary.R\n')
