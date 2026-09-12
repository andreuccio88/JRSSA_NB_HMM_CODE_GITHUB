rm(list = ls(all.names = TRUE)); gc()
ROOT <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
source(file.path(ROOT, "R", "project_loader.R"), local = FALSE)
source_project_files(ROOT); assert_project_api(); load_project_packages()
cfg <- read_project_config(file.path(ROOT, "config", "default.yml"), ROOT)
configure_runtime(cfg); make_output_directories(cfg); check_cmdstan_installation()

if (!file.exists(cfg$data$dx_file) || !file.exists(cfg$data$ex_file)) {
  stop("HMD input objects are missing from data/.", call. = FALSE)
}

manifest <- country_screening_manifest(cfg, write = TRUE)
if (sum(manifest$included == 1L) != cfg$data$expected_selected_populations)
  stop("Population screen did not return the expected 21 populations.", call. = FALSE)

dat <- prepare_mortality_data(cfg)
full_grid <- seq(min(dat$Cohort), max(dat$Cohort), by = cfg$data$cohort_step)
cp_pack <- build_changepoint_data(dat, cfg, cohorts_global = full_grid)
sm_pack <- build_smooth_data(dat, cfg)
if (cp_pack$stan_data$I != 21L || sm_pack$stan_data$I != 21L) stop("Prepared model data do not contain 21 populations.")
if (any(cp_pack$stan_data$include_likelihood != 1L) || any(sm_pack$stan_data$include_likelihood != 1L)) stop("Unmasked preflight unexpectedly contains masked cells.")

cat("Compiling change-point model...\n")
compile_project_stan_model(cfg$stan_changepoint_file, cfg)
cat("Compiling smooth-cohort comparator...\n")
compile_project_stan_model(cfg$stan_smooth_file, cfg)

out <- file.path(cfg$paths$analysis, "preflight")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
writeLines(c(
  paste0("version=", cfg$project$version),
  paste0("n_populations=", cp_pack$stan_data$I),
  paste0("age_groups=", cp_pack$stan_data$A),
  paste0("period_groups=", cp_pack$stan_data$P),
  paste0("global_cohort_groups=", cp_pack$stan_data$C),
  paste0("switch_grid_points=", cp_pack$stan_data$J_switch),
  "stan_compile=PASS"
), file.path(out, "preflight_summary.txt"))
cat("\nPREFLIGHT PASS.\n")
