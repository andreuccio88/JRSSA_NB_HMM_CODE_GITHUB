rm(list = ls(all.names = TRUE)); gc()
ROOT <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
source(file.path(ROOT, "R", "project_loader.R"), local = FALSE)
source_project_files(ROOT); load_project_packages()
cfg <- read_project_config(file.path(ROOT, "config", "default.yml"), ROOT)

cp_file <- file.path(cfg$paths$analysis, "changepoint", "masked_1925_1935", "tables", "masked_predictive_density_by_country.csv")
sm_file <- file.path(cfg$paths$analysis, "smooth", "masked_1925_1935", "tables", "masked_predictive_density_by_country.csv")
if (!file.exists(cp_file) || !file.exists(sm_file)) stop("Run both masked models first.")

cp <- readr::read_csv(cp_file, show_col_types = FALSE) |> dplyr::rename(lpd_cp = log_predictive_density)
sm <- readr::read_csv(sm_file, show_col_types = FALSE) |> dplyr::rename(lpd_smooth = log_predictive_density)
by_country <- dplyr::inner_join(cp, sm, by = "country") |>
  dplyr::mutate(delta_lpd_cp_minus_smooth = lpd_cp - lpd_smooth)

outdir <- file.path(cfg$paths$analysis, "comparison")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
readr::write_csv(by_country, file.path(outdir, "masked_target_comparison_by_country.csv"))

pack <- readRDS(file.path(cfg$paths$analysis, "changepoint", "masked_1925_1935", "stan_pack.rds"))
n_masked <- sum(pack$stan_data$include_likelihood == 0L)
summary <- data.frame(
  target_start = cfg$analysis$target_mask_start,
  target_end = cfg$analysis$target_mask_end,
  masked_cells = n_masked,
  lpd_changepoint = sum(by_country$lpd_cp),
  lpd_smooth = sum(by_country$lpd_smooth),
  delta_lpd_cp_minus_smooth = sum(by_country$delta_lpd_cp_minus_smooth),
  delta_lpd_per_masked_cell = sum(by_country$delta_lpd_cp_minus_smooth) / n_masked
)
readr::write_csv(summary, file.path(outdir, "masked_target_model_comparison.csv"))
print(summary)
