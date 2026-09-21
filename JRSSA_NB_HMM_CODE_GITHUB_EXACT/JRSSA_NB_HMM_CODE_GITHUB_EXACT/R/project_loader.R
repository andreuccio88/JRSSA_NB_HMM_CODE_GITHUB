source_project_files <- function() {
  files <- c(
    'R/00_packages.R',
    'R/01_config.R',
    'R/20_country_selection.R',
    'R/02_data.R',
    'R/03_fit.R',
    'R/15_changepoint_v7_light.R',
    'R/10_smooth_benchmark.R'
  )
  for (f in files) {
    if (!file.exists(f)) stop('Missing project file: ', f)
    source(f, local = FALSE)
  }
  invisible(files)
}
