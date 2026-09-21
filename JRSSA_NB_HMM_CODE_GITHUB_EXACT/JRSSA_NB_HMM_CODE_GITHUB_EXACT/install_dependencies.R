required <- c(
  'posterior','yaml','dplyr','tidyr','readr','ggplot2','Rcpp','splines',
  'future','future.apply','scales','stringr','matrixStats','jsonlite','zip'
)
missing <- required[!vapply(required, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing)) install.packages(missing, repos='https://cloud.r-project.org')
if (!requireNamespace('cmdstanr', quietly=TRUE)) {
  install.packages('cmdstanr', repos=c('https://stan-dev.r-universe.dev','https://cloud.r-project.org'))
}
if (!requireNamespace('cmdstanr', quietly=TRUE)) stop('cmdstanr installation failed.')
if (is.null(cmdstanr::cmdstan_version(error_on_NA=FALSE))) {
  message('CmdStan is not installed. Installing current CmdStan...')
  cores <- parallel::detectCores(logical=FALSE)
  if (!is.finite(cores) || cores < 2) cores <- 2
  cmdstanr::install_cmdstan(cores=max(1L, cores-1L))
}
message('Dependencies ready. Next run: source("00_preflight.R")')
