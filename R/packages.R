required_packages <- c(
  "cmdstanr", "posterior", "yaml", "dplyr", "readr", "ggplot2", "splines"
)

ensure_project_packages <- function(pkgs = required_packages) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop(
      "Missing R packages: ", paste(missing, collapse = ", "), "\n",
      "Run source('install_dependencies.R') before fitting.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

load_project_packages <- function() {
  ensure_project_packages()
  suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(ggplot2)
    library(splines)
    library(posterior)
  })
  invisible(TRUE)
}

check_cmdstan_installation <- function() {
  ensure_project_packages("cmdstanr")
  path <- tryCatch(cmdstanr::cmdstan_path(), error = function(e) "")
  if (!nzchar(path) || !dir.exists(path)) {
    stop(
      "CmdStan is not installed or is not visible to CmdStanR. ",
      "Install it with cmdstanr::install_cmdstan() or set cmdstanr::set_cmdstan_path().",
      call. = FALSE
    )
  }
  invisible(path)
}
