# ============================================================================
# Objective study-scope + HMD support selection - v6.1 STAN JRSS A
#
# STUDY SCOPE (fixed before fitting):
#   Countries that entered the OECD before 1 January 1990.
#   This defines the long-standing industrialized-economy comparison and
#   excludes the post-socialist/post-Cold-War OECD enlargement by construction.
#   Scope is institutional and independent of mortality outcomes.
#
# HMD SUPPORT GATES (also fixed before fitting):
#   within calendar years 1930-2019 and ages 20-89,
#   1. one primary total-population HMD series per country;
#   2. >= 60 usable calendar years;
#   3. >= 12 quinquennial period groups;
#   4. >= 12 quinquennial cohort groups observed at >= 4 adult age groups;
#   5. >= 95% complete joint death/exposure cells for BOTH sexes.
#
# No fitted mortality quantity, state probability, transition cohort, ELPD,
# regime RR or agreement with earlier results enters population selection.
# ============================================================================

extract_data_frame_from_rdata <- function(path, preferred_object) {
  env <- new.env(parent=emptyenv())
  loaded <- load(path, envir=env)
  candidates <- unique(c(preferred_object, loaded))
  for (nm in candidates) {
    if (!exists(nm, envir=env, inherits=FALSE)) next
    obj <- get(nm, envir=env)
    if (is.data.frame(obj)) return(obj)
    if (is.list(obj) && "data" %in% names(obj) && is.data.frame(obj$data)) return(obj$data)
  }
  stop("No usable data frame found in ", path, call.=FALSE)
}

standardize_hmd_source <- function(x, kind=c("Dx","Ex")) {
  kind <- match.arg(kind)
  needed <- c("country","Year","Age","Female","Male")
  miss <- setdiff(needed, names(x))
  if (length(miss)) stop(kind, " input missing columns: ", paste(miss,collapse=", "), call.=FALSE)
  z <- x |>
    dplyr::transmute(
      country=trimws(as.character(country)),
      Year=suppressWarnings(as.integer(Year)),
      Age=suppressWarnings(as.integer(Age)),
      Female=as.numeric(Female),
      Male=as.numeric(Male)
    ) |>
    dplyr::filter(!is.na(country), nzchar(country), is.finite(Year), is.finite(Age))
  if (kind=="Dx") dplyr::rename(z, Deaths_F=Female, Deaths_M=Male)
  else dplyr::rename(z, Exp_F=Female, Exp_M=Male)
}

read_hmd_sources_separately <- function(cfg) {
  dx <- extract_data_frame_from_rdata(cfg$data$dx_file, "HMD_Dx") |> standardize_hmd_source("Dx")
  ex <- extract_data_frame_from_rdata(cfg$data$ex_file, "HMD_Ex") |> standardize_hmd_source("Ex")
  list(dx=dx, ex=ex)
}

read_series_metadata <- function(cfg) {
  x <- readr::read_csv(cfg$data$series_metadata, show_col_types=FALSE)
  needed <- c("hmd_code","population","primary_series","series_note")
  miss <- setdiff(needed,names(x))
  if (length(miss)) stop("series_metadata.csv missing: ",paste(miss,collapse=", "),call.=FALSE)
  x |>
    dplyr::mutate(
      hmd_code=as.character(hmd_code),
      primary_series=as.integer(primary_series)
    )
}

read_longstanding_oecd_universe <- function(cfg) {
  x <- readr::read_csv(cfg$data$oecd_universe, show_col_types=FALSE, na=c("","NA"))
  needed <- c("country","oecd_accession_date","hmd_code","scope_note")
  miss <- setdiff(needed,names(x))
  if (length(miss)) stop("longstanding_oecd_universe.csv missing: ",paste(miss,collapse=", "),call.=FALSE)
  x <- x |>
    dplyr::mutate(
      country=as.character(country),
      oecd_accession_date=as.Date(oecd_accession_date),
      hmd_code=as.character(hmd_code)
    )
  if (any(is.na(x$oecd_accession_date))) stop("OECD accession dates must be complete.",call.=FALSE)
  if (any(x$oecd_accession_date >= as.Date("1990-01-01"))) {
    stop("Study-scope file contains a country entering OECD in/after 1990.",call.=FALSE)
  }
  if (anyDuplicated(x$country)) stop("Duplicate country in OECD study-scope file.",call.=FALSE)
  x
}

country_screening_manifest <- function(cfg, write=TRUE) {
  src <- read_hmd_sources_separately(cfg)
  meta <- read_series_metadata(cfg)
  scope <- read_longstanding_oecd_universe(cfg)

  dx_codes <- unique(src$dx$country)
  ex_codes <- unique(src$ex$country)
  all_codes <- sort(unique(c(dx_codes, ex_codes)))

  unknown <- setdiff(all_codes, meta$hmd_code)
  if (length(unknown)) {
    stop(
      "The supplied RData contain HMD codes absent from config/series_metadata.csv: ",
      paste(unknown,collapse=", "),
      ". Add their population-level classification before fitting.",
      call.=FALSE
    )
  }

  scope_codes <- stats::na.omit(scope$hmd_code)
  scope_codes <- scope_codes[nzchar(scope_codes)]
  missing_scope_metadata <- setdiff(scope_codes, meta$hmd_code)
  if (length(missing_scope_metadata)) {
    stop("Study-scope HMD codes absent from series_metadata.csv: ",
         paste(missing_scope_metadata,collapse=", "),call.=FALSE)
  }

  joined <- dplyr::full_join(src$dx, src$ex, by=c("country","Year","Age")) |>
    dplyr::filter(
      Year >= as.integer(cfg$data$year_min),
      Year <= as.integer(cfg$data$year_max),
      Age >= as.integer(cfg$data$age_min),
      Age <= as.integer(cfg$data$age_max)
    ) |>
    dplyr::mutate(
      complete_F = is.finite(Deaths_F) & is.finite(Exp_F) & Deaths_F >= 0 & Exp_F > 0,
      complete_M = is.finite(Deaths_M) & is.finite(Exp_M) & Deaths_M >= 0 & Exp_M > 0,
      complete_both = complete_F & complete_M
    )

  one_country <- function(code) {
    m <- meta |> dplyr::filter(hmd_code == code)
    d <- joined |> dplyr::filter(country == code)

    in_dx <- code %in% dx_codes
    in_ex <- code %in% ex_codes
    in_scope <- code %in% scope_codes
    srow <- scope |> dplyr::filter(hmd_code == code)
    accession <- if (nrow(srow)) as.character(srow$oecd_accession_date[1]) else NA_character_

    # Compute support metrics for every primary series for transparent auditing,
    # but only in-scope series can ever be included.
    if (!nrow(d)) {
      first_joint <- NA_integer_; last_joint <- NA_integer_
      complete_fraction <- 0; usable_years <- integer(0)
      n_qperiods <- 0L; n_supported_cohorts <- 0L
    } else {
      yd <- src$dx |> dplyr::filter(country==code, Year>=cfg$data$year_min, Year<=cfg$data$year_max,
                                   Age>=cfg$data$age_min, Age<=cfg$data$age_max) |> dplyr::pull(Year)
      ye <- src$ex |> dplyr::filter(country==code, Year>=cfg$data$year_min, Year<=cfg$data$year_max,
                                   Age>=cfg$data$age_min, Age<=cfg$data$age_max) |> dplyr::pull(Year)
      if (!length(yd) || !length(ye)) {
        first_joint <- NA_integer_; last_joint <- NA_integer_
      } else {
        first_joint <- max(min(yd,na.rm=TRUE),min(ye,na.rm=TRUE),as.integer(cfg$data$year_min))
        last_joint <- min(max(yd,na.rm=TRUE),max(ye,na.rm=TRUE),as.integer(cfg$data$year_max))
      }

      if (!is.finite(first_joint) || !is.finite(last_joint) || first_joint > last_joint) {
        complete_fraction <- 0; usable_years <- integer(0)
        n_qperiods <- 0L; n_supported_cohorts <- 0L
      } else {
        dspan <- d |> dplyr::filter(Year>=first_joint,Year<=last_joint)
        expected_ages <- as.integer(cfg$data$age_max)-as.integer(cfg$data$age_min)+1L
        year_qc <- dspan |>
          dplyr::group_by(Year) |>
          dplyr::summarise(
            complete_cells=dplyr::n_distinct(Age[complete_both]),
            complete_fraction=complete_cells/expected_ages,
            .groups="drop"
          )
        usable_years <- year_qc |>
          dplyr::filter(complete_fraction >= cfg$data$min_complete_cell_fraction) |>
          dplyr::pull(Year)

        expected_n <- (last_joint-first_joint+1L) * expected_ages
        complete_n <- dspan |>
          dplyr::filter(complete_both) |>
          dplyr::distinct(Year,Age) |>
          nrow()
        complete_fraction <- if (expected_n>0) complete_n/expected_n else 0

        qdat <- dspan |>
          dplyr::filter(complete_both) |>
          dplyr::mutate(
            Age5=5L*floor(Age/5L)+2L,
            Period5=5L*floor(Year/5L)+2L,
            Cohort5=5L*floor((Period5-Age5)/5L)
          )
        n_qperiods <- dplyr::n_distinct(qdat$Period5)
        supported <- qdat |>
          dplyr::distinct(Cohort5,Age5) |>
          dplyr::count(Cohort5,name="n_age_groups") |>
          dplyr::filter(n_age_groups >= as.integer(cfg$data$min_age_groups_per_supported_cohort))
        n_supported_cohorts <- nrow(supported)
      }
    }

    reasons <- character(0)
    if (m$primary_series[1] != 1L) {
      reasons <- "overlapping/non-primary HMD series"
    } else if (!in_scope) {
      reasons <- "outside pre-1990 OECD study universe"
    } else {
      if (!in_dx) reasons <- c(reasons,"absent from deaths file")
      if (!in_ex) reasons <- c(reasons,"absent from exposures file")
      if (length(usable_years) < as.integer(cfg$data$min_calendar_years)) reasons <- c(reasons,"< minimum usable calendar years")
      if (n_qperiods < as.integer(cfg$data$min_quinquennial_periods)) reasons <- c(reasons,"< minimum quinquennial periods")
      if (n_supported_cohorts < as.integer(cfg$data$min_supported_cohorts)) reasons <- c(reasons,"< minimum supported cohort groups")
      if (complete_fraction < as.numeric(cfg$data$min_complete_cell_fraction)) reasons <- c(reasons,"insufficient joint cell completeness")
    }

    included <- as.integer(length(reasons)==0)

    data.frame(
      hmd_code=code,
      population=m$population[1],
      primary_series=m$primary_series[1],
      in_longstanding_oecd_scope=as.integer(in_scope),
      oecd_accession_date=accession,
      in_deaths=in_dx,
      in_exposures=in_ex,
      first_year=ifelse(is.finite(first_joint),as.integer(first_joint),NA_integer_),
      last_year=ifelse(is.finite(last_joint),as.integer(last_joint),NA_integer_),
      n_usable_years=length(usable_years),
      n_quinquennial_periods=as.integer(n_qperiods),
      n_supported_cohorts=as.integer(n_supported_cohorts),
      complete_cell_fraction=as.numeric(complete_fraction),
      included=included,
      exclusion_reason=if(included==1L) "included" else paste(unique(reasons),collapse="; "),
      series_note=m$series_note[1],
      stringsAsFactors=FALSE
    )
  }

  manifest <- dplyr::bind_rows(lapply(meta$hmd_code, one_country)) |>
    dplyr::arrange(dplyr::desc(included), dplyr::desc(in_longstanding_oecd_scope), population)

  if (write) {
    out <- file.path(cfg$output_dir,"data")
    dir.create(out,recursive=TRUE,showWarnings=FALSE)
    readr::write_csv(manifest,file.path(out,"country_selection_manifest.csv"))
    selected <- manifest |> dplyr::filter(included==1L) |> dplyr::pull(hmd_code)
    writeLines(selected,file.path(out,"selected_hmd_codes.txt"))

    # Full institutional study universe, including members not represented by
    # a usable HMD total-population series in the supplied files (e.g. Türkiye).
    universe_manifest <- scope |>
      dplyr::mutate(
        present_in_series_metadata = !is.na(hmd_code) & hmd_code %in% meta$hmd_code,
        present_in_deaths = !is.na(hmd_code) & hmd_code %in% dx_codes,
        present_in_exposures = !is.na(hmd_code) & hmd_code %in% ex_codes
      ) |>
      dplyr::left_join(
        manifest |>
          dplyr::select(hmd_code,included,exclusion_reason,n_usable_years,
                        n_quinquennial_periods,n_supported_cohorts,complete_cell_fraction),
        by="hmd_code"
      ) |>
      dplyr::mutate(
        included=dplyr::coalesce(included,0L),
        exclusion_reason=dplyr::case_when(
          included==1L ~ "included",
          is.na(hmd_code) | !present_in_series_metadata ~ "no primary HMD series in supplied data",
          TRUE ~ exclusion_reason
        )
      )
    readr::write_csv(universe_manifest,file.path(out,"study_universe_manifest.csv"))

    summary <- data.frame(
      longstanding_oecd_universe=nrow(scope),
      scope_countries_represented_in_metadata=sum(!is.na(scope$hmd_code) & scope$hmd_code %in% meta$hmd_code),
      selected_after_hmd_support=sum(manifest$included==1L),
      out_of_scope_primary_hmd_series=sum(manifest$primary_series==1L & manifest$in_longstanding_oecd_scope==0L),
      analysis_year_min=cfg$data$year_min,
      analysis_year_max=cfg$data$year_max,
      min_calendar_years=cfg$data$min_calendar_years,
      min_complete_cell_fraction=cfg$data$min_complete_cell_fraction
    )
    readr::write_csv(summary,file.path(out,"country_selection_summary.csv"))
  }

  manifest
}

selected_country_codes <- function(cfg, ensure_manifest=TRUE) {
  path <- file.path(cfg$output_dir,"data","country_selection_manifest.csv")
  man <- if (file.exists(path)) readr::read_csv(path,show_col_types=FALSE) else country_screening_manifest(cfg,write=ensure_manifest)
  x <- man |> dplyr::filter(included==1L) |> dplyr::pull(hmd_code) |> as.character()
  if (length(x)<2L) stop("Fewer than two populations satisfy the pre-1990 OECD + HMD support rule.",call.=FALSE)
  x
}
