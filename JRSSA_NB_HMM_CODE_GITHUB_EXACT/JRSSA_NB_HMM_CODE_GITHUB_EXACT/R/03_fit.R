fit_age_baseline <- function(dat_i,ages,X,death_col,exposure_col) {
  agg <- dat_i |> dplyr::group_by(Age) |> dplyr::summarise(deaths=sum(.data[[death_col]],na.rm=TRUE),exposure=sum(.data[[exposure_col]],na.rm=TRUE),.groups='drop') |> dplyr::arrange(Age)
  y <- rep(0,length(ages)); m <- match(ages,agg$Age); valid <- !is.na(m) & agg$exposure[m]>0
  y[valid] <- log((agg$deaths[m[valid]]+.5)/(agg$exposure[m[valid]]+1)); fallback <- if(any(valid)) mean(y[valid]) else -6
  if(sum(valid)<ncol(X)) return(rep(fallback,ncol(X)))
  tryCatch(as.numeric(qr.solve(X[valid,,drop=FALSE],y[valid])),error=function(e) rep(fallback,ncol(X)))
}

make_independent_initialization <- function(pack,cfg) {
  sd <- pack$stan_data; I <- sd$I; K <- sd$K; meta <- pack$metadata
  betaF <- betaM <- matrix(0,I,sd$B_age)
  for(i in seq_len(I)) {
    di <- dplyr::filter(meta$dat_long,country_id==i)
    betaF[i,] <- fit_age_baseline(di,meta$ages,sd$X_age_baseline,'Deaths_F','Exp_F')
    betaM[i,] <- fit_age_baseline(di,meta$ages,sd$X_age_baseline,'Deaths_M','Exp_M')
  }
  globalF <- colMeans(betaF); globalM <- colMeans(betaM)
  function(chain_id=1L) {
    set.seed(cfg$sampling$seed+100L*chain_id)
    list(
      baseline_age_global_F=globalF+rnorm(sd$B_age,0,.01), baseline_age_global_M=globalM+rnorm(sd$B_age,0,.01),
      baseline_global_rw2_sd_F=.08, baseline_global_rw2_sd_M=.08,
      baseline_country_level_raw_F=rnorm(I,0,.03), baseline_country_level_raw_M=rnorm(I,0,.03),
      baseline_country_shape_raw_F=matrix(rnorm(I*sd$B_age_dev,0,.02),I,sd$B_age_dev),
      baseline_country_shape_raw_M=matrix(rnorm(I*sd$B_age_dev,0,.02),I,sd$B_age_dev),
      drift_global_F=-.08, drift_global_M=-.07,
      drift_country_sd_F=.025, drift_country_sd_M=.025,
      drift_country_raw_F=rnorm(I,0,.02), drift_country_raw_M=rnorm(I,0,.02),
      period_global_coef_F=rep(0,sd$B_period), period_global_coef_M=rep(0,sd$B_period),
      period_country_raw_F=matrix(0,I,sd$B_period), period_country_raw_M=matrix(0,I,sd$B_period),
      cohort_global_coef_F=rep(0,sd$B_cohort), cohort_global_coef_M=rep(0,sd$B_cohort),
      cohort_country_raw_F=matrix(0,I,sd$B_cohort), cohort_country_raw_M=matrix(0,I,sd$B_cohort),
      regime_range_F=.20, regime_range_M=.20,
      regime_gap_weights=rep(1/(K-1L),K-1L),
      global_regime_shape_free_F=matrix(0,K-1L,sd$B_regime), global_regime_shape_free_M=matrix(0,K-1L,sd$B_regime),
      country_regime_level_raw_F=matrix(0,I,K-1L), country_regime_level_raw_M=matrix(0,I,K-1L),
      country_regime_level_sd_F=.04, country_regime_level_sd_M=.04,
      country_regime_shape_raw_F=matrix(0,I*(K-1L),sd$B_regime_dev), country_regime_shape_raw_M=matrix(0,I*(K-1L),sd$B_regime_dev),
      country_regime_shape_sd_F=.025, country_regime_shape_sd_M=.025,
      phi_F=8, phi_M=8,
      initial_occupancy_global_logits=rep(0,K-1L),
      initial_occupancy_country_raw=matrix(0,I,K-1L), initial_occupancy_country_sd=.15,
      transition_offdiag=array(1/(K-1L),c(I,K,K-1L)),
      log_mu_duration_global=rep(log(cfg$model$duration_prior_mean_intervals),K),
      log_mu_duration_raw=matrix(0,I,K), duration_country_sd=.08, duration_dispersion=rep(2,K)
    )
  }
}


windows_rtools_runtime_candidates <- function() {
  if (.Platform$OS.type != 'windows') return(character())
  split_path <- function(x) {
    if (is.null(x) || !length(x) || !nzchar(x)) return(character())
    unlist(strsplit(x, .Platform$path.sep, fixed=TRUE), use.names=FALSE)
  }
  env_paths <- unique(c(
    split_path(Sys.getenv('R_RTOOLS44_PATH', unset='')),
    split_path(Sys.getenv('R_RTOOLS45_PATH', unset='')),
    split_path(Sys.getenv('R_RTOOLS43_PATH', unset=''))
  ))
  homes <- unique(Filter(nzchar, c(
    Sys.getenv('RTOOLS44_HOME', unset=''), Sys.getenv('RTOOLS45_HOME', unset=''),
    Sys.getenv('RTOOLS43_HOME', unset=''), Sys.getenv('RTOOLS_HOME', unset='')
  )))
  home_paths <- unlist(lapply(homes, function(h) c(
    file.path(h,'x86_64-w64-mingw32.static.posix','bin'),
    file.path(h,'ucrt64','bin'), file.path(h,'usr','bin')
  )), use.names=FALSE)
  compiler_paths <- unname(Sys.which(c('g++','gcc','gfortran')))
  compiler_paths <- dirname(compiler_paths[nzchar(compiler_paths)])
  common_roots <- c('C:/rtools45','C:/rtools44','C:/rtools43',
                    'C:/RBuildTools/4.5','C:/RBuildTools/4.4','C:/RBuildTools/4.3')
  common_paths <- unlist(lapply(common_roots, function(h) c(
    file.path(h,'x86_64-w64-mingw32.static.posix','bin'),
    file.path(h,'ucrt64','bin'), file.path(h,'usr','bin')
  )), use.names=FALSE)
  out <- unique(c(env_paths, compiler_paths, home_paths, common_paths))
  out <- out[nzchar(out) & dir.exists(out)]
  if (length(out)) {
    out <- vapply(out, function(x) normalizePath(x, winslash='/', mustWork=TRUE), character(1))
    out <- out[!duplicated(tolower(out))]
  }
  out
}

ensure_cmdstanr_model_methods_runtime <- function(quiet=FALSE, diagnostic_file=NULL) {
  if (.Platform$OS.type != 'windows') {
    info <- list(os=.Platform$OS.type, adjusted=FALSE, candidates=character(),
                 runtime_dirs=character(), compiler=unname(Sys.which('g++')))
    if (!is.null(diagnostic_file)) { dir.create(dirname(diagnostic_file),recursive=TRUE,showWarnings=FALSE); saveRDS(info, diagnostic_file) }
    return(invisible(info))
  }
  cand <- windows_rtools_runtime_candidates()
  runtime_names <- c('libstdc++-6.dll','libwinpthread-1.dll','libgcc_s_seh-1.dll')
  has_runtime <- vapply(cand, function(d) any(file.exists(file.path(d, runtime_names))), logical(1))
  runtime_dirs <- cand[has_runtime]
  # Prepend every detected Rtools bin directory. This is deliberate: model-method
  # wrappers are loaded by Windows LoadLibrary and may depend on GCC/OpenMP runtime
  # DLLs even when ordinary CmdStan compilation/sampling already works.
  current <- strsplit(Sys.getenv('PATH'), .Platform$path.sep, fixed=TRUE)[[1]]
  prepend <- unique(c(runtime_dirs, cand))
  merged <- c(prepend, current)
  merged <- merged[nzchar(merged)]
  merged <- merged[!duplicated(tolower(gsub('\\\\','/',merged)))]
  if (length(prepend)) Sys.setenv(PATH=paste(merged, collapse=.Platform$path.sep))
  info <- list(
    os=.Platform$OS.type, adjusted=length(prepend)>0L,
    candidates=cand, runtime_dirs=runtime_dirs,
    compiler=unname(Sys.which('g++')),
    path_head=head(strsplit(Sys.getenv('PATH'), .Platform$path.sep, fixed=TRUE)[[1]], 12L)
  )
  if (!is.null(diagnostic_file)) {
    dir.create(dirname(diagnostic_file), recursive=TRUE, showWarnings=FALSE)
    saveRDS(info, diagnostic_file)
    txt <- c(
      paste0('compiler=', info$compiler),
      paste0('adjusted=', info$adjusted),
      paste0('runtime_dirs=', paste(info$runtime_dirs, collapse=';')),
      paste0('candidates=', paste(info$candidates, collapse=';')),
      paste0('path_head=', paste(info$path_head, collapse=';'))
    )
    writeLines(txt, sub('\\.rds$','.txt',diagnostic_file))
  }
  if (!quiet) {
    if (length(runtime_dirs)) {
      message('Windows CmdStanR model-method runtime PATH prepared: ', paste(runtime_dirs, collapse='; '))
    } else if (length(cand)) {
      message('Windows Rtools bin directories added to PATH for CmdStanR model methods: ', paste(cand, collapse='; '))
    } else {
      warning('Could not detect an Rtools runtime bin directory. init_model_methods() may fail with LoadLibrary. ',
              'For Rtools44 the required PATH normally includes x86_64-w64-mingw32.static.posix/bin and usr/bin.',
              call.=FALSE)
    }
  }
  invisible(info)
}

init_cmdstanr_model_methods_safe <- function(fit, seed=1L, verbose=FALSE, context='model methods', diagnostic_file=NULL) {
  info <- ensure_cmdstanr_model_methods_runtime(quiet=TRUE, diagnostic_file=diagnostic_file)
  tryCatch({
    fit$init_model_methods(seed=seed, verbose=verbose)
    invisible(fit)
  }, error=function(e) {
    msg <- conditionMessage(e)
    if (.Platform$OS.type=='windows' && grepl('LoadLibrary|shared object|DLL', msg, ignore.case=TRUE)) {
      rt <- if(length(info$runtime_dirs)) paste(info$runtime_dirs,collapse='; ') else '<none detected>'
      stop(context, ' failed while Windows was loading the Rcpp/CmdStanR model-method DLL. ',
           'This is a toolchain runtime PATH problem, not evidence of a Stan density mismatch. ',
           'Detected Rtools runtime directories: ', rt, '. ',
           'Restart R and rerun this script from the v7-light project. If it persists, ensure the Rtools ',
           'x86_64-w64-mingw32.static.posix/bin directory is on PATH. Original error: ', msg,
           call.=FALSE)
    }
    stop(e)
  })
}

clean_directory <- function(path) {
  dir.create(path,recursive=TRUE,showWarnings=FALSE)
  x <- list.files(path,full.names=TRUE,all.files=TRUE,no..=TRUE); if(length(x)) unlink(x,recursive=TRUE,force=TRUE)
  invisible(path)
}

normalise_diagnostic_column_names <- function(x) {
  nm <- names(x)
  for(target in c('rhat','ess_bulk','ess_tail')) {
    exact <- which(nm==target); if(length(exact)) { names(x)[exact[1]] <- target; next }
    suffix <- which(grepl(paste0('(^|::)',target,'$'),nm)); if(length(suffix)) names(x)[suffix[1]] <- target
  }
  x
}

core_parameter_mask <- function(parameters) {
  grepl(paste(c(
    '^regime_range_[FM]$','^regime_level_[FM]\\[','^country_regime_level_[FM]\\[',
    '^country_regime_level_sd_[FM]$','^country_regime_shape_sd_[FM]$',
    '^drift_global_[FM]$','^drift_country_sd_[FM]$','^drift_country_[FM]\\[',
    '^mu_duration_country\\[','^duration_dispersion\\[',
    '^initial_occupancy_matrix\\[','^initial_occupancy_country_sd$',
    '^transition_offdiag\\[','^phi_[FM]$','^duration_country_sd$'
  ),collapse='|'),parameters)
}


parameter_group <- function(parameters) {
  out <- rep('nuisance_or_basis',length(parameters)); out[core_parameter_mask(parameters)] <- 'scientific_core'; out
}

summarise_fit_diagnostics <- function(fit,cfg,pack,prefix='independent') {
  q025 <- function(x) stats::quantile(x,.025); q975 <- function(x) stats::quantile(x,.975)
  diagnostics <- posterior::summarise_draws(fit$draws(),mean,median,sd,q025,q975,posterior::rhat,posterior::ess_bulk,posterior::ess_tail) |>
    as.data.frame() |> normalise_diagnostic_column_names()
  names(diagnostics)[names(diagnostics)=='variable'] <- 'parameter'
  core <- core_parameter_mask(diagnostics$parameter); diagnostics$parameter_group <- parameter_group(diagnostics$parameter); diagnostics$core_supported_parameter <- core
  qr <- cfg$quality_control
  diagnostics$report_bad_rhat <- is.finite(diagnostics$rhat)&diagnostics$rhat>qr$reporting_rhat_threshold
  diagnostics$report_low_bulk_ess <- is.finite(diagnostics$ess_bulk)&diagnostics$ess_bulk<qr$reporting_min_bulk_ess
  diagnostics$report_low_tail_ess <- is.finite(diagnostics$ess_tail)&diagnostics$ess_tail<qr$reporting_min_tail_ess
  diagnostics$hard_bad_rhat <- is.finite(diagnostics$rhat)&diagnostics$rhat>qr$hard_rhat_threshold
  diagnostics$hard_low_bulk_ess <- is.finite(diagnostics$ess_bulk)&diagnostics$ess_bulk<qr$hard_min_bulk_ess
  diagnostics$hard_low_tail_ess <- is.finite(diagnostics$ess_tail)&diagnostics$ess_tail<qr$hard_min_tail_ess
  sampler <- posterior::as_draws_matrix(fit$sampler_diagnostics())
  divergent <- if('divergent__'%in%colnames(sampler)) sum(sampler[,'divergent__']) else NA_real_
  treedepth <- if('treedepth__'%in%colnames(sampler)) sum(sampler[,'treedepth__']>=cfg$sampling$max_treedepth) else NA_real_
  cmd <- suppressWarnings(fit$diagnostic_summary()) |> as.data.frame(); cmd$chain <- seq_len(nrow(cmd))
  bfmi_name <- intersect(c('ebfmi','e_bfmi','E-BFMI','bfmi'),names(cmd)); low_bfmi <- if(length(bfmi_name)) sum(cmd[[bfmi_name[1]]]<qr$bfmi_threshold,na.rm=TRUE) else NA_integer_
  report_fail <- diagnostics$report_bad_rhat|diagnostics$report_low_bulk_ess|diagnostics$report_low_tail_ess
  hard_fail <- diagnostics$hard_bad_rhat|diagnostics$hard_low_bulk_ess|diagnostics$hard_low_tail_ess
  sampler_pass <- (is.na(divergent)||isTRUE(qr$allow_divergences)||divergent==0L) && (is.na(treedepth)||treedepth<=qr$max_treedepth_allowed) && (is.na(low_bfmi)||low_bfmi==0L)
  eligible <- sampler_pass && !any(hard_fail[core],na.rm=TRUE); reporting_pass <- sampler_pass && !any(report_fail,na.rm=TRUE)
  status <- if(!eligible) 'FAIL' else if(reporting_pass) 'PASS' else 'PASS_WITH_WARNINGS'
  qc <- data.frame(
    model=paste0('K',pack$stan_data$K),diagnostic_status=status,scientifically_eligible=eligible,
    strict_reporting_targets_passed=reporting_pass,
    full_parameters_above_reporting_rhat=sum(diagnostics$report_bad_rhat,na.rm=TRUE),
    full_parameters_below_reporting_bulk_ess=sum(diagnostics$report_low_bulk_ess,na.rm=TRUE),
    full_parameters_below_reporting_tail_ess=sum(diagnostics$report_low_tail_ess,na.rm=TRUE),
    core_parameters_above_hard_rhat=sum(diagnostics$hard_bad_rhat[core],na.rm=TRUE),
    core_parameters_below_hard_bulk_ess=sum(diagnostics$hard_low_bulk_ess[core],na.rm=TRUE),
    core_parameters_below_hard_tail_ess=sum(diagnostics$hard_low_tail_ess[core],na.rm=TRUE),
    divergent_transitions=divergent,max_treedepth_hits=treedepth,chains_below_bfmi_threshold=low_bfmi,passed=eligible
  )
  readr::write_csv(diagnostics,file.path(cfg$output_dir,'tables',paste0(prefix,'_mcmc_parameter_diagnostics.csv')))
  readr::write_csv(cmd,file.path(cfg$output_dir,'tables',paste0(prefix,'_sampler_diagnostics_by_chain.csv')))
  readr::write_csv(qc,file.path(cfg$output_dir,'tables',paste0(prefix,'_mcmc_quality_control.csv')))
  readr::write_csv(diagnostics[report_fail,,drop=FALSE],file.path(cfg$output_dir,'tables',paste0(prefix,'_parameters_warning_reporting_targets.csv')))
  readr::write_csv(diagnostics[core & hard_fail,,drop=FALSE],file.path(cfg$output_dir,'tables',paste0(prefix,'_parameters_failing_hard_core_qc.csv')))
  capture.output(print(fit$diagnostic_summary()),file=file.path(cfg$output_dir,'logs',paste0(prefix,'_cmdstan_diagnostic_summary.txt')))
  list(parameter_summary=diagnostics,sampler_by_chain=cmd,quality_control=qc,status=status,eligible=eligible,passed=eligible)
}

run_k3_pathfinder_init <- function(model, pack, cfg) {
  # v6: Pathfinder is deliberately not used as a default initializer. The v5 K3
  # approximation showed severe importance-resampling instability; deterministic
  # chain-specific initializations are safer and make K3 an optional sensitivity.
  NULL
}

fit_independent_emission_hsmm <- function(pack,cfg) {
  model <- compile_project_stan_model(cfg$stan_independent_file,cfg,force_recompile=FALSE,allow_no_range_checks=TRUE)
  final_csv <- file.path(cfg$output_dir,'fit','csv'); clean_directory(final_csv)
  init_arg <- make_independent_initialization(pack,cfg)
  if(pack$stan_data$K==3L) {
    pf <- run_k3_pathfinder_init(model,pack,cfg)
    if(!is.null(pf)) init_arg <- pf
  }
  message('Hierarchical HSMM K=',pack$stan_data$K,': ',cfg$sampling$chains,' chains; ',cfg$sampling$iter_warmup,' warmup + ',cfg$sampling$iter_sampling,' sampling; adapt_delta=',cfg$sampling$adapt_delta,'.')
  fit <- model$sample(
    data=pack$stan_data, seed=cfg$sampling$seed, chains=cfg$sampling$chains,
    parallel_chains=min(cfg$sampling$parallel_chains,cfg$sampling$chains),
    threads_per_chain=as.integer(cfg$sampling$threads_per_chain),
    iter_warmup=cfg$sampling$iter_warmup, iter_sampling=cfg$sampling$iter_sampling,
    adapt_delta=cfg$sampling$adapt_delta, step_size=cfg$sampling$initial_step_size,
    max_treedepth=cfg$sampling$max_treedepth, refresh=cfg$sampling$refresh,
    init=init_arg, output_basename=paste0('k',pack$stan_data$K), save_metric=TRUE
  )
  if(any(fit$return_codes()!=0L)) stop('One or more Stan chains failed for K=',pack$stan_data$K,'.',call.=FALSE)
  if(isTRUE(cfg$sampling$save_cmdstan_csv)) {
    dir.create(final_csv,recursive=TRUE,showWarnings=FALSE)
    fit$save_output_files(dir=final_csv,basename=paste0('k',pack$stan_data$K),timestamp=FALSE,random=FALSE)
  }
  fit$save_object(file.path(cfg$output_dir,'fit','fit_cmdstanr.rds'))
  fit
}
