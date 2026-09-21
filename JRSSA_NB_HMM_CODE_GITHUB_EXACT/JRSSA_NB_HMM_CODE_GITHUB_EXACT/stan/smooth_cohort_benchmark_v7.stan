functions {
  real clamp_log_mean(real x) { return fmin(30, fmax(-30, x)); }
  real rw2_anchor_lpdf(vector x, real sigma) {
    int n = num_elements(x);
    real lp = normal_lpdf(x[1] | 0, 2);
    if (n >= 2) lp += normal_lpdf(x[2] | 0, 2);
    if (n > 2) for (j in 3:n) lp += normal_lpdf(x[j] - 2 * x[j - 1] + x[j - 2] | 0, sigma);
    return lp;
  }
  real smooth_obs_partial_sum(
    array[] int obs_slice, int slice_start, int slice_end,
    matrix X_age_baseline, matrix X_age_deviation, matrix X_period, matrix X_cohort,
    vector drift_index,
    array[] int country_idx, array[] int age_idx, array[] int period_idx, array[] int cohort_idx, array[] int include_likelihood,
    array[] int has_F, array[] int has_M, array[] int D_F, array[] int D_M,
    vector logE_F, vector logE_M,
    vector baseline_age_global_F, vector baseline_age_global_M,
    vector baseline_country_level_F, vector baseline_country_level_M,
    matrix baseline_country_shape_F, matrix baseline_country_shape_M,
    vector drift_country_F, vector drift_country_M,
    matrix period_country_coef_F, matrix period_country_coef_M,
    matrix cohort_country_coef_F, matrix cohort_country_coef_M,
    real phi_F, real phi_M
  ) {
    real lp = 0;
    for (ii in 1:size(obs_slice)) {
      int n = obs_slice[ii];
      if (include_likelihood[n] == 1) {
        int i = country_idx[n];
        int a = age_idx[n];
        int p = period_idx[n];
        int c = cohort_idx[n];
        real baselineF = dot_product(row(X_age_baseline, a), baseline_age_global_F) + baseline_country_level_F[i] + dot_product(row(X_age_deviation, a), baseline_country_shape_F[i]');
        real baselineM = dot_product(row(X_age_baseline, a), baseline_age_global_M) + baseline_country_level_M[i] + dot_product(row(X_age_deviation, a), baseline_country_shape_M[i]');
        real periodF = dot_product(row(X_period, p), period_country_coef_F[i]');
        real periodM = dot_product(row(X_period, p), period_country_coef_M[i]');
        real cohortF = dot_product(row(X_cohort, c), cohort_country_coef_F[i]');
        real cohortM = dot_product(row(X_cohort, c), cohort_country_coef_M[i]');
        if (has_F[n] == 1) {
          real etaF = baselineF + drift_country_F[i] * drift_index[p] + periodF + cohortF;
          lp += neg_binomial_2_log_lpmf(D_F[n] | clamp_log_mean(logE_F[n] + etaF), phi_F);
        }
        if (has_M[n] == 1) {
          real etaM = baselineM + drift_country_M[i] * drift_index[p] + periodM + cohortM;
          lp += neg_binomial_2_log_lpmf(D_M[n] | clamp_log_mean(logE_M[n] + etaM), phi_M);
        }
      }
    }
    return lp;
  }
}

data {
  int<lower=2> I;
  int<lower=1> reduce_grainsize;
  int<lower=4> A;
  int<lower=4> B_age;
  int<lower=2> B_age_dev;
  matrix[A, B_age] X_age_baseline;
  matrix[A, B_age_dev] X_age_deviation;
  int<lower=5> P;
  int<lower=3> B_period;
  matrix[P, B_period] X_period;
  vector[P] drift_index;
  int<lower=5> C;
  int<lower=3> B_cohort;
  matrix[C, B_cohort] X_cohort;
  matrix[C - 2, C] D2_cohort;
  int<lower=1> N;
  array[N] int<lower=1, upper=I> country_idx;
  array[N] int<lower=1, upper=A> age_idx;
  array[N] int<lower=1, upper=P> period_idx;
  array[N] int<lower=1, upper=C> cohort_idx;
  array[N] int<lower=0, upper=1> include_likelihood;
  array[N] int<lower=0, upper=1> has_F;
  array[N] int<lower=0, upper=1> has_M;
  array[N] int<lower=0> D_F;
  array[N] int<lower=0> D_M;
  vector[N] logE_F;
  vector[N] logE_M;
  real<lower=0> baseline_country_level_prior_sd;
  real<lower=0> baseline_country_shape_prior_sd;
  real<lower=0> baseline_global_rw2_prior_sd;
  real<lower=0> drift_global_prior_sd;
  real<lower=0> drift_country_sd_prior;
  real<lower=0> period_global_prior_sd;
  real<lower=0> period_country_deviation_prior_sd;
  real<lower=0> cohort_global_prior_sd;
  real<lower=0> cohort_country_deviation_prior_sd;
  real<lower=0> cohort_smoothness_sigma;
  real<lower=0> numerical_floor_phi;
  real<lower=0> numerical_floor_baseline_rw2_sd;
}

transformed data {
  array[N] int obs_ids;
  matrix[C - 2, B_cohort] D2X_cohort = D2_cohort * X_cohort;
  for (n in 1:N) obs_ids[n] = n;
}

parameters {
  vector[B_age] baseline_age_global_F;
  vector[B_age] baseline_age_global_M;
  real<lower=numerical_floor_baseline_rw2_sd> baseline_global_rw2_sd_F;
  real<lower=numerical_floor_baseline_rw2_sd> baseline_global_rw2_sd_M;
  vector[I] baseline_country_level_raw_F;
  vector[I] baseline_country_level_raw_M;
  matrix[I, B_age_dev] baseline_country_shape_raw_F;
  matrix[I, B_age_dev] baseline_country_shape_raw_M;
  real drift_global_F;
  real drift_global_M;
  real<lower=0> drift_country_sd_F;
  real<lower=0> drift_country_sd_M;
  vector[I] drift_country_raw_F;
  vector[I] drift_country_raw_M;
  vector[B_period] period_global_coef_F;
  vector[B_period] period_global_coef_M;
  matrix[I, B_period] period_country_raw_F;
  matrix[I, B_period] period_country_raw_M;
  vector[B_cohort] cohort_global_coef_F;
  vector[B_cohort] cohort_global_coef_M;
  matrix[I, B_cohort] cohort_country_raw_F;
  matrix[I, B_cohort] cohort_country_raw_M;
  real<lower=numerical_floor_phi> phi_F;
  real<lower=numerical_floor_phi> phi_M;
}

transformed parameters {
  vector[I] baseline_country_level_F;
  vector[I] baseline_country_level_M;
  matrix[I, B_age_dev] baseline_country_shape_F;
  matrix[I, B_age_dev] baseline_country_shape_M;
  vector[I] drift_country_F;
  vector[I] drift_country_M;
  matrix[I, B_period] period_country_coef_F;
  matrix[I, B_period] period_country_coef_M;
  matrix[I, B_cohort] cohort_country_coef_F;
  matrix[I, B_cohort] cohort_country_coef_M;

  {
    real mF = mean(baseline_country_level_raw_F);
    real mM = mean(baseline_country_level_raw_M);
    for (i in 1:I) {
      baseline_country_level_F[i] = baseline_country_level_prior_sd * (baseline_country_level_raw_F[i] - mF);
      baseline_country_level_M[i] = baseline_country_level_prior_sd * (baseline_country_level_raw_M[i] - mM);
    }
  }
  for (b in 1:B_age_dev) {
    real mF = mean(col(baseline_country_shape_raw_F, b));
    real mM = mean(col(baseline_country_shape_raw_M, b));
    for (i in 1:I) {
      baseline_country_shape_F[i, b] = baseline_country_shape_prior_sd * (baseline_country_shape_raw_F[i, b] - mF);
      baseline_country_shape_M[i, b] = baseline_country_shape_prior_sd * (baseline_country_shape_raw_M[i, b] - mM);
    }
  }
  {
    real mF = mean(drift_country_raw_F);
    real mM = mean(drift_country_raw_M);
    for (i in 1:I) {
      drift_country_F[i] = drift_global_F + drift_country_sd_F * (drift_country_raw_F[i] - mF);
      drift_country_M[i] = drift_global_M + drift_country_sd_M * (drift_country_raw_M[i] - mM);
    }
  }
  for (b in 1:B_period) {
    real mF = mean(col(period_country_raw_F, b));
    real mM = mean(col(period_country_raw_M, b));
    for (i in 1:I) {
      period_country_coef_F[i, b] = period_global_coef_F[b] + period_country_deviation_prior_sd * (period_country_raw_F[i, b] - mF);
      period_country_coef_M[i, b] = period_global_coef_M[b] + period_country_deviation_prior_sd * (period_country_raw_M[i, b] - mM);
    }
  }
  for (b in 1:B_cohort) {
    real mF = mean(col(cohort_country_raw_F, b));
    real mM = mean(col(cohort_country_raw_M, b));
    for (i in 1:I) {
      cohort_country_coef_F[i, b] = cohort_global_coef_F[b] + cohort_country_deviation_prior_sd * (cohort_country_raw_F[i, b] - mF);
      cohort_country_coef_M[i, b] = cohort_global_coef_M[b] + cohort_country_deviation_prior_sd * (cohort_country_raw_M[i, b] - mM);
    }
  }
}

model {
  baseline_global_rw2_sd_F ~ normal(0, baseline_global_rw2_prior_sd);
  baseline_global_rw2_sd_M ~ normal(0, baseline_global_rw2_prior_sd);
  target += rw2_anchor_lpdf(baseline_age_global_F | baseline_global_rw2_sd_F);
  target += rw2_anchor_lpdf(baseline_age_global_M | baseline_global_rw2_sd_M);
  baseline_country_level_raw_F ~ std_normal();
  baseline_country_level_raw_M ~ std_normal();
  to_vector(baseline_country_shape_raw_F) ~ std_normal();
  to_vector(baseline_country_shape_raw_M) ~ std_normal();

  drift_global_F ~ normal(0, drift_global_prior_sd);
  drift_global_M ~ normal(0, drift_global_prior_sd);
  drift_country_sd_F ~ normal(0, drift_country_sd_prior);
  drift_country_sd_M ~ normal(0, drift_country_sd_prior);
  drift_country_raw_F ~ std_normal();
  drift_country_raw_M ~ std_normal();

  period_global_coef_F ~ normal(0, period_global_prior_sd);
  period_global_coef_M ~ normal(0, period_global_prior_sd);
  to_vector(period_country_raw_F) ~ std_normal();
  to_vector(period_country_raw_M) ~ std_normal();

  cohort_global_coef_F ~ normal(0, cohort_global_prior_sd);
  cohort_global_coef_M ~ normal(0, cohort_global_prior_sd);
  to_vector(cohort_country_raw_F) ~ std_normal();
  to_vector(cohort_country_raw_M) ~ std_normal();
  target += normal_lpdf(D2X_cohort * cohort_global_coef_F | 0, cohort_smoothness_sigma);
  target += normal_lpdf(D2X_cohort * cohort_global_coef_M | 0, cohort_smoothness_sigma);
  // The same pre-specified functional roughness scale is imposed on country deviations.
  // This keeps a common global basis (and therefore meaningful coefficient pooling) while
  // preventing national smooth terms from freely reconstructing step-like cohort changes.
  for (i in 1:I) {
    target += normal_lpdf(D2X_cohort * (cohort_country_coef_F[i]' - cohort_global_coef_F) | 0, cohort_smoothness_sigma);
    target += normal_lpdf(D2X_cohort * (cohort_country_coef_M[i]' - cohort_global_coef_M) | 0, cohort_smoothness_sigma);
  }

  phi_F ~ lognormal(1.5, 0.7);
  phi_M ~ lognormal(1.5, 0.7);

  target += reduce_sum(
    smooth_obs_partial_sum, obs_ids, reduce_grainsize,
    X_age_baseline, X_age_deviation, X_period, X_cohort, drift_index,
    country_idx, age_idx, period_idx, cohort_idx, include_likelihood,
    has_F, has_M, D_F, D_M, logE_F, logE_M,
    baseline_age_global_F, baseline_age_global_M,
    baseline_country_level_F, baseline_country_level_M,
    baseline_country_shape_F, baseline_country_shape_M,
    drift_country_F, drift_country_M,
    period_country_coef_F, period_country_coef_M,
    cohort_country_coef_F, cohort_country_coef_M, phi_F, phi_M
  );
}

generated quantities {
  vector[I] log_masked_predictive_kernel_country = rep_vector(0, I);
  for (n in 1:N) {
    if (include_likelihood[n] == 0) {
      int i = country_idx[n];
      int a = age_idx[n];
      int p = period_idx[n];
      int c = cohort_idx[n];
      real eta_F = dot_product(row(X_age_baseline, a), baseline_age_global_F)
        + baseline_country_level_F[i]
        + dot_product(row(X_age_deviation, a), baseline_country_shape_F[i]')
        + drift_country_F[i] * drift_index[p]
        + dot_product(row(X_period, p), period_country_coef_F[i]')
        + dot_product(row(X_cohort, c), cohort_country_coef_F[i]');
      real eta_M = dot_product(row(X_age_baseline, a), baseline_age_global_M)
        + baseline_country_level_M[i]
        + dot_product(row(X_age_deviation, a), baseline_country_shape_M[i]')
        + drift_country_M[i] * drift_index[p]
        + dot_product(row(X_period, p), period_country_coef_M[i]')
        + dot_product(row(X_cohort, c), cohort_country_coef_M[i]');
      if (has_F[n] == 1)
        log_masked_predictive_kernel_country[i] += neg_binomial_2_log_lpmf(
          D_F[n] | clamp_log_mean(logE_F[n] + eta_F), phi_F
        );
      if (has_M[n] == 1)
        log_masked_predictive_kernel_country[i] += neg_binomial_2_log_lpmf(
          D_M[n] | clamp_log_mean(logE_M[n] + eta_M), phi_M
        );
    }
  }
}
