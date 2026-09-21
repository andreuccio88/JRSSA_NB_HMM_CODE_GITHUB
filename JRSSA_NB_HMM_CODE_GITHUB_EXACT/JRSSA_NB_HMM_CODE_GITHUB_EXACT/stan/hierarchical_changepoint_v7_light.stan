functions {
  real clamp_log_mean(real x) {
    return fmin(30, fmax(-30, x));
  }

  real rw2_anchor_lpdf(vector x, real sigma) {
    int n = num_elements(x);
    real lp = normal_lpdf(x[1] | 0, 2);
    if (n >= 2) lp += normal_lpdf(x[2] | 0, 2);
    if (n > 2) {
      for (j in 3:n)
        lp += normal_lpdf(x[j] - 2 * x[j - 1] + x[j - 2] | 0, sigma);
    }
    return lp;
  }

  // Columns of log_emission are:
  //   1 = post-switch / lower-female-mortality state
  //   2 = pre-switch  / higher-female-mortality state
  // A global switch-grid index g is the first cohort assigned to state 1.
  vector changepoint_log_terms(
    int T,
    int J_switch,
    int active_start,
    matrix log_emission,
    vector log_switch_prior
  ) {
    vector[T + 1] cum_post;
    vector[T + 1] cum_pre;
    vector[J_switch] out;

    cum_post[1] = 0;
    cum_pre[1] = 0;
    for (t in 1:T) {
      cum_post[t + 1] = cum_post[t] + log_emission[t, 1];
      cum_pre[t + 1] = cum_pre[t] + log_emission[t, 2];
    }

    for (g in 1:J_switch) {
      // first local cohort in the post-switch state
      int j_local = g - active_start + 1;
      real ll;

      if (j_local <= 1) {
        // switch occurred at/before the first supported cohort
        ll = cum_post[T + 1];
      } else if (j_local > T + 1) {
        // switch occurs after the last supported cohort
        ll = cum_pre[T + 1];
      } else {
        // cohorts 1,...,j_local-1 are pre; j_local,...,T are post
        ll = cum_pre[j_local]
             + (cum_post[T + 1] - cum_post[j_local]);
      }
      out[g] = log_switch_prior[g] + ll;
    }
    return out;
  }

  real changepoint_loglik(
    int T,
    int J_switch,
    int active_start,
    matrix log_emission,
    vector log_switch_prior
  ) {
    return log_sum_exp(
      changepoint_log_terms(T, J_switch, active_start,
                            log_emission, log_switch_prior)
    );
  }

  real country_changepoint_partial_sum(
    array[] int country_slice,
    int slice_start,
    int slice_end,
    int A,
    int P,
    int J_switch,
    matrix X_age_baseline,
    matrix X_age_deviation,
    matrix X_period,
    matrix X_age_regime,
    vector drift_index,
    array[] int T,
    array[] int active_start,
    array[,] int n_obs,
    array[,] int start_idx,
    array[,] int end_idx,
    array[] int age_idx,
    array[] int period_idx,
    array[] int include_likelihood,
    array[] int has_F,
    array[] int has_M,
    array[] int D_F,
    array[] int D_M,
    vector logE_F,
    vector logE_M,
    vector baseline_age_global_F,
    vector baseline_age_global_M,
    vector baseline_country_level_F,
    vector baseline_country_level_M,
    matrix baseline_country_shape_F,
    matrix baseline_country_shape_M,
    vector drift_country_F,
    vector drift_country_M,
    matrix period_country_coef_F,
    matrix period_country_coef_M,
    vector regime_contrast_F,
    vector regime_contrast_M,
    vector country_regime_contrast_F,
    vector country_regime_contrast_M,
    real phi_F,
    real phi_M,
    vector log_switch_prior
  ) {
    real lp = 0;

    for (ii in 1:size(country_slice)) {
      int i = country_slice[ii];
      vector[A] baseline_F =
        X_age_baseline * baseline_age_global_F
        + rep_vector(baseline_country_level_F[i], A)
        + X_age_deviation * baseline_country_shape_F[i]';
      vector[A] baseline_M =
        X_age_baseline * baseline_age_global_M
        + rep_vector(baseline_country_level_M[i], A)
        + X_age_deviation * baseline_country_shape_M[i]';
      vector[P] period_F = X_period * period_country_coef_F[i]';
      vector[P] period_M = X_period * period_country_coef_M[i]';
      vector[A] contrast_F = regime_contrast_F
        + rep_vector(country_regime_contrast_F[i], A);
      vector[A] contrast_M = regime_contrast_M
        + rep_vector(country_regime_contrast_M[i], A);
      matrix[T[i], 2] log_emission = rep_matrix(0, T[i], 2);

      for (t in 1:T[i]) {
        if (n_obs[i, t] > 0) {
          for (n in start_idx[i, t]:end_idx[i, t]) {
            if (include_likelihood[n] == 1) {
              int a = age_idx[n];
              int p = period_idx[n];
              real common_F = baseline_F[a]
                + drift_country_F[i] * drift_index[p]
                + period_F[p];
              real common_M = baseline_M[a]
                + drift_country_M[i] * drift_index[p]
                + period_M[p];

              // State 1 = post-switch, centered at -1/2 contrast.
              // State 2 = pre-switch,  centered at +1/2 contrast.
              if (has_F[n] == 1) {
                log_emission[t, 1] += neg_binomial_2_log_lpmf(
                  D_F[n] |
                  clamp_log_mean(logE_F[n] + common_F - 0.5 * contrast_F[a]),
                  phi_F
                );
                log_emission[t, 2] += neg_binomial_2_log_lpmf(
                  D_F[n] |
                  clamp_log_mean(logE_F[n] + common_F + 0.5 * contrast_F[a]),
                  phi_F
                );
              }
              if (has_M[n] == 1) {
                log_emission[t, 1] += neg_binomial_2_log_lpmf(
                  D_M[n] |
                  clamp_log_mean(logE_M[n] + common_M - 0.5 * contrast_M[a]),
                  phi_M
                );
                log_emission[t, 2] += neg_binomial_2_log_lpmf(
                  D_M[n] |
                  clamp_log_mean(logE_M[n] + common_M + 0.5 * contrast_M[a]),
                  phi_M
                );
              }
            }
          }
        }
      }

      lp += changepoint_loglik(
        T[i], J_switch, active_start[i], log_emission, log_switch_prior
      );
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
  int<lower=3> B_regime;
  matrix[A, B_age] X_age_baseline;
  matrix[A, B_age_dev] X_age_deviation;
  matrix[A, B_regime] X_age_regime;

  int<lower=5> P;
  int<lower=3> B_period;
  matrix[P, B_period] X_period;
  vector[P] drift_index;

  int<lower=5> C;
  int<lower=6> J_switch; // normally C+1
  vector[J_switch] switch_grid_std;
  real switch_reference_year;
  real<lower=1> cohort_step_years;

  int<lower=3> maxT;
  array[I] int<lower=3, upper=maxT> T;
  array[I] int<lower=1, upper=C> active_start;

  int<lower=1> N;
  array[N] int<lower=1, upper=A> age_idx;
  array[N] int<lower=1, upper=P> period_idx;
  array[N] int<lower=0, upper=1> include_likelihood;
  array[N] int<lower=0, upper=1> has_F;
  array[N] int<lower=0, upper=1> has_M;
  array[N] int<lower=0> D_F;
  array[N] int<lower=0> D_M;
  vector[N] logE_F;
  vector[N] logE_M;
  array[I, maxT] int<lower=0> n_obs;
  array[I, maxT] int<lower=1, upper=N> start_idx;
  array[I, maxT] int<lower=1, upper=N> end_idx;

  real<lower=0> baseline_country_level_prior_sd;
  real<lower=0> baseline_country_shape_prior_sd;
  real<lower=0> baseline_global_rw2_prior_sd;

  real<lower=0> drift_global_prior_sd;
  real<lower=0> drift_country_deviation_prior_sd;

  real<lower=0> period_global_prior_sd;
  real<lower=0> period_country_deviation_prior_sd;

  real<lower=0> regime_range_prior_sd;
  real<lower=0> regime_global_shape_prior_sd;
  real<lower=0> country_regime_level_prior_sd;

  real<lower=0> switch_center_prior_sd_intervals;
  real<lower=0> switch_sd_prior_intervals;
  real<lower=0> switch_sd_floor_intervals;

  real<lower=0> numerical_floor_phi;
  real<lower=0> numerical_floor_baseline_rw2_sd;
}

transformed data {
  array[I] int country_ids;
  for (i in 1:I) country_ids[i] = i;
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
  vector[I] drift_country_raw_F;
  vector[I] drift_country_raw_M;

  vector[B_period] period_global_coef_F;
  vector[B_period] period_global_coef_M;
  matrix[I, B_period] period_country_raw_F;
  matrix[I, B_period] period_country_raw_M;

  // Female range anchors labels. Male range is signed.
  real<lower=0> regime_range_F;
  real regime_range_M;
  vector[B_regime] regime_shape_F;
  vector[B_regime] regime_shape_M;
  vector[I] country_regime_level_raw_F;
  vector[I] country_regime_level_raw_M;

  real<lower=numerical_floor_phi> phi_F;
  real<lower=numerical_floor_phi> phi_M;

  // Hierarchical distribution of country switch years, in 5-year-grid units.
  real switch_center_std;
  real<lower=0> switch_sd_raw;
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
  vector[A] regime_contrast_F;
  vector[A] regime_contrast_M;
  vector[I] country_regime_contrast_F;
  vector[I] country_regime_contrast_M;
  real switch_sd_intervals;
  real switch_center_year;
  real switch_sd_years;

  {
    real mF = mean(baseline_country_level_raw_F);
    real mM = mean(baseline_country_level_raw_M);
    for (i in 1:I) {
      baseline_country_level_F[i] = baseline_country_level_prior_sd
        * (baseline_country_level_raw_F[i] - mF);
      baseline_country_level_M[i] = baseline_country_level_prior_sd
        * (baseline_country_level_raw_M[i] - mM);
    }
  }

  for (b in 1:B_age_dev) {
    real mF = mean(col(baseline_country_shape_raw_F, b));
    real mM = mean(col(baseline_country_shape_raw_M, b));
    for (i in 1:I) {
      baseline_country_shape_F[i, b] = baseline_country_shape_prior_sd
        * (baseline_country_shape_raw_F[i, b] - mF);
      baseline_country_shape_M[i, b] = baseline_country_shape_prior_sd
        * (baseline_country_shape_raw_M[i, b] - mM);
    }
  }

  {
    real mF = mean(drift_country_raw_F);
    real mM = mean(drift_country_raw_M);
    for (i in 1:I) {
      drift_country_F[i] = drift_global_F
        + drift_country_deviation_prior_sd * (drift_country_raw_F[i] - mF);
      drift_country_M[i] = drift_global_M
        + drift_country_deviation_prior_sd * (drift_country_raw_M[i] - mM);
    }
  }

  for (b in 1:B_period) {
    real mF = mean(col(period_country_raw_F, b));
    real mM = mean(col(period_country_raw_M, b));
    for (i in 1:I) {
      period_country_coef_F[i, b] = period_global_coef_F[b]
        + period_country_deviation_prior_sd * (period_country_raw_F[i, b] - mF);
      period_country_coef_M[i, b] = period_global_coef_M[b]
        + period_country_deviation_prior_sd * (period_country_raw_M[i, b] - mM);
    }
  }

  regime_contrast_F = rep_vector(regime_range_F, A)
    + X_age_regime * regime_shape_F;
  regime_contrast_M = rep_vector(regime_range_M, A)
    + X_age_regime * regime_shape_M;

  {
    real mF = mean(country_regime_level_raw_F);
    real mM = mean(country_regime_level_raw_M);
    for (i in 1:I) {
      country_regime_contrast_F[i] = country_regime_level_prior_sd
        * (country_regime_level_raw_F[i] - mF);
      country_regime_contrast_M[i] = country_regime_level_prior_sd
        * (country_regime_level_raw_M[i] - mM);
    }
  }

  switch_sd_intervals = switch_sd_floor_intervals + switch_sd_raw;
  switch_center_year = switch_reference_year
    + cohort_step_years * switch_center_std;
  switch_sd_years = cohort_step_years * switch_sd_intervals;
}

model {
  vector[J_switch] log_switch_prior;

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
  drift_country_raw_F ~ std_normal();
  drift_country_raw_M ~ std_normal();

  period_global_coef_F ~ normal(0, period_global_prior_sd);
  period_global_coef_M ~ normal(0, period_global_prior_sd);
  to_vector(period_country_raw_F) ~ std_normal();
  to_vector(period_country_raw_M) ~ std_normal();

  regime_range_F ~ normal(0, regime_range_prior_sd);
  regime_range_M ~ normal(0, regime_range_prior_sd);
  regime_shape_F ~ normal(0, regime_global_shape_prior_sd);
  regime_shape_M ~ normal(0, regime_global_shape_prior_sd);
  country_regime_level_raw_F ~ std_normal();
  country_regime_level_raw_M ~ std_normal();

  phi_F ~ lognormal(1.5, 0.7);
  phi_M ~ lognormal(1.5, 0.7);

  // Very broad location prior; numerical origin is arbitrary and fixed in config.
  switch_center_std ~ normal(0, switch_center_prior_sd_intervals);
  switch_sd_raw ~ normal(0, switch_sd_prior_intervals);

  for (g in 1:J_switch) {
    log_switch_prior[g] = normal_lpdf(
      switch_grid_std[g] | switch_center_std, switch_sd_intervals
    );
  }
  log_switch_prior -= log_sum_exp(log_switch_prior);

  target += reduce_sum(
    country_changepoint_partial_sum,
    country_ids,
    reduce_grainsize,
    A,
    P,
    J_switch,
    X_age_baseline,
    X_age_deviation,
    X_period,
    X_age_regime,
    drift_index,
    T,
    active_start,
    n_obs,
    start_idx,
    end_idx,
    age_idx,
    period_idx,
    include_likelihood,
    has_F,
    has_M,
    D_F,
    D_M,
    logE_F,
    logE_M,
    baseline_age_global_F,
    baseline_age_global_M,
    baseline_country_level_F,
    baseline_country_level_M,
    baseline_country_shape_F,
    baseline_country_shape_M,
    drift_country_F,
    drift_country_M,
    period_country_coef_F,
    period_country_coef_M,
    regime_contrast_F,
    regime_contrast_M,
    country_regime_contrast_F,
    country_regime_contrast_M,
    phi_F,
    phi_M,
    log_switch_prior
  );
}

generated quantities {
  vector[J_switch] switch_prior_probability;
  matrix[I, J_switch] switch_probability;
  vector[I] switch_mean_year;
  vector[I] switch_posterior_sd_years;
  vector[I] log_marginal_lik_country;
  vector[I] log_masked_predictive_kernel_country;

  vector[J_switch] log_switch_prior;
  vector[J_switch] switch_grid_year;

  for (g in 1:J_switch) {
    switch_grid_year[g] = switch_reference_year
      + cohort_step_years * switch_grid_std[g];
    log_switch_prior[g] = normal_lpdf(
      switch_grid_std[g] | switch_center_std, switch_sd_intervals
    );
  }
  log_switch_prior -= log_sum_exp(log_switch_prior);
  switch_prior_probability = softmax(log_switch_prior);

  for (i in 1:I) {
    vector[A] baseline_F =
      X_age_baseline * baseline_age_global_F
      + rep_vector(baseline_country_level_F[i], A)
      + X_age_deviation * baseline_country_shape_F[i]';
    vector[A] baseline_M =
      X_age_baseline * baseline_age_global_M
      + rep_vector(baseline_country_level_M[i], A)
      + X_age_deviation * baseline_country_shape_M[i]';
    vector[P] period_F = X_period * period_country_coef_F[i]';
    vector[P] period_M = X_period * period_country_coef_M[i]';
    vector[A] contrast_F = regime_contrast_F
      + rep_vector(country_regime_contrast_F[i], A);
    vector[A] contrast_M = regime_contrast_M
      + rep_vector(country_regime_contrast_M[i], A);
    matrix[T[i], 2] log_emission_fit = rep_matrix(0, T[i], 2);
    matrix[T[i], 2] log_emission_all = rep_matrix(0, T[i], 2);

    for (t in 1:T[i]) {
      if (n_obs[i, t] > 0) {
        for (n in start_idx[i, t]:end_idx[i, t]) {
          int a = age_idx[n];
          int p = period_idx[n];
          real common_F = baseline_F[a]
            + drift_country_F[i] * drift_index[p]
            + period_F[p];
          real common_M = baseline_M[a]
            + drift_country_M[i] * drift_index[p]
            + period_M[p];
          real llF_post = 0;
          real llF_pre = 0;
          real llM_post = 0;
          real llM_pre = 0;

          if (has_F[n] == 1) {
            llF_post = neg_binomial_2_log_lpmf(
              D_F[n] |
              clamp_log_mean(logE_F[n] + common_F - 0.5 * contrast_F[a]),
              phi_F
            );
            llF_pre = neg_binomial_2_log_lpmf(
              D_F[n] |
              clamp_log_mean(logE_F[n] + common_F + 0.5 * contrast_F[a]),
              phi_F
            );
          }
          if (has_M[n] == 1) {
            llM_post = neg_binomial_2_log_lpmf(
              D_M[n] |
              clamp_log_mean(logE_M[n] + common_M - 0.5 * contrast_M[a]),
              phi_M
            );
            llM_pre = neg_binomial_2_log_lpmf(
              D_M[n] |
              clamp_log_mean(logE_M[n] + common_M + 0.5 * contrast_M[a]),
              phi_M
            );
          }

          log_emission_all[t, 1] += llF_post + llM_post;
          log_emission_all[t, 2] += llF_pre + llM_pre;
          if (include_likelihood[n] == 1) {
            log_emission_fit[t, 1] += llF_post + llM_post;
            log_emission_fit[t, 2] += llF_pre + llM_pre;
          }
        }
      }
    }

    {
      vector[J_switch] terms_fit = changepoint_log_terms(
        T[i], J_switch, active_start[i], log_emission_fit, log_switch_prior
      );
      vector[J_switch] terms_all = changepoint_log_terms(
        T[i], J_switch, active_start[i], log_emission_all, log_switch_prior
      );
      vector[J_switch] p = softmax(terms_fit);
      real m = dot_product(p, switch_grid_year);
      real v = dot_product(p, square(switch_grid_year - rep_vector(m, J_switch)));

      switch_probability[i] = p';
      switch_mean_year[i] = m;
      switch_posterior_sd_years[i] = sqrt(fmax(v, 0));
      log_marginal_lik_country[i] = log_sum_exp(terms_fit);
      // For a masked fit this is log p(y_mask | y_fit, theta).
      // Integrate across posterior draws in R with log-mean-exp.
      log_masked_predictive_kernel_country[i] =
        log_sum_exp(terms_all) - log_sum_exp(terms_fit);
    }
  }
}
