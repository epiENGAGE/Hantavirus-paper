
data {
  // ==========================================================
  // Published offspring-count studies
  // ==========================================================
  int<lower=1> J_old;

  int<lower=0> N_ind;
  array[N_ind] int<lower=0> offspring;
  array[N_ind] int<lower=1, upper=J_old> study_ind;

  int<lower=0> N_agg;
  array[N_agg] int<lower=1> primary_count;
  array[N_agg] int<lower=0> secondary_count;
  array[N_agg] int<lower=1, upper=J_old> study_agg;

  // ==========================================================
  // 2026 cruise symptom-onset curve
  // ==========================================================
  int<lower=1> T_curve;
  array[T_curve] int<lower=0> y_curve_desc;

  // Latent infection-time grid begins PRE days before
  // the first observed onset day.
  int<lower=1> PRE;

  // Generation-interval PMF:
  // w_gen[g] = P(GI = g days), g = 1,...,G
  int<lower=1> G;
  vector<lower=0>[G] w_gen;

  // Incubation-period PMF:
  // w_inc[d+1] = P(IP = d days), d = 0,...,D-1
  int<lower=1> D;
  vector<lower=0>[D] w_inc;

  // Distribution of the index infection date conditional
  // on the first observed onset being day 1.
  // Element p corresponds to latent day p-(PRE+1),
  // i.e. -PRE,...,-1.
  vector<lower=0>[PRE] seed_profile;
}

parameters {
  // J_old published studies plus one cruise-specific R
  real mu_log_R;
  real<lower=0> tau_log_R;
  vector[J_old + 1] z_R;

  // Shared person-level offspring dispersion for studies
  // that actually contain offspring information.
  real log_k;
}

transformed parameters {
  vector[J_old + 1] log_R;
  vector<lower=0>[J_old + 1] R;
  real<lower=0> k;

  // Work entirely on the log scale for the cruise renewal model.
  // This prevents numerical overflow when R is large.
  vector[PRE + T_curve] log_seed;
  vector[PRE + T_curve] log_inf_total;
  vector[PRE + T_curve] log_inf_desc;
  vector[T_curve] log_onset_desc;
  simplex[T_curve] onset_prob;

  log_R = mu_log_R + tau_log_R * z_R;
  R = exp(log_R);
  k = exp(log_k);

  // ----------------------------------------------------------
  // Seed distribution for the conditioned-on index infection
  // ----------------------------------------------------------
  for (t in 1:(PRE + T_curve)) {
    log_seed[t] = negative_infinity();
  }

  for (p in 1:PRE) {
    if (seed_profile[p] > 0) {
      log_seed[p] = log(seed_profile[p]);
    }
  }

  // ----------------------------------------------------------
  // Log-scale renewal recursion
  //
  // log_inf_total[t] = index seed contribution
  //                    + descendant contribution
  //
  // log_inf_desc[t] is the descendant-only contribution.
  // ----------------------------------------------------------
  for (t in 1:(PRE + T_curve)) {
    vector[G] tmp;
    int n_valid = 0;

    for (g in 1:G) {
      tmp[g] = negative_infinity();

      if (t - g >= 1 && w_gen[g] > 0) {
        tmp[g] = log_inf_total[t - g] + log(w_gen[g]);
        n_valid += 1;
      }
    }

    if (n_valid > 0) {
      log_inf_desc[t] = log_R[J_old + 1] + log_sum_exp(tmp);
    } else {
      log_inf_desc[t] = negative_infinity();
    }

    log_inf_total[t] = log_sum_exp(log_seed[t], log_inf_desc[t]);
  }

  // ----------------------------------------------------------
  // Descendant infection incidence -> symptom-onset timing
  // via the incubation-period PMF, still on log scale.
  // ----------------------------------------------------------
  for (tt in 1:T_curve) {
    vector[D] tmp_inc;
    int t_lat = PRE + tt;

    for (d1 in 1:D) {
      int d0 = d1 - 1;
      int inf_idx = t_lat - d0;

      tmp_inc[d1] = negative_infinity();

      if (inf_idx >= 1 &&
          inf_idx <= PRE + T_curve &&
          w_inc[d1] > 0) {
        tmp_inc[d1] = log_inf_desc[inf_idx] + log(w_inc[d1]);
      }
    }

    log_onset_desc[tt] = log_sum_exp(tmp_inc);
  }

  // Condition on the observed number of descendant cases.
  // The epi curve informs the *shape* of onset timing, not
  // an unconstrained epidemic-size scale.
  onset_prob = softmax(log_onset_desc);
}

model {
  // ==========================================================
  // Hierarchical priors
  // ==========================================================
  mu_log_R ~ normal(log(0.5), 1.0);
  tau_log_R ~ normal(0, 0.75);
  z_R ~ std_normal();

  log_k ~ normal(log(0.3), 1.25);

  // ==========================================================
  // Individual offspring-count likelihood
  // ==========================================================
  for (n in 1:N_ind) {
    offspring[n] ~ neg_binomial_2(
      R[study_ind[n]],
      k
    );
  }

  // ==========================================================
  // Aggregate offspring-count likelihood
  //
  // If X_i ~ NB2(R_j, k), independent,
  // sum_i X_i ~ NB2(N_j R_j, N_j k).
  // ==========================================================
  for (a in 1:N_agg) {
    secondary_count[a] ~ neg_binomial_2(
      primary_count[a] * R[study_agg[a]],
      primary_count[a] * k
    );
  }

  // ==========================================================
  // Cruise epi-curve timing likelihood
  //
  // Conditional on the total number of descendant cases.
  // This avoids the previous R -> huge / phi -> tiny pathology.
  // ==========================================================
  y_curve_desc ~ multinomial(onset_prob);
}

generated quantities {
  real<lower=0> pooled_R = exp(mu_log_R);
  real<lower=0> pooled_k = k;
  real<lower=0> R_cruise = R[J_old + 1];

  array[N_ind] int offspring_rep;
  array[N_agg] int secondary_rep;
  array[T_curve] int y_curve_rep;

  vector[N_ind] log_lik_ind;
  vector[N_agg] log_lik_agg;
  real log_lik_curve;

  int N_curve = sum(y_curve_desc);

  for (n in 1:N_ind) {
    offspring_rep[n] = neg_binomial_2_rng(
      R[study_ind[n]],
      k
    );

    log_lik_ind[n] = neg_binomial_2_lpmf(
      offspring[n] |
      R[study_ind[n]],
      k
    );
  }

  for (a in 1:N_agg) {
    real mu_a = primary_count[a] * R[study_agg[a]];
    real shape_a = primary_count[a] * k;

    secondary_rep[a] = neg_binomial_2_rng(
      mu_a,
      shape_a
    );

    log_lik_agg[a] = neg_binomial_2_lpmf(
      secondary_count[a] |
      mu_a,
      shape_a
    );
  }

  y_curve_rep = multinomial_rng(onset_prob, N_curve);
  log_lik_curve = multinomial_lpmf(y_curve_desc | onset_prob);
}
