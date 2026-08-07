
data {
  int<lower=1> J;

  int<lower=0> N_ind;
  array[N_ind] int<lower=0> offspring;
  array[N_ind] int<lower=1, upper=J> study_ind;

  int<lower=0> N_agg;
  array[N_agg] int<lower=1> primary_count;
  array[N_agg] int<lower=0> secondary_count;
  array[N_agg] int<lower=1, upper=J> study_agg;
}

parameters {
  // Population mean of log R
  real mu_log_R;

  // Between-study heterogeneity in log R
  real<lower=0> tau_log_R;

  // Standard-normal non-centred study effects
  vector[J] z_R;

  // Shared individual-level NB dispersion
  real log_k;
}

transformed parameters {
  vector[J] log_R;
  vector<lower=0>[J] R;
  real<lower=0> k;

  log_R = mu_log_R + tau_log_R * z_R;
  R = exp(log_R);
  k = exp(log_k);
}

model {
  // Priors
  mu_log_R ~ normal(log(0.5), 1.0);
  tau_log_R ~ normal(0, 0.75);
  z_R ~ std_normal();

  // Broad prior for dispersion
  log_k ~ normal(log(0.3), 1.25);

  // Individual offspring likelihood
  for (n in 1:N_ind) {
    offspring[n] ~ neg_binomial_2(
      R[study_ind[n]],
      k
    );
  }

  // Aggregate likelihood
  //
  // If each individual offspring count is NB2(R_j, k),
  // then the sum over N individuals is:
  //
  // NB2(mean = N * R_j, shape = N * k)
  for (a in 1:N_agg) {
    secondary_count[a] ~ neg_binomial_2(
      primary_count[a] * R[study_agg[a]],
      primary_count[a] * k
    );
  }
}

generated quantities {
  real<lower=0> pooled_R0;
  real<lower=0> pooled_k;

  array[N_ind] int offspring_rep;
  array[N_agg] int secondary_rep;

  pooled_R0 = exp(mu_log_R);
  pooled_k = k;

  for (n in 1:N_ind) {
    offspring_rep[n] = neg_binomial_2_rng(
      R[study_ind[n]],
      k
    );
  }

  for (a in 1:N_agg) {
    secondary_rep[a] = neg_binomial_2_rng(
      primary_count[a] * R[study_agg[a]],
      primary_count[a] * k
    );
  }
}

