functions {
  real log_vload(real t, real log_r, real log_d, real m, real log_kappa) {
    real r = exp(log_r);
    real d = exp(log_d);
    real kappa = exp(log_kappa);

    // Paper-ready pathogen-load trajectory:
    // V(t) proportional to exp(r*t) * {1 + (r/d) * exp[kappa*(t-m)]}^(-(r+d)/kappa).
    // With this parameterisation, m is the exact peak time.
    return r * t
           - ((r + d) / kappa) * log1p_exp(log(r / d) + kappa * (t - m));
  }

  real log_mech_hazard(real t,
                       real log_lambda,
                       real log_r,
                       real log_d,
                       real m,
                       real log_kappa,
                       int use_gate,
                       real gate_onset,
                       real gate_slope) {
    real out = log_lambda + log_vload(t, log_r, log_d, m, log_kappa);
    if (use_gate == 1) {
      out += log_inv_logit((t - gate_onset) / gate_slope);
    }
    return out;
  }

  // Build cumulative hazard once per leapfrog evaluation, not once per case.
  // H[1] is H(0); H[k+1] is H(k * dt_int).
  vector build_cumhaz_grid(int K_int,
                           real dt_int,
                           real log_lambda,
                           real log_r,
                           real log_d,
                           real m,
                           real log_kappa,
                           int use_gate,
                           real gate_onset,
                           real gate_slope) {
    vector[K_int + 1] H;
    H[1] = 0;
    for (k in 1:K_int) {
      real mid = (k - 0.5) * dt_int;
      H[k + 1] = H[k] + exp(log_mech_hazard(mid, log_lambda, log_r, log_d, m, log_kappa,
                                             use_gate, gate_onset, gate_slope)) * dt_int;
    }
    return H;
  }

  // Linear interpolation of cumulative hazard from precomputed grid.
  // This avoids recomputing K_int exponentials for every case. The lookup loop is
  // only comparisons/arithmetic, so it is much cheaper than the old cumhaz code.
  real cumhaz_lookup(real x, int K_int, real dt_int, vector H) {
    if (x <= 0) return 0;
    if (x >= K_int * dt_int) return H[K_int + 1];

    for (k in 1:K_int) {
      real right = k * dt_int;
      if (x <= right) {
        real left = (k - 1) * dt_int;
        real frac = (x - left) / dt_int;
        return H[k] + frac * (H[k + 1] - H[k]);
      }
    }
    return H[K_int + 1];
  }

  real mech_first_event_logdens_fast(real x,
                                     int K_int,
                                     real dt_int,
                                     real t_max,
                                     vector H,
                                     real log_lambda,
                                     real log_r,
                                     real log_d,
                                     real m,
                                     real log_kappa,
                                     int use_gate,
                                     real gate_onset,
                                     real gate_slope) {
    real H_x;
    real H_max;
    if (x <= 0 || x >= t_max) return negative_infinity();
    H_x = cumhaz_lookup(x, K_int, dt_int, H);
    H_max = H[K_int + 1];
    if (H_max <= 1e-10) return negative_infinity();
    return log_mech_hazard(x, log_lambda, log_r, log_d, m, log_kappa,
                           use_gate, gate_onset, gate_slope)
           - H_x - log1m_exp(-H_max);
  }

  real log_R_at(real t, int B, array[] real knot_day, vector log_R) {
    if (t <= knot_day[1]) return log_R[1];
    if (t >= knot_day[B]) return log_R[B];
    for (b in 1:(B - 1)) {
      if (t >= knot_day[b] && t < knot_day[b + 1]) {
        real w = (t - knot_day[b]) / (knot_day[b + 1] - knot_day[b]);
        return (1 - w) * log_R[b] + w * log_R[b + 1];
      }
    }
    return log_R[B];
  }

  vector mech_grid_pmf_fast(int K_grid,
                            array[] real grid_day,
                            real dt_grid,
                            int K_int,
                            real dt_int,
                            real t_max,
                            vector H,
                            real log_lambda,
                            real log_r,
                            real log_d,
                            real m,
                            real log_kappa,
                            int use_gate,
                            real gate_onset,
                            real gate_slope) {
    vector[K_grid] p;
    for (k in 1:K_grid) {
      real mid = grid_day[k] + 0.5 * dt_grid;
      p[k] = exp(mech_first_event_logdens_fast(mid, K_int, dt_int, t_max, H,
                                                log_lambda, log_r, log_d, m, log_kappa,
                                                use_gate, gate_onset, gate_slope)) * dt_grid;
    }
    return p / sum(p);
  }

  real mech_grid_mean_fast(int K_grid,
                           array[] real grid_day,
                           real dt_grid,
                           int K_int,
                           real dt_int,
                           real t_max,
                           vector H,
                           real log_lambda,
                           real log_r,
                           real log_d,
                           real m,
                           real log_kappa,
                           int use_gate,
                           real gate_onset,
                           real gate_slope) {
    vector[K_grid] p;
    real denom;
    real out = 0;
    for (k in 1:K_grid) {
      real mid = grid_day[k] + 0.5 * dt_grid;
      p[k] = exp(mech_first_event_logdens_fast(mid, K_int, dt_int, t_max, H,
                                                log_lambda, log_r, log_d, m, log_kappa,
                                                use_gate, gate_onset, gate_slope)) * dt_grid;
    }
    denom = sum(p);
    for (k in 1:K_grid) {
      out += (grid_day[k] + 0.5 * dt_grid) * p[k] / denom;
    }
    return out;
  }

  real control_log_weight(real delta,
                          real iso_delay_days,
                          real residual_contact_prob) {
    if (delta <= iso_delay_days) {
      return 0;
    } else {
      return log(residual_contact_prob);
    }
  }

  real mech_first_event_cdf_fast(real x,
                                 int K_int,
                                 real dt_int,
                                 real t_max,
                                 vector H) {
    real H_x;
    real H_max = H[K_int + 1];

    if (x <= 0) return 0;
    if (x >= t_max) return 1;
    if (H_max <= 1e-10) return 0;

    H_x = cumhaz_lookup(x, K_int, dt_int, H);
    return (1 - exp(-H_x)) / (1 - exp(-H_max));
  }



  // CDF for a normalized intensity shape over finite support.
  // Used for the viral-load shape V(t) when normalizing the post-infectiousness
  // transmission delay U conditional on an observed transmission event.
  real normalized_intensity_cdf_fast(real x,
                                     int K_int,
                                     real dt_int,
                                     real t_max,
                                     vector H_shape) {
    real H_x;
    real H_max = H_shape[K_int + 1];

    if (x <= 0) return 0;
    if (x >= t_max) return 1;
    if (H_max <= 1e-10) return 0;

    H_x = cumhaz_lookup(x, K_int, dt_int, H_shape);
    return H_x / H_max;
  }

  // Observed-pair generation-time density under the paper event-time model:
  //   LP is a first-event infectiousness-onset time with hazard beta_I V(t)
  //   U | LP = l, obs is proportional to V(l + u), normalized over u >= 0
  //   GT = LP + U
  // This implements f_GT(g) = int_0^g f_LP(l) f_U(g-l | l, obs) dl.
  real observed_pair_gen_logdens_fast(real g,
                                      int K_int,
                                      real dt_int,
                                      real t_max,
                                      vector H_lp,
                                      vector H_shape,
                                      real log_lambda_inf,
                                      real log_r,
                                      real log_d,
                                      real m,
                                      real log_kappa) {
    real acc = 0;
    real log_V_g;

    if (g <= 0 || g >= t_max) return negative_infinity();

    log_V_g = log_vload(g, log_r, log_d, m, log_kappa);

    for (k in 1:K_int) {
      real l = (k - 0.5) * dt_int;
      if (l < g) {
        real tail_shape = H_shape[K_int + 1] - cumhaz_lookup(l, K_int, dt_int, H_shape);
        if (tail_shape > 1e-12) {
          real lp_l = mech_first_event_logdens_fast(l, K_int, dt_int, t_max, H_lp,
                                                    log_lambda_inf, log_r, log_d, m, log_kappa,
                                                    0, 0, 1);
          acc += exp(lp_l + log_V_g - log(tail_shape)) * dt_int;
        }
      }
    }

    if (acc <= 1e-300) return log(1e-300);
    return log(acc);
  }

  vector observed_pair_gen_grid_pmf_fast(int K_grid,
                                         array[] real grid_day,
                                         real dt_grid,
                                         int K_int,
                                         real dt_int,
                                         real t_max,
                                         vector H_lp,
                                         vector H_shape,
                                         real log_lambda_inf,
                                         real log_r,
                                         real log_d,
                                         real m,
                                         real log_kappa) {
    vector[K_grid] p;
    real z;

    for (k in 1:K_grid) {
      real mid = grid_day[k] + 0.5 * dt_grid;
      p[k] = exp(observed_pair_gen_logdens_fast(mid, K_int, dt_int, t_max,
                                                H_lp, H_shape, log_lambda_inf,
                                                log_r, log_d, m, log_kappa)) * dt_grid;
    }

    z = sum(p);
    if (z <= 1e-12) return rep_vector(1.0 / K_grid, K_grid);
    return p / z;
  }

  real observed_pair_gen_grid_mean_fast(int K_grid,
                                        array[] real grid_day,
                                        real dt_grid,
                                        int K_int,
                                        real dt_int,
                                        real t_max,
                                        vector H_lp,
                                        vector H_shape,
                                        real log_lambda_inf,
                                        real log_r,
                                        real log_d,
                                        real m,
                                        real log_kappa) {
    vector[K_grid] p = observed_pair_gen_grid_pmf_fast(K_grid, grid_day, dt_grid,
                                                       K_int, dt_int, t_max,
                                                       H_lp, H_shape, log_lambda_inf,
                                                       log_r, log_d, m, log_kappa);
    real out = 0;
    for (k in 1:K_grid) {
      out += (grid_day[k] + 0.5 * dt_grid) * p[k];
    }
    return out;
  }


  // Much faster gridded implementation of the paper event-time generation-time PMF.
  // It avoids nested calls to cumhaz_lookup inside the likelihood.
  vector observed_pair_gen_grid_pmf_direct(int K_grid,
                                           array[] real grid_day,
                                           real dt_grid,
                                           int K_int,
                                           real dt_int,
                                           real t_max,
                                           vector H_lp,
                                           vector H_shape,
                                           real log_lambda_inf,
                                           real log_r,
                                           real log_d,
                                           real m,
                                           real log_kappa) {
    vector[K_grid] p;
    vector[K_int] log_f_lp;
    vector[K_int] log_tail_shape;
    real H_lp_max = H_lp[K_int + 1];
    real H_shape_max = H_shape[K_int + 1];

    for (j in 1:K_int) {
      real l = (j - 0.5) * dt_int;
      real H_lp_mid = 0.5 * (H_lp[j] + H_lp[j + 1]);
      real H_shape_mid = 0.5 * (H_shape[j] + H_shape[j + 1]);
      real tail = H_shape_max - H_shape_mid;

      if (H_lp_max <= 1e-10 || tail <= 1e-12) {
        // Set log_f_lp to -inf so this grid point contributes zero.
        // Keep log_tail_shape finite to avoid NaN from -inf - -inf below.
        log_f_lp[j] = negative_infinity();
        log_tail_shape[j] = 0;
      } else {
        log_f_lp[j] = log_lambda_inf + log_vload(l, log_r, log_d, m, log_kappa)
                      - H_lp_mid - log1m_exp(-H_lp_max);
        log_tail_shape[j] = log(tail);
      }
    }

    for (k in 1:K_grid) {
      real g = grid_day[k] + 0.5 * dt_grid;
      real log_V_g;
      real acc = 0;

      if (g <= 0 || g >= t_max) {
        p[k] = 0;
      } else {
        log_V_g = log_vload(g, log_r, log_d, m, log_kappa);
        for (j in 1:K_int) {
          real l = (j - 0.5) * dt_int;
          if (l < g) {
            acc += exp(log_f_lp[j] + log_V_g - log_tail_shape[j]) * dt_int;
          }
        }
        p[k] = acc * dt_grid;
      }
    }

    if (sum(p) <= 1e-12) return rep_vector(1.0 / K_grid, K_grid);
    return p / sum(p);
  }

  real log_grid_density_at(real x,
                           int K_grid,
                           array[] real grid_day,
                           real dt_grid,
                           vector p_mass) {
    if (x <= grid_day[1] || x >= grid_day[K_grid] + dt_grid) return negative_infinity();

    for (k in 1:K_grid) {
      if (x >= grid_day[k] && x < grid_day[k] + dt_grid) {
        return log(fmax(p_mass[k], 1e-300)) - log(dt_grid);
      }
    }
    return negative_infinity();
  }

  real control_log_normalizer_from_grid(real source_incubation,
                                        real iso_delay_days,
                                        real residual_contact_prob,
                                        int K_grid,
                                        array[] real grid_day,
                                        real dt_grid,
                                        vector p_gen) {
    real cutoff = source_incubation + iso_delay_days;
    real F = 0;

    if (cutoff <= 0) {
      F = 0;
    } else {
      for (k in 1:K_grid) {
        real mid = grid_day[k] + 0.5 * dt_grid;
        if (mid <= cutoff) F += p_gen[k];
      }
      if (F > 1) F = 1;
    }

    return log(residual_contact_prob + (1 - residual_contact_prob) * F);
  }


}

data {
  int<lower=1> N;
  array[N] int<lower=0, upper=N> source_idx;  // 0 means index / exogenous infection
  array[N] int<lower=0> Zobs;

  array[N] real onset_lo_day;
  array[N] real onset_hi_day;
  array[N] real exp_lo_day;
  array[N] real exp_hi_day;
  real<lower=0> index_pre_onset_max;

  int<lower=2> B;
  array[B] real knot_day;

  real<lower=0> t_max;
  int<lower=1> K_int;
  real<lower=0> dt_int;

  int<lower=1> K_grid;
  array[K_grid] real grid_day;
  real<lower=0> dt_grid;

  int<lower=0, upper=1> prior_only;
  // Retained for compatibility with old R data; ignored in the exact paper model.
  int<lower=0, upper=1> use_transmission_gate;

  // Optional prior on the mechanistic population mean delta:
  //   delta_mean_mech = mean(generation interval) - mean(source incubation).
  // Negative values favor presymptomatic realized transmission on average.
  int<lower=0, upper=1> use_delta_mean_prior;
  real prior_delta_mean;
  real<lower=0> prior_delta_sd;

  // Optional selection/truncation correction for post-control observed pairs.
  // For cases whose source was post-control, observed transmission times are
  // modelled as natural generation times weighted by residual contact after isolation.
  int<lower=0, upper=1> use_control_selection;
  array[N] int<lower=0, upper=1> post_control_source;
  real<lower=0> iso_delay_days;
  real<lower=0, upper=1> residual_contact_prob;

  real prior_log_lambda_sym_mean;
  real<lower=0> prior_log_lambda_sym_sd;
  real prior_log_lambda_inf_mean;
  real<lower=0> prior_log_lambda_inf_sd;
  real prior_log_r_mean;
  real<lower=0> prior_log_r_sd;
  real prior_log_d_mean;
  real<lower=0> prior_log_d_sd;
  real prior_m_mean;
  real<lower=0> prior_m_sd;
  real prior_log_kappa_mean;
  real<lower=0> prior_log_kappa_sd;

  real prior_gate_onset_mean;
  real<lower=0> prior_gate_onset_sd;
  real prior_log_gate_slope_mean;
  real<lower=0> prior_log_gate_slope_sd;

  real prior_log_R_init_mean;
  real<lower=0> prior_log_R_init_sd;
  real<lower=0> prior_sigma_rw_sd;
  real<lower=0> prior_phi_inv_sqrt_sd;
}

parameters {
  vector<lower=0, upper=1>[N] onset_raw;
  vector<lower=0, upper=1>[N] inf_raw;

  real log_lambda_sym;
  real log_lambda_inf;
  real log_r;
  real log_d;
  real m;
  real log_kappa;

  real gate_onset;
  real log_gate_slope;

  real log_R_init;
  vector[B - 1] eps_R;
  real<lower=0> sigma_rw;
  real<lower=0> phi_inv_sqrt;
}

transformed parameters {
  vector[N] T_onset;
  vector[N] T_inf;
  vector[B] log_R;
  vector[K_int + 1] H_sym;
  vector[K_int + 1] H_inf;
  vector[K_int + 1] H_shape;
  vector[K_grid] p_gen_mech;
  real inc_mean_mech;
  real gen_mean_mech;
  real delta_mean_mech;
  real<lower=0> gate_slope = exp(log_gate_slope);
  real<lower=0> k_nb = inv_square(phi_inv_sqrt);

  log_R[1] = log_R_init;
  for (b in 2:B) {
    log_R[b] = log_R[b - 1] + sigma_rw * eps_R[b - 1];
  }

  H_sym = build_cumhaz_grid(K_int, dt_int,
                              log_lambda_sym, log_r, log_d, m, log_kappa,
                              0, 0, 1);
  // Infectiousness onset uses the paper first-event hazard beta_I V(t).
  H_inf = build_cumhaz_grid(K_int, dt_int,
                            log_lambda_inf, log_r, log_d, m, log_kappa,
                            0, 0, 1);

  // Shape integral of V(t), without an event scale. This is used to normalize
  // U | LP = l, obs proportional to V(l + u); the transmission scale cancels.
  H_shape = build_cumhaz_grid(K_int, dt_int,
                              0, log_r, log_d, m, log_kappa,
                              0, 0, 1);

  p_gen_mech = observed_pair_gen_grid_pmf_direct(K_grid, grid_day, dt_grid,
                                                 K_int, dt_int, t_max,
                                                 H_inf, H_shape, log_lambda_inf,
                                                 log_r, log_d, m, log_kappa);

  inc_mean_mech = mech_grid_mean_fast(K_grid, grid_day, dt_grid, K_int, dt_int, t_max, H_sym,
                                      log_lambda_sym, log_r, log_d, m, log_kappa,
                                      0, 0, 1);
  gen_mean_mech = 0;
  for (k in 1:K_grid) {
    gen_mean_mech += (grid_day[k] + 0.5 * dt_grid) * p_gen_mech[k];
  }
  delta_mean_mech = gen_mean_mech - inc_mean_mech;

  for (i in 1:N) {
    T_onset[i] = onset_lo_day[i] + onset_raw[i] * (onset_hi_day[i] - onset_lo_day[i]);

    if (source_idx[i] == 0) {
      real lo = T_onset[i] - index_pre_onset_max;
      real hi = T_onset[i] - 1e-6;
      T_inf[i] = lo + inf_raw[i] * (hi - lo);
    } else {
      T_inf[i] = exp_lo_day[i] + inf_raw[i] * (exp_hi_day[i] - exp_lo_day[i]);
    }
  }
}

model {
  log_lambda_sym ~ normal(prior_log_lambda_sym_mean, prior_log_lambda_sym_sd);
  log_lambda_inf ~ normal(prior_log_lambda_inf_mean, prior_log_lambda_inf_sd);
  log_r ~ normal(prior_log_r_mean, prior_log_r_sd);
  log_d ~ normal(prior_log_d_mean, prior_log_d_sd);
  m ~ normal(prior_m_mean, prior_m_sd);
  log_kappa ~ normal(prior_log_kappa_mean, prior_log_kappa_sd);

  gate_onset ~ normal(prior_gate_onset_mean, prior_gate_onset_sd);
  log_gate_slope ~ normal(prior_log_gate_slope_mean, prior_log_gate_slope_sd);

  log_R_init ~ normal(prior_log_R_init_mean, prior_log_R_init_sd);
  eps_R ~ std_normal();
  sigma_rw ~ normal(0, prior_sigma_rw_sd);
  phi_inv_sqrt ~ normal(0, prior_phi_inv_sqrt_sd);

  if (use_delta_mean_prior == 1) {
    target += normal_lpdf(delta_mean_mech | prior_delta_mean, prior_delta_sd);
  }

  if (prior_only == 0) {
    for (i in 1:N) {
      real P_i = T_onset[i] - T_inf[i];
      real log_R_i = log_R_at(T_onset[i], B, knot_day, log_R);
      real R_i = exp(log_R_i);

      target += mech_first_event_logdens_fast(P_i, K_int, dt_int, t_max, H_sym,
                                           log_lambda_sym, log_r, log_d, m, log_kappa,
                                           0, 0, 1);

      if (source_idx[i] != 0) {
        int src = source_idx[i];
        real G_i = T_inf[i] - T_inf[src];
        real delta_i = T_inf[i] - T_onset[src];
        real P_src = T_onset[src] - T_inf[src];
        if (G_i <= 0) {
          target += negative_infinity();
        } else {
          real lp_G = log_grid_density_at(G_i, K_grid, grid_day, dt_grid, p_gen_mech);
          if (use_control_selection == 1 && post_control_source[i] == 1) {
            real lp_sel = control_log_weight(delta_i, iso_delay_days, residual_contact_prob);
            real log_Z_sel = control_log_normalizer_from_grid(P_src, iso_delay_days, residual_contact_prob,
                                                             K_grid, grid_day, dt_grid, p_gen_mech);
            target += lp_G + lp_sel - log_Z_sel;
          } else {
            target += lp_G;
          }
        }
      }

      Zobs[i] ~ neg_binomial_2(R_i, k_nb);
    }
  }
}

generated quantities {
  vector[N] log_lik_inc;
  vector[N] log_lik_trans;
  vector[N] log_lik_Z;
  vector[N] inc_obs;
  vector[N] gen_obs;
  vector[N] delta_obs;
  vector[N] si_obs;
  vector[K_grid] p_inc;
  vector[K_grid] p_lp;
  vector[K_grid] p_gen;
  vector[K_grid] p_delta;
  vector[K_grid] p_si;
  real inc_mean = 0;
  real gen_mean = 0;
  real delta_mean = 0;
  real si_mean = 0;
  real p_delta_presym = 0;
  real p_delta_by_day0 = 0;
  real p_delta_by_day1 = 0;
  real p_delta_by_day2 = 0;
  real p_gen_positive = 1;
  real vload_peak_time;

  vload_peak_time = m;

  p_inc = mech_grid_pmf_fast(K_grid, grid_day, dt_grid, K_int, dt_int, t_max, H_sym,
                           log_lambda_sym, log_r, log_d, m, log_kappa,
                           0, 0, 1);
  p_lp = mech_grid_pmf_fast(K_grid, grid_day, dt_grid, K_int, dt_int, t_max, H_inf,
                            log_lambda_inf, log_r, log_d, m, log_kappa,
                            0, 0, 1);
  p_gen = p_gen_mech;

  p_delta = rep_vector(0, K_grid);
  p_si = rep_vector(0, K_grid);

  // grid_day starts at 0. delta = generation interval - source incubation.
  // Store only delta values that fall on the same grid range after shifting by +delta_offset.
  // Here p_delta[k] corresponds to delta = grid_day[k] - 0.5 * t_max, approximately.
  // For stable export, use pair-level delta_obs as the primary delta estimate.
  for (a in 1:K_grid) {
    inc_mean += (grid_day[a] + 0.5 * dt_grid) * p_inc[a];
    gen_mean += (grid_day[a] + 0.5 * dt_grid) * p_gen[a];
  }

  // Serial interval distribution from independent delta and recipient incubation is easier
  // to summarise outside Stan from posterior draws of inc_obs/delta_obs. This placeholder
  // keeps dimensions valid for CSV extraction.
  for (k in 1:K_grid) {
    p_delta[k] = 1.0 / K_grid;
    p_si[k] = 1.0 / K_grid;
  }

  for (i in 1:N) {
    real P_i = T_onset[i] - T_inf[i];
    real log_R_i = log_R_at(T_onset[i], B, knot_day, log_R);
    real R_i = exp(log_R_i);

    inc_obs[i] = P_i;
    gen_obs[i] = -999;
    delta_obs[i] = -999;
    si_obs[i] = -999;

    log_lik_inc[i] = mech_first_event_logdens_fast(P_i, K_int, dt_int, t_max, H_sym,
                                           log_lambda_sym, log_r, log_d, m, log_kappa,
                                           0, 0, 1);

    if (source_idx[i] == 0) {
      log_lik_trans[i] = 0;
    } else {
      int src = source_idx[i];
      real G_i = T_inf[i] - T_inf[src];
      gen_obs[i] = G_i;
      delta_obs[i] = T_inf[i] - T_onset[src];
      si_obs[i] = T_onset[i] - T_onset[src];
      {
        real lp_G = log_grid_density_at(G_i, K_grid, grid_day, dt_grid, p_gen_mech);
        if (use_control_selection == 1 && post_control_source[i] == 1) {
          real P_src = T_onset[src] - T_inf[src];
          real lp_sel = control_log_weight(delta_obs[i], iso_delay_days, residual_contact_prob);
          real log_Z_sel = control_log_normalizer_from_grid(P_src, iso_delay_days, residual_contact_prob,
                                                         K_grid, grid_day, dt_grid, p_gen_mech);
          log_lik_trans[i] = lp_G + lp_sel - log_Z_sel;
        } else {
          log_lik_trans[i] = lp_G;
        }
      }
    }

    log_lik_Z[i] = neg_binomial_2_lpmf(Zobs[i] | R_i, k_nb);
  }

  // Pair-level summaries, excluding index cases.
  {
    int E = 0;
    real dsum = 0;
    real sisum = 0;
    for (i in 1:N) {
      if (source_idx[i] != 0) {
        E += 1;
        dsum += delta_obs[i];
        sisum += si_obs[i];
        if (delta_obs[i] < 0) p_delta_presym += 1;
        if (delta_obs[i] <= 0) p_delta_by_day0 += 1;
        if (delta_obs[i] <= 1) p_delta_by_day1 += 1;
        if (delta_obs[i] <= 2) p_delta_by_day2 += 1;
      }
    }
    if (E > 0) {
      delta_mean = dsum / E;
      si_mean = sisum / E;
      p_delta_presym /= E;
      p_delta_by_day0 /= E;
      p_delta_by_day1 /= E;
      p_delta_by_day2 /= E;
    }
  }
}
