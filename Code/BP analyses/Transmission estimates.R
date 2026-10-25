# Hierarchical negative-binomial model for ANDV transmission

library(cmdstanr)
library(posterior)
library(bayesplot)

# Individual-level offspring distributions
epuyen_pre <- c(
  5, 6, 0, 0, 1, 0, 0, 1, 10,
  3, 0, 0, 1, 1, 0, 0, 0, 0
)

lazaro_2007 <- c(
  3,
  rep(1, 7),
  rep(0, 43)
)

martinez_2005 <- c(
  rep(1, 3),
  rep(0, 10)
)

martinez_valdebenito_2014 <- c(
  3, 1, 0, 0, 0
)

# Primary reconstruction
alonso_2020 <- c(
  1, 1, 0
)

# Alternative sensitivity reconstruction:
alonso_2020_strict <- c(
  0, 1, 0
)

toro_1998 <- c(
  3,
  rep(1, 2),
  rep(0, 13)
)

individual_list <- list(
  "Epuyen pre-control"          = epuyen_pre,
  "Lazaro 2007"                 = lazaro_2007,
  "Martinez 2005"               = martinez_2005,
  "Martinez-Valdebenito 2014"   = martinez_valdebenito_2014,
  "Alonso 2020"                 = alonso_2020,
  "Toro 1998"                   = toro_1998
)

# Convert list to long-format data
individual_data <- do.call(
  rbind,
  lapply(seq_along(individual_list), function(j) {
    data.frame(
      study = names(individual_list)[j],
      offspring = individual_list[[j]]
    )
  })
)

rownames(individual_data) <- NULL

# Aggregate study data
aggregate_data <- data.frame(
  study = c(
    "Barrera 2025",
    "Bellomo 2015",
    "Bellomo 2025",
    "Castillo 2000",
    "Ferres 2007",
    "Ferres 2020",
    "Ferres 2024",
    "Iglesias 2022",
    "Kofman 2018"
  ),
  confirmed_cases = c(
    2, 73, 1, 20, 92, 2, 131, 154, 1
  ),
  secondary_cases = c(
    2, 0, 0, 0, 16, 1, 21, 38, 0
  )
)

# Create common study identifiers
study_names <- c(
  names(individual_list),
  aggregate_data$study
)

stopifnot(!anyDuplicated(study_names))

J <- length(study_names)

individual_data$study_id <- match(
  individual_data$study,
  study_names
)

aggregate_data$study_id <- match(
  aggregate_data$study,
  study_names
)

# Data supplied to Stan
stan_data <- list(
  J = J,
  
  # Individual-level observations
  N_ind = nrow(individual_data),
  offspring = individual_data$offspring,
  study_ind = individual_data$study_id,
  
  # Aggregate observations
  N_agg = nrow(aggregate_data),
  primary_count = aggregate_data$confirmed_cases,
  secondary_count = aggregate_data$secondary_cases,
  study_agg = aggregate_data$study_id
)

# Stan model
stan_code <- "
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
"

writeLines(
  stan_code,
  con = "andv_hierarchical_nb.stan"
)

# Compile and fit
model <- cmdstan_model(
  "andv_hierarchical_nb.stan"
)

fit <- model$sample(
  data = stan_data,
  seed = 2026,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 2000,
  iter_sampling = 4000,
  adapt_delta = 0.99,
  max_treedepth = 15,
  refresh = 500
)

# Robust posterior summaries
library(posterior)
library(bayesplot)

# Helper function: summarise selected parameters directly
summarise_parameters <- function(fit_object, variables) {
  
  draws_df <- posterior::as_draws_df(
    fit_object$draws(variables = variables)
  )
  
  out <- lapply(variables, function(v) {
    
    x <- draws_df[[v]]
    
    data.frame(
      variable = v,
      mean = mean(x),
      sd = sd(x),
      median = median(x),
      q2.5 = unname(stats::quantile(x, 0.025)),
      q97.5 = unname(stats::quantile(x, 0.975))
    )
  })
  
  do.call(rbind, out)
}


# Pooled estimates
pooled_summary <- summarise_parameters(
  fit_object = fit,
  variables = c(
    "pooled_R0",
    "pooled_k",
    "tau_log_R"
  )
)

# Add convergence diagnostics separately
pooled_diagnostics <- fit$summary(
  variables = c(
    "pooled_R0",
    "pooled_k",
    "tau_log_R"
  )
)

pooled_summary <- merge(
  pooled_summary,
  pooled_diagnostics[
    ,
    c(
      "variable",
      "rhat",
      "ess_bulk",
      "ess_tail"
    )
  ],
  by = "variable",
  all.x = TRUE,
  sort = FALSE
)

print(pooled_summary)


# Study-specific R estimates
R_variables <- paste0("R[", seq_len(J), "]")

study_R_summary <- summarise_parameters(
  fit_object = fit,
  variables = R_variables
)

study_R_diagnostics <- fit$summary(
  variables = R_variables
)

study_R_summary <- merge(
  study_R_summary,
  study_R_diagnostics[
    ,
    c(
      "variable",
      "rhat",
      "ess_bulk",
      "ess_tail"
    )
  ],
  by = "variable",
  all.x = TRUE,
  sort = FALSE
)

# Restore intended R[1], R[2], ... ordering after merge
study_R_summary$study_id <- as.integer(
  sub(
    pattern = "R\\[([0-9]+)\\]",
    replacement = "\\1",
    x = study_R_summary$variable
  )
)

study_R_summary <- study_R_summary[
  order(study_R_summary$study_id),
]

study_R_summary$study <- study_names[
  study_R_summary$study_id
]

study_R_summary <- study_R_summary[
  ,
  c(
    "study",
    "mean",
    "sd",
    "median",
    "q2.5",
    "q97.5",
    "rhat",
    "ess_bulk",
    "ess_tail"
  )
]

rownames(study_R_summary) <- NULL

print(study_R_summary)


# Sampling diagnostics
diagnostic_summary <- fit$diagnostic_summary()

print(diagnostic_summary)

if ("num_divergent" %in% names(diagnostic_summary)) {
  cat(
    "Total divergent transitions:",
    sum(diagnostic_summary$num_divergent),
    "\n"
  )
}

if ("num_max_treedepth" %in% names(diagnostic_summary)) {
  cat(
    "Total maximum-treedepth transitions:",
    sum(diagnostic_summary$num_max_treedepth),
    "\n"
  )
}


# Trace plots
draws_main <- fit$draws(
  variables = c(
    "pooled_R0",
    "pooled_k",
    "tau_log_R"
  )
)

bayesplot::mcmc_trace(draws_main)

bayesplot::mcmc_pairs(
  draws_main,
  pars = c(
    "pooled_R0",
    "pooled_k",
    "tau_log_R"
  )
)


# Posterior probabilities
main_draws <- posterior::as_draws_df(
  fit$draws(
    variables = c(
      "pooled_R0",
      "pooled_k",
      "tau_log_R"
    )
  )
)

cat(
  "P(pooled R0 > 1) =",
  mean(main_draws$pooled_R0 > 1),
  "\n"
)

cat(
  "P(k < 0.2) =",
  mean(main_draws$pooled_k < 0.2),
  "\n"
)

cat(
  "P(k < 1) =",
  mean(main_draws$pooled_k < 1),
  "\n"
)


# Posterior predictive checks: individual-level data
offspring_rep_draws <- posterior::as_draws_matrix(
  fit$draws("offspring_rep")
)

set.seed(2026)

n_ppc_ind <- min(
  100,
  nrow(offspring_rep_draws)
)

ppc_rows_ind <- sample(
  seq_len(nrow(offspring_rep_draws)),
  size = n_ppc_ind
)

bayesplot::ppc_bars(
  y = stan_data$offspring,
  yrep = offspring_rep_draws[
    ppc_rows_ind,
    ,
    drop = FALSE
  ]
)


# Posterior predictive checks: aggregate data
secondary_rep_draws <- posterior::as_draws_matrix(
  fit$draws("secondary_rep")
)

n_ppc_agg <- min(
  100,
  nrow(secondary_rep_draws)
)

ppc_rows_agg <- sample(
  seq_len(nrow(secondary_rep_draws)),
  size = n_ppc_agg
)

bayesplot::ppc_intervals(
  y = stan_data$secondary_count,
  yrep = secondary_rep_draws[
    ppc_rows_agg,
    ,
    drop = FALSE
  ],
  x = seq_len(stan_data$N_agg)
)


# Sensitivity analysis: strict Alonso reconstruction
individual_list_strict <- individual_list

individual_list_strict[["Alonso 2020"]] <-
  alonso_2020_strict

individual_data_strict <- do.call(
  rbind,
  lapply(
    seq_along(individual_list_strict),
    function(j) {
      data.frame(
        study = names(individual_list_strict)[j],
        offspring = individual_list_strict[[j]]
      )
    }
  )
)

rownames(individual_data_strict) <- NULL

individual_data_strict$study_id <- match(
  individual_data_strict$study,
  study_names
)

stan_data_strict <- stan_data

stan_data_strict$N_ind <-
  nrow(individual_data_strict)

stan_data_strict$offspring <-
  individual_data_strict$offspring

stan_data_strict$study_ind <-
  individual_data_strict$study_id

fit_strict <- model$sample(
  data = stan_data_strict,
  seed = 2027,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 2000,
  iter_sampling = 4000,
  adapt_delta = 0.99,
  max_treedepth = 15,
  refresh = 500
)


# Summarise strict Alonso model
strict_summary <- summarise_parameters(
  fit_object = fit_strict,
  variables = c(
    "pooled_R0",
    "pooled_k",
    "tau_log_R"
  )
)

strict_diagnostics <- fit_strict$summary(
  variables = c(
    "pooled_R0",
    "pooled_k",
    "tau_log_R"
  )
)

strict_summary <- merge(
  strict_summary,
  strict_diagnostics[
    ,
    c(
      "variable",
      "rhat",
      "ess_bulk",
      "ess_tail"
    )
  ],
  by = "variable",
  all.x = TRUE,
  sort = FALSE
)

print(strict_summary)


# Compare primary and strict analyses
primary_compare <- pooled_summary[
  pooled_summary$variable %in%
    c("pooled_R0", "pooled_k"),
  c(
    "variable",
    "mean",
    "median",
    "q2.5",
    "q97.5"
  )
]

primary_compare$model <-
  "Primary Alonso reconstruction"

strict_compare <- strict_summary[
  strict_summary$variable %in%
    c("pooled_R0", "pooled_k"),
  c(
    "variable",
    "mean",
    "median",
    "q2.5",
    "q97.5"
  )
]

strict_compare$model <-
  "Strict Alonso reconstruction"

comparison <- rbind(
  primary_compare,
  strict_compare
)

comparison <- comparison[
  ,
  c(
    "model",
    "variable",
    "mean",
    "median",
    "q2.5",
    "q97.5"
  )
]

rownames(comparison) <- NULL

print(comparison)

fit_nb <- fitdistrplus::fitdist(
  epuyen_pre,
  "nbinom",
  method = "mle",
  fix.arg = list(mu = 2.12)
)

summary(fit_nb)