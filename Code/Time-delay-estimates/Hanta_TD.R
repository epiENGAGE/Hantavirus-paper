# Estimate the time delays
setwd("~/Hantavirus-paper/Code/Time-delay-estimates")

suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
})

# Define directories
this_file_dir <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) getwd())
stan_file <- file.path(this_file_dir, "Hanta_TD.stan")
if (!file.exists(stan_file)) stan_file <- "Hanta_TD.stan"
if (!file.exists(stan_file)) stop("Could not find Hanta_TD.stan")

out_dir <- file.path(getwd(), "andes_wh_event_time_line_list_paper_model_FAST_COMPILES_outputs")
out_dir <- "~/Hantavirus-paper/TD results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

linelist_file <- NULL
linelist_url <- "https://raw.githubusercontent.com/epiforecasts/andv-linelist-analysis/main/data/linelist.csv"

multi_source_policy <- "first"   # "first" or "drop"
index_pre_onset_max <- 80.0

# Mechanistic integration support for first-event LP/IP and observed-pair GT.
t_max <- 120.0
K_int <- 240L
# dt_int = t_max / K_int = 0.5 days by default. Increase K_int for accuracy.
dt_int <- t_max / K_int

# Export grid for posterior PMFs.
grid_day <- seq(0, t_max - 1, by = 1)
dt_grid <- 1.0

# Retained only so old data lists remain compatible with the Stan file.
# The paper event-time model below ignores this gate.
use_transmission_gate <- 0L

# Optional weak prior that shifts the *mechanistic population mean* delta negative:
#   delta_mean_mech = mean(generation interval) - mean(source incubation)
# This does NOT force every pair to be presymptomatic. It only says the natural
# mechanism probably has a slightly negative mean delta.
use_delta_mean_prior <- 1L
prior_delta_mean <- -1.0
prior_delta_sd <- 1.5

use_control_selection <- 1L
control_date <- as.Date("2018-12-31")
post_control_iso_delay_days <- 2.0
residual_contact_prob <- 0.05

prior <- list(
  log_lambda_sym_mean = log(0.03),
  log_lambda_sym_sd = 0.35,

  log_lambda_inf_mean = log(0.04),
  log_lambda_inf_sd = 0.35,

  log_r_mean = log(0.65),
  log_r_sd = 0.25,

  log_d_mean = log(0.13),
  log_d_sd = 0.25,

  m_mean = 14.0,
  m_sd = 2.0,

  log_kappa_mean = log(0.35),
  log_kappa_sd = 0.50,

  gate_onset_mean = 15.0,
  gate_onset_sd = 7.0,
  log_gate_slope_mean = log(2.0),
  log_gate_slope_sd = 0.5,

  log_R_init_mean = log(1.5),
  log_R_init_sd = 1.0,
  sigma_rw_sd = 0.5,
  phi_inv_sqrt_sd = 1.0
)

chains <- 4
parallel_chains <- min(4, chains)
iter_warmup <- 500
iter_sampling <- 500
adapt_delta <- 0.98
max_treedepth <- 12
seed <- 20260519

# Fast debugging controls. Turn run_prior_predictive back on only for final checks.
run_prior_predictive <- FALSE


# Data loading
if (is.null(linelist_file)) {
  linelist_file <- file.path(out_dir, "linelist.csv")
  if (!file.exists(linelist_file)) {
    download.file(linelist_url, linelist_file, mode = "wb", quiet = FALSE)
  }
}

ll_raw <- read.csv(linelist_file, stringsAsFactors = FALSE, na.strings = c("NA", ""))

required_cols <- c("patient_id", "source_case", "onset_date")
missing_cols <- setdiff(required_cols, names(ll_raw))
if (length(missing_cols) > 0) stop("Missing required columns: ", paste(missing_cols, collapse = ", "))

ll <- ll_raw %>%
  filter(!grepl("_alt$", as.character(patient_id))) %>%
  mutate(
    patient_id = as.character(patient_id),
    source_case = as.character(source_case),
    onset_date = as.Date(onset_date)
  )

if (all(c("onset_lower", "onset_upper") %in% names(ll))) {
  ll <- ll %>%
    mutate(
      onset_lower = as.Date(onset_lower),
      onset_upper = as.Date(onset_upper)
    )
} else {
  ll <- ll %>%
    mutate(
      onset_lower = onset_date,
      onset_upper = onset_date + 1
    )
}

# Exposure windows are required for sourced cases. If absent for a sourced case,
# use a broad fallback relative to the source onset. The fallback is diagnostic only;
# if many cases use it, fix the data instead of trusting the fit.
if ("exposure_lower" %in% names(ll)) ll$exposure_lower <- as.Date(ll$exposure_lower)
if ("exposure_upper" %in% names(ll)) ll$exposure_upper <- as.Date(ll$exposure_upper)

parse_source <- function(x, policy = c("first", "drop")) {
  policy <- match.arg(policy)
  if (length(x) != 1 || is.na(x) || x == "" || x == "index") return(NA_character_)
  x <- as.character(x)
  if (grepl("/", x, fixed = TRUE)) {
    if (policy == "drop") return(NA_character_)
    return(strsplit(x, "/", fixed = TRUE)[[1]][1])
  }
  x
}

ll$source_id <- vapply(ll$source_case, parse_source, character(1), policy = multi_source_policy)

id_levels <- ll$patient_id
id_to_idx <- setNames(seq_along(id_levels), id_levels)
ll$source_idx <- ifelse(is.na(ll$source_id), 0L, as.integer(id_to_idx[ll$source_id]))
if (any(is.na(ll$source_idx))) {
  bad <- ll$source_id[is.na(ll$source_idx)]
  stop("Some source IDs are not present as patient_id: ", paste(unique(bad), collapse = ", "))
}

# Offspring counts from the directed source list.
Zobs <- tabulate(ll$source_idx[ll$source_idx > 0], nbins = nrow(ll))
ll$Zobs <- as.integer(Zobs)

# Control-selection indicator. For each recipient, this is 1 when its source's
# onset was after implementation of isolation/quarantine, so observed late
# post-onset transmission opportunities are treated as selected/truncated.
source_onset_by_idx <- ll$onset_date
ll$post_control_source <- ifelse(
  ll$source_idx > 0,
  as.integer(source_onset_by_idx[ll$source_idx] >= control_date),
  0L
)
message("Post-control sourced pairs: ", sum(ll$post_control_source), " / ", sum(ll$source_idx > 0))

# Day origin.
t0 <- min(ll$onset_lower, ll$exposure_lower, na.rm = TRUE)
ll <- ll %>%
  mutate(
    onset_lo_day = as.numeric(onset_lower - t0),
    onset_hi_day = as.numeric(onset_upper - t0),
    onset_hi_day = if_else(onset_hi_day <= onset_lo_day, onset_lo_day + 1, onset_hi_day)
  )

# Fill exposure windows for index and any missing sourced cases.
source_onset_lookup <- ll %>% transmute(patient_id, source_onset_date = onset_date)
ll <- ll %>% left_join(source_onset_lookup, by = c("source_id" = "patient_id"))

if (!all(c("exposure_lower", "exposure_upper") %in% names(ll))) {
  ll$exposure_lower <- as.Date(NA)
  ll$exposure_upper <- as.Date(NA)
}

ll <- ll %>%
  mutate(
    fallback_exp_lower = if_else(
      source_idx == 0L,
      onset_date - as.integer(index_pre_onset_max),
      source_onset_date - 40L
    ),
    fallback_exp_upper = if_else(
      source_idx == 0L,
      onset_date - 1L,
      source_onset_date + 80L
    ),
    exposure_lower_filled = coalesce(exposure_lower, fallback_exp_lower),
    exposure_upper_filled = coalesce(exposure_upper, fallback_exp_upper),
    used_fallback_exposure = is.na(exposure_lower) | is.na(exposure_upper),
    exp_lo_day = as.numeric(exposure_lower_filled - t0),
    exp_hi_day = as.numeric(exposure_upper_filled - t0),
    exp_hi_day = if_else(exp_hi_day <= exp_lo_day, exp_lo_day + 1, exp_hi_day)
  )

write.csv(ll, file.path(out_dir, "linelist_used.csv"), row.names = FALSE)

# Weekly R(t) knots spanning onset times.
knot_day <- seq(floor(min(ll$onset_lo_day)), ceiling(max(ll$onset_hi_day)) + 7, by = 7)
if (length(knot_day) < 2) knot_day <- c(min(ll$onset_lo_day), max(ll$onset_hi_day) + 7)

stan_data <- list(
  N = nrow(ll),
  source_idx = as.integer(ll$source_idx),
  Zobs = as.integer(ll$Zobs),
  onset_lo_day = as.numeric(ll$onset_lo_day),
  onset_hi_day = as.numeric(ll$onset_hi_day),
  exp_lo_day = as.numeric(ll$exp_lo_day),
  exp_hi_day = as.numeric(ll$exp_hi_day),
  index_pre_onset_max = index_pre_onset_max,
  B = length(knot_day),
  knot_day = as.numeric(knot_day),
  t_max = t_max,
  K_int = as.integer(K_int),
  dt_int = dt_int,
  K_grid = length(grid_day),
  grid_day = as.numeric(grid_day),
  dt_grid = dt_grid,
  prior_only = 0L,
  use_transmission_gate = as.integer(use_transmission_gate),
  use_delta_mean_prior = as.integer(use_delta_mean_prior),
  prior_delta_mean = prior_delta_mean,
  prior_delta_sd = prior_delta_sd,
  use_control_selection = as.integer(use_control_selection),
  post_control_source = as.integer(ll$post_control_source),
  iso_delay_days = post_control_iso_delay_days,
  residual_contact_prob = residual_contact_prob,
  prior_log_lambda_sym_mean = prior$log_lambda_sym_mean,
  prior_log_lambda_sym_sd = prior$log_lambda_sym_sd,
  prior_log_lambda_inf_mean = prior$log_lambda_inf_mean,
  prior_log_lambda_inf_sd = prior$log_lambda_inf_sd,
  prior_log_r_mean = prior$log_r_mean,
  prior_log_r_sd = prior$log_r_sd,
  prior_log_d_mean = prior$log_d_mean,
  prior_log_d_sd = prior$log_d_sd,
  prior_m_mean = prior$m_mean,
  prior_m_sd = prior$m_sd,
  prior_log_kappa_mean = prior$log_kappa_mean,
  prior_log_kappa_sd = prior$log_kappa_sd,
  prior_gate_onset_mean = prior$gate_onset_mean,
  prior_gate_onset_sd = prior$gate_onset_sd,
  prior_log_gate_slope_mean = prior$log_gate_slope_mean,
  prior_log_gate_slope_sd = prior$log_gate_slope_sd,
  prior_log_R_init_mean = prior$log_R_init_mean,
  prior_log_R_init_sd = prior$log_R_init_sd,
  prior_sigma_rw_sd = prior$sigma_rw_sd,
  prior_phi_inv_sqrt_sd = prior$phi_inv_sqrt_sd
)

prior_data <- stan_data
prior_data$prior_only <- 1L
saveRDS(stan_data, file.path(out_dir, "stan_data.rds"))

# Initial values that keep latent infection before onset for early adaptation.
make_inits <- function() {
  function() {
    onset_raw <- rep(0.5, stan_data$N)
    inf_raw <- rep(0.5, stan_data$N)
    list(
      onset_raw = onset_raw,
      inf_raw = inf_raw,
      log_lambda_sym = prior$log_lambda_sym_mean,
      log_lambda_inf = prior$log_lambda_inf_mean,
      log_r = prior$log_r_mean,
      log_d = prior$log_d_mean,
      m = prior$m_mean,
      log_kappa = prior$log_kappa_mean,
      gate_onset = prior$gate_onset_mean,
      log_gate_slope = prior$log_gate_slope_mean,
      log_R_init = prior$log_R_init_mean,
      eps_R = rep(0, stan_data$B - 1),
      sigma_rw = 0.1,
      phi_inv_sqrt = 0.5
    )
  }
}

# Compile and fit
mod <- cmdstan_model(stan_file, force_recompile = TRUE)

if (isTRUE(run_prior_predictive)) {
  message("Running prior predictive fit...")
  fit_prior <- mod$sample(
    data = prior_data,
    chains = chains,
    parallel_chains = parallel_chains,
    iter_warmup = 250,
    iter_sampling = 250,
    seed = seed + 1,
    init = make_inits(),
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth,
    refresh = 25
  )
  fit_prior$save_object(file.path(out_dir, "fit_prior.rds"))
} else {
  message("Skipping prior predictive fit for speed. Set run_prior_predictive <- TRUE for final checks.")
}

message("Running posterior fit...")
fit <- mod$sample(
  data = stan_data,
  chains = chains,
  parallel_chains = parallel_chains,
  iter_warmup = iter_warmup,
  iter_sampling = iter_sampling,
  seed = seed,
  init = make_inits(),
  adapt_delta = adapt_delta,
  max_treedepth = max_treedepth,
  refresh = 25
)
fit$save_object(file.path(out_dir, "fit_posterior.rds"))

# Summaries and diagnostics
scalar_vars <- c(
  "log_lambda_sym", "log_lambda_inf", "log_r", "log_d", "m", "log_kappa",
  "gate_onset", "gate_slope", "inc_mean_mech", "gen_mean_mech", "delta_mean_mech",
  "inc_mean", "gen_mean", "delta_mean", "si_mean",
  "p_delta_presym", "p_delta_by_day0", "p_delta_by_day1", "p_delta_by_day2",
  "sigma_rw", "phi_inv_sqrt", "k_nb", "vload_peak_time"
)

summ <- fit$summary(variables = scalar_vars)
write.csv(summ, file.path(out_dir, "posterior_scalar_summary.csv"), row.names = FALSE)
print(summ)

diag_df <- fit$sampler_diagnostics(format = "df")
write.csv(diag_df, file.path(out_dir, "sampler_diagnostics.csv"), row.names = FALSE)
message("Divergences: ", sum(diag_df$divergent__), " / ", nrow(diag_df))
message("Max treedepth hits: ", sum(diag_df$treedepth__ >= max_treedepth), " / ", nrow(diag_df))

extract_draw_matrix <- function(fit_obj, variable, n) {
  mat <- as.matrix(fit_obj$draws(variables = variable, format = "draws_matrix"))
  cols <- paste0(variable, "[", seq_len(n), "]")
  missing <- setdiff(cols, colnames(mat))
  if (length(missing) > 0) stop("Missing draw columns: ", paste(head(missing), collapse = ", "))
  mat[, cols, drop = FALSE]
}

summarise_matrix <- function(mat, x_name, x_values) {
  tibble(
    !!x_name := x_values,
    mean = colMeans(mat),
    q05 = apply(mat, 2, quantile, 0.05),
    q50 = apply(mat, 2, quantile, 0.50),
    q95 = apply(mat, 2, quantile, 0.95)
  )
}

p_inc_mat <- extract_draw_matrix(fit, "p_inc", length(grid_day))
p_lp_mat <- extract_draw_matrix(fit, "p_lp", length(grid_day))
p_gen_mat <- extract_draw_matrix(fit, "p_gen", length(grid_day))
inc_obs_mat <- extract_draw_matrix(fit, "inc_obs", nrow(ll))
delta_obs_mat <- extract_draw_matrix(fit, "delta_obs", nrow(ll))
gen_obs_mat <- extract_draw_matrix(fit, "gen_obs", nrow(ll))
si_obs_mat <- extract_draw_matrix(fit, "si_obs", nrow(ll))

p_inc_df <- summarise_matrix(p_inc_mat, "day_since_infection", grid_day)
p_lp_df <- summarise_matrix(p_lp_mat, "day_since_infection", grid_day)
p_gen_df <- summarise_matrix(p_gen_mat, "day_since_infection", grid_day)

write.csv(p_inc_df, file.path(out_dir, "posterior_mechanistic_incubation_pmf.csv"), row.names = FALSE)
write.csv(p_lp_df, file.path(out_dir, "posterior_mechanistic_latent_period_pmf.csv"), row.names = FALSE)
write.csv(p_gen_df, file.path(out_dir, "posterior_mechanistic_generation_pmf.csv"), row.names = FALSE)
write.csv(as.data.frame(p_inc_mat), file.path(out_dir, "bp_draws_incubation_pmf_matrix.csv"), row.names = FALSE)
write.csv(as.data.frame(p_lp_mat), file.path(out_dir, "bp_draws_latent_period_pmf_matrix.csv"), row.names = FALSE)
write.csv(as.data.frame(p_gen_mat), file.path(out_dir, "bp_draws_generation_pmf_matrix.csv"), row.names = FALSE)
write.csv(data.frame(day_since_infection = grid_day), file.path(out_dir, "bp_day_grid.csv"), row.names = FALSE)

# Pair-level posterior distributions, in the epiforecasts style.
sourced_idx <- which(ll$source_idx > 0)
pair_draws <- bind_rows(lapply(seq_along(sourced_idx), function(j) {
  i <- sourced_idx[j]
  
  tibble(
    draw = seq_len(nrow(as.matrix(delta_obs_mat))),
    recipient_id = ll$patient_id[i],
    source_id = ll$patient_id[ll$source_idx[i]],
    incubation = as.numeric(as.matrix(inc_obs_mat)[, i]),
    generation_interval = as.numeric(as.matrix(gen_obs_mat)[, i]),
    delta = as.numeric(as.matrix(delta_obs_mat)[, i]),
    serial_interval = as.numeric(as.matrix(si_obs_mat)[, i])
  )
}))
write.csv(pair_draws, file.path(out_dir, "posterior_pair_level_draws.csv"), row.names = FALSE)

pair_summary <- pair_draws %>%
  group_by(recipient_id, source_id) %>%
  summarise(
    incubation_q50 = median(incubation),
    incubation_q05 = quantile(incubation, 0.05),
    incubation_q95 = quantile(incubation, 0.95),
    generation_q50 = median(generation_interval),
    generation_q05 = quantile(generation_interval, 0.05),
    generation_q95 = quantile(generation_interval, 0.95),
    delta_q50 = median(delta),
    delta_q05 = quantile(delta, 0.05),
    delta_q95 = quantile(delta, 0.95),
    serial_q50 = median(serial_interval),
    serial_q05 = quantile(serial_interval, 0.05),
    serial_q95 = quantile(serial_interval, 0.95),
    .groups = "drop"
  )
write.csv(pair_summary, file.path(out_dir, "posterior_pair_level_summary.csv"), row.names = FALSE)

# Plot for transmission timing: pair-specific posterior intervals.
p_delta_pair_intervals <- pair_summary %>%
  mutate(pair_label = paste0(source_id, " -> ", recipient_id)) %>%
  arrange(delta_q50) %>%
  mutate(pair_label = factor(pair_label, levels = pair_label)) %>%
  ggplot(aes(y = pair_label, x = delta_q50)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_errorbarh(aes(xmin = delta_q05, xmax = delta_q95), height = 0.2, alpha = 0.8) +
  geom_point(size = 1.8) +
  labs(
    x = "Transmission timing relative to source symptom onset, days",
    y = "Linked transmission pair",
    title = "Pair-specific posterior transmission timing",
    subtitle = "Points are posterior medians; intervals are 90% credible intervals. Negative values imply pre-symptomatic transmission."
  ) +
  theme_minimal()

ggsave(file.path(out_dir, "posterior_delta_pair_specific_intervals.png"),
       p_delta_pair_intervals, width = 8, height = max(4.8, 0.22 * nrow(pair_summary)), dpi = 300)

# Observed source-onset serial intervals for comparison.
obs_pairs <- ll %>%
  filter(source_idx > 0) %>%
  transmute(
    recipient_id = patient_id,
    source_id = patient_id[source_idx],
    observed_si = as.numeric(onset_date - onset_date[source_idx])
  )
write.csv(obs_pairs, file.path(out_dir, "observed_pairs.csv"), row.names = FALSE)

# Plot incubation PMF.
p_inc_plot <- ggplot(p_inc_df, aes(x = day_since_infection)) +
  geom_ribbon(aes(ymin = q05, ymax = q95), alpha = 0.25) +
  geom_line(aes(y = q50), linewidth = 0.9) +
  labs(x = "Days since infection", y = "Probability", title = "Mechanistic incubation distribution") +
  theme_minimal()
ggsave(file.path(out_dir, "posterior_mechanistic_incubation_pmf.png"), p_inc_plot, width = 7, height = 4.5, dpi = 300)

# Plot latent-period PMF.
p_lp_plot <- ggplot(p_lp_df, aes(x = day_since_infection)) +
  geom_ribbon(aes(ymin = q05, ymax = q95), alpha = 0.25) +
  geom_line(aes(y = q50), linewidth = 0.9) +
  labs(x = "Days since infection", y = "Probability", title = "Mechanistic latent-period distribution") +
  theme_minimal()
ggsave(file.path(out_dir, "posterior_mechanistic_latent_period_pmf.png"), p_lp_plot, width = 7, height = 4.5, dpi = 300)

# Plot generation PMF.
p_gen_plot <- ggplot(p_gen_df, aes(x = day_since_infection)) +
  geom_ribbon(aes(ymin = q05, ymax = q95), alpha = 0.25) +
  geom_line(aes(y = q50), linewidth = 0.9) +
  labs(x = "Days since source infection", y = "Probability", title = "Observed-pair generation-interval distribution (GT = LP + U)") +
  theme_minimal()
ggsave(file.path(out_dir, "posterior_mechanistic_generation_pmf.png"), p_gen_plot, width = 7, height = 4.5, dpi = 300)

# Plot pair-level delta distribution relative to source onset.
delta_df <- pair_draws %>%
  group_by(draw) %>%
  summarise(delta_mean = mean(delta), .groups = "drop")

p_delta_plot <- ggplot(pair_draws, aes(x = delta)) +
  geom_histogram(aes(y = after_stat(density)), bins = 40, alpha = 0.35) +
  geom_vline(xintercept = 0, linetype = 2) +
  labs(
    x = "Transmission timing relative to source symptom onset, days",
    y = "Pooled posterior density across pair-level latent times",
    title = "Diagnostic pooled pair-level transmission timing"
  ) +
  theme_minimal()
ggsave(file.path(out_dir, "posterior_delta_pair_level.png"), p_delta_plot, width = 7, height = 4.5, dpi = 300)

# PPC-style plot: posterior latent serial intervals against observed source-onset intervals.
si_ppc_df <- pair_draws %>%
  mutate(si_bin = round(serial_interval)) %>%
  count(draw, si_bin, name = "n") %>%
  group_by(si_bin) %>%
  summarise(
    q05 = quantile(n, 0.05),
    q50 = quantile(n, 0.50),
    q95 = quantile(n, 0.95),
    .groups = "drop"
  ) %>%
  complete(si_bin = seq(floor(min(pair_draws$serial_interval)), ceiling(max(pair_draws$serial_interval))),
           fill = list(q05 = 0, q50 = 0, q95 = 0))

obs_counts <- obs_pairs %>%
  count(observed_si, name = "observed")

p_ppc <- ggplot(si_ppc_df, aes(x = si_bin)) +
  geom_ribbon(aes(ymin = q05, ymax = q95), alpha = 0.25) +
  geom_line(aes(y = q50), linewidth = 0.9) +
  geom_col(data = obs_counts, aes(x = observed_si, y = observed), alpha = 0.45, width = 0.8, inherit.aes = FALSE) +
  labs(
    x = "Serial interval, days",
    y = "Count",
    title = "Posterior latent serial intervals vs observed onset-to-onset intervals",
    subtitle = "Bars are observed onset-date differences; line/ribbon summarise latent T_onset differences"
  ) +
  theme_minimal()
ggsave(file.path(out_dir, "posterior_latent_serial_interval_ppc.png"), p_ppc, width = 8, height = 4.8, dpi = 300)
print(p_ppc)

message("Done. Outputs written to: ", out_dir)

library(dplyr)
library(tibble)
library(purrr)

pair_draws_joint <- bind_rows(lapply(seq_along(sourced_idx), function(j) {
  i <- sourced_idx[j]
  s <- ll$source_idx[i]
  
  tibble(
    draw = seq_len(nrow(gen_obs_mat)),
    recipient_id = ll$patient_id[i],
    source_id = ll$patient_id[s],
    
    generation_interval = as.numeric(gen_obs_mat[, i]),
    delta = as.numeric(delta_obs_mat[, i]),
    serial_interval = as.numeric(si_obs_mat[, i]),
    
    recipient_incubation = as.numeric(inc_obs_mat[, i]),
    source_incubation = as.numeric(inc_obs_mat[, s])
  )
}))

library(ggplot2)

ggplot(pair_draws_joint,
       aes(x = source_incubation, y = generation_interval)) +
  stat_density_2d(
    aes(fill = after_stat(level)),
    geom = "polygon",
    contour = TRUE,
    alpha = 0.5
  ) +
  geom_density_2d(color = "black", linewidth = 0.6) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  labs(
    x = "Source incubation period (days since source infection)",
    y = "Generation interval (days since source infection)",
    title = "Joint posterior of generation interval and source incubation",
    subtitle = "Dashed line is generation interval = source incubation; points below imply presymptomatic transmission"
  ) +
  theme_minimal()

# Plot 1: GT vs source incubation (joint posterior)

# If not already built:
pair_draws_joint <- bind_rows(lapply(seq_along(sourced_idx), function(j) {
  i <- sourced_idx[j]
  s <- ll$source_idx[i]
  
  tibble(
    draw = seq_len(nrow(gen_obs_mat)),
    recipient_id = ll$patient_id[i],
    source_id = ll$patient_id[s],
    generation_interval = as.numeric(gen_obs_mat[, i]),
    delta = as.numeric(delta_obs_mat[, i]),
    serial_interval = as.numeric(si_obs_mat[, i]),
    recipient_incubation = as.numeric(inc_obs_mat[, i]),
    source_incubation = as.numeric(inc_obs_mat[, s])
  )
}))

# Optional: thin draws for plotting speed/clarity
set.seed(1)
plot_df <- pair_draws_joint %>%
  sample_n(min(20000, n()))

p_gt_vs_ip <- ggplot(plot_df, aes(x = source_incubation, y = generation_interval)) +
  geom_point(alpha = 0.03, size = 0.5) +
  geom_density_2d(linewidth = 0.7) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  coord_equal() +
  labs(
    x = "Source incubation period (days)",
    y = "Generation interval (days)",
    title = "Generation interval vs source incubation period",
    subtitle = "Dashed line: GT = source IP; below the line implies pre-symptomatic transmission"
  ) +
  theme_minimal()

ggsave(
  file.path(out_dir, "gt_vs_source_ip_joint.png"),
  p_gt_vs_ip,
  width = 7,
  height = 6,
  dpi = 300
)

print(p_gt_vs_ip)

gt_ip_summary <- pair_draws_joint %>%
  summarise(
    p_gt_lt_ip = mean(generation_interval < source_incubation),
    p_gt_eq_or_lt_ip = mean(generation_interval <= source_incubation),
    mean_delta = mean(delta),
    median_delta = median(delta)
  )

print(gt_ip_summary)

# Plot 2: observed serial intervals with posterior predictive overlay

# Observed onset-to-onset intervals
obs_si <- obs_pairs$observed_si
obs_si <- obs_si[is.finite(obs_si)]

# Build posterior predictive density for serial intervals
# We estimate a density within each draw, then summarise across draws.

si_grid <- seq(
  floor(min(c(obs_si, pair_draws$serial_interval), na.rm = TRUE)),
  ceiling(max(c(obs_si, pair_draws$serial_interval), na.rm = TRUE)),
  by = 0.25
)

dens_by_draw <- lapply(split(pair_draws$serial_interval, pair_draws$draw), function(x) {
  d <- density(x, from = min(si_grid), to = max(si_grid), n = length(si_grid), bw = "nrd0")
  tibble(
    x = d$x,
    y = d$y
  )
})

dens_df <- bind_rows(lapply(seq_along(dens_by_draw), function(i) {
  dens_by_draw[[i]] %>%
    mutate(draw = i)
}))

dens_sum <- dens_df %>%
  group_by(x) %>%
  summarise(
    q05 = quantile(y, 0.05),
    q50 = quantile(y, 0.50),
    q95 = quantile(y, 0.95),
    .groups = "drop"
  )

p_si_overlay <- ggplot() +
  geom_histogram(
    aes(x = obs_si, y = after_stat(density)),
    bins = 18,
    fill = "grey70",
    color = "white",
    alpha = 0.7
  ) +
  geom_ribbon(
    data = dens_sum,
    aes(x = x, ymin = q05, ymax = q95),
    alpha = 0.25
  ) +
  geom_line(
    data = dens_sum,
    aes(x = x, y = q50),
    linewidth = 1
  ) +
  labs(
    x = "Serial interval (days)",
    y = "Density",
    title = "Observed serial intervals with posterior predictive overlay",
    subtitle = "Histogram = observed onset-to-onset intervals; line/ribbon = fitted model posterior predictive density"
  ) +
  theme_minimal()

ggsave(
  file.path(out_dir, "serial_interval_overlay_density.png"),
  p_si_overlay,
  width = 8,
  height = 5,
  dpi = 300
)

print(p_si_overlay)

# Delay-distribution reporting estimates
# Outputs actual estimates and uncertainty for the reporting checklist.
report_dir <- file.path(out_dir, "delay_reporting_estimates")
dir.create(report_dir, showWarnings = FALSE, recursive = TRUE)

clean_x <- function(x) {
  x <- as.numeric(x)
  x[is.finite(x) & x > -900]
}

per_draw_distribution_stats <- function(mat) {
  bind_rows(lapply(seq_len(nrow(mat)), function(i) {
    x <- clean_x(mat[i, ])
    tibble(
      draw = i,
      mean = mean(x),
      median = median(x),
      q025 = unname(quantile(x, 0.025)),
      q25 = unname(quantile(x, 0.25)),
      q75 = unname(quantile(x, 0.75)),
      q975 = unname(quantile(x, 0.975))
    )
  }))
}

summarise_delay <- function(mat, label) {
  s <- per_draw_distribution_stats(mat)
  tibble(
    quantity = label,
    mean_estimate = median(s$mean),
    mean_cri_low = unname(quantile(s$mean, 0.025)),
    mean_cri_high = unname(quantile(s$mean, 0.975)),
    median_estimate = median(s$median),
    median_cri_low = unname(quantile(s$median, 0.025)),
    median_cri_high = unname(quantile(s$median, 0.975)),
    distribution_q025 = median(s$q025),
    distribution_q25 = median(s$q25),
    distribution_q75 = median(s$q75),
    distribution_q975 = median(s$q975)
  )
}

fmt <- function(x, digits = 2) sprintf(paste0("%.", digits, "f"), x)
fmt_cri <- function(est, lo, hi, digits = 2) {
  paste0(fmt(est, digits), " (95% CrI: ", fmt(lo, digits), "--", fmt(hi, digits), ")")
}
fmt_range <- function(lo, hi, digits = 2) paste0(fmt(lo, digits), "--", fmt(hi, digits))

# Distributional estimates from fitted population PMFs and pair-level latent draws.
gen_pair_mat <- gen_obs_mat[, sourced_idx, drop = FALSE]
delta_pair_mat <- delta_obs_mat[, sourced_idx, drop = FALSE]
si_pair_mat <- si_obs_mat[, sourced_idx, drop = FALSE]

weighted_quantile <- function(x, w, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)) {
  ok <- is.finite(x) & is.finite(w) & w >= 0
  x <- x[ok]; w <- w[ok]
  if (length(x) == 0 || sum(w) <= 0) return(rep(NA_real_, length(probs)))
  ord <- order(x); x <- x[ord]; w <- w[ord] / sum(w)
  cw <- cumsum(w)
  sapply(probs, function(p) x[which(cw >= p)[1]])
}

summarise_pmf_delay <- function(pmat, label) {
  grid_mid <- grid_day + 0.5 * dt_grid
  s <- bind_rows(lapply(seq_len(nrow(pmat)), function(i) {
    p <- as.numeric(pmat[i, ]); p[!is.finite(p) | p < 0] <- 0
    p <- if (sum(p) > 0) p / sum(p) else rep(1 / length(p), length(p))
    mu <- sum(grid_mid * p)
    qs <- weighted_quantile(grid_mid, p)
    tibble(
      draw = i,
      mean = mu,
      median = qs[3],
      q025 = qs[1],
      q25 = qs[2],
      q75 = qs[4],
      q975 = qs[5]
    )
  }))
  tibble(
    quantity = label,
    mean_estimate = median(s$mean),
    mean_cri_low = unname(quantile(s$mean, 0.025)),
    mean_cri_high = unname(quantile(s$mean, 0.975)),
    median_estimate = median(s$median),
    median_cri_low = unname(quantile(s$median, 0.025)),
    median_cri_high = unname(quantile(s$median, 0.975)),
    distribution_q025 = median(s$q025),
    distribution_q25 = median(s$q25),
    distribution_q75 = median(s$q75),
    distribution_q975 = median(s$q975)
  )
}

delay_estimates <- bind_rows(
  summarise_pmf_delay(p_inc_mat, "Incubation period"),
  summarise_pmf_delay(p_lp_mat, "Latent period"),
  summarise_pmf_delay(p_gen_mat, "Generation interval"),
  summarise_delay(si_pair_mat, "Serial interval"),
  summarise_delay(delta_pair_mat, "Transmission timing relative to source onset")
) %>%
  mutate(
    mean_with_cri = fmt_cri(mean_estimate, mean_cri_low, mean_cri_high),
    median_with_cri = fmt_cri(median_estimate, median_cri_low, median_cri_high),
    central_95_distributional_interval = fmt_range(distribution_q025, distribution_q975),
    iqr = fmt_range(distribution_q25, distribution_q75)
  )

presym_by_draw <- apply(delta_pair_mat, 1, function(x) mean(clean_x(x) < 0))
delta_le_1_by_draw <- apply(delta_pair_mat, 1, function(x) mean(clean_x(x) <= 1))
delta_le_2_by_draw <- apply(delta_pair_mat, 1, function(x) mean(clean_x(x) <= 2))

prob_estimates <- tibble(
  quantity = c(
    "Transmission before source symptom onset, Pr(Delta < 0)",
    "Transmission by 1 day after source symptom onset, Pr(Delta <= 1)",
    "Transmission by 2 days after source symptom onset, Pr(Delta <= 2)"
  ),
  estimate = c(median(presym_by_draw), median(delta_le_1_by_draw), median(delta_le_2_by_draw)),
  cri_low = c(unname(quantile(presym_by_draw, 0.025)), unname(quantile(delta_le_1_by_draw, 0.025)), unname(quantile(delta_le_2_by_draw, 0.025))),
  cri_high = c(unname(quantile(presym_by_draw, 0.975)), unname(quantile(delta_le_1_by_draw, 0.975)), unname(quantile(delta_le_2_by_draw, 0.975)))
) %>%
  mutate(
    estimate_percent = 100 * estimate,
    cri_low_percent = 100 * cri_low,
    cri_high_percent = 100 * cri_high,
    estimate_with_cri = paste0(fmt(estimate_percent, 1), "% (95% CrI: ", fmt(cri_low_percent, 1), "--", fmt(cri_high_percent, 1), "%)")
  )

write_csv(delay_estimates, file.path(report_dir, "delay_reporting_estimates.csv"))
write_csv(prob_estimates, file.path(report_dir, "presymptomatic_transmission_probabilities.csv"))

latex_delay_rows <- delay_estimates %>%
  transmute(row = paste0(quantity, " & ", median_with_cri, " & ", mean_with_cri, " & ", central_95_distributional_interval, " \\\\")) %>%
  pull(row)

latex_delay_table <- c(
  "\\begin{table}[ht]",
  "\\centering",
  "\\caption{Estimated epidemiological delay distributions from the within-host-informed event-time model. Median and mean estimates are posterior medians, with 95\\% credible intervals. The final column gives the posterior median central 95\\% distributional interval.}",
  "\\begin{tabular}{llll}",
  "\\hline",
  "Quantity & Median estimate & Mean estimate & Central 95\\% distributional interval \\\\",
  "\\hline",
  latex_delay_rows,
  "\\hline",
  "\\end{tabular}",
  "\\end{table}"
)
writeLines(latex_delay_table, file.path(report_dir, "delay_reporting_estimates.tex"))

latex_prob_rows <- prob_estimates %>%
  transmute(row = paste0(quantity, " & ", estimate_with_cri, " \\\\")) %>%
  pull(row)

latex_prob_table <- c(
  "\\begin{table}[ht]",
  "\\centering",
  "\\caption{Posterior estimates of transmission timing relative to source symptom onset.}",
  "\\begin{tabular}{ll}",
  "\\hline",
  "Quantity & Estimate \\\\",
  "\\hline",
  latex_prob_rows,
  "\\hline",
  "\\end{tabular}",
  "\\end{table}"
)
writeLines(latex_prob_table, file.path(report_dir, "presymptomatic_transmission_probabilities.tex"))

message("Paper-ready delay reporting files written to: ", report_dir)
print(delay_estimates %>% select(quantity, median_with_cri, mean_with_cri, central_95_distributional_interval))
print(prob_estimates %>% select(quantity, estimate_with_cri))

# 4-panel figure
library(readr)

final_dir <- file.path(out_dir, "final_four_panel")
dir.create(final_dir, showWarnings = FALSE, recursive = TRUE)

grid_mid <- grid_day + 0.5 * dt_grid
sourced_idx <- which(ll$source_idx > 0)
source_idx_for_pair <- ll$source_idx[sourced_idx]

normalize_p <- function(p) {
  p <- as.numeric(p)
  p[!is.finite(p) | p < 0] <- 0
  if (sum(p) <= 0) rep(1 / length(p), length(p)) else p / sum(p)
}

sample_from_pmf <- function(x, p, n) {
  sample(x, size = n, replace = TRUE, prob = normalize_p(p))
}

make_density_summary <- function(df, bw = "nrd0", n_grid = 512) {
  x_lo <- quantile(df$value, 0.001, na.rm = TRUE)
  x_hi <- quantile(df$value, 0.999, na.rm = TRUE)
  
  dens_by_draw <- df %>%
    group_by(draw) %>%
    group_split() %>%
    lapply(function(dd) {
      d <- density(
        dd$value,
        from = x_lo,
        to = x_hi,
        n = n_grid,
        bw = bw
      )
      
      tibble(
        draw = unique(dd$draw),
        x = d$x,
        density = d$y
      )
    }) %>%
    bind_rows()
  
  dens_by_draw %>%
    group_by(x) %>%
    summarise(
      q025 = quantile(density, 0.025, na.rm = TRUE),
      q05  = quantile(density, 0.05,  na.rm = TRUE),
      q50  = quantile(density, 0.50,  na.rm = TRUE),
      q95  = quantile(density, 0.95,  na.rm = TRUE),
      q975 = quantile(density, 0.975, na.rm = TRUE),
      .groups = "drop"
    )
}

set.seed(20260519)
n_samp_per_draw <- 5000

# Smooth fitted IP, GT, and SI from model-implied PMFs
smooth_samples <- bind_rows(lapply(seq_len(nrow(p_inc_mat)), function(i) {
  ip_source <- sample_from_pmf(grid_mid, p_inc_mat[i, ], n_samp_per_draw)
  ip_recipient <- sample_from_pmf(grid_mid, p_inc_mat[i, ], n_samp_per_draw)
  gt <- sample_from_pmf(grid_mid, p_gen_mat[i, ], n_samp_per_draw)
  
  # Model-implied serial interval:
  # SI = GT + IP_recipient - IP_source
  si <- gt + ip_recipient - ip_source
  
  bind_rows(
    tibble(draw = i, quantity = "Generation interval", value = gt),
    tibble(draw = i, quantity = "Incubation period", value = ip_recipient),
    tibble(draw = i, quantity = "Serial interval", value = si)
  )
}))

gt_sum <- make_density_summary(
  smooth_samples %>% filter(quantity == "Generation interval"),
  bw = "nrd0"
) %>% mutate(quantity = "Generation interval")

ip_sum <- make_density_summary(
  smooth_samples %>% filter(quantity == "Incubation period"),
  bw = "nrd0"
) %>% mutate(quantity = "Incubation period")

si_sum <- make_density_summary(
  smooth_samples %>% filter(quantity == "Serial interval"),
  bw = "nrd0"
) %>% mutate(quantity = "Serial interval")


# Delta: joint posterior from same draw and same pair
source_ip_mat <- inc_obs_mat[, source_idx_for_pair, drop = FALSE]
gt_pair_mat <- gen_obs_mat[, sourced_idx, drop = FALSE]

delta_samples <- bind_rows(lapply(seq_len(nrow(gt_pair_mat)), function(i) {
  pair_sample <- sample(seq_along(sourced_idx), size = n_samp_per_draw, replace = TRUE)
  
  source_ip <- as.numeric(source_ip_mat[i, pair_sample])
  gt_pair <- as.numeric(gt_pair_mat[i, pair_sample])
  
  tibble(
    draw = i,
    quantity = "Transmission timing \u03b4",
    value = gt_pair - source_ip
  )
})) %>%
  filter(is.finite(value), value > -20, value < 20)

delta_sum <- make_density_summary(
  delta_samples,
  bw = 0.35
) %>% mutate(quantity = "Transmission timing \u03b4")


# Combine and plot
panel_order <- c(
  "Generation interval",
  "Incubation period",
  "Serial interval",
  "Transmission timing \u03b4"
)

dens_summary_all <- bind_rows(gt_sum, ip_sum, si_sum, delta_sum) %>%
  mutate(quantity = factor(quantity, levels = panel_order))

vline_df <- tibble(
  quantity = factor(panel_order, levels = panel_order),
  xint = c(0, 0, 0, 0)
)

p4_final <- ggplot(dens_summary_all, aes(x = x)) +
  geom_ribbon(aes(ymin = q025, ymax = q975), fill = "grey60", alpha = 0.22) +
  geom_ribbon(aes(ymin = q05, ymax = q95), fill = "grey35", alpha = 0.28) +
  geom_line(aes(y = q50), linewidth = 1.2, color = "black") +
  geom_vline(
    data = vline_df,
    aes(xintercept = xint),
    linetype = 2,
    linewidth = 0.7,
    inherit.aes = FALSE
  ) +
  facet_wrap(~ quantity, scales = "free", ncol = 2) +
  labs(
    x = "days",
    y = "density"
  ) +
  theme_bw(base_size = 15) +
  theme(
    strip.text = element_text(face = "bold", size = 18),
    panel.grid.minor = element_line(linewidth = 0.25),
    panel.grid.major = element_line(linewidth = 0.35)
  )

ggsave(
  file.path(final_dir, "final_four_panel_clean_model_implied_si.png"),
  p4_final,
  width = 11,
  height = 8.5,
  dpi = 300
)

library(dplyr)
library(ggplot2)
library(patchwork)

make_panel <- function(quantity_name) {
  ggplot(
    filter(dens_summary_all, quantity == quantity_name),
    aes(x = x)
  ) +
    geom_ribbon(
      aes(ymin = q025, ymax = q975),
      fill = "grey60",
      alpha = 0.22
    ) +
    geom_ribbon(
      aes(ymin = q05, ymax = q95),
      fill = "grey35",
      alpha = 0.28
    ) +
    geom_line(
      aes(y = q50),
      linewidth = 1.2,
      colour = "black"
    ) +
    geom_vline(
      data = filter(vline_df, quantity == quantity_name),
      aes(xintercept = xint),
      linetype = 2,
      linewidth = 0.7,
      inherit.aes = FALSE
    ) +
    labs(
      x = "Days",
      y = "Density"
    ) +
    theme_bw(base_size = 15) +
    theme(
      panel.grid.minor = element_line(linewidth = 0.25),
      panel.grid.major = element_line(linewidth = 0.35)
    )
}

pA <- make_panel("Transmission timing δ")
pB <- make_panel("Incubation period")
pC <- make_panel("Generation interval")
pD <- make_panel("Serial interval")

p4_final <-
  (pA | pB) /
  (pC | pD) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(
      face = "bold",
      size = 18
    ),
    plot.tag.position = c(0, 1)
  )

ggsave(
  file.path(final_dir, "final_four_panel_clean_model_implied_si.png"),
  p4_final,
  width = 11,
  height = 8.5,
  dpi = 300
)

print(p4_final)

write_csv(
  dens_summary_all,
  file.path(final_dir, "final_four_panel_clean_density_summary.csv")
)

message("Saved final figure to: ",
        file.path(final_dir, "final_four_panel_clean_model_implied_si.png"))

# Charniga-style delay reporting outputs
report_dir <- file.path(out_dir, "charniga_delay_reporting")
dir.create(report_dir, showWarnings = FALSE, recursive = TRUE)

# Helpers
clean_x <- function(x) {
  x <- as.numeric(x)
  x[is.finite(x)]
}

fmt <- function(x, digits = 2) {
  sprintf(paste0("%.", digits, "f"), x)
}

fmt_cri <- function(est, lo, hi, digits = 2) {
  paste0(fmt(est, digits), " (95% CrI: ", fmt(lo, digits), "--", fmt(hi, digits), ")")
}

weighted_quantile <- function(x, w, probs = c(0.025, 0.25, 0.5, 0.75, 0.95, 0.975, 0.99)) {
  ok <- is.finite(x) & is.finite(w) & w >= 0
  x <- x[ok]
  w <- w[ok]
  
  if (length(x) == 0 || sum(w) <= 0) {
    return(rep(NA_real_, length(probs)))
  }
  
  ord <- order(x)
  x <- x[ord]
  w <- w[ord] / sum(w)
  cw <- cumsum(w)
  
  sapply(probs, function(p) x[which(cw >= p)[1]])
}

summarise_delay_draws <- function(per_draw, label) {
  tibble(
    quantity = label,
    
    mean = median(per_draw$mean, na.rm = TRUE),
    mean_l95 = unname(quantile(per_draw$mean, 0.025, na.rm = TRUE)),
    mean_u95 = unname(quantile(per_draw$mean, 0.975, na.rm = TRUE)),
    
    sd = median(per_draw$sd, na.rm = TRUE),
    sd_l95 = unname(quantile(per_draw$sd, 0.025, na.rm = TRUE)),
    sd_u95 = unname(quantile(per_draw$sd, 0.975, na.rm = TRUE)),
    
    median = median(per_draw$median, na.rm = TRUE),
    median_l95 = unname(quantile(per_draw$median, 0.025, na.rm = TRUE)),
    median_u95 = unname(quantile(per_draw$median, 0.975, na.rm = TRUE)),
    
    q025 = median(per_draw$q025, na.rm = TRUE),
    q25 = median(per_draw$q25, na.rm = TRUE),
    q75 = median(per_draw$q75, na.rm = TRUE),
    q95 = median(per_draw$q95, na.rm = TRUE),
    q975 = median(per_draw$q975, na.rm = TRUE),
    q99 = median(per_draw$q99, na.rm = TRUE)
  )
}

summarise_pmf_delay <- function(pmat, grid_day, dt_grid, label) {
  grid_mid <- grid_day + 0.5 * dt_grid
  
  per_draw <- bind_rows(lapply(seq_len(nrow(pmat)), function(i) {
    p <- as.numeric(pmat[i, ])
    p[!is.finite(p) | p < 0] <- 0
    
    if (sum(p) > 0) {
      p <- p / sum(p)
    } else {
      p <- rep(1 / length(p), length(p))
    }
    
    qs <- weighted_quantile(grid_mid, p)
    mu <- sum(grid_mid * p)
    sig <- sqrt(sum((grid_mid - mu)^2 * p))
    
    tibble(
      draw = i,
      mean = mu,
      sd = sig,
      q025 = qs[1],
      q25 = qs[2],
      median = qs[3],
      q75 = qs[4],
      q95 = qs[5],
      q975 = qs[6],
      q99 = qs[7]
    )
  }))
  
  summarise_delay_draws(per_draw, label)
}

summarise_empirical_delay <- function(mat, label) {
  per_draw <- bind_rows(lapply(seq_len(nrow(mat)), function(i) {
    x <- clean_x(mat[i, ])
    
    if (length(x) == 0) {
      return(tibble(
        draw = i,
        mean = NA_real_,
        sd = NA_real_,
        q025 = NA_real_,
        q25 = NA_real_,
        median = NA_real_,
        q75 = NA_real_,
        q95 = NA_real_,
        q975 = NA_real_,
        q99 = NA_real_
      ))
    }
    
    tibble(
      draw = i,
      mean = mean(x),
      sd = sd(x),
      q025 = unname(quantile(x, 0.025)),
      q25 = unname(quantile(x, 0.25)),
      median = median(x),
      q75 = unname(quantile(x, 0.75)),
      q95 = unname(quantile(x, 0.95)),
      q975 = unname(quantile(x, 0.975)),
      q99 = unname(quantile(x, 0.99))
    )
  }))
  
  summarise_delay_draws(per_draw, label)
}

# Main delay estimates
sourced_idx <- which(ll$source_idx > 0)

delay_estimates <- bind_rows(
  summarise_pmf_delay(
    p_inc_mat,
    grid_day,
    dt_grid,
    "Incubation period"
  ),
  summarise_pmf_delay(
    p_lp_mat,
    grid_day,
    dt_grid,
    "Latent period"
  ),
  summarise_pmf_delay(
    p_gen_mat,
    grid_day,
    dt_grid,
    "Generation interval"
  ),
  summarise_empirical_delay(
    si_obs_mat[, sourced_idx, drop = FALSE],
    "Serial interval"
  ),
  summarise_empirical_delay(
    delta_obs_mat[, sourced_idx, drop = FALSE],
    "Transmission timing relative to source symptom onset"
  )
) %>%
  mutate(
    mean_report = fmt_cri(mean, mean_l95, mean_u95),
    sd_report = fmt_cri(sd, sd_l95, sd_u95),
    median_report = fmt_cri(median, median_l95, median_u95),
    iqr = paste0(fmt(q25), "--", fmt(q75)),
    central_95_distribution = paste0(fmt(q025), "--", fmt(q975)),
    q95_report = fmt(q95),
    q99_report = fmt(q99)
  )

write_csv(
  delay_estimates,
  file.path(report_dir, "delay_distribution_estimates.csv")
)

# Transmission timing probabilities
delta_mat <- delta_obs_mat[, sourced_idx, drop = FALSE]

prob_draws <- bind_rows(
  tibble(
    draw = seq_len(nrow(delta_mat)),
    quantity = "Pr(delta < 0): transmission before source symptom onset",
    value = apply(delta_mat, 1, function(x) mean(clean_x(x) < 0))
  ),
  tibble(
    draw = seq_len(nrow(delta_mat)),
    quantity = "Pr(delta <= 1): transmission by 1 day after source symptom onset",
    value = apply(delta_mat, 1, function(x) mean(clean_x(x) <= 1))
  ),
  tibble(
    draw = seq_len(nrow(delta_mat)),
    quantity = "Pr(delta <= 2): transmission by 2 days after source symptom onset",
    value = apply(delta_mat, 1, function(x) mean(clean_x(x) <= 2))
  )
)

prob_estimates <- prob_draws %>%
  group_by(quantity) %>%
  summarise(
    estimate = median(value, na.rm = TRUE),
    l95 = unname(quantile(value, 0.025, na.rm = TRUE)),
    u95 = unname(quantile(value, 0.975, na.rm = TRUE)),
    estimate_report = paste0(
      fmt(100 * estimate, 1), "% (95% CrI: ",
      fmt(100 * l95, 1), "--", fmt(100 * u95, 1), "%)"
    ),
    .groups = "drop"
  )

write_csv(
  prob_draws,
  file.path(report_dir, "transmission_timing_probability_draws.csv")
)

write_csv(
  prob_estimates,
  file.path(report_dir, "transmission_timing_probabilities.csv")
)

# Diagnostics
available_vars <- colnames(as.matrix(fit$draws()))

wanted_scalar_vars <- c(
  "inc_mean",
  "gen_mean",
  "delta_mean",
  "si_mean",
  "inc_mean_mech",
  "gen_mean_mech",
  "delta_mean_mech",
  "p_delta_presym",
  "p_delta_by_day0",
  "p_delta_by_day1",
  "p_delta_by_day2",
  "sigma_rw",
  "k_nb",
  "log_lambda_sym",
  "log_lambda_inf",
  "log_r",
  "log_d",
  "m",
  "log_kappa"
)

scalar_vars <- wanted_scalar_vars[wanted_scalar_vars %in% available_vars]

if (length(scalar_vars) > 0) {
  diagnostics <- fit$summary(variables = scalar_vars)
  write_csv(
    diagnostics,
    file.path(report_dir, "parameter_diagnostics.csv")
  )
}

sampler_diag <- fit$sampler_diagnostics(format = "df")

diag_summary <- tibble(
  n_sampler_rows = nrow(sampler_diag),
  divergences = sum(sampler_diag$divergent__),
  max_treedepth = max(sampler_diag$treedepth__),
  mean_stepsize = mean(sampler_diag$stepsize__, na.rm = TRUE)
)

if (exists("diagnostics")) {
  diag_summary <- diag_summary %>%
    mutate(
      max_rhat = max(diagnostics$rhat, na.rm = TRUE),
      min_bulk_ess = min(diagnostics$ess_bulk, na.rm = TRUE),
      min_tail_ess = min(diagnostics$ess_tail, na.rm = TRUE)
    )
}

write_csv(
  diag_summary,
  file.path(report_dir, "diagnostics_summary.csv")
)

# Data/context summary
context_summary <- tibble(
  item = c(
    "Outbreak",
    "Cases included",
    "Linked transmission pairs",
    "Interval censoring",
    "Model type",
    "Delay quantities reported"
  ),
  value = c(
    "2018--2019 Epuyen Andes virus outbreak",
    as.character(nrow(ll)),
    as.character(length(sourced_idx)),
    "Latent infection and onset times estimated within exposure/onset windows",
    "Within-host-informed event-time model",
    "Incubation period; latent period; generation interval; serial interval; transmission timing relative to source symptom onset"
  )
)

write_csv(
  context_summary,
  file.path(report_dir, "data_context_summary.csv")
)

# Latex-ready table
latex_rows <- delay_estimates %>%
  transmute(
    row = paste0(
      quantity, " & ",
      median_report, " & ",
      mean_report, " & ",
      sd_report, " & ",
      central_95_distribution, " & ",
      q95_report, " & ",
      q99_report, " \\\\"
    )
  ) %>%
  pull(row)

latex_table <- c(
  "\\begin{table}[ht]",
  "\\centering",
  "\\caption{Estimated epidemiological delay distributions from the within-host-informed event-time model. Estimates are posterior medians with 95\\% credible intervals. Distributional intervals and quantiles summarize fitted individual-level variation.}",
  "\\begin{tabular}{lllllll}",
  "\\hline",
  "Quantity & Median & Mean & SD & Central 95\\% distribution & 95th percentile & 99th percentile \\\\",
  "\\hline",
  latex_rows,
  "\\hline",
  "\\end{tabular}",
  "\\end{table}"
)

writeLines(
  latex_table,
  file.path(report_dir, "delay_distribution_estimates.tex")
)

message("Done. Reporting outputs written to: ", report_dir)

print(delay_estimates %>%
        select(quantity, median_report, mean_report, sd_report, central_95_distribution))

print(prob_estimates %>%
        select(quantity, estimate_report))

