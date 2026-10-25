
# ============================================================
# ANDV hierarchical R/k model + 2026 cruise onset curve
# ROBUST VERSION
#
# Cruise curve is fitted with a multinomial timing likelihood
# conditional on its total number of descendant cases.
#
# This avoids the previous pathological solution:
# huge R_cruise + tiny daily-count overdispersion.
# ============================================================

library(cmdstanr)
library(posterior)
library(readr)
library(dplyr)
library(tibble)

# ------------------------------------------------------------
# 1. Published offspring data
# ------------------------------------------------------------

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

alonso_2020 <- c(
  1, 1, 0
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
  primary_count = c(
    2, 73, 1, 20, 92, 2, 131, 154, 1
  ),
  secondary_count = c(
    2, 0, 0, 0, 16, 1, 21, 38, 0
  )
)

study_names <- c(
  names(individual_list),
  aggregate_data$study
)

stopifnot(!anyDuplicated(study_names))

J_old <- length(study_names)

individual_data$study_id <- match(
  individual_data$study,
  study_names
)

aggregate_data$study_id <- match(
  aggregate_data$study,
  study_names
)

# ------------------------------------------------------------
# 2. Read time-delay posterior PMFs from Hanta_TD.R outputs
# ------------------------------------------------------------

td_dir <- path.expand("~/Hantavirus-paper/TD results")

grid_day <- read_csv(
  file.path(td_dir, "bp_day_grid.csv"),
  show_col_types = FALSE
)$day_since_infection

inc_draws <- as.matrix(
  read_csv(
    file.path(td_dir, "bp_draws_incubation_pmf_matrix.csv"),
    show_col_types = FALSE
  )
)

gen_draws <- as.matrix(
  read_csv(
    file.path(td_dir, "bp_draws_generation_pmf_matrix.csv"),
    show_col_types = FALSE
  )
)

stopifnot(
  ncol(inc_draws) == length(grid_day),
  ncol(gen_draws) == length(grid_day)
)

# Plug-in posterior median PMFs
w_inc0 <- apply(inc_draws, 2, median, na.rm = TRUE)
w_gen0 <- apply(gen_draws, 2, median, na.rm = TRUE)

w_inc0[!is.finite(w_inc0) | w_inc0 < 0] <- 0
w_gen0[!is.finite(w_gen0) | w_gen0 < 0] <- 0

stopifnot(sum(w_inc0) > 0, sum(w_gen0) > 0)

w_inc0 <- w_inc0 / sum(w_inc0)
w_gen0 <- w_gen0 / sum(w_gen0)

# ------------------------------------------------------------
# 3. Map PMFs to integer-day support used by Stan
# ------------------------------------------------------------

# Original Hanta_TD PMF bins have midpoint grid_day + 0.5.
inc_delay <- pmax(
  0L,
  as.integer(round(grid_day + 0.5))
)

D <- max(inc_delay) + 1L
w_inc <- numeric(D)

for (j in seq_along(w_inc0)) {
  w_inc[inc_delay[j] + 1L] <-
    w_inc[inc_delay[j] + 1L] + w_inc0[j]
}

w_inc <- w_inc / sum(w_inc)

# Generation interval must be positive in the renewal model.
gen_delay <- pmax(
  1L,
  as.integer(round(grid_day + 0.5))
)

G <- max(gen_delay)
w_gen <- numeric(G)

for (j in seq_along(w_gen0)) {
  w_gen[gen_delay[j]] <-
    w_gen[gen_delay[j]] + w_gen0[j]
}

w_gen <- w_gen / sum(w_gen)

# Trim only extreme tails for computational efficiency.
G_keep <- which(cumsum(w_gen) >= 0.999)[1]
w_gen <- w_gen[seq_len(G_keep)]
w_gen <- w_gen / sum(w_gen)
G <- length(w_gen)

D_keep <- which(cumsum(w_inc) >= 0.999)[1]
w_inc <- w_inc[seq_len(D_keep)]
w_inc <- w_inc / sum(w_inc)
D <- length(w_inc)

cat("Generation support: 1 to", G, "days\n")
cat("Incubation support: 0 to", D - 1L, "days\n")

# ------------------------------------------------------------
# 4. WHO 2026 cruise onset curve
# ------------------------------------------------------------

cruise_events <- tibble(
  date = as.Date(c(
    "2026-04-03",
    "2026-04-21",
    "2026-04-22",
    "2026-04-23",
    "2026-04-28",
    "2026-04-30",
    "2026-05-01",
    "2026-05-08",
    "2026-05-11",
    "2026-05-14",
    "2026-05-20",
    "2026-05-30"
  )),
  cases = c(
    1L,
    1L,
    1L,
    1L,
    2L,
    1L,
    1L,
    1L,
    1L,
    1L,
    1L,
    1L
  )
)

# Complete daily series, including all zero-case days.
all_dates <- tibble(
  date = seq(
    min(cruise_events$date),
    max(cruise_events$date),
    by = "day"
  )
)

cruise <- left_join(
  all_dates,
  cruise_events,
  by = "date"
)

cruise$cases[is.na(cruise$cases)] <- 0L
cruise$cases <- as.integer(cruise$cases)

stopifnot(
  sum(cruise$cases) == 13L,
  all(diff(cruise$date) == 1)
)

# ------------------------------------------------------------
# 5. Condition on first onset as index case
# ------------------------------------------------------------

first_case_day <- which(cruise$cases > 0)[1]

cruise <- cruise[
  first_case_day:nrow(cruise),
  ,
  drop = FALSE
]

y_curve_desc <- as.integer(cruise$cases)

# Remove exactly ONE conditioned-on index case.
stopifnot(y_curve_desc[1] >= 1L)
y_curve_desc[1] <- y_curve_desc[1] - 1L

T_curve <- length(y_curve_desc)

cat(
  "Cruise descendant cases in likelihood:",
  sum(y_curve_desc),
  "\n"
)

stopifnot(sum(y_curve_desc) == 12L)

# ------------------------------------------------------------
# 6. Distribution of index infection time given first onset
# ------------------------------------------------------------

PRE <- max(1L, D - 1L)

seed_profile <- numeric(PRE)

# If onset is day 0 and incubation is d,
# infection occurred on day -d.
for (d in 1:(D - 1L)) {
  latent_day <- -d
  p <- PRE + latent_day + 1L

  if (p >= 1L && p <= PRE) {
    seed_profile[p] <-
      seed_profile[p] + w_inc[d + 1L]
  }
}

stopifnot(sum(seed_profile) > 0)

seed_profile <- seed_profile / sum(seed_profile)

# ------------------------------------------------------------
# 7. Stan data
# ------------------------------------------------------------

stan_data_curve <- list(
  J_old = J_old,

  N_ind = nrow(individual_data),
  offspring = as.integer(individual_data$offspring),
  study_ind = as.integer(individual_data$study_id),

  N_agg = nrow(aggregate_data),
  primary_count = as.integer(aggregate_data$primary_count),
  secondary_count = as.integer(aggregate_data$secondary_count),
  study_agg = as.integer(aggregate_data$study_id),

  T_curve = T_curve,
  y_curve_desc = as.integer(y_curve_desc),

  PRE = PRE,

  G = G,
  w_gen = as.numeric(w_gen),

  D = D,
  w_inc = as.numeric(w_inc),

  seed_profile = as.numeric(seed_profile)
)

# ------------------------------------------------------------
# 8. Compile and fit
# ------------------------------------------------------------

stan_file <- "andv_hierarchical_nb_plus_cruise_curve_MULTINOMIAL.stan"

mod_curve <- cmdstan_model(
  stan_file,
  force_recompile = TRUE
)

fit_curve <- mod_curve$sample(
  data = stan_data_curve,
  seed = 2026,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 2000,
  iter_sampling = 4000,
  adapt_delta = 0.99,
  max_treedepth = 15,
  refresh = 500
)

# ------------------------------------------------------------
# 9. Main summaries
# ------------------------------------------------------------

main_summary <- fit_curve$summary(
  variables = c(
    "pooled_R",
    "pooled_k",
    "R_cruise",
    "tau_log_R"
  )
)

print(main_summary)

# Study-specific R values
R_vars <- paste0(
  "R[",
  seq_len(J_old + 1L),
  "]"
)

R_sum <- fit_curve$summary(
  variables = R_vars
)

R_sum$study <- c(
  study_names,
  "2026 cruise onset curve"
)

print(
  R_sum[, c(
    "study",
    "mean",
    "median",
    "q5",
    "q95",
    "rhat",
    "ess_bulk",
    "ess_tail"
  )]
)

# ------------------------------------------------------------
# 10. Numerical sanity checks
# ------------------------------------------------------------

check_draws <- posterior::as_draws_df(
  fit_curve$draws(
    variables = c(
      "pooled_R",
      "pooled_k",
      "R_cruise",
      "tau_log_R"
    )
  )
)

finite_check <- sapply(
  check_draws,
  function(x) sum(!is.finite(x))
)

print(finite_check)

if (any(finite_check > 0)) {
  warning(
    "Non-finite posterior draws remain. Inspect before interpreting."
  )
}

# ------------------------------------------------------------
# 11. Posterior predictive check for cruise curve
# ------------------------------------------------------------

yrep <- posterior::as_draws_matrix(
  fit_curve$draws("y_curve_rep")
)

prob_draws <- posterior::as_draws_matrix(
  fit_curve$draws("onset_prob")
)

stopifnot(
  all(is.finite(yrep)),
  all(is.finite(prob_draws))
)

curve_ppc <- tibble(
  date = cruise$date,
  observed_desc = y_curve_desc,

  prob_q025 = apply(
    prob_draws,
    2,
    quantile,
    probs = 0.025,
    na.rm = TRUE
  ),

  prob_q50 = apply(
    prob_draws,
    2,
    quantile,
    probs = 0.50,
    na.rm = TRUE
  ),

  prob_q975 = apply(
    prob_draws,
    2,
    quantile,
    probs = 0.975,
    na.rm = TRUE
  ),

  yrep_q025 = apply(
    yrep,
    2,
    quantile,
    probs = 0.025,
    na.rm = TRUE
  ),

  yrep_q50 = apply(
    yrep,
    2,
    quantile,
    probs = 0.50,
    na.rm = TRUE
  ),

  yrep_q975 = apply(
    yrep,
    2,
    quantile,
    probs = 0.975,
    na.rm = TRUE
  )
)

print(curve_ppc)

write_csv(
  curve_ppc,
  "cruise_curve_multinomial_ppc.csv"
)

# ------------------------------------------------------------
# 12. Save fitted object
# ------------------------------------------------------------

fit_curve$save_object(
  "fit_andv_hierarchical_plus_cruise_multinomial.rds"
)

cat("\nDONE\n")
