# ============================================================
# 2026 CRUISE OUTBREAK
# COMPATIBILITY WITH HIGH / LOW (R0, k) PARAMETERISATIONS
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

set.seed(20260826)

# ------------------------------------------------------------
# SETTINGS
# ------------------------------------------------------------

input_dir <- "~/Hantavirus-paper/TD results"

n_sims <- 10000L
max_cases <- 5000L
max_days <- 150

# ------------------------------------------------------------
# OBSERVED CRUISE ONSET DATES
# ------------------------------------------------------------

observed_dates <- as.Date(c(
  "2026-04-03",
  "2026-04-21",
  "2026-04-22",
  "2026-04-23",
  "2026-04-28",
  "2026-04-28",
  "2026-04-30",
  "2026-05-01",
  "2026-05-08",
  "2026-05-11",
  "2026-05-14",
  "2026-05-20",
  "2026-05-30"
))

stopifnot(length(observed_dates) == 13)

observed_onsets <- as.numeric(
  observed_dates - min(observed_dates)
)

# ------------------------------------------------------------
# TRANSMISSION PARAMETER SETS
# ------------------------------------------------------------

regimes <- tibble(
  regime = c(
    "High transmission",
    "Low transmission"
  ),
  R0 = c(
    2.12,
    0.2416027
  ),
  k = c(
    0.3157164,
    0.3998142
  )
)

# ------------------------------------------------------------
# FIND TIMING DRAWS
# ------------------------------------------------------------

delay_candidates <- c(
  file.path(
    input_dir,
    "tables",
    "joint_delay_draws_for_branching.csv"
  ),
  file.path(
    input_dir,
    "joint_delay_draws_for_branching.csv"
  ),
  file.path(
    input_dir,
    "posterior_pair_level_draws.csv"
  ),
  file.path(
    input_dir,
    "tables",
    "posterior_pair_level_draws.csv"
  )
)

delay_file <- delay_candidates[
  file.exists(delay_candidates)
][1]

if (length(delay_file) == 0 || is.na(delay_file)) {
  
  hits <- list.files(
    input_dir,
    pattern =
      "^(joint_delay_draws_for_branching|posterior_pair_level_draws)\\.csv$",
    recursive = TRUE,
    full.names = TRUE
  )
  
  if (length(hits) > 0) {
    delay_file <- hits[1]
  }
}

if (length(delay_file) == 0 || is.na(delay_file)) {
  stop(
    "Could not find timing draws under: ",
    input_dir,
    call. = FALSE
  )
}

message("Using timing file: ", delay_file)

timing <- read_csv(
  delay_file,
  show_col_types = FALSE
)

# ------------------------------------------------------------
# STANDARDISE TIMING COLUMNS
# ------------------------------------------------------------

if ("generation_time" %in% names(timing)) {
  
  timing$generation_time <-
    as.numeric(timing$generation_time)
  
} else if ("generation_interval" %in% names(timing)) {
  
  timing$generation_time <-
    as.numeric(timing$generation_interval)
  
} else {
  
  stop(
    "Need generation_time or generation_interval.",
    call. = FALSE
  )
}

if ("recipient_incubation" %in% names(timing)) {
  
  timing$incubation <-
    as.numeric(timing$recipient_incubation)
  
} else if ("incubation" %in% names(timing)) {
  
  timing$incubation <-
    as.numeric(timing$incubation)
  
} else {
  
  stop(
    "Need recipient_incubation or incubation.",
    call. = FALSE
  )
}

timing <- timing %>%
  filter(
    is.finite(generation_time),
    generation_time > 0,
    is.finite(incubation),
    incubation > 0
  )

generation_draws <- timing$generation_time
incubation_draws <- timing$incubation

# ------------------------------------------------------------
# SIMULATE ONE OUTBREAK
# ------------------------------------------------------------

simulate_outbreak <- function(
    R0,
    k
) {
  
  infection_time <- 0
  
  onset_time <- sample(
    incubation_draws,
    1
  )
  
  processed <- FALSE
  
  i <- 1L
  
  while (
    i <= length(infection_time) &&
    length(infection_time) < max_cases
  ) {
    
    if (!processed[i]) {
      
      z <- rnbinom(
        1,
        size = k,
        mu = R0
      )
      
      if (z > 0) {
        
        generation_times <- sample(
          generation_draws,
          z,
          replace = TRUE
        )
        
        child_infection_times <-
          infection_time[i] +
          generation_times
        
        keep <-
          is.finite(child_infection_times) &
          child_infection_times <= max_days
        
        child_infection_times <-
          child_infection_times[keep]
        
        if (length(child_infection_times) > 0) {
          
          child_incubations <- sample(
            incubation_draws,
            length(child_infection_times),
            replace = TRUE
          )
          
          child_onsets <-
            child_infection_times +
            child_incubations
          
          infection_time <- c(
            infection_time,
            child_infection_times
          )
          
          onset_time <- c(
            onset_time,
            child_onsets
          )
          
          processed <- c(
            processed,
            rep(
              FALSE,
              length(child_infection_times)
            )
          )
        }
      }
      
      processed[i] <- TRUE
    }
    
    i <- i + 1L
  }
  
  onset_time <- sort(onset_time)
  
  # Align first symptom onset to day zero
  onset_time - min(onset_time)
}

# ------------------------------------------------------------
# WEEKLY CURVE FUNCTION
#
# IMPORTANT:
# No hist(), so simulated values cannot fall outside breaks.
# ------------------------------------------------------------

n_weeks <- ceiling(max_days / 7) + 1L

make_weekly_curve <- function(onsets) {
  
  out <- integer(n_weeks)
  
  if (length(onsets) == 0) {
    return(out)
  }
  
  week_id <- floor(onsets / 7) + 1L
  
  # Anything beyond the chosen horizon is ignored
  week_id <- week_id[
    week_id >= 1 &
      week_id <= n_weeks
  ]
  
  if (length(week_id) > 0) {
    
    counts <- tabulate(
      week_id,
      nbins = n_weeks
    )
    
    out[] <- counts
  }
  
  out
}

# Observed curve on EXACTLY the same scale
observed_weekly <- make_weekly_curve(
  observed_onsets
)

message("\nObserved weekly epidemic curve:")
print(observed_weekly)

message(
  "\nObserved total cases in weekly vector: ",
  sum(observed_weekly)
)

# ------------------------------------------------------------
# COMPATIBILITY TEST
# ------------------------------------------------------------

compatibility_test <- function(
    R0,
    k,
    n_sims
) {
  
  message(
    "\nSimulating R0 = ",
    R0,
    ", k = ",
    k
  )
  
  X <- matrix(
    0,
    nrow = n_sims,
    ncol = n_weeks
  )
  
  total_cases <- integer(n_sims)
  
  for (s in seq_len(n_sims)) {
    
    onsets <- simulate_outbreak(
      R0 = R0,
      k = k
    )
    
    total_cases[s] <- length(onsets)
    
    X[s, ] <- make_weekly_curve(
      onsets
    )
  }
  
  # ----------------------------------------------------------
  # Remove weeks with zero variance in every simulation
  # ----------------------------------------------------------
  
  variances <- apply(
    X,
    2,
    var
  )
  
  keep_cols <- (
    is.finite(variances) &
      variances > 0
  )
  
  X_use <- X[, keep_cols, drop = FALSE]
  
  obs_use <- observed_weekly[
    keep_cols
  ]
  
  # ----------------------------------------------------------
  # Mean and covariance under this model
  # ----------------------------------------------------------
  
  mu <- colMeans(
    X_use
  )
  
  S <- cov(
    X_use
  )
  
  # Numerical stabilisation
  ridge <- max(
    1e-6,
    mean(diag(S)) * 1e-4
  )
  
  S <- S +
    diag(
      ridge,
      ncol(S)
    )
  
  Sinv <- solve(S)
  
  # ----------------------------------------------------------
  # Observed Mahalanobis distance
  # ----------------------------------------------------------
  
  obs_diff <- obs_use - mu
  
  d_obs <- as.numeric(
    t(obs_diff) %*%
      Sinv %*%
      obs_diff
  )
  
  # ----------------------------------------------------------
  # Simulated Mahalanobis distances
  # ----------------------------------------------------------
  
  d_sim <- apply(
    X_use,
    1,
    function(x) {
      
      dx <- x - mu
      
      as.numeric(
        t(dx) %*%
          Sinv %*%
          dx
      )
    }
  )
  
  # ----------------------------------------------------------
  # Empirical compatibility p-value
  # ----------------------------------------------------------
  
  p_compat <- mean(
    d_sim >= d_obs
  )
  
  tibble(
    R0 = R0,
    k = k,
    
    compatibility_p = p_compat,
    
    observed_distance = d_obs,
    
    simulated_distance_median =
      median(d_sim),
    
    simulated_distance_q95 =
      as.numeric(
        quantile(
          d_sim,
          0.95
        )
      ),
    
    p_at_least_13_cases =
      mean(
        total_cases >= 13
      )
  )
}

# ------------------------------------------------------------
# RUN BOTH PARAMETERISATIONS
# ------------------------------------------------------------

high_result <- compatibility_test(
  R0 = regimes$R0[1],
  k = regimes$k[1],
  n_sims = n_sims
) %>%
  mutate(
    regime = regimes$regime[1]
  )

low_result <- compatibility_test(
  R0 = regimes$R0[2],
  k = regimes$k[2],
  n_sims = n_sims
) %>%
  mutate(
    regime = regimes$regime[2]
  )

results <- bind_rows(
  high_result,
  low_result
) %>%
  select(
    regime,
    R0,
    k,
    compatibility_p,
    observed_distance,
    simulated_distance_median,
    simulated_distance_q95,
    p_at_least_13_cases
  ) %>%
  mutate(
    compatibility = case_when(
      compatibility_p < 0.05 ~
        "Poorly compatible",
      
      TRUE ~
        "Compatible"
    )
  )

# ------------------------------------------------------------
# RESULTS
# ------------------------------------------------------------

message(
  "\n========================================"
)

message(
  "2026 CRUISE OUTBREAK COMPATIBILITY"
)

message(
  "========================================\n"
)

print(
  results,
  width = Inf
)

# ------------------------------------------------------------
# SAVE
# ------------------------------------------------------------

write_csv(
  results,
  "cruise_epi_curve_compatibility.csv"
)

print(results, width = Inf)

cat("\n\n--- OBSERVED WEEKLY CURVE ---\n")
print(observed_weekly)

cat("\n\n--- FULL RESULTS AS DATA.FRAME ---\n")
print(as.data.frame(results), row.names = FALSE)

cat("\n\n--- EXACT VALUES ---\n")
for (i in seq_len(nrow(results))) {
  cat("\n", results$regime[i], "\n", sep = "")
  cat("R0 =", format(results$R0[i], digits = 10), "\n")
  cat("k =", format(results$k[i], digits = 10), "\n")
  cat("compatibility_p =", format(results$compatibility_p[i], digits = 10), "\n")
  cat("observed_distance =", format(results$observed_distance[i], digits = 10), "\n")
  cat("simulated_distance_median =", format(results$simulated_distance_median[i], digits = 10), "\n")
  cat("simulated_distance_q95 =", format(results$simulated_distance_q95[i], digits = 10), "\n")
  cat("p_at_least_13_cases =", format(results$p_at_least_13_cases[i], digits = 10), "\n")
  cat("compatibility =", results$compatibility[i], "\n")
}
