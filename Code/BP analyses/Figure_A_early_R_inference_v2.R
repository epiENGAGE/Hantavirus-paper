#!/usr/bin/env Rscript

# ============================================================
# Figure A: How quickly can the underlying transmission regime
# be inferred from accumulating DETECTED cases?
#
# This is designed to sit immediately before the current delayed-
# intervention Figure 3. It uses the same two transmission regimes,
# posterior-derived timing distributions, case-detection probability,
# and reporting delay assumptions as Figure 3.
#
# Key idea
# --------
# For each simulated outbreak, pause just BEFORE the response would
# begin at 1, 2, 3, 5, 10, 15, or 20 detected cases. At that time,
# use only detected cases/transmission links available by then.
#
# For a detected source i, y_i is the number of its direct offspring
# that have themselves been detected by the observation time t.
# Because recently infected offspring may not yet have been reported,
# the likelihood corrects for right-truncation/reporting delay using
#
#   a_i(t) = p_detect * P(child detection time <= t | source onset),
#
# where child detection time relative to source onset is sampled from
# serial interval + reporting delay. Under binomial thinning of an
# NB2 offspring distribution,
#
#   y_i ~ NegBin2(mean = R0 * a_i, shape = k).
#
# R0 and k are estimated with weakly regularising priors on log scale.
# For every run x detection threshold, uncertainty is approximated from
# the local posterior curvature (Laplace approximation) on log(R0), log(k).
# This produces a within-outbreak posterior interval for R0 and P(R0 > 1),
# rather than using the across-simulation spread of MAP point estimates as
# if it were inferential uncertainty. A supplementary plot does the same
# for k.
#
# IMPORTANT INTERPRETATION
# ------------------------
# Estimates are conditional on an outbreak having reached the specified
# number of detected cases. This is deliberate: rare low-R0 outbreaks
# that reach later thresholds are selected transmission-rich realizations.
# That selection is part of the public-health identifiability problem.
# ============================================================

suppressPackageStartupMessages({
  library(Rcpp)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(tibble)
  library(scales)
})

parse_args <- function() {
  x <- list(
    input = "~/Hantavirus-paper/TD results",
    out = "~/Hantavirus-paper/BP results/Fig A early inference",
    n_runs = 10000L,
    seed = 20260826L,
    reported_case_triggers = c(1L, 2L, 3L, 5L, 10L, 15L, 20L),
    max_days = 260,
    max_cases = 5000L,
    initial_infected = 1L,

    # Same detection assumptions as current Figure 3
    p_detect = 0.95,
    detection_delay_mean = 2.0,
    detection_delay_sd = 1.0,

    # Number of Monte Carlo draws used to approximate the probability
    # that an offspring case would have been detected by observation time.
    n_report_delay_draws = 100000L,

    # Weak priors for real-time inference. These deliberately span both
    # the subcritical and supercritical scenarios.
    prior_log_R_mean = log(1.0),
    prior_log_R_sd = 1.5,
    prior_log_k_mean = log(0.3),
    prior_log_k_sd = 1.25,

    # Do not draw inferential summaries when only a handful of simulations
    # reached a threshold. This prevents a line based on 1--3 extreme
    # low-transmission outbreaks from looking like a stable estimate.
    min_reached_for_plot = 20L,
    make_plots = TRUE
  )

  a <- commandArgs(trailingOnly = TRUE)
  if (length(a) == 0) return(x)

  i <- 1L
  while (i <= length(a)) {
    key <- a[[i]]

    if (key %in% c("--no_plots", "--no-plots")) {
      x$make_plots <- FALSE
      i <- i + 1L
      next
    }

    if (i == length(a)) stop("Missing value for ", key, call. = FALSE)
    val <- a[[i + 1L]]

    if (key == "--input") x$input <- val
    else if (key %in% c("--out", "--outdir")) x$out <- val
    else if (key %in% c("--n_runs", "--n-runs")) x$n_runs <- as.integer(val)
    else if (key == "--seed") x$seed <- as.integer(val)
    else if (key %in% c("--reported_case_triggers", "--reported-case-triggers")) {
      x$reported_case_triggers <- as.integer(strsplit(val, ",", fixed = TRUE)[[1]])
    }
    else if (key %in% c("--max_days", "--max-days")) x$max_days <- as.numeric(val)
    else if (key %in% c("--max_cases", "--max-cases")) x$max_cases <- as.integer(val)
    else if (key %in% c("--initial_infected", "--initial-infected")) x$initial_infected <- as.integer(val)
    else if (key == "--p_detect") x$p_detect <- as.numeric(val)
    else if (key %in% c("--detection_delay_mean", "--detection-delay-mean")) x$detection_delay_mean <- as.numeric(val)
    else if (key %in% c("--detection_delay_sd", "--detection-delay-sd")) x$detection_delay_sd <- as.numeric(val)
    else if (key %in% c("--n_report_delay_draws", "--n-report-delay-draws")) x$n_report_delay_draws <- as.integer(val)
    else if (key %in% c("--prior_log_R_mean", "--prior-log-R-mean")) x$prior_log_R_mean <- as.numeric(val)
    else if (key %in% c("--prior_log_R_sd", "--prior-log-R-sd")) x$prior_log_R_sd <- as.numeric(val)
    else if (key %in% c("--prior_log_k_mean", "--prior-log-k-mean")) x$prior_log_k_mean <- as.numeric(val)
    else if (key %in% c("--prior_log_k_sd", "--prior-log-k-sd")) x$prior_log_k_sd <- as.numeric(val)
    else if (key %in% c("--min_reached_for_plot", "--min-reached-for-plot")) x$min_reached_for_plot <- as.integer(val)
    else stop("Unknown argument: ", key, call. = FALSE)

    i <- i + 2L
  }

  x$reported_case_triggers <- sort(unique(
    x$reported_case_triggers[is.finite(x$reported_case_triggers) & x$reported_case_triggers >= 1]
  ))

  if (length(x$reported_case_triggers) == 0) stop("No valid reported-case triggers.", call. = FALSE)
  if (!is.finite(x$p_detect) || x$p_detect <= 0 || x$p_detect > 1) stop("p_detect must be in (0,1].", call. = FALSE)
  x
}

args <- parse_args()
set.seed(args$seed)

dir.create(args$out, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(args$out, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(args$out, "figures"), recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Read the same posterior timing draws used by Figure 3
# ------------------------------------------------------------
input_dir <- normalizePath(path.expand(args$input), winslash = "/", mustWork = FALSE)
parent_dir <- normalizePath(dirname(input_dir), winslash = "/", mustWork = FALSE)

delay_candidates <- unique(c(
  file.path(input_dir, "tables", "joint_delay_draws_for_branching.csv"),
  file.path(input_dir, "joint_delay_draws_for_branching.csv"),
  file.path(input_dir, "posterior_pair_level_draws.csv"),
  file.path(input_dir, "tables", "posterior_pair_level_draws.csv"),
  file.path(input_dir, "andes_wh_event_time_line_list_paper_model_FAST_COMPILES_outputs", "posterior_pair_level_draws.csv"),
  file.path(input_dir, "andes_wh_event_time_line_list_paper_model_FAST_COMPILES_outputs", "tables", "posterior_pair_level_draws.csv"),
  file.path(parent_dir, "tables", "joint_delay_draws_for_branching.csv"),
  file.path(parent_dir, "joint_delay_draws_for_branching.csv"),
  file.path(parent_dir, "posterior_pair_level_draws.csv"),
  file.path(parent_dir, "tables", "posterior_pair_level_draws.csv")
))

delay_file <- delay_candidates[file.exists(delay_candidates)][1]

if (is.na(delay_file) || length(delay_file) == 0) {
  recursive_hits <- unique(c(
    list.files(input_dir,
      pattern = "^(joint_delay_draws_for_branching|posterior_pair_level_draws)\\.csv$",
      recursive = TRUE, full.names = TRUE
    ),
    list.files(parent_dir,
      pattern = "^(joint_delay_draws_for_branching|posterior_pair_level_draws)\\.csv$",
      recursive = TRUE, full.names = TRUE
    )
  ))
  delay_file <- recursive_hits[file.exists(recursive_hits)][1]
}

if (is.na(delay_file) || length(delay_file) == 0) {
  stop(
    "Could not find joint_delay_draws_for_branching.csv or posterior_pair_level_draws.csv under: ",
    input_dir,
    "\nPass the Figure-3 timing-output directory with --input.",
    call. = FALSE
  )
}

message("Using posterior timing file: ", normalizePath(delay_file, winslash = "/"))
pair_draws_raw <- readr::read_csv(delay_file, show_col_types = FALSE)

if (!"generation_time" %in% names(pair_draws_raw)) {
  if ("generation_interval" %in% names(pair_draws_raw)) {
    pair_draws_raw <- pair_draws_raw %>% mutate(generation_time = as.numeric(generation_interval))
  } else {
    stop("Timing file must contain generation_time or generation_interval.", call. = FALSE)
  }
}

if (!"recipient_incubation" %in% names(pair_draws_raw)) {
  if ("incubation" %in% names(pair_draws_raw)) {
    pair_draws_raw <- pair_draws_raw %>% mutate(recipient_incubation = as.numeric(incubation))
  } else {
    stop("Timing file must contain recipient_incubation or incubation.", call. = FALSE)
  }
}

if (!"delta" %in% names(pair_draws_raw)) {
  stop("Timing file must contain delta (transmission timing relative to source onset).", call. = FALSE)
}

pair_draws_raw <- pair_draws_raw %>%
  mutate(
    generation_time = as.numeric(generation_time),
    recipient_incubation = as.numeric(recipient_incubation),
    delta = as.numeric(delta),
    source_incubation = if ("source_incubation" %in% names(.)) {
      as.numeric(source_incubation)
    } else {
      generation_time - delta
    },
    serial_interval = if ("serial_interval" %in% names(.)) {
      as.numeric(serial_interval)
    } else {
      delta + recipient_incubation
    }
  ) %>%
  filter(
    is.finite(generation_time), generation_time > 0,
    is.finite(recipient_incubation), recipient_incubation > 0,
    is.finite(source_incubation), source_incubation > 0,
    is.finite(delta),
    is.finite(serial_interval)
  )

if (nrow(pair_draws_raw) == 0) stop("No usable joint timing draws found.", call. = FALSE)

readr::write_csv(
  pair_draws_raw,
  file.path(args$out, "tables", "joint_delay_draws_used.csv")
)

generation_draws <- pair_draws_raw$generation_time
recipient_incubation_draws <- pair_draws_raw$recipient_incubation
source_incubation_draws <- pair_draws_raw$source_incubation
delta_draws <- pair_draws_raw$delta
serial_interval_draws <- pair_draws_raw$serial_interval

# Child DETECTION time relative to source symptom onset is approximately
# serial interval + child reporting delay. This is the right-truncation
# distribution used in the offspring likelihood.
gamma_from_mean_sd <- function(n, mean, sd) {
  if (!is.finite(mean) || mean <= 0) return(rep(0, n))
  if (!is.finite(sd) || sd <= 0) return(rep(mean, n))
  shape <- mean^2 / sd^2
  scale <- sd^2 / mean
  rgamma(n, shape = shape, scale = scale)
}

set.seed(args$seed + 71L)
si_sample <- sample(serial_interval_draws, args$n_report_delay_draws, replace = TRUE)
report_delay_sample <- gamma_from_mean_sd(
  args$n_report_delay_draws,
  args$detection_delay_mean,
  args$detection_delay_sd
)
child_detection_lag_draws <- sort(si_sample + report_delay_sample)
child_detection_lag_draws <- child_detection_lag_draws[is.finite(child_detection_lag_draws)]

if (length(child_detection_lag_draws) < 1000) {
  stop("Too few valid child-detection-lag draws.", call. = FALSE)
}

# ------------------------------------------------------------
# Transmission regimes: match the current Figure 3 script
# ------------------------------------------------------------
transmission_regimes <- tibble::tribble(
  ~transmission_regime, ~transmission_label, ~R0,       ~k,
  "supercritical",      "High transmission", 2.12,      0.3157164,
  "subcritical",        "Low transmission",  0.2416027, 0.3998142
)

# ------------------------------------------------------------
# C++: simulate uncontrolled potential trees and extract the
# information available at each detected-case threshold.
#
# No intervention acts BEFORE the threshold. This exactly matches
# the pre-response portion of the current Figure 3 analysis; the
# current Figure 3 then turns intervention on at that threshold.
# ------------------------------------------------------------
Rcpp::sourceCpp(code = '
#include <Rcpp.h>
#include <algorithm>
#include <numeric>
#include <cmath>
using namespace Rcpp;

inline double runif01() { return R::runif(0.0, 1.0); }

inline double gamma_delay_one(double mean, double sd, double min_value = 0.0) {
  if (!R_finite(mean) || mean <= 0.0) return min_value;
  if (!R_finite(sd) || sd <= 0.0) return std::max(mean, min_value);
  double shape = (mean * mean) / (sd * sd);
  double scale = (sd * sd) / mean;
  double y = R::rgamma(shape, scale);
  if (!R_finite(y) || y < min_value) y = min_value;
  return y;
}

inline int sample_index(int n) {
  int idx = (int)std::floor(R::runif(0.0, (double)n));
  if (idx < 0) idx = 0;
  if (idx >= n) idx = n - 1;
  return idx;
}

inline double sample_from_vec(const NumericVector& x) {
  return x[sample_index(x.size())];
}

inline bool finite_num(double x) { return R_finite(x); }

struct JointDelayDraw {
  double generation_time;
  double recipient_incubation;
};

inline int sample_joint_index_conditioned_on_source_ip(
    const NumericVector& source_incubation_draws,
    double source_incubation) {
  int n = source_incubation_draws.size();
  if (n <= 1 || !R_finite(source_incubation)) return sample_index(n);

  int n_try = std::min(64, n);
  int best_idx = sample_index(n);
  double best_dist = R_PosInf;
  for (int tt = 0; tt < n_try; ++tt) {
    int idx = sample_index(n);
    double x = source_incubation_draws[idx];
    if (!R_finite(x)) continue;
    double dist = std::fabs(x - source_incubation);
    if (dist < best_dist) {
      best_dist = dist;
      best_idx = idx;
    }
  }
  return best_idx;
}

inline JointDelayDraw sample_joint_delay(
    const NumericVector& generation_draws,
    const NumericVector& recipient_incubation_draws,
    const NumericVector& source_incubation_draws,
    const NumericVector& delta_draws,
    double source_incubation,
    double max_days) {

  int n = generation_draws.size();

  for (int guard = 0; guard < 200; ++guard) {
    int idx = sample_joint_index_conditioned_on_source_ip(
      source_incubation_draws, source_incubation
    );
    double inc = recipient_incubation_draws[idx];
    double delta = delta_draws[idx];
    double G = source_incubation + delta;

    if (R_finite(G) && G > 0.0 && G < max_days && R_finite(inc) && inc > 0.0) {
      JointDelayDraw out;
      out.generation_time = G;
      out.recipient_incubation = inc;
      return out;
    }
  }

  for (int guard = 0; guard < 200; ++guard) {
    int idx = sample_index(n);
    double G = generation_draws[idx];
    double inc = recipient_incubation_draws[idx];
    if (R_finite(G) && G > 0.0 && G < max_days && R_finite(inc) && inc > 0.0) {
      JointDelayDraw out;
      out.generation_time = G;
      out.recipient_incubation = inc;
      return out;
    }
  }

  JointDelayDraw out;
  out.generation_time = 1.0;
  out.recipient_incubation = 21.0;
  return out;
}

inline double empirical_cdf_sorted(const NumericVector& sorted_x, double x) {
  if (!R_finite(x)) return 0.0;
  NumericVector::const_iterator it = std::upper_bound(sorted_x.begin(), sorted_x.end(), x);
  double n_le = (double)std::distance(sorted_x.begin(), it);
  return n_le / (double)sorted_x.size();
}

// [[Rcpp::export]]
DataFrame simulate_detection_snapshots_cpp(
    int n_runs,
    IntegerVector reported_case_triggers,
    NumericVector generation_draws,
    NumericVector recipient_incubation_draws,
    NumericVector source_incubation_draws,
    NumericVector delta_draws,
    NumericVector child_detection_lag_draws,
    int max_cases,
    double max_days,
    int initial_infected,
    double R0,
    double k,
    double p_detect,
    double detection_delay_mean,
    double detection_delay_sd,
    int seed) {

  Environment base("package:base");
  Function set_seed = base["set.seed"];
  set_seed(seed);

  std::vector<int> out_run_id;
  std::vector<int> out_trigger;
  std::vector<double> out_observation_time;
  std::vector<int> out_source_rank;
  std::vector<int> out_observed_offspring;
  std::vector<double> out_ascertainment;
  std::vector<double> out_source_onset;
  std::vector<double> out_source_detection;
  std::vector<int> out_source_generation;

  int max_trigger = 0;
  for (int j = 0; j < reported_case_triggers.size(); ++j) {
    max_trigger = std::max(max_trigger, reported_case_triggers[j]);
  }

  out_run_id.reserve(n_runs * max_trigger);
  out_trigger.reserve(n_runs * max_trigger);
  out_observation_time.reserve(n_runs * max_trigger);
  out_source_rank.reserve(n_runs * max_trigger);
  out_observed_offspring.reserve(n_runs * max_trigger);
  out_ascertainment.reserve(n_runs * max_trigger);
  out_source_onset.reserve(n_runs * max_trigger);
  out_source_detection.reserve(n_runs * max_trigger);
  out_source_generation.reserve(n_runs * max_trigger);

  for (int rr = 0; rr < n_runs; ++rr) {
    std::vector<int> infector; infector.reserve(max_cases);
    std::vector<int> generation; generation.reserve(max_cases);
    std::vector<double> infection_time; infection_time.reserve(max_cases);
    std::vector<double> incubation_period; incubation_period.reserve(max_cases);
    std::vector<double> symptom_onset; symptom_onset.reserve(max_cases);
    std::vector<int> detected; detected.reserve(max_cases);
    std::vector<double> detection_time; detection_time.reserve(max_cases);
    std::vector<int> processed; processed.reserve(max_cases);

    auto add_case = [&](int inf, int gen, double itime, double inc) {
      double onset = itime + inc;
      int det = (runif01() < p_detect) ? 1 : 0;
      double dtime = R_PosInf;
      if (det == 1) {
        dtime = onset + gamma_delay_one(detection_delay_mean, detection_delay_sd, 0.0);
      }

      infector.push_back(inf);
      generation.push_back(gen);
      infection_time.push_back(itime);
      incubation_period.push_back(inc);
      symptom_onset.push_back(onset);
      detected.push_back(det);
      detection_time.push_back(dtime);
      processed.push_back(0);

      return (int)infection_time.size() - 1;
    };

    for (int ii = 0; ii < initial_infected; ++ii) {
      double inc0 = sample_from_vec(recipient_incubation_draws);
      add_case(NA_INTEGER, 0, 0.0, inc0);
    }

    double nb_prob = k / (k + R0);

    while ((int)infection_time.size() < max_cases) {
      int todo = -1;
      double best_t = R_PosInf;

      for (int ii = 0; ii < (int)infection_time.size(); ++ii) {
        if (processed[ii] == 0 && infection_time[ii] < best_t) {
          best_t = infection_time[ii];
          todo = ii;
        }
      }

      if (todo < 0) break;
      processed[todo] = 1;

      int z0 = (int)R::rnbinom(k, nb_prob);

      for (int zz = 0; zz < z0; ++zz) {
        if ((int)infection_time.size() >= max_cases) break;

        JointDelayDraw dd = sample_joint_delay(
          generation_draws,
          recipient_incubation_draws,
          source_incubation_draws,
          delta_draws,
          incubation_period[todo],
          max_days
        );

        double child_time = infection_time[todo] + dd.generation_time;
        if (!finite_num(child_time) || child_time > max_days) continue;

        add_case(todo, generation[todo] + 1, child_time, dd.recipient_incubation);
      }
    }

    int n_cases = (int)infection_time.size();

    std::vector< std::pair<double,int> > detected_order;
    detected_order.reserve(n_cases);
    for (int ii = 0; ii < n_cases; ++ii) {
      if (detected[ii] == 1 && finite_num(detection_time[ii])) {
        detected_order.push_back(std::make_pair(detection_time[ii], ii));
      }
    }
    std::sort(detected_order.begin(), detected_order.end());

    for (int tt = 0; tt < reported_case_triggers.size(); ++tt) {
      int trigger = reported_case_triggers[tt];
      if ((int)detected_order.size() < trigger) continue;

      double obs_time = detected_order[trigger - 1].first;

      // All cases detected by the threshold time (continuous delays make ties rare).
      std::vector<int> observed_sources;
      observed_sources.reserve(trigger);
      for (int jj = 0; jj < (int)detected_order.size(); ++jj) {
        if (detected_order[jj].first <= obs_time + 1e-10) {
          observed_sources.push_back(detected_order[jj].second);
        } else {
          break;
        }
      }

      for (int ss = 0; ss < (int)observed_sources.size(); ++ss) {
        int source = observed_sources[ss];
        int y = 0;

        for (int cc = 0; cc < n_cases; ++cc) {
          if (infector[cc] != source) continue;
          if (detected[cc] == 1 && finite_num(detection_time[cc]) && detection_time[cc] <= obs_time + 1e-10) {
            y++;
          }
        }

        double followup_from_source_onset = obs_time - symptom_onset[source];
        double ascertainment = p_detect * empirical_cdf_sorted(
          child_detection_lag_draws,
          followup_from_source_onset
        );

        // Avoid exactly zero mean in numerical fitting. A source with essentially
        // no opportunity for an offspring to have been detected contributes
        // almost no information anyway.
        ascertainment = std::max(ascertainment, 1e-10);

        out_run_id.push_back(rr + 1);
        out_trigger.push_back(trigger);
        out_observation_time.push_back(obs_time);
        out_source_rank.push_back(ss + 1);
        out_observed_offspring.push_back(y);
        out_ascertainment.push_back(ascertainment);
        out_source_onset.push_back(symptom_onset[source]);
        out_source_detection.push_back(detection_time[source]);
        out_source_generation.push_back(generation[source]);
      }
    }
  }

  return DataFrame::create(
    Named("run_id") = out_run_id,
    Named("reported_case_trigger") = out_trigger,
    Named("observation_time") = out_observation_time,
    Named("source_rank") = out_source_rank,
    Named("observed_offspring") = out_observed_offspring,
    Named("ascertainment_probability") = out_ascertainment,
    Named("source_onset") = out_source_onset,
    Named("source_detection") = out_source_detection,
    Named("source_generation") = out_source_generation
  );
}
')

# ------------------------------------------------------------
# Simulate snapshots under each true regime
# ------------------------------------------------------------
message("Simulating early-outbreak observation snapshots...")
t0 <- Sys.time()

snapshots_list <- vector("list", nrow(transmission_regimes))
for (ii in seq_len(nrow(transmission_regimes))) {
  reg <- transmission_regimes[ii, ]
  regime_seed <- as.integer(args$seed + ii * 1000003L)

  message(
    "  ", reg$transmission_label,
    " (true R0=", reg$R0, ", true k=", reg$k, ")"
  )

  snapshots_list[[ii]] <- simulate_detection_snapshots_cpp(
    n_runs = args$n_runs,
    reported_case_triggers = args$reported_case_triggers,
    generation_draws = generation_draws,
    recipient_incubation_draws = recipient_incubation_draws,
    source_incubation_draws = source_incubation_draws,
    delta_draws = delta_draws,
    child_detection_lag_draws = child_detection_lag_draws,
    max_cases = args$max_cases,
    max_days = args$max_days,
    initial_infected = args$initial_infected,
    R0 = reg$R0,
    k = reg$k,
    p_detect = args$p_detect,
    detection_delay_mean = args$detection_delay_mean,
    detection_delay_sd = args$detection_delay_sd,
    seed = regime_seed
  ) %>%
    as_tibble() %>%
    mutate(
      transmission_regime = reg$transmission_regime,
      transmission_label = reg$transmission_label,
      true_R0 = reg$R0,
      true_k = reg$k,
      regime_seed = regime_seed
    )
}

snapshots <- bind_rows(snapshots_list)

message(
  "Simulation time: ",
  round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 2),
  " seconds"
)

readr::write_csv(
  snapshots,
  file.path(args$out, "tables", "early_inference_source_level_snapshots.csv")
)

# Number of simulations reaching each detected-case threshold.
# Each reached run contributes one row for source_rank == 1, so count there.
reached_summary <- snapshots %>%
  filter(source_rank == 1) %>%
  count(
    transmission_regime, transmission_label, true_R0, true_k,
    reported_case_trigger,
    name = "n_reached"
  ) %>%
  tidyr::complete(
    nesting(transmission_regime, transmission_label, true_R0, true_k),
    reported_case_trigger = args$reported_case_triggers,
    fill = list(n_reached = 0L)
  ) %>%
  mutate(
    n_simulated = args$n_runs,
    proportion_reached = n_reached / n_simulated
  )

readr::write_csv(
  reached_summary,
  file.path(args$out, "tables", "early_inference_threshold_reached.csv")
)

# ------------------------------------------------------------
# Real-time NB2 posterior approximation at each run x trigger
# ------------------------------------------------------------
# We retain the same right-truncation-aware likelihood as v1:
#   y_i ~ NB2(mean = R0 * a_i, shape = k)
# where a_i is the probability that an offspring of source i would already
# have been detected by the observation time.
#
# The key change is that we now retain WITHIN-OUTBREAK uncertainty. We find
# the posterior mode on log(R0), log(k), calculate the Hessian there, and use
# a Laplace approximation. This gives an approximate posterior SD for log R,
# a 95% credible interval for R0, and P(R0 > 1) for every run x threshold.
#
# This directly fixes the main problem in v1: its ribbon was the 2.5--97.5%
# spread of MAP estimates ACROSS selected simulations. It was not a credible
# interval for R0 and therefore could misleadingly sit entirely above 1.

safe_inverse_2x2 <- function(H) {
  if (any(!is.finite(H)) || !all(dim(H) == c(2, 2))) return(NULL)
  H <- (H + t(H)) / 2
  ee <- eigen(H, symmetric = TRUE)
  # Numerical floor. Priors should make H positive definite, but tiny early
  # samples can still produce nearly flat directions, especially for k.
  vals <- pmax(ee$values, 1e-6)
  ee$vectors %*% diag(1 / vals, 2, 2) %*% t(ee$vectors)
}

fit_one_snapshot <- function(dat) {
  y <- dat$observed_offspring
  a <- pmax(dat$ascertainment_probability, 1e-10)

  R_start <- sum(y) / sum(a)
  if (!is.finite(R_start) || R_start <= 0) R_start <- exp(args$prior_log_R_mean)
  R_start <- min(max(R_start, 0.02), 8)

  objective <- function(theta) {
    log_R <- theta[1]
    log_k <- theta[2]
    R0 <- exp(log_R)
    k <- exp(log_k)
    mu <- pmax(R0 * a, 1e-12)

    loglik <- sum(stats::dnbinom(y, size = k, mu = mu, log = TRUE))
    logprior <- stats::dnorm(
      log_R, mean = args$prior_log_R_mean,
      sd = args$prior_log_R_sd, log = TRUE
    ) + stats::dnorm(
      log_k, mean = args$prior_log_k_mean,
      sd = args$prior_log_k_sd, log = TRUE
    )
    -(loglik + logprior)
  }

  fit <- tryCatch(
    stats::optim(
      par = c(log(R_start), args$prior_log_k_mean),
      fn = objective,
      method = "L-BFGS-B",
      lower = c(log(0.005), log(0.005)),
      upper = c(log(10), log(10)),
      hessian = TRUE,
      control = list(maxit = 250)
    ),
    error = function(e) NULL
  )

  if (is.null(fit) || fit$convergence != 0 || any(!is.finite(fit$par))) {
    return(tibble(
      R_map = NA_real_, k_map = NA_real_,
      R_post_median = NA_real_, R_post_lower = NA_real_, R_post_upper = NA_real_,
      prob_R_gt1 = NA_real_,
      k_post_median = NA_real_, k_post_lower = NA_real_, k_post_upper = NA_real_,
      logR_post_sd = NA_real_, logk_post_sd = NA_real_,
      fit_converged = FALSE, posterior_ok = FALSE,
      n_sources = length(y),
      observed_transmission_links = sum(y),
      effective_source_exposure = sum(a)
    ))
  }

  V <- safe_inverse_2x2(fit$hessian)
  posterior_ok <- !is.null(V) && all(is.finite(V)) && V[1,1] > 0 && V[2,2] > 0

  if (!posterior_ok) {
    return(tibble(
      R_map = exp(fit$par[1]), k_map = exp(fit$par[2]),
      R_post_median = NA_real_, R_post_lower = NA_real_, R_post_upper = NA_real_,
      prob_R_gt1 = NA_real_,
      k_post_median = NA_real_, k_post_lower = NA_real_, k_post_upper = NA_real_,
      logR_post_sd = NA_real_, logk_post_sd = NA_real_,
      fit_converged = TRUE, posterior_ok = FALSE,
      n_sources = length(y),
      observed_transmission_links = sum(y),
      effective_source_exposure = sum(a)
    ))
  }

  sd_logR <- sqrt(V[1,1])
  sd_logk <- sqrt(V[2,2])
  z975 <- stats::qnorm(0.975)

  # Under the Gaussian approximation on log scale, exp(mode_logR) is the
  # approximate posterior median. Importantly, the interval is free to cross 1.
  R_med <- exp(fit$par[1])
  R_lo <- exp(fit$par[1] - z975 * sd_logR)
  R_hi <- exp(fit$par[1] + z975 * sd_logR)
  p_gt1 <- 1 - stats::pnorm((0 - fit$par[1]) / sd_logR)

  k_med <- exp(fit$par[2])
  k_lo <- exp(fit$par[2] - z975 * sd_logk)
  k_hi <- exp(fit$par[2] + z975 * sd_logk)

  tibble(
    R_map = exp(fit$par[1]), k_map = exp(fit$par[2]),
    R_post_median = R_med,
    R_post_lower = R_lo,
    R_post_upper = R_hi,
    prob_R_gt1 = p_gt1,
    k_post_median = k_med,
    k_post_lower = k_lo,
    k_post_upper = k_hi,
    logR_post_sd = sd_logR,
    logk_post_sd = sd_logk,
    fit_converged = TRUE,
    posterior_ok = TRUE,
    n_sources = length(y),
    observed_transmission_links = sum(y),
    effective_source_exposure = sum(a)
  )
}

message("Fitting R0 and k with within-outbreak posterior uncertainty...")
t1 <- Sys.time()

fits <- snapshots %>%
  group_by(
    transmission_regime, transmission_label, true_R0, true_k,
    regime_seed, run_id, reported_case_trigger, observation_time
  ) %>%
  group_modify(~ fit_one_snapshot(.x)) %>%
  ungroup()

message(
  "Inference time: ",
  round(as.numeric(difftime(Sys.time(), t1, units = "secs")), 2),
  " seconds"
)

readr::write_csv(
  fits,
  file.path(args$out, "tables", "early_inference_run_level_posterior.csv")
)

# ------------------------------------------------------------
# Summaries
# ------------------------------------------------------------
qfun <- function(x, p) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  as.numeric(stats::quantile(x, probs = p, names = FALSE, type = 8))
}

summary_estimates <- fits %>%
  filter(fit_converged, posterior_ok) %>%
  group_by(
    transmission_regime, transmission_label, true_R0, true_k,
    reported_case_trigger
  ) %>%
  summarise(
    n_reached = n(),

    # Central curve: median of each outbreak's posterior median.
    R_median = median(R_post_median, na.rm = TRUE),

    # Ribbon: a TYPICAL within-outbreak 95% credible interval, summarized by
    # the median lower and upper credible limit across reached simulations.
    # This is deliberately not the spread of point estimates across runs.
    R_typical_lower = median(R_post_lower, na.rm = TRUE),
    R_typical_upper = median(R_post_upper, na.rm = TRUE),

    # Across-outbreak variability in the posterior medians is retained in the
    # table as a separate quantity rather than being mislabeled uncertainty.
    R_median_q25 = qfun(R_post_median, 0.25),
    R_median_q75 = qfun(R_post_median, 0.75),

    # Operationally interpretable inferential quantity.
    median_prob_R_gt1 = median(prob_R_gt1, na.rm = TRUE),
    prob_R_gt1_q25 = qfun(prob_R_gt1, 0.25),
    prob_R_gt1_q75 = qfun(prob_R_gt1, 0.75),

    k_median = median(k_post_median, na.rm = TRUE),
    k_typical_lower = median(k_post_lower, na.rm = TRUE),
    k_typical_upper = median(k_post_upper, na.rm = TRUE),

    median_observed_links = median(observed_transmission_links, na.rm = TRUE),
    median_effective_source_exposure = median(effective_source_exposure, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    reached_summary %>%
      select(transmission_regime, reported_case_trigger, proportion_reached),
    by = c("transmission_regime", "reported_case_trigger")
  )

readr::write_csv(
  summary_estimates,
  file.path(args$out, "tables", "early_inference_summary_v2.csv")
)

# Threshold-attainment labels make it explicit why the low-transmission line
# stops: few or no low-R0 outbreaks reach later case counts.
attainment_labels <- reached_summary %>%
  mutate(
    label = paste0(n_reached, "/", n_simulated),
    transmission_label = factor(
      transmission_label,
      levels = c("High transmission", "Low transmission")
    )
  )

readr::write_csv(
  attainment_labels,
  file.path(args$out, "tables", "early_inference_threshold_attainment_labels.csv")
)

# Fixed-input audit table
fixed_inputs <- bind_rows(
  tibble(
    parameter = c(
      "input", "delay_file_used", "out", "n_runs_per_regime", "seed",
      "reported_case_triggers", "max_days", "max_cases", "initial_infected",
      "p_detect", "detection_delay_mean", "detection_delay_sd",
      "n_report_delay_draws", "prior_log_R_mean", "prior_log_R_sd",
      "prior_log_k_mean", "prior_log_k_sd", "min_reached_for_plot",
      "posterior_approximation"
    ),
    value = as.character(c(
      args$input, delay_file, args$out, args$n_runs, args$seed,
      paste(args$reported_case_triggers, collapse = ","),
      args$max_days, args$max_cases, args$initial_infected,
      args$p_detect, args$detection_delay_mean, args$detection_delay_sd,
      args$n_report_delay_draws,
      args$prior_log_R_mean, args$prior_log_R_sd,
      args$prior_log_k_mean, args$prior_log_k_sd,
      args$min_reached_for_plot,
      "Laplace approximation on log(R0), log(k) at posterior mode"
    ))
  ),
  transmission_regimes %>%
    transmute(
      parameter = paste0(transmission_regime, "_truth"),
      value = paste0("R0=", R0, "; k=", k)
    )
)

readr::write_csv(
  fixed_inputs,
  file.path(args$out, "tables", "early_inference_fixed_inputs_v2.csv")
)

# ------------------------------------------------------------
# Figure A
# ------------------------------------------------------------
# Panel A: posterior R0 estimate + typical 95% CrI.
# Panel B: posterior probability that R0 > 1. This makes the decision problem
# explicit and avoids asking the reader to infer it from overlapping ribbons.
#
# Trigger 1 is retained, because "one detected case" was part of the proposed
# analysis, but it is drawn as an open point and not connected to trigger 2.
# At one detected case the offspring likelihood is generally prior-dominated;
# the wide credible interval, rather than the center, is the meaningful result.

if (isTRUE(args$make_plots)) {
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required for the two-panel Figure A.", call. = FALSE)
  }

  plot_dat <- summary_estimates %>%
    filter(n_reached >= args$min_reached_for_plot) %>%
    mutate(
      transmission_label = factor(
        transmission_label,
        levels = c("High transmission", "Low transmission")
      )
    )

  plot_main <- plot_dat %>% filter(reported_case_trigger >= 2)
  plot_one <- plot_dat %>% filter(reported_case_trigger == 1)

  p_R <- ggplot(
    plot_main,
    aes(
      x = reported_case_trigger,
      y = R_median,
      group = transmission_label,
      linetype = transmission_label,
      shape = transmission_label
    )
  ) +
    geom_ribbon(
      aes(
        ymin = R_typical_lower,
        ymax = R_typical_upper,
        fill = transmission_label
      ),
      alpha = 0.18,
      colour = NA,
      show.legend = FALSE
    ) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 2.4) +
    geom_errorbar(
      data = plot_one,
      aes(
        x = reported_case_trigger,
        ymin = R_typical_lower,
        ymax = R_typical_upper,
        linetype = transmission_label
      ),
      width = 0.12,
      inherit.aes = FALSE,
      show.legend = FALSE
    ) +
    geom_point(
      data = plot_one,
      aes(
        x = reported_case_trigger,
        y = R_median,
        shape = transmission_label
      ),
      size = 2.6,
      fill = "white",
      inherit.aes = FALSE
    ) +
    geom_hline(yintercept = 1, linewidth = 0.45, linetype = "dotted") +
    scale_x_continuous(
      breaks = args$reported_case_triggers,
      minor_breaks = NULL
    ) +
    scale_y_log10(
      breaks = c(0.05, 0.1, 0.25, 0.5, 1, 2, 4, 8, 16),
      labels = scales::number_format(accuracy = 0.01)
    ) +
    coord_cartesian(ylim = c(0.03, 20)) +
    scale_linetype_manual(
      values = c("High transmission" = "solid", "Low transmission" = "dashed")
    ) +
    scale_shape_manual(
      values = c("High transmission" = 16, "Low transmission" = 15)
    ) +
    scale_fill_manual(
      values = c("High transmission" = "grey55", "Low transmission" = "grey20")
    ) +
    labs(
      x = "Detected cases available for inference",
      y = expression("Estimated " * R[0] * " (typical 95% CrI)"),
      linetype = NULL,
      shape = NULL,
      caption = paste0(
        "Open symbols at 1 detected case are prior-dominated. ",
        "Regime-threshold summaries are omitted when fewer than ",
        args$min_reached_for_plot, "/", args$n_runs, " simulations reached the threshold."
      )
    ) +
    annotate(
      "text", x = max(args$reported_case_triggers), y = 1.08,
      label = "R[0] == 1", parse = TRUE,
      hjust = 1, vjust = -0.2, size = 3.0
    ) +
    annotate(
      "text", x = 1, y = 0.04,
      label = "prior-dominated", angle = 90,
      hjust = 0, vjust = -0.4, size = 2.6
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.text = element_text(colour = "black"),
      axis.title.x = element_text(size = 12.5, margin = margin(t = 8)),
      axis.title.y = element_text(size = 12.5, margin = margin(r = 8)),
      legend.position = "bottom",
      legend.margin = margin(t = -3),
      plot.margin = margin(6, 8, 6, 6)
    )

  p_prob <- ggplot(
    plot_main,
    aes(
      x = reported_case_trigger,
      y = median_prob_R_gt1,
      group = transmission_label,
      linetype = transmission_label,
      shape = transmission_label
    )
  ) +
    geom_ribbon(
      aes(
        ymin = prob_R_gt1_q25,
        ymax = prob_R_gt1_q75,
        fill = transmission_label
      ),
      alpha = 0.18,
      colour = NA,
      show.legend = FALSE
    ) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 2.4) +
    geom_point(
      data = plot_one,
      aes(
        x = reported_case_trigger,
        y = median_prob_R_gt1,
        shape = transmission_label
      ),
      size = 2.6,
      fill = "white",
      inherit.aes = FALSE
    ) +
    scale_x_continuous(
      breaks = args$reported_case_triggers,
      minor_breaks = NULL
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.25),
      labels = scales::number_format(accuracy = 0.01),
      expand = expansion(mult = c(0, 0.02))
    ) +
    scale_linetype_manual(
      values = c("High transmission" = "solid", "Low transmission" = "dashed")
    ) +
    scale_shape_manual(
      values = c("High transmission" = 16, "Low transmission" = 15)
    ) +
    scale_fill_manual(
      values = c("High transmission" = "grey55", "Low transmission" = "grey20")
    ) +
    labs(
      x = "Detected cases available for inference",
      y = expression("Posterior probability " * P(R[0] > 1)),
      linetype = NULL,
      shape = NULL
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.text = element_text(colour = "black"),
      axis.title.x = element_text(size = 12.5, margin = margin(t = 8)),
      axis.title.y = element_text(size = 12.5, margin = margin(r = 8)),
      legend.position = "none",
      plot.margin = margin(6, 8, 6, 6)
    )

  p_combined <- patchwork::wrap_plots(
    p_R, p_prob,
    ncol = 1,
    guides = "collect"
  ) +
    patchwork::plot_annotation(tag_levels = "A") &
    theme(legend.position = "bottom")

  ggsave(
    file.path(args$out, "figures", "early_R_inference_Figure_A_v2.pdf"),
    p_combined,
    width = 6.3, height = 8.0, units = "in",
    device = cairo_pdf
  )
  ggsave(
    file.path(args$out, "figures", "early_R_inference_Figure_A_v2.png"),
    p_combined,
    width = 6.3, height = 8.0, units = "in", dpi = 600
  )

  # Supplementary dispersion plot. Keep k out of the main figure because it
  # is expected to be very weakly identified with only a few detected cases.
  p_k <- ggplot(
    plot_main,
    aes(
      x = reported_case_trigger,
      y = k_median,
      group = transmission_label,
      linetype = transmission_label,
      shape = transmission_label
    )
  ) +
    geom_ribbon(
      aes(ymin = k_typical_lower, ymax = k_typical_upper, fill = transmission_label),
      alpha = 0.18, colour = NA, show.legend = FALSE
    ) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 2.4) +
    scale_x_continuous(breaks = args$reported_case_triggers, minor_breaks = NULL) +
    scale_y_log10() +
    scale_linetype_manual(
      values = c("High transmission" = "solid", "Low transmission" = "dashed")
    ) +
    scale_shape_manual(
      values = c("High transmission" = 16, "Low transmission" = 15)
    ) +
    scale_fill_manual(
      values = c("High transmission" = "grey55", "Low transmission" = "grey20")
    ) +
    labs(
      x = "Detected cases available for inference",
      y = "Estimated dispersion, k (typical 95% CrI)",
      linetype = NULL, shape = NULL
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.text = element_text(colour = "black"),
      axis.title.x = element_text(size = 13, margin = margin(t = 8)),
      axis.title.y = element_text(size = 13, margin = margin(r = 8)),
      legend.position = "bottom"
    )

  ggsave(
    file.path(args$out, "figures", "early_k_inference_supplement_v2.pdf"),
    p_k,
    width = 5.8, height = 4.6, units = "in",
    device = cairo_pdf
  )
  ggsave(
    file.path(args$out, "figures", "early_k_inference_supplement_v2.png"),
    p_k,
    width = 5.8, height = 4.6, units = "in", dpi = 600
  )
}

message("Done.")
message("Main figure: ", file.path(args$out, "figures", "early_R_inference_Figure_A_v2.pdf"))
message("Summary table: ", file.path(args$out, "tables", "early_inference_summary_v2.csv"))

# ============================================================
# IMPROVED ONE-PANEL R0 PLOT
# ============================================================

plot_simple <- summary_estimates %>%
  filter(
    reported_case_trigger >= 2,
    n_reached >= args$min_reached_for_plot
  ) %>%
  mutate(
    transmission_label = factor(
      transmission_label,
      levels = c("High transmission", "Low transmission")
    )
  )

p_simple <- ggplot(
  plot_simple,
  aes(
    x = reported_case_trigger,
    y = R_median,
    group = transmission_label
  )
) +
  
  # Uncertainty ribbons
  geom_ribbon(
    data = plot_simple %>%
      filter(transmission_label == "High transmission"),
    aes(
      ymin = R_typical_lower,
      ymax = R_typical_upper
    ),
    fill = "grey75",
    alpha = 0.35,
    colour = NA
  ) +
  
  geom_ribbon(
    data = plot_simple %>%
      filter(transmission_label == "Low transmission"),
    aes(
      ymin = R_typical_lower,
      ymax = R_typical_upper
    ),
    fill = "grey45",
    alpha = 0.22,
    colour = NA
  ) +
  
  # True generating R0 values
  geom_hline(
    yintercept = 2.12,
    linetype = "longdash",
    linewidth = 0.4,
    colour = "grey55"
  ) +
  
  geom_hline(
    yintercept = 0.2416027,
    linetype = "longdash",
    linewidth = 0.4,
    colour = "grey55"
  ) +
  
  # Supercritical boundary
  geom_hline(
    yintercept = 1,
    linetype = "dotted",
    linewidth = 0.65,
    colour = "black"
  ) +
  
  # High-transmission estimates
  geom_line(
    data = plot_simple %>%
      filter(transmission_label == "High transmission"),
    aes(
      linetype = transmission_label
    ),
    linewidth = 0.9
  ) +
  
  geom_point(
    data = plot_simple %>%
      filter(transmission_label == "High transmission"),
    aes(
      shape = transmission_label
    ),
    size = 3
  ) +
  
  # Low-transmission estimates
  geom_line(
    data = plot_simple %>%
      filter(transmission_label == "Low transmission"),
    aes(
      linetype = transmission_label
    ),
    linewidth = 0.9
  ) +
  
  geom_point(
    data = plot_simple %>%
      filter(transmission_label == "Low transmission"),
    aes(
      shape = transmission_label
    ),
    size = 3
  ) +
  
  scale_x_continuous(
    breaks = c(2, 3, 5, 10, 15, 20),
    minor_breaks = NULL
  ) +
  
  scale_y_log10(
    breaks = c(
      0.1, 0.25, 0.5, 1,
      2, 4, 8, 16
    ),
    labels = scales::number_format(
      accuracy = 0.01
    )
  ) +
  
  coord_cartesian(
    ylim = c(0.08, 20)
  ) +
  
  scale_linetype_manual(
    values = c(
      "High transmission" = "solid",
      "Low transmission" = "dashed"
    )
  ) +
  
  scale_shape_manual(
    values = c(
      "High transmission" = 16,
      "Low transmission" = 15
    )
  ) +
  
  labs(
    x = "Detected cases",
    y = expression("Real-time estimate of " * R[0]),
    linetype = NULL,
    shape = NULL
  ) +
  
  theme_classic(
    base_size = 13
  ) +
  
  theme(
    axis.text = element_text(
      colour = "black"
    ),
    
    axis.title.x = element_text(
      size = 14,
      margin = margin(t = 8)
    ),
    
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 8)
    ),
    
    legend.position = "bottom",
    
    legend.box.spacing = unit(
      0.2,
      "cm"
    ),
    
    legend.margin = margin(
      t = 2,
      r = 0,
      b = 0,
      l = 0
    )
  )

ggsave(
  file.path(
    args$out,
    "figures",
    "R0_one_panel_improved.pdf"
  ),
  p_simple,
  width = 6.2,
  height = 4.6
)

ggsave(
  file.path(
    args$out,
    "figures",
    "R0_one_panel_improved.png"
  ),
  p_simple,
  width = 6.2,
  height = 4.6,
  dpi = 600
)

print(p_simple)