# Final main text analysis
input_dir <- "~/Hantavirus-paper/TD results"
out_dir   <- "~/Hantavirus-paper/BP results/Regime discrimination time"

n_runs_discrimination <- 50000L
n_runs_intervention   <- 20000L

days_after_first_report <- c(0, 7, 14, 21, 28, 35, 42, 56)

seed <- 20260813L

max_days  <- 260
max_cases <- 5000L

# Regimes used elsewhere in the paper
R_high <- 2.12
k_high <- 0.3157164

R_low  <- 0.2416027
k_low  <- 0.3998142

# Observation process
p_detect <- 0.95
detection_delay_mean <- 2.0
detection_delay_sd   <- 1.0

# Intervention process: match Figure 3 as closely as possible
q_system <- 0.95

isolation_delay_mean <- 0.5
isolation_delay_sd   <- 0.5
isolation_effectiveness <- 0.99

quarantine_duration <- 42.0
quarantine_effectiveness <- 1.0
quarantine_adherence <- 1.0

backward_tracing_window <- 14.0
daily_contact_rate <- 13.50
contact_count_dispersion <- 2.11
contact_multiplier <- 1.0

# Bootstrap CIs for Panel A AUC
n_boot_auc <- 500L


# Packages
required <- c(
  "Rcpp", "dplyr", "tidyr", "readr", "ggplot2",
  "tibble", "scales", "patchwork", "pROC"
)

missing <- required[
  !vapply(required, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing) > 0) {
  stop(
    "Missing packages: ",
    paste(missing, collapse = ", "),
    "\nInstall with:\ninstall.packages(c(",
    paste(sprintf('"%s"', missing), collapse = ", "),
    "))",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(Rcpp)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(tibble)
  library(scales)
  library(patchwork)
})

set.seed(seed)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)


# Load posterior joint timing draws
delay_candidates <- c(
  file.path(input_dir, "posterior_pair_level_draws.csv"),
  file.path(input_dir, "tables", "posterior_pair_level_draws.csv"),
  file.path(input_dir, "joint_delay_draws_for_branching.csv"),
  file.path(input_dir, "tables", "joint_delay_draws_for_branching.csv")
)

delay_file <- delay_candidates[file.exists(delay_candidates)][1]

if (length(delay_file) == 0 || is.na(delay_file)) {
  hits <- list.files(
    input_dir,
    pattern = "^(posterior_pair_level_draws|joint_delay_draws_for_branching)\\.csv$",
    recursive = TRUE,
    full.names = TRUE
  )
  
  if (length(hits) > 0) {
    posterior_hits <- hits[basename(hits) == "posterior_pair_level_draws.csv"]
    delay_file <- if (length(posterior_hits) > 0) posterior_hits[1] else hits[1]
  }
}

if (length(delay_file) == 0 || is.na(delay_file) || !file.exists(delay_file)) {
  stop(
    "Could not find posterior_pair_level_draws.csv or ",
    "joint_delay_draws_for_branching.csv under:\n",
    input_dir,
    call. = FALSE
  )
}

pair_draws <- read_csv(delay_file, show_col_types = FALSE)

if ("generation_time" %in% names(pair_draws)) {
  pair_draws$generation_time <- as.numeric(pair_draws$generation_time)
} else if ("generation_interval" %in% names(pair_draws)) {
  pair_draws$generation_time <- as.numeric(pair_draws$generation_interval)
} else {
  stop("No generation_time/generation_interval column found.")
}

if ("recipient_incubation" %in% names(pair_draws)) {
  pair_draws$recipient_incubation <- as.numeric(pair_draws$recipient_incubation)
} else if ("incubation" %in% names(pair_draws)) {
  pair_draws$recipient_incubation <- as.numeric(pair_draws$incubation)
} else {
  stop("No recipient_incubation/incubation column found.")
}

if (!"delta" %in% names(pair_draws)) {
  stop("No delta column found; joint timing sampling requires delta.")
}
pair_draws$delta <- as.numeric(pair_draws$delta)

if ("source_incubation" %in% names(pair_draws)) {
  pair_draws$source_incubation <- as.numeric(pair_draws$source_incubation)
} else {
  pair_draws$source_incubation <-
    pair_draws$generation_time - pair_draws$delta
}

pair_draws <- pair_draws %>%
  filter(
    is.finite(generation_time), generation_time > 0,
    is.finite(recipient_incubation), recipient_incubation > 0,
    is.finite(source_incubation), source_incubation > 0,
    is.finite(delta)
  )

if (nrow(pair_draws) == 0) {
  stop("No usable joint timing draws remained after filtering.")
}

generation_draws           <- pair_draws$generation_time
recipient_incubation_draws <- pair_draws$recipient_incubation
source_incubation_draws    <- pair_draws$source_incubation
delta_draws                <- pair_draws$delta


# Rcpp simulation core
Rcpp::sourceCpp(code = '
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

inline double runif01() {
  return R::runif(0.0, 1.0);
}

inline bool finite_num(double x) {
  return R_finite(x);
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

inline double gamma_delay_one(double mean, double sd, double min_value = 0.0) {
  if (!R_finite(mean) || mean <= 0.0) return min_value;
  if (!R_finite(sd) || sd <= 0.0) return std::max(mean, min_value);

  double shape = (mean * mean) / (sd * sd);
  double scale = (sd * sd) / mean;
  double y = R::rgamma(shape, scale);

  if (!R_finite(y) || y < min_value) y = min_value;
  return y;
}

inline int conditioned_joint_index(
    const NumericVector& source_incubation_draws,
    double source_incubation) {

  int n = source_incubation_draws.size();

  if (n <= 1 || !R_finite(source_incubation)) {
    return sample_index(n);
  }

  int n_try = std::min(64, n);
  int best_idx = sample_index(n);
  double best_dist = R_PosInf;

  for (int tt = 0; tt < n_try; ++tt) {
    int idx = sample_index(n);
    double x = source_incubation_draws[idx];

    if (!R_finite(x)) continue;

    double d = std::fabs(x - source_incubation);

    if (d < best_dist) {
      best_dist = d;
      best_idx = idx;
    }
  }

  return best_idx;
}

struct DelayDraw {
  double G;
  double recipient_inc;
};

inline DelayDraw sample_joint_delay(
    const NumericVector& generation_draws,
    const NumericVector& recipient_incubation_draws,
    const NumericVector& source_incubation_draws,
    const NumericVector& delta_draws,
    double source_incubation,
    double max_days) {

  int n = generation_draws.size();

  // Preserve dependence between source incubation and transmission timing.
  for (int guard = 0; guard < 200; ++guard) {

    int idx = conditioned_joint_index(
      source_incubation_draws,
      source_incubation
    );

    double G = source_incubation + delta_draws[idx];
    double inc = recipient_incubation_draws[idx];

    if (R_finite(G) && G > 0.0 && G < max_days &&
        R_finite(inc) && inc > 0.0) {

      DelayDraw out;
      out.G = G;
      out.recipient_inc = inc;
      return out;
    }
  }

  // Fallback to direct sampled joint row.
  for (int guard = 0; guard < 200; ++guard) {

    int idx = sample_index(n);

    double G = generation_draws[idx];
    double inc = recipient_incubation_draws[idx];

    if (R_finite(G) && G > 0.0 && G < max_days &&
        R_finite(inc) && inc > 0.0) {

      DelayDraw out;
      out.G = G;
      out.recipient_inc = inc;
      return out;
    }
  }

  DelayDraw out;
  out.G = 22.5;
  out.recipient_inc = 18.5;
  return out;
}

struct PotentialTx {
  int source;
  int child;
  double time;
};


// ============================================================
// A. NATURAL OUTBREAKS FOR DISCRIMINATION
// ============================================================
//
// Every returned trajectory has at least one reported case.
// The clock is then reset to the first report.
// At each requested elapsed day, return the cumulative number
// of reported cases observed by then.
//
// [[Rcpp::export]]
DataFrame simulate_observed_counts_cpp(
    int n_runs,
    NumericVector snapshot_days,
    NumericVector generation_draws,
    NumericVector recipient_incubation_draws,
    NumericVector source_incubation_draws,
    NumericVector delta_draws,
    double Rmean,
    double k,
    double p_detect,
    double detection_delay_mean,
    double detection_delay_sd,
    int max_cases,
    double max_days,
    int seed) {

  RNGScope scope;

  Environment base_env("package:base");
  Function set_seed = base_env["set.seed"];
  set_seed(seed);

  int S = snapshot_days.size();
  int total_rows = n_runs * S;

  IntegerVector out_run_id(total_rows);
  NumericVector out_day(total_rows);
  IntegerVector out_reported_count(total_rows);
  IntegerVector out_infections_by_time(total_rows);
  IntegerVector out_total_infections(total_rows);
  IntegerVector out_total_reported(total_rows);

  int row = 0;
  int accepted = 0;

  while (accepted < n_runs) {

    std::vector<double> infection_time;
    std::vector<double> incubation;
    std::vector<double> detection_time;
    std::vector<int> processed;

    infection_time.reserve(max_cases);
    incubation.reserve(max_cases);
    detection_time.reserve(max_cases);
    processed.reserve(max_cases);

    auto add_case = [&](double itime, double inc) {

      double onset = itime + inc;
      double dtime = R_PosInf;

      if (runif01() < p_detect) {
        dtime = onset + gamma_delay_one(
          detection_delay_mean,
          detection_delay_sd,
          0.0
        );
      }

      infection_time.push_back(itime);
      incubation.push_back(inc);
      detection_time.push_back(dtime);
      processed.push_back(0);
    };

    add_case(
      0.0,
      sample_from_vec(recipient_incubation_draws)
    );

    double nb_prob = k / (k + Rmean);

    while ((int)infection_time.size() < max_cases) {

      int todo = -1;
      double earliest = R_PosInf;

      for (int ii = 0; ii < (int)infection_time.size(); ++ii) {
        if (processed[ii] == 0 && infection_time[ii] < earliest) {
          earliest = infection_time[ii];
          todo = ii;
        }
      }

      if (todo < 0) break;

      processed[todo] = 1;

      int z = (int)R::rnbinom(k, nb_prob);

      for (int zz = 0; zz < z; ++zz) {

        if ((int)infection_time.size() >= max_cases) break;

        DelayDraw dd = sample_joint_delay(
          generation_draws,
          recipient_incubation_draws,
          source_incubation_draws,
          delta_draws,
          incubation[todo],
          max_days
        );

        double child_time = infection_time[todo] + dd.G;

        if (!finite_num(child_time) || child_time > max_days) continue;

        add_case(child_time, dd.recipient_inc);
      }
    }

    int N = infection_time.size();

    std::vector<double> reports;
    reports.reserve(N);

    for (int ii = 0; ii < N; ++ii) {
      if (finite_num(detection_time[ii]) &&
          detection_time[ii] <= max_days) {
        reports.push_back(detection_time[ii]);
      }
    }

    std::sort(reports.begin(), reports.end());

    // The analysis is explicitly anchored to the first observed report.
    if (reports.size() == 0) continue;

    accepted++;

    double first_report = reports[0];

    for (int ss = 0; ss < S; ++ss) {

      double target = first_report + snapshot_days[ss];

      int nrep = 0;
      for (int jj = 0; jj < (int)reports.size(); ++jj) {
        if (reports[jj] <= target) nrep++;
        else break;
      }

      int ninf = 0;
      for (int ii = 0; ii < N; ++ii) {
        if (infection_time[ii] <= target) ninf++;
      }

      out_run_id[row] = accepted;
      out_day[row] = snapshot_days[ss];
      out_reported_count[row] = nrep;
      out_infections_by_time[row] = ninf;
      out_total_infections[row] = N;
      out_total_reported[row] = reports.size();

      row++;
    }
  }

  return DataFrame::create(
    Named("run_id") = out_run_id,
    Named("day_after_first_report") = out_day,
    Named("reported_count") = out_reported_count,
    Named("infections_by_time") = out_infections_by_time,
    Named("total_infections_natural") = out_total_infections,
    Named("total_reported_natural") = out_total_reported
  );
}


// ============================================================
// B. FIXED-TIME INTERVENTION AFTER FIRST REPORT
// ============================================================
//
// True regime is supplied by Rmean/k (main analysis = high).
//
// First generate the potential transmission tree.
// Define response time = first reported case + requested delay.
// Then replay transmission chronologically with:
//   - quarantine at response with q_system
//   - case isolation
//   - tracing/quarantine of downstream infected contacts
//   - reconstructed uninfected contact burden
//
// This follows the event logic of the existing reported-case
// trigger analysis, but the response trigger is elapsed TIME,
// not an observed case-count threshold.
//
// [[Rcpp::export]]
DataFrame simulate_fixed_delay_intervention_cpp(
    int n_runs,
    NumericVector response_delays,
    NumericVector generation_draws,
    NumericVector recipient_incubation_draws,
    NumericVector source_incubation_draws,
    NumericVector delta_draws,
    int max_cases,
    double max_days,
    double Rmean,
    double k,
    double p_detect,
    double detection_delay_mean,
    double detection_delay_sd,
    double isolation_delay_mean,
    double isolation_delay_sd,
    double isolation_effectiveness,
    double q_system,
    double quarantine_duration,
    double quarantine_effectiveness,
    double quarantine_adherence,
    double backward_tracing_window,
    double daily_contact_rate,
    double contact_count_dispersion,
    double contact_multiplier,
    int seed) {

  RNGScope scope;

  Environment base_env("package:base");
  Function set_seed = base_env["set.seed"];
  set_seed(seed);

  int D = response_delays.size();
  int total_rows = n_runs * D;

  IntegerVector out_run_id(total_rows);
  NumericVector out_delay(total_rows);
  NumericVector out_response_time(total_rows);

  IntegerVector out_reported_at_response(total_rows);
  IntegerVector out_infections_at_response(total_rows);

  IntegerVector out_total_infections(total_rows);
  IntegerVector out_post_response_infections(total_rows);

  IntegerVector out_total_quarantined(total_rows);
  IntegerVector out_uninfected_quarantined(total_rows);

  NumericVector out_response_duration(total_rows);

  int row = 0;
  int accepted = 0;

  while (accepted < n_runs) {

    // --------------------------------------------------------
    // Build ONE potential natural epidemic.
    // We replay the same potential tree under each delay.
    // --------------------------------------------------------

    std::vector<int> infector;
    std::vector<double> infection_time;
    std::vector<double> incubation;
    std::vector<double> symptom_onset;
    std::vector<int> detected;
    std::vector<double> detection_time;
    std::vector<double> isolation_time;
    std::vector<int> processed;
    std::vector<PotentialTx> potential_tx;

    infector.reserve(max_cases);
    infection_time.reserve(max_cases);
    incubation.reserve(max_cases);
    symptom_onset.reserve(max_cases);
    detected.reserve(max_cases);
    detection_time.reserve(max_cases);
    isolation_time.reserve(max_cases);
    processed.reserve(max_cases);
    potential_tx.reserve(max_cases * 2);

    auto add_case = [&](int src, double itime, double inc) {

      double onset = itime + inc;

      int det = (runif01() < p_detect) ? 1 : 0;

      double dtime = R_PosInf;
      double isotime = R_PosInf;

      if (det == 1) {
        dtime = onset + gamma_delay_one(
          detection_delay_mean,
          detection_delay_sd,
          0.0
        );

        isotime = dtime + gamma_delay_one(
          isolation_delay_mean,
          isolation_delay_sd,
          0.0
        );
      }

      infector.push_back(src);
      infection_time.push_back(itime);
      incubation.push_back(inc);
      symptom_onset.push_back(onset);
      detected.push_back(det);
      detection_time.push_back(dtime);
      isolation_time.push_back(isotime);
      processed.push_back(0);

      return (int)infection_time.size() - 1;
    };

    add_case(
      NA_INTEGER,
      0.0,
      sample_from_vec(recipient_incubation_draws)
    );

    double nb_prob = k / (k + Rmean);

    while ((int)infection_time.size() < max_cases) {

      int todo = -1;
      double earliest = R_PosInf;

      for (int ii = 0; ii < (int)infection_time.size(); ++ii) {
        if (processed[ii] == 0 && infection_time[ii] < earliest) {
          earliest = infection_time[ii];
          todo = ii;
        }
      }

      if (todo < 0) break;

      processed[todo] = 1;

      int z = (int)R::rnbinom(k, nb_prob);

      for (int zz = 0; zz < z; ++zz) {

        if ((int)infection_time.size() >= max_cases) break;

        DelayDraw dd = sample_joint_delay(
          generation_draws,
          recipient_incubation_draws,
          source_incubation_draws,
          delta_draws,
          incubation[todo],
          max_days
        );

        double child_time = infection_time[todo] + dd.G;

        if (!finite_num(child_time) || child_time > max_days) continue;

        int child_idx = add_case(
          todo,
          child_time,
          dd.recipient_inc
        );

        PotentialTx tx;
        tx.source = todo;
        tx.child = child_idx;
        tx.time = child_time;
        potential_tx.push_back(tx);
      }
    }

    int N = infection_time.size();

    std::vector<double> reports;

    for (int ii = 0; ii < N; ++ii) {
      if (detected[ii] == 1 &&
          finite_num(detection_time[ii]) &&
          detection_time[ii] <= max_days) {
        reports.push_back(detection_time[ii]);
      }
    }

    std::sort(reports.begin(), reports.end());

    // Need an observable first report to define day zero.
    if (reports.size() == 0) continue;

    accepted++;

    double first_report = reports[0];

    // Fixed order for potential transmissions.
    std::vector<int> tx_order(potential_tx.size());
    std::iota(tx_order.begin(), tx_order.end(), 0);

    std::sort(
      tx_order.begin(),
      tx_order.end(),
      [&](int a, int b) {
        return potential_tx[a].time < potential_tx[b].time;
      }
    );

    // --------------------------------------------------------
    // Replay independently for every intervention delay.
    // --------------------------------------------------------

    for (int dd_idx = 0; dd_idx < D; ++dd_idx) {

      double response_delay = response_delays[dd_idx];
      double response_time = first_report + response_delay;

      std::vector<int> retained(N, 0);
      std::vector<int> quarantined(N, 0);
      std::vector<double> q_start(N, R_PosInf);
      std::vector<double> q_end(N, R_PosInf);

      retained[0] = 1;

      // Number of infections/reports that had already occurred
      // by the response time in the potential natural epidemic.
      int infections_at_response = 0;
      int reported_at_response = 0;

      for (int ii = 0; ii < N; ++ii) {
        if (infection_time[ii] <= response_time) {
          infections_at_response++;
        }

        if (detected[ii] == 1 &&
            finite_num(detection_time[ii]) &&
            detection_time[ii] <= response_time) {
          reported_at_response++;
        }
      }

      // Quarantine infected members of the known exposed cohort
      // at the moment the response begins, with completeness q.
      int infected_quarantined_at_response = 0;

      for (int ii = 0; ii < N; ++ii) {
        if (infection_time[ii] <= response_time &&
            runif01() < q_system) {

          quarantined[ii] = 1;
          q_start[ii] = response_time;
          q_end[ii] = response_time + quarantine_duration;
          infected_quarantined_at_response++;
        }
      }

      int future_infected_quarantined = 0;

      // Replay the potential tree chronologically.
      for (int oo = 0; oo < (int)tx_order.size(); ++oo) {

        PotentialTx tx = potential_tx[tx_order[oo]];

        int s = tx.source;
        int c = tx.child;
        double t = tx.time;

        if (s < 0 || s >= N || c < 0 || c >= N) continue;
        if (retained[s] == 0) continue;

        bool after_response = (t >= response_time);

        bool source_quarantined =
          after_response &&
          quarantined[s] == 1 &&
          finite_num(q_start[s]) &&
          t >= q_start[s] &&
          t <= q_end[s];

        bool source_isolated =
          after_response &&
          detected[s] == 1 &&
          finite_num(isolation_time[s]) &&
          t >= isolation_time[s];

        double retain_prob = 1.0;

        if (after_response) {

          if (source_quarantined &&
              runif01() < quarantine_adherence) {
            retain_prob *=
              (1.0 - quarantine_effectiveness);
          }

          if (source_isolated) {
            retain_prob *=
              (1.0 - isolation_effectiveness);
          }
        }

        if (runif01() > retain_prob) {
          continue;
        }

        retained[c] = 1;

        // Successfully traced infected contact after response.
        if (after_response && runif01() < q_system) {

          // Contact notification cannot predate infection.
          // If child is later detected, detection time is used as
          // a conservative proxy exactly as in the prior trigger code.
          double trace_time = detection_time[c];

          if (!finite_num(trace_time) || trace_time < t) {
            trace_time = t;
          }

          if (quarantined[c] == 0) {
            future_infected_quarantined++;
          }

          quarantined[c] = 1;
          q_start[c] = trace_time;
          q_end[c] = trace_time + quarantine_duration;
        }
      }

      // Count retained infections.
      int total_retained = 0;
      int post_response_infections = 0;

      for (int ii = 0; ii < N; ++ii) {

        if (retained[ii] == 0) continue;

        total_retained++;

        if (infection_time[ii] >= response_time) {
          post_response_infections++;
        }
      }

      // ------------------------------------------------------
      // Reconstruct uninfected contact quarantine burden
      // for detected retained cases after response.
      // ------------------------------------------------------

      int future_uninfected_contacts_quarantined = 0;
      double max_end = response_time + quarantine_duration;

      for (int jj = 0; jj < N; ++jj) {

        if (retained[jj] == 0) continue;

        if (detected[jj] != 1 ||
            !finite_num(detection_time[jj]) ||
            detection_time[jj] < response_time) {
          continue;
        }

        bool already_quarantined_at_detection =
          quarantined[jj] == 1 &&
          finite_num(q_start[jj]) &&
          detection_time[jj] >= q_start[jj] &&
          detection_time[jj] <= q_end[jj];

        if (already_quarantined_at_detection) {
          if (finite_num(q_end[jj]) && q_end[jj] > max_end) {
            max_end = q_end[jj];
          }
          continue;
        }

        double trace_start = std::max(
          symptom_onset[jj] - backward_tracing_window,
          infection_time[jj]
        );

        double trace_end = detection_time[jj];

        double trace_duration =
          std::max(0.0, trace_end - trace_start);

        int infected_contacts_in_pool = 0;

        for (int cc = 0; cc < N; ++cc) {

          if (retained[cc] == 0) continue;
          if (infector[cc] != jj) continue;

          if (infection_time[cc] >= trace_start &&
              infection_time[cc] <= trace_end) {
            infected_contacts_in_pool++;
          }
        }

        double mean_contacts =
          daily_contact_rate *
          contact_multiplier *
          trace_duration;

        double size_contacts =
          contact_count_dispersion *
          trace_duration;

        int total_traceable_contacts =
          infected_contacts_in_pool;

        if (mean_contacts > 0.0) {

          if (size_contacts > 0.0) {

            double prob_contacts =
              size_contacts /
              (size_contacts + mean_contacts);

            total_traceable_contacts =
              (int)R::rnbinom(
                size_contacts,
                prob_contacts
              );

          } else {

            total_traceable_contacts =
              (int)R::rpois(mean_contacts);
          }

          if (total_traceable_contacts <
              infected_contacts_in_pool) {

            total_traceable_contacts =
              infected_contacts_in_pool;
          }
        }

        int uninfected_contacts =
          total_traceable_contacts -
          infected_contacts_in_pool;

        int traced_uninfected = 0;

        if (uninfected_contacts > 0) {
          traced_uninfected =
            (int)R::rbinom(
              (double)uninfected_contacts,
              q_system
            );
        }

        future_uninfected_contacts_quarantined +=
          traced_uninfected;

        double this_end =
          detection_time[jj] +
          quarantine_duration;

        if (this_end > max_end) {
          max_end = this_end;
        }
      }

      int total_quarantined =
        infected_quarantined_at_response +
        future_infected_quarantined +
        future_uninfected_contacts_quarantined;

      double response_duration =
        std::max(0.0, max_end - response_time);

      out_run_id[row] = accepted;
      out_delay[row] = response_delay;
      out_response_time[row] = response_time;

      out_reported_at_response[row] =
        reported_at_response;

      out_infections_at_response[row] =
        infections_at_response;

      out_total_infections[row] =
        total_retained;

      out_post_response_infections[row] =
        post_response_infections;

      out_total_quarantined[row] =
        total_quarantined;

      out_uninfected_quarantined[row] =
        future_uninfected_contacts_quarantined;

      out_response_duration[row] =
        response_duration;

      row++;
    }
  }

  return DataFrame::create(
    Named("run_id") = out_run_id,
    Named("days_after_first_report") = out_delay,
    Named("response_time") = out_response_time,
    Named("reported_cases_at_response") = out_reported_at_response,
    Named("infections_at_response") = out_infections_at_response,
    Named("total_infections") = out_total_infections,
    Named("post_response_infections") = out_post_response_infections,
    Named("total_quarantined") = out_total_quarantined,
    Named("uninfected_contacts_quarantined") = out_uninfected_quarantined,
    Named("response_duration_days") = out_response_duration
  );
}
')


# PANEL A
obs_high <- simulate_observed_counts_cpp(
  n_runs = n_runs_discrimination,
  snapshot_days = days_after_first_report,
  generation_draws = generation_draws,
  recipient_incubation_draws = recipient_incubation_draws,
  source_incubation_draws = source_incubation_draws,
  delta_draws = delta_draws,
  Rmean = R_high,
  k = k_high,
  p_detect = p_detect,
  detection_delay_mean = detection_delay_mean,
  detection_delay_sd = detection_delay_sd,
  max_cases = max_cases,
  max_days = max_days,
  seed = seed
) %>%
  as_tibble() %>%
  mutate(
    regime = "high",
    high = 1L
  )

obs_low <- simulate_observed_counts_cpp(
  n_runs = n_runs_discrimination,
  snapshot_days = days_after_first_report,
  generation_draws = generation_draws,
  recipient_incubation_draws = recipient_incubation_draws,
  source_incubation_draws = source_incubation_draws,
  delta_draws = delta_draws,
  Rmean = R_low,
  k = k_low,
  p_detect = p_detect,
  detection_delay_mean = detection_delay_mean,
  detection_delay_sd = detection_delay_sd,
  max_cases = max_cases,
  max_days = max_days,
  seed = seed + 1000003L
) %>%
  as_tibble() %>%
  mutate(
    regime = "low",
    high = 0L
  )

observed_runs <- bind_rows(
  obs_high,
  obs_low
)

write_csv(
  observed_runs,
  file.path(
    out_dir, "tables",
    "observed_reported_counts_by_time.csv"
  )
)


# AUC at each time point
auc_one_day <- function(dat, n_boot = 500, seed_local = 1) {
  
  roc_obj <- pROC::roc(
    response = dat$high,
    predictor = dat$reported_count,
    levels = c(0, 1),
    direction = "<",
    quiet = TRUE
  )
  
  auc_est <- as.numeric(
    pROC::auc(roc_obj)
  )
  
  set.seed(seed_local)
  
  boot_auc <- replicate(
    n_boot,
    {
      ids_high <- sample(
        which(dat$high == 1),
        replace = TRUE
      )
      
      ids_low <- sample(
        which(dat$high == 0),
        replace = TRUE
      )
      
      z <- bind_rows(
        dat[ids_high, , drop = FALSE],
        dat[ids_low, , drop = FALSE]
      )
      
      r <- pROC::roc(
        response = z$high,
        predictor = z$reported_count,
        levels = c(0, 1),
        direction = "<",
        quiet = TRUE
      )
      
      as.numeric(
        pROC::auc(r)
      )
    }
  )
  
  tibble(
    auc = auc_est,
    lower = quantile(
      boot_auc, 0.025,
      na.rm = TRUE,
      names = FALSE
    ),
    upper = quantile(
      boot_auc, 0.975,
      na.rm = TRUE,
      names = FALSE
    )
  )
}

auc_summary <- lapply(
  seq_along(days_after_first_report),
  function(ii) {
    
    dd <- days_after_first_report[ii]
    
    z <- observed_runs %>%
      filter(
        day_after_first_report == dd
      )
    
    bind_cols(
      tibble(
        days_after_first_report = dd
      ),
      auc_one_day(
        z,
        n_boot = n_boot_auc,
        seed_local = seed + 100L + ii
      )
    )
  }
) %>%
  bind_rows()

write_csv(
  auc_summary,
  file.path(
    out_dir, "tables",
    "regime_discrimination_auc.csv"
  )
)


# Also save the observable case-count distributions
count_summary <- observed_runs %>%
  group_by(
    regime,
    day_after_first_report
  ) %>%
  summarise(
    median_reported =
      median(reported_count),
    
    q25_reported =
      quantile(reported_count, 0.25),
    
    q75_reported =
      quantile(reported_count, 0.75),
    
    q90_reported =
      quantile(reported_count, 0.90),
    
    q95_reported =
      quantile(reported_count, 0.95),
    
    p_still_one_reported =
      mean(reported_count == 1),
    
    p_at_least_2 =
      mean(reported_count >= 2),
    
    p_at_least_3 =
      mean(reported_count >= 3),
    
    p_at_least_5 =
      mean(reported_count >= 5),
    
    .groups = "drop"
  )

write_csv(
  count_summary,
  file.path(
    out_dir, "tables",
    "observable_case_count_summary.csv"
  )
)


# PANEL B:
intervention_runs <- simulate_fixed_delay_intervention_cpp(
  n_runs = n_runs_intervention,
  response_delays = days_after_first_report,
  generation_draws = generation_draws,
  recipient_incubation_draws = recipient_incubation_draws,
  source_incubation_draws = source_incubation_draws,
  delta_draws = delta_draws,
  max_cases = max_cases,
  max_days = max_days,
  Rmean = R_high,
  k = k_high,
  p_detect = p_detect,
  detection_delay_mean = detection_delay_mean,
  detection_delay_sd = detection_delay_sd,
  isolation_delay_mean = isolation_delay_mean,
  isolation_delay_sd = isolation_delay_sd,
  isolation_effectiveness = isolation_effectiveness,
  q_system = q_system,
  quarantine_duration = quarantine_duration,
  quarantine_effectiveness = quarantine_effectiveness,
  quarantine_adherence = quarantine_adherence,
  backward_tracing_window = backward_tracing_window,
  daily_contact_rate = daily_contact_rate,
  contact_count_dispersion = contact_count_dispersion,
  contact_multiplier = contact_multiplier,
  seed = seed + 2000003L
) %>%
  as_tibble()

write_csv(
  intervention_runs,
  file.path(
    out_dir, "tables",
    "fixed_time_intervention_run_level_results.csv"
  )
)


# Intervention summaries + exact binomial CIs
binom_ci <- function(events, n) {
  if (n <= 0) return(c(NA_real_, NA_real_))
  as.numeric(
    stats::binom.test(events, n)$conf.int
  )
}

intervention_summary <- intervention_runs %>%
  group_by(
    days_after_first_report
  ) %>%
  summarise(
    n = n(),
    
    infections_ge_100 =
      sum(total_infections >= 100),
    
    infections_ge_500 =
      sum(total_infections >= 500),
    
    p_infections_ge_100 =
      mean(total_infections >= 100),
    
    p_infections_ge_500 =
      mean(total_infections >= 500),
    
    median_total_infections =
      median(total_infections),
    
    q25_total_infections =
      quantile(total_infections, 0.25),
    
    q75_total_infections =
      quantile(total_infections, 0.75),
    
    q95_total_infections =
      quantile(total_infections, 0.95),
    
    median_infections_at_response =
      median(infections_at_response),
    
    q25_infections_at_response =
      quantile(infections_at_response, 0.25),
    
    q75_infections_at_response =
      quantile(infections_at_response, 0.75),
    
    median_reported_at_response =
      median(reported_cases_at_response),
    
    median_total_quarantined =
      median(total_quarantined),
    
    q25_total_quarantined =
      quantile(total_quarantined, 0.25),
    
    q75_total_quarantined =
      quantile(total_quarantined, 0.75),
    
    q95_total_quarantined =
      quantile(total_quarantined, 0.95),
    
    median_uninfected_quarantined =
      median(uninfected_contacts_quarantined),
    
    q25_uninfected_quarantined =
      quantile(uninfected_contacts_quarantined, 0.25),
    
    q75_uninfected_quarantined =
      quantile(uninfected_contacts_quarantined, 0.75),
    
    median_response_duration_weeks =
      median(response_duration_days) / 7,
    
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    p100_lower =
      binom_ci(
        infections_ge_100,
        n
      )[1],
    
    p100_upper =
      binom_ci(
        infections_ge_100,
        n
      )[2],
    
    p500_lower =
      binom_ci(
        infections_ge_500,
        n
      )[1],
    
    p500_upper =
      binom_ci(
        infections_ge_500,
        n
      )[2]
  ) %>%
  ungroup()

write_csv(
  intervention_summary,
  file.path(
    out_dir, "tables",
    "fixed_time_intervention_summary.csv"
  )
)


# ============================================================
# SELF-CONTAINED ADD-ON
#
# X axis:
#   Probability outbreak is the higher-transmission scenario
#
# Y axis:
#   Probability of:
#      >10 infections
#      >=50 infections
#      >=100 infections
#
# Uses EXISTING Figure 3 simulation outputs.
# NO RESIMULATION.
# NO CASE-COUNT NUMBERS PRINTED ON GRAPH.
# ============================================================

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(scales)


# ============================================================
# PATHS
# ============================================================

fig3_dir <-
  "~/Hantavirus-paper/BP results/Fig 3"

out_dir <-
  "~/Hantavirus-paper/BP results/Regime discrimination time"

figure_dir <-
  file.path(out_dir, "figures")

table_dir <-
  file.path(out_dir, "tables")

dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  table_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# SETTINGS
# ============================================================

target_triggers <- c(
  1, 2, 3, 5, 10, 15, 20
)

set.seed(20260813)


# ============================================================
# PART 1
# PROBABILITY OF HIGHER-TRANSMISSION SCENARIO
#
# Uses probability of reaching each reported-case milestone
# under HIGH versus LOW transmission.
#
# Assumes equal prior probability (50:50) of the two modeled
# transmission scenarios.
# ============================================================

attainment <- read_csv(
  file.path(
    fig3_dir,
    "tables",
    "reported_case_trigger_attainment.csv"
  ),
  show_col_types = FALSE
)


# ------------------------------------------------------------
# Keep one q value to remove duplicate attainment rows.
#
# q does not affect attainment because the trigger is reached
# before the intervention is applied.
# ------------------------------------------------------------

attain_use <- attainment %>%
  filter(
    reported_case_trigger %in% target_triggers,
    abs(q_system - 0.95) < 1e-8
  ) %>%
  mutate(
    regime_simple = case_when(
      transmission_regime == "supercritical" ~ "high",
      transmission_regime == "subcritical" ~ "low",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(
    !is.na(regime_simple)
  ) %>%
  select(
    regime_simple,
    reported_case_trigger,
    n_generated,
    n_reached_trigger,
    probability_reached_trigger
  )


# ------------------------------------------------------------
# Put HIGH and LOW side-by-side
# ------------------------------------------------------------

attain_wide <- attain_use %>%
  pivot_wider(
    names_from = regime_simple,
    values_from = c(
      n_generated,
      n_reached_trigger,
      probability_reached_trigger
    )
  )


# ------------------------------------------------------------
# Probability HIGH assuming equal prior probabilities
#
# P(HIGH | reached m)
# =
# P(reached m | HIGH)
# -----------------------------------
# P(reached m | HIGH) + P(reached m | LOW)
# ------------------------------------------------------------

panelA_summary <- attain_wide %>%
  mutate(
    
    p_reach_high =
      probability_reached_trigger_high,
    
    p_reach_low =
      probability_reached_trigger_low,
    
    prob_high =
      p_reach_high /
      (
        p_reach_high +
          p_reach_low
      ),
    
    prob_high_percent =
      100 * prob_high
  )


# ============================================================
# UNCERTAINTY IN PROBABILITY HIGH
#
# Draw uncertainty in each attainment probability using
# Jeffreys Beta distributions, then propagate into P(HIGH).
# ============================================================

B <- 10000

high_probability_intervals <- lapply(
  seq_len(nrow(panelA_summary)),
  function(i) {
    
    p_high_draw <- rbeta(
      B,
      panelA_summary$n_reached_trigger_high[i] + 0.5,
      panelA_summary$n_generated_high[i] -
        panelA_summary$n_reached_trigger_high[i] + 0.5
    )
    
    p_low_draw <- rbeta(
      B,
      panelA_summary$n_reached_trigger_low[i] + 0.5,
      panelA_summary$n_generated_low[i] -
        panelA_summary$n_reached_trigger_low[i] + 0.5
    )
    
    prob_high_draw <-
      p_high_draw /
      (
        p_high_draw +
          p_low_draw
      )
    
    tibble(
      x_lower =
        100 * quantile(
          prob_high_draw,
          0.025,
          names = FALSE,
          na.rm = TRUE
        ),
      
      x_upper =
        100 * quantile(
          prob_high_draw,
          0.975,
          names = FALSE,
          na.rm = TRUE
        )
    )
  }
) %>%
  bind_rows()


panelA_summary <- bind_cols(
  panelA_summary,
  high_probability_intervals
)


print(
  panelA_summary %>%
    select(
      reported_case_trigger,
      p_reach_high,
      p_reach_low,
      prob_high_percent,
      x_lower,
      x_upper
    ),
  n = Inf
)


# ============================================================
# PART 2
# OUTBREAK CONSEQUENCES IF INTERVENTION WAITS UNTIL EACH
# REPORTED-CASE TRIGGER
#
# HIGH-transmission scenario only
# q = 0.95
# conditional on reaching the requested trigger
# ============================================================

fig3_runs <- read_csv(
  file.path(
    fig3_dir,
    "tables",
    "reported_case_trigger_run_level_results.csv"
  ),
  show_col_types = FALSE
)


high_trigger_runs <- fig3_runs %>%
  filter(
    transmission_regime == "supercritical",
    abs(q_system - 0.95) < 1e-8,
    reached_trigger == 1,
    reported_case_trigger %in% target_triggers
  )


# ------------------------------------------------------------
# Exact binomial CI helper
# ------------------------------------------------------------

binom_ci <- function(events, n) {
  
  if (n == 0) {
    return(
      c(
        lower = NA_real_,
        upper = NA_real_
      )
    )
  }
  
  ci <- binom.test(
    events,
    n
  )$conf.int
  
  c(
    lower = unname(ci[1]),
    upper = unname(ci[2])
  )
}


# ------------------------------------------------------------
# Calculate >10, >=50 and >=100 infection probabilities
# ------------------------------------------------------------

panelB_summary <- high_trigger_runs %>%
  group_by(
    reported_case_trigger
  ) %>%
  summarise(
    
    n = n(),
    
    events_gt10 =
      sum(
        total_retained_infections > 10,
        na.rm = TRUE
      ),
    
    events_50 =
      sum(
        total_retained_infections >= 50,
        na.rm = TRUE
      ),
    
    events_100 =
      sum(
        total_retained_infections >= 100,
        na.rm = TRUE
      ),
    
    p_gt10 =
      mean(
        total_retained_infections > 10,
        na.rm = TRUE
      ),
    
    p_50 =
      mean(
        total_retained_infections >= 50,
        na.rm = TRUE
      ),
    
    p_100 =
      mean(
        total_retained_infections >= 100,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  ) %>%
  
  rowwise() %>%
  
  mutate(
    
    pgt10_lower =
      binom_ci(
        events_gt10,
        n
      )[["lower"]],
    
    pgt10_upper =
      binom_ci(
        events_gt10,
        n
      )[["upper"]],
    
    p50_lower =
      binom_ci(
        events_50,
        n
      )[["lower"]],
    
    p50_upper =
      binom_ci(
        events_50,
        n
      )[["upper"]],
    
    p100_lower =
      binom_ci(
        events_100,
        n
      )[["lower"]],
    
    p100_upper =
      binom_ci(
        events_100,
        n
      )[["upper"]]
  ) %>%
  
  ungroup()


print(
  panelB_summary,
  n = Inf
)


# ============================================================
# PART 3
# COMBINE THE TWO ANALYSES
#
# X = probability higher-transmission scenario
# Y = probability of outbreak threshold
# ============================================================

tradeoff_data <- panelA_summary %>%
  
  select(
    reported_case_trigger,
    prob_high_percent,
    x_lower,
    x_upper
  ) %>%
  
  left_join(
    panelB_summary %>%
      select(
        reported_case_trigger,
        
        p_gt10,
        pgt10_lower,
        pgt10_upper,
        
        p_50,
        p50_lower,
        p50_upper,
        
        p_100,
        p100_lower,
        p100_upper
      ),
    by = "reported_case_trigger"
  ) %>%
  
  pivot_longer(
    cols = c(
      p_gt10,
      p_50,
      p_100
    ),
    names_to = "outcome",
    values_to = "probability"
  ) %>%
  
  mutate(
    
    y_lower = case_when(
      outcome == "p_gt10" ~
        100 * pgt10_lower,
      
      outcome == "p_50" ~
        100 * p50_lower,
      
      outcome == "p_100" ~
        100 * p100_lower
    ),
    
    y_upper = case_when(
      outcome == "p_gt10" ~
        100 * pgt10_upper,
      
      outcome == "p_50" ~
        100 * p50_upper,
      
      outcome == "p_100" ~
        100 * p100_upper
    ),
    
    outcome = case_when(
      outcome == "p_gt10" ~
        ">10 infections",
      
      outcome == "p_50" ~
        "\u226550 infections",
      
      outcome == "p_100" ~
        "\u2265100 infections"
    ),
    
    x =
      prob_high_percent,
    
    y =
      100 * probability
  ) %>%
  
  mutate(
    outcome = factor(
      outcome,
      levels = c(
        ">10 infections",
        "\u226550 infections",
        "\u2265100 infections"
      )
    )
  )


print(
  tradeoff_data,
  n = Inf
)


# ============================================================
# SAVE NUMBERS
# ============================================================

write_csv(
  tradeoff_data,
  file.path(
    table_dir,
    "certainty_vs_outbreak_risk_tradeoff.csv"
  )
)


# ============================================================
# PART 4
# CLEAN TRADE-OFF PLOT
#
# NO REPORTED-CASE NUMBERS PRINTED ON GRAPH
# ============================================================

p_tradeoff_clean <- ggplot(
  tradeoff_data,
  aes(
    x = x,
    y = y,
    group = outcome
  )
) +
  
  # ----------------------------------------------------------
# Horizontal uncertainty:
# uncertainty in probability of HIGH transmission
# ----------------------------------------------------------

geom_errorbarh(
  aes(
    xmin = x_lower,
    xmax = x_upper,
    colour = outcome
  ),
  height = 0,
  linewidth = 0.45,
  alpha = 0.75
) +
  
  # ----------------------------------------------------------
# Vertical uncertainty:
# uncertainty in outbreak-size probability
# ----------------------------------------------------------

geom_errorbar(
  aes(
    ymin = y_lower,
    ymax = y_upper,
    colour = outcome
  ),
  width = 0,
  linewidth = 0.45,
  alpha = 0.75
) +
  
  # ----------------------------------------------------------
# Curves
# ----------------------------------------------------------

geom_line(
  aes(
    colour = outcome,
    linetype = outcome
  ),
  linewidth = 1.1
) +
  
  geom_point(
    aes(
      colour = outcome,
      shape = outcome
    ),
    size = 2.8
  ) +
  
  # ----------------------------------------------------------
# Colours
# ----------------------------------------------------------

scale_colour_manual(
  values = c(
    ">10 infections" =
      "#009E73",
    
    "\u226550 infections" =
      "#D55E00",
    
    "\u2265100 infections" =
      "#CC79A7"
  ),
  name = NULL
) +
  
  scale_linetype_manual(
    values = c(
      ">10 infections" =
        "solid",
      
      "\u226550 infections" =
        "longdash",
      
      "\u2265100 infections" =
        "dotted"
    ),
    name = NULL
  ) +
  
  scale_shape_manual(
    values = c(
      ">10 infections" = 16,
      "\u226550 infections" = 17,
      "\u2265100 infections" = 15
    ),
    name = NULL
  ) +
  
  # ----------------------------------------------------------
# X axis
# ----------------------------------------------------------

scale_x_continuous(
  limits = c(
    45,
    100
  ),
  breaks = seq(
    50,
    100,
    by = 10
  ),
  labels = function(x) {
    paste0(
      x,
      "%"
    )
  }
) +
  
  # ----------------------------------------------------------
# Y axis
# ----------------------------------------------------------

scale_y_continuous(
  limits = c(
    0,
    100
  ),
  breaks = seq(
    0,
    100,
    by = 20
  ),
  labels = function(x) {
    paste0(
      x,
      "%"
    )
  }
) +
  
  # ----------------------------------------------------------
# Labels
# ----------------------------------------------------------

labs(
  x =
    "Probability outbreak is the higher-transmission scenario",
  
  y =
    "Probability of outbreak size threshold",
  
  title =
    "Increasing certainty about high transmission versus outbreak risk",
  
  subtitle =
    paste0(
      "Outbreak risk when intervention is delayed as evidence for higher transmission accumulates; ",
      "probability of the higher-transmission scenario assumes equal prior probability of the two modeled scenarios."
    )
) +
  
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    legend.position =
      "bottom",
    
    legend.title =
      element_blank()
  )


# ============================================================
# DISPLAY
# ============================================================

print(
  p_tradeoff_clean
)


# ============================================================
# SAVE FIGURE
# ============================================================

ggsave(
  file.path(
    figure_dir,
    "Figure_certainty_vs_waiting_tradeoff_CLEAN.png"
  ),
  p_tradeoff_clean,
  width = 8.8,
  height = 6.2,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(
    figure_dir,
    "Figure_certainty_vs_waiting_tradeoff_CLEAN.pdf"
  ),
  p_tradeoff_clean,
  width = 8.8,
  height = 6.2,
  bg = "white"
)


message(
  "DONE -- self-contained add-on completed; no simulations rerun."
)