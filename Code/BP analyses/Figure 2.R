# Figure 2 code
setwd("~/Hantavirus-paper/TD results")

parse_args <- function() {
  x <- list(
    input = getwd(),
    out = "~/Hantavirus-paper/BP results/Fig 2",
    n_runs = 10000L,
    chunk_size = 1000L,
    seed = 20260519L,
    q_values = c(1.00, 0.95, 0.90, 0.75, 0.50),
    current_generation = 2L,
    
    response_day = 35.0,
    response_lag_from_first_reported_onset = 26.0,
    align_response_to_first_reported_onset = TRUE,
    
    exposed_cohort_size = 100L,
    
    backward_tracing_window = 14.0,
    daily_contact_rate = 13.50,
    contact_count_dispersion = 2.11,
    contact_multiplier = 1.0,
    notification_delay_mean = 1.0,
    notification_delay_sd = 0.5,
    
    max_days = 1000.0,
    max_cases = 100000L,
    initial_infected = 1L,
    
    p_detect = 0.95,
    detection_delay_mean = 2.0,
    detection_delay_sd = 1.0,
    isolation_delay_mean = 0.5,
    isolation_delay_sd = 0.5,
    isolation_effectiveness = 0.99,
    
    quarantine_duration = 42.0,
    quarantine_effectiveness = 1.0,
    quarantine_adherence = 1.0,
    
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
    if (key %in% c("--no_align_response", "--no-align-response")) {
      x$align_response_to_first_reported_onset <- FALSE
      i <- i + 1L
      next
    }
    
    if (i == length(a)) stop("Missing value for ", key, call. = FALSE)
    val <- a[[i + 1L]]
    
    if (key == "--input") x$input <- val
    else if (key %in% c("--out", "--outdir")) x$out <- val
    else if (key %in% c("--n_runs", "--n-runs")) x$n_runs <- as.integer(val)
    else if (key %in% c("--chunk_size", "--chunk-size")) {
      x$chunk_size <- as.integer(val)
    }
    else if (key == "--seed") x$seed <- as.integer(val)
    else if (key %in% c("--q_values", "--q-values")) {
      x$q_values <- as.numeric(strsplit(val, ",", fixed = TRUE)[[1]])
    }
    else if (key %in% c("--current_generation", "--current-generation")) {
      x$current_generation <- as.integer(val)
    }
    else if (key %in% c("--response_day", "--response-day")) {
      x$response_day <- as.numeric(val)
    }
    else if (key %in% c("--response_lag_from_first_reported_onset",
                        "--response-lag-from-first-reported-onset")) {
      x$response_lag_from_first_reported_onset <- as.numeric(val)
    }
    else if (key %in% c("--exposed_cohort_size", "--exposed-cohort-size")) {
      x$exposed_cohort_size <- as.integer(val)
    }
    else if (key %in% c("--backward_tracing_window", "--backward-tracing-window")) {
      x$backward_tracing_window <- as.numeric(val)
    }
    else if (key %in% c("--daily_contact_rate", "--daily-contact-rate")) {
      x$daily_contact_rate <- as.numeric(val)
    }
    else if (key %in% c("--contact_count_dispersion", "--contact-count-dispersion")) {
      x$contact_count_dispersion <- as.numeric(val)
    }
    else if (key %in% c("--contact_multiplier", "--contact-multiplier")) {
      x$contact_multiplier <- as.numeric(val)
    }
    else if (key %in% c("--notification_delay_mean", "--notification-delay-mean")) {
      x$notification_delay_mean <- as.numeric(val)
    }
    else if (key %in% c("--notification_delay_sd", "--notification-delay-sd")) {
      x$notification_delay_sd <- as.numeric(val)
    }
    else if (key %in% c("--max_days", "--max-days")) x$max_days <- as.numeric(val)
    else if (key %in% c("--max_cases", "--max-cases")) x$max_cases <- as.integer(val)
    else if (key %in% c("--initial_infected", "--initial-infected")) {
      x$initial_infected <- as.integer(val)
    }
    else if (key == "--p_detect") x$p_detect <- as.numeric(val)
    else if (key %in% c("--detection_delay_mean", "--detection-delay-mean")) {
      x$detection_delay_mean <- as.numeric(val)
    }
    else if (key %in% c("--detection_delay_sd", "--detection-delay-sd")) {
      x$detection_delay_sd <- as.numeric(val)
    }
    else if (key %in% c("--isolation_delay_mean", "--isolation-delay-mean")) {
      x$isolation_delay_mean <- as.numeric(val)
    }
    else if (key %in% c("--isolation_delay_sd", "--isolation-delay-sd")) {
      x$isolation_delay_sd <- as.numeric(val)
    }
    else if (key %in% c("--isolation_effectiveness", "--isolation-effectiveness")) {
      x$isolation_effectiveness <- as.numeric(val)
    }
    else if (key %in% c("--quarantine_duration", "--quarantine-duration")) {
      x$quarantine_duration <- as.numeric(val)
    }
    else if (key %in% c("--quarantine_effectiveness", "--quarantine-effectiveness")) {
      x$quarantine_effectiveness <- as.numeric(val)
    }
    else if (key %in% c("--quarantine_adherence", "--quarantine-adherence")) {
      x$quarantine_adherence <- as.numeric(val)
    }
    else stop("Unknown argument: ", key, call. = FALSE)
    
    i <- i + 2L
  }
  
  x$q_values <- x$q_values[
    is.finite(x$q_values) & x$q_values >= 0 & x$q_values <= 1
  ]
  if (length(x$q_values) == 0) stop("No valid q_values supplied.", call. = FALSE)
  
  x
}

required <- c(
  "Rcpp", "dplyr", "tidyr", "ggplot2", "readr",
  "tibble", "scales", "patchwork"
)
missing <- required[
  !vapply(required, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing) > 0) {
  stop(
    "Missing packages: ", paste(missing, collapse = ", "),
    "\nInstall with install.packages(c(",
    paste(sprintf('"%s"', missing), collapse = ", "),
    "))",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(Rcpp)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(tibble)
  library(scales)
  library(patchwork)
})

args <- parse_args()
set.seed(args$seed)

dir.create(args$out, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(args$out, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(args$out, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(args$out, "chunks"), recursive = TRUE, showWarnings = FALSE)

binom_ci <- function(events, n) {
  if (!is.finite(events) || !is.finite(n) || n <= 0) {
    return(c(lower = NA_real_, upper = NA_real_))
  }
  ci <- stats::binom.test(events, n)$conf.int
  c(lower = unname(ci[1]), upper = unname(ci[2]))
}

quant <- function(x, p) {
  as.numeric(stats::quantile(x, p, na.rm = TRUE, names = FALSE))
}

make_q_label <- function(q) paste0(round(100 * q), "%")
make_missed_label <- function(q) paste0(round(100 * (1 - q)), "%")

# Read joint posterior timing draws
delay_candidates <- c(
  file.path(args$input, "tables", "joint_delay_draws_for_branching.csv"),
  file.path(args$input, "joint_delay_draws_for_branching.csv"),
  file.path(args$input, "posterior_pair_level_draws.csv"),
  file.path(args$input, "tables", "posterior_pair_level_draws.csv")
)

delay_file <- delay_candidates[file.exists(delay_candidates)][1]
if (is.na(delay_file) || length(delay_file) == 0) {
  stop(
    "Could not find a supported delay-draw file under --input: ",
    args$input,
    call. = FALSE
  )
}

message("Reading delay draws from: ",
        normalizePath(delay_file, winslash = "/", mustWork = FALSE))

pair_draws_raw <- readr::read_csv(delay_file, show_col_types = FALSE)

if ("generation_time" %in% names(pair_draws_raw)) {
  pair_draws_raw$generation_time <- as.numeric(pair_draws_raw$generation_time)
} else if ("generation_interval" %in% names(pair_draws_raw)) {
  pair_draws_raw$generation_time <- as.numeric(
    pair_draws_raw$generation_interval
  )
} else {
  stop("No generation-time column found.", call. = FALSE)
}

if ("recipient_incubation" %in% names(pair_draws_raw)) {
  pair_draws_raw$recipient_incubation <- as.numeric(
    pair_draws_raw$recipient_incubation
  )
} else if ("incubation" %in% names(pair_draws_raw)) {
  pair_draws_raw$recipient_incubation <- as.numeric(
    pair_draws_raw$incubation
  )
} else {
  stop("No recipient-incubation column found.", call. = FALSE)
}

pair_draws_joint <- pair_draws_raw %>%
  filter(
    is.finite(generation_time),
    generation_time > 0,
    is.finite(recipient_incubation),
    recipient_incubation > 0
  ) %>%
  transmute(
    generation_time,
    recipient_incubation
  )

if (nrow(pair_draws_joint) == 0) {
  stop("No usable joint timing draws found.", call. = FALSE)
}

readr::write_csv(
  pair_draws_joint,
  file.path(args$out, "tables", "joint_timing_draws_used.csv")
)

generation_draws <- pair_draws_joint$generation_time
incubation_draws <- pair_draws_joint$recipient_incubation

# C++ simulation core
# Detection, isolation, and contact-notification delays are gamma distributed.
# Gamma shape and scale are derived internally from the supplied mean and SD.
Rcpp::sourceCpp(code = '
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <queue>
#include <string>
#include <vector>
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

inline double runif01() {
  return R::runif(0.0, 1.0);
}

inline bool finite_num(double x) {
  return R_finite(x);
}

inline double gamma_delay_one(
    double mean,
    double sd,
    double min_value = 0.0
) {
  if (!R_finite(mean) || mean <= 0.0) return min_value;
  if (!R_finite(sd) || sd <= 0.0) return std::max(mean, min_value);

  double shape = (mean / sd) * (mean / sd);
  double scale = (sd * sd) / mean;
  double y = R::rgamma(shape, scale);

  if (!R_finite(y) || y < min_value) {
    y = std::max(mean, min_value);
  }
  return y;
}

inline int sample_row(int n) {
  int idx = (int)std::floor(R::runif(0.0, (double)n));
  if (idx < 0) idx = 0;
  if (idx >= n) idx = n - 1;
  return idx;
}

struct SimEvent {
  double time;
  int type;       // 1=detection, 2=potential transmission
  int source;     // case id
  double child_incubation;
};

struct EarlierEvent {
  bool operator()(const SimEvent& a, const SimEvent& b) const {
    if (a.time != b.time) return a.time > b.time;
    if (a.type != b.type) return a.type > b.type;
    return a.source > b.source;
  }
};

// [[Rcpp::export]]
DataFrame run_dynamic_tracing_cpp(
    NumericVector q_values,
    int n_runs,
    int run_offset,
    NumericVector generation_draws,
    NumericVector incubation_draws,
    int current_generation,
    double response_day,
    double response_lag_from_first_reported_onset,
    bool align_response_to_first_reported_onset,
    int exposed_cohort_size,
    double backward_tracing_window,
    double daily_contact_rate,
    double contact_count_dispersion,
    double contact_multiplier,
    double notification_delay_mean,
    double notification_delay_sd,
    int max_cases,
    double max_days,
    int initial_infected,
    double R0,
    double k,
    double p_detect,
    double detection_delay_mean,
    double detection_delay_sd,
    double isolation_delay_mean,
    double isolation_delay_sd,
    double isolation_effectiveness,
    double quarantine_duration,
    double quarantine_effectiveness,
    double quarantine_adherence,
    int seed
) {
  RNGScope scope;
  Environment base_env("package:base");
  Function set_seed = base_env["set.seed"];
  set_seed(seed);

  CharacterVector scenario_names = CharacterVector::create(
    "perfect_intervention",
    "imperfect_intervention"
  );

  int n_scenarios = scenario_names.size();
  int n_q = q_values.size();
  int total_rows = n_scenarios * n_q * n_runs;
  int n_joint = generation_draws.size();

  IntegerVector out_run_id(total_rows);
  CharacterVector out_scenario(total_rows);
  NumericVector out_q(total_rows);
  NumericVector out_contacts_missed(total_rows);
  IntegerVector out_reached(total_rows);
  NumericVector out_response_time(total_rows);
  NumericVector out_calendar_shift(total_rows);
  NumericVector out_first_reported_onset(total_rows);

  IntegerVector out_total_potential(total_rows);
  IntegerVector out_total_retained(total_rows);
  IntegerVector out_max_generation_retained(total_rows);

  IntegerVector out_reached_case_limit(total_rows);
  IntegerVector out_reached_time_limit(total_rows);
  IntegerVector out_proposed_after_time_limit(total_rows);

  IntegerVector out_known_infected_at_response(total_rows);
  IntegerVector out_quarantined_known_infected(total_rows);
  IntegerVector out_missed_known_infected(total_rows);
  IntegerVector out_uninfected_initial_quarantined(total_rows);
  IntegerVector out_exposed_cohort_too_small(total_rows);

  IntegerVector out_detected_cases_traced(total_rows);
  IntegerVector out_infected_contacts_traced(total_rows);
  IntegerVector out_uninfected_contacts_traced(total_rows);

  IntegerVector out_post_response_transmissions(total_rows);
  IntegerVector out_post_response_transmissions_outside(total_rows);
  IntegerVector out_any_leakage(total_rows);
  IntegerVector out_post_response_infections(total_rows);
  IntegerVector out_additional_after_current_generation(total_rows);
  IntegerVector out_additional_after_response(total_rows);

  IntegerVector out_prevented_by_quarantine(total_rows);
  IntegerVector out_prevented_by_isolation(total_rows);
  IntegerVector out_max_secondary(total_rows);

  NumericVector out_response_duration(total_rows);
  NumericVector out_days_to_first_leakage(total_rows);

  int row = 0;
  double nb_prob = k / (k + R0);

  for (int sc = 0; sc < n_scenarios; ++sc) {
    std::string scenario = as<std::string>(scenario_names[sc]);

    double pdet_sc = p_detect;
    double det_mean_sc = detection_delay_mean;
    double det_sd_sc = detection_delay_sd;
    double iso_delay_mean_sc = isolation_delay_mean;
    double iso_delay_sd_sc = isolation_delay_sd;
    double iso_eff_sc = isolation_effectiveness;
    double q_eff_sc = quarantine_effectiveness;
    double q_adh_sc = quarantine_adherence;

    if (scenario == "perfect_intervention") {
      pdet_sc = 1.0;
      det_mean_sc = 0.0;
      det_sd_sc = 0.0;
      iso_delay_mean_sc = 0.0;
      iso_delay_sd_sc = 0.0;
      iso_eff_sc = 1.0;
      q_eff_sc = 1.0;
      q_adh_sc = 1.0;
    }

    for (int qidx = 0; qidx < n_q; ++qidx) {
      double q_system = q_values[qidx];

      for (int rr = 0; rr < n_runs; ++rr) {
        bool reached_case_limit = false;
        bool reached_time_limit = false;
        int proposed_after_time_limit = 0;

        std::vector<int> infector;
        std::vector<int> generation;
        std::vector<double> infection_time;
        std::vector<double> incubation_period;
        std::vector<double> symptom_onset;
        std::vector<int> detected;
        std::vector<double> detection_time;
        std::vector<double> isolation_time;
        std::vector<int> quarantined;
        std::vector<double> q_start;
        std::vector<double> q_end;
        std::vector<int> retained_offspring;
        std::vector< std::vector<int> > children;

        int reserve_n = std::min(max_cases, 10000);
        infector.reserve(reserve_n);
        generation.reserve(reserve_n);
        infection_time.reserve(reserve_n);
        incubation_period.reserve(reserve_n);
        symptom_onset.reserve(reserve_n);
        detected.reserve(reserve_n);
        detection_time.reserve(reserve_n);
        isolation_time.reserve(reserve_n);
        quarantined.reserve(reserve_n);
        q_start.reserve(reserve_n);
        q_end.reserve(reserve_n);
        retained_offspring.reserve(reserve_n);
        children.reserve(reserve_n);

        std::priority_queue<
          SimEvent,
          std::vector<SimEvent>,
          EarlierEvent
        > events;

        double first_reported_onset_base = R_PosInf;
        double response_time_base = R_PosInf;
        bool response_started = false;
        bool reached = false;

        int known_inf = 0;
        int q_known_inf = 0;
        int missed_known_inf = 0;
        int uninfected_initial_quarantined = 0;
        int exposed_cohort_too_small = 0;

        int detected_cases_traced = 0;
        int infected_contacts_traced = 0;
        int uninfected_contacts_traced = 0;

        int post_response_transmissions = 0;
        int post_response_transmissions_outside = 0;
        int prevented_q_total = 0;
        int prevented_iso_total = 0;
        double first_leakage_time_base = R_PosInf;

        auto schedule_case_events = [&](int case_id) {
          if (detected[case_id] == 1 &&
              finite_num(detection_time[case_id])) {
            SimEvent de;
            de.time = detection_time[case_id];
            de.type = 1;
            de.source = case_id;
            de.child_incubation = NA_REAL;
            events.push(de);
          }

          int z0 = (int)R::rnbinom(k, nb_prob);
          for (int zz = 0; zz < z0; ++zz) {
            int draw_idx = sample_row(n_joint);
            double G = generation_draws[draw_idx];
            double child_inc = incubation_draws[draw_idx];
            double tx_time = infection_time[case_id] + G;

            if (!finite_num(tx_time)) continue;

            SimEvent te;
            te.time = tx_time;
            te.type = 2;
            te.source = case_id;
            te.child_incubation = child_inc;
            events.push(te);
          }
        };

        auto add_case = [&](int inf, int gen, double itime,
                            double incubation_value) {
          int id = (int)infection_time.size();
          double onset = itime + incubation_value;
          int det = (runif01() < pdet_sc) ? 1 : 0;

          double dtime = R_PosInf;
          double isotime = R_PosInf;

          if (det == 1) {
            dtime = onset + gamma_delay_one(
              det_mean_sc, det_sd_sc, 0.0
            );
            isotime = dtime + gamma_delay_one(
              iso_delay_mean_sc, iso_delay_sd_sc, 0.0
            );

            if (onset < first_reported_onset_base) {
              first_reported_onset_base = onset;
              response_time_base =
                onset + response_lag_from_first_reported_onset;
            }
          }

          infector.push_back(inf);
          generation.push_back(gen);
          infection_time.push_back(itime);
          incubation_period.push_back(incubation_value);
          symptom_onset.push_back(onset);
          detected.push_back(det);
          detection_time.push_back(dtime);
          isolation_time.push_back(isotime);
          quarantined.push_back(0);
          q_start.push_back(R_PosInf);
          q_end.push_back(R_PosInf);
          retained_offspring.push_back(0);
          children.push_back(std::vector<int>());

          schedule_case_events(id);
          return id;
        };

        for (int ii = 0; ii < initial_infected; ++ii) {
          int draw_idx = sample_row(n_joint);
          add_case(
            NA_INTEGER,
            0,
            0.0,
            incubation_draws[draw_idx]
          );
        }

        while (true) {
          double next_event_time =
            events.empty() ? R_PosInf : events.top().time;

          // Response is inserted dynamically before the next later event.
          if (!response_started &&
              finite_num(response_time_base) &&
              response_time_base <= next_event_time) {
            response_started = true;

            for (int ii = 0; ii < (int)infection_time.size(); ++ii) {
              if (infection_time[ii] <= response_time_base) {
                known_inf++;

                if (generation[ii] == current_generation) {
                  reached = true;
                }

                if (runif01() <= q_system) {
                  quarantined[ii] = 1;
                  q_start[ii] = response_time_base;
                  q_end[ii] =
                    response_time_base + quarantine_duration;
                  q_known_inf++;
                } else {
                  missed_known_inf++;
                }
              }
            }

            if (known_inf > exposed_cohort_size) {
              exposed_cohort_too_small = 1;
            }

            int uninfected_exposed =
              exposed_cohort_size - known_inf;
            if (uninfected_exposed < 0) {
              uninfected_exposed = 0;
            }

            if (uninfected_exposed > 0) {
              uninfected_initial_quarantined =
                (int)R::rbinom(
                  (double)uninfected_exposed,
                  q_system
                );
            }

            continue;
          }

          if (events.empty()) break;

          SimEvent ev = events.top();
          events.pop();

          if (ev.time > max_days) {
            reached_time_limit = true;
            proposed_after_time_limit++;
            break;
          }

          if (ev.type == 1) {
            int s = ev.source;

            if (!response_started) continue;
            if (!reached) continue;
            if (s < 0 || s >= (int)infection_time.size()) continue;
            if (ev.time < response_time_base) continue;

            detected_cases_traced++;

            double trace_start = std::max(
              symptom_onset[s] - backward_tracing_window,
              infection_time[s]
            );
            double trace_end = ev.time;

            std::vector<int> infected_traceable_contacts;
            infected_traceable_contacts.reserve(children[s].size());

            for (int kk = 0; kk < (int)children[s].size(); ++kk) {
              int c = children[s][kk];

              if (infection_time[c] >= trace_start &&
                  infection_time[c] <= trace_end) {
                infected_traceable_contacts.push_back(c);
              }
            }

            int n_infected_traceable =
              (int)infected_traceable_contacts.size();

            for (int kk = 0;
                 kk < n_infected_traceable;
                 ++kk) {
              int c = infected_traceable_contacts[kk];

              if (runif01() <= q_system) {
                double notification_delay =
                  gamma_delay_one(
                    notification_delay_mean,
                    notification_delay_sd,
                    0.0
                  );

                double proposed_q_start =
                  ev.time + notification_delay;

                if (!finite_num(q_start[c]) ||
                    proposed_q_start < q_start[c]) {
                  quarantined[c] = 1;
                  q_start[c] = proposed_q_start;
                  q_end[c] =
                    proposed_q_start + quarantine_duration;
                }

                infected_contacts_traced++;
              }
            }

            double trace_duration =
              std::max(0.0, trace_end - trace_start);

            double mean_contacts =
              daily_contact_rate *
              contact_multiplier *
              trace_duration;

            double size_contacts =
              contact_count_dispersion *
              trace_duration;

            int total_traceable_contacts =
              n_infected_traceable;

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
                  n_infected_traceable) {
                total_traceable_contacts =
                  n_infected_traceable;
              }
            }

            int uninfected_contacts =
              total_traceable_contacts -
              n_infected_traceable;

            if (uninfected_contacts > 0) {
              uninfected_contacts_traced +=
                (int)R::rbinom(
                  (double)uninfected_contacts,
                  q_system
                );
            }

            continue;
          }

          if (ev.type == 2) {
            int s = ev.source;
            if (s < 0 || s >= (int)infection_time.size()) continue;

            bool after_response =
              response_started &&
              finite_num(response_time_base) &&
              ev.time >= response_time_base;

            bool source_quarantined =
              after_response &&
              quarantined[s] == 1 &&
              finite_num(q_start[s]) &&
              ev.time >= q_start[s] &&
              ev.time <= q_end[s];

            bool source_isolated =
              after_response &&
              detected[s] == 1 &&
              finite_num(isolation_time[s]) &&
              ev.time >= isolation_time[s];

            double retain_prob = 1.0;

            if (after_response) {
              if (source_quarantined &&
                  runif01() < q_adh_sc) {
                retain_prob *= (1.0 - q_eff_sc);
              }

              if (source_isolated) {
                retain_prob *= (1.0 - iso_eff_sc);
              }
            }

            if (runif01() > retain_prob) {
              if (source_quarantined) {
                prevented_q_total++;
              } else if (source_isolated) {
                prevented_iso_total++;
              }
              continue;
            }

            if ((int)infection_time.size() >= max_cases) {
              reached_case_limit = true;
              break;
            }

            int child = add_case(
              s,
              generation[s] + 1,
              ev.time,
              ev.child_incubation
            );

            children[s].push_back(child);
            retained_offspring[s] += 1;

            if (after_response) {
              post_response_transmissions++;

              bool recipient_quarantined_at_infection =
                quarantined[child] == 1 &&
                finite_num(q_start[child]) &&
                ev.time >= q_start[child] &&
                ev.time <= q_end[child];

              if (!source_quarantined &&
                  !recipient_quarantined_at_infection) {
                post_response_transmissions_outside++;

                if (ev.time < first_leakage_time_base) {
                  first_leakage_time_base = ev.time;
                }
              }
            }
          }
        }

        int total_retained = (int)infection_time.size();
        int max_generation_retained = 0;
        int post_response_infections = 0;
        int additional_after_cg = 0;
        int additional_after_response = 0;
        int max_secondary = 0;

        double response_end_base =
          finite_num(response_time_base) ?
          response_time_base + quarantine_duration :
          R_PosInf;

        for (int ii = 0; ii < total_retained; ++ii) {
          if (generation[ii] > max_generation_retained) {
            max_generation_retained = generation[ii];
          }

          if (generation[ii] > current_generation) {
            additional_after_cg++;
          }

          if (response_started &&
              infection_time[ii] >= response_time_base) {
            post_response_infections++;
            additional_after_response++;
          }

          if (retained_offspring[ii] > max_secondary) {
            max_secondary = retained_offspring[ii];
          }

          if (response_started) {
            if (infection_time[ii] > response_end_base) {
              response_end_base = infection_time[ii];
            }
            if (symptom_onset[ii] > response_end_base) {
              response_end_base = symptom_onset[ii];
            }
            if (finite_num(detection_time[ii]) &&
                detection_time[ii] > response_end_base) {
              response_end_base = detection_time[ii];
            }
            if (finite_num(q_end[ii]) &&
                q_end[ii] > response_end_base) {
              response_end_base = q_end[ii];
            }
          }
        }

        double calendar_shift = 0.0;
        double first_reported_onset = R_PosInf;
        double reported_response_time = response_day;

        if (finite_num(first_reported_onset_base)) {
          if (align_response_to_first_reported_onset) {
            calendar_shift =
              response_day -
              response_lag_from_first_reported_onset -
              first_reported_onset_base;
            first_reported_onset =
              first_reported_onset_base + calendar_shift;
            reported_response_time = response_day;
          } else {
            first_reported_onset = first_reported_onset_base;
            reported_response_time = response_time_base;
          }
        }

        double response_duration = NA_REAL;
        if (response_started && reached) {
          response_duration =
            std::max(
              0.0,
              response_end_base - response_time_base
            );
        }

        out_run_id[row] = run_offset + rr + 1;
        out_scenario[row] = scenario;
        out_q[row] = q_system;
        out_contacts_missed[row] = 1.0 - q_system;
        out_reached[row] = reached ? 1 : 0;
        out_response_time[row] = reported_response_time;
        out_calendar_shift[row] = calendar_shift;
        out_first_reported_onset[row] =
          finite_num(first_reported_onset) ?
          first_reported_onset :
          NA_REAL;

        out_total_potential[row] = total_retained;
        out_total_retained[row] = total_retained;
        out_max_generation_retained[row] =
          max_generation_retained;

        out_reached_case_limit[row] =
          reached_case_limit ? 1 : 0;
        out_reached_time_limit[row] =
          reached_time_limit ? 1 : 0;
        out_proposed_after_time_limit[row] =
          proposed_after_time_limit;

        out_known_infected_at_response[row] = known_inf;
        out_quarantined_known_infected[row] =
          q_known_inf;
        out_missed_known_infected[row] =
          missed_known_inf;
        out_uninfected_initial_quarantined[row] =
          uninfected_initial_quarantined;
        out_exposed_cohort_too_small[row] =
          exposed_cohort_too_small;

        out_detected_cases_traced[row] =
          detected_cases_traced;
        out_infected_contacts_traced[row] =
          infected_contacts_traced;
        out_uninfected_contacts_traced[row] =
          uninfected_contacts_traced;

        out_post_response_transmissions[row] =
          post_response_transmissions;
        out_post_response_transmissions_outside[row] =
          post_response_transmissions_outside;
        out_any_leakage[row] =
          post_response_transmissions_outside > 0 ? 1 : 0;
        out_post_response_infections[row] =
          post_response_infections;
        out_additional_after_current_generation[row] =
          additional_after_cg;
        out_additional_after_response[row] =
          additional_after_response;

        out_prevented_by_quarantine[row] =
          prevented_q_total;
        out_prevented_by_isolation[row] =
          prevented_iso_total;
        out_max_secondary[row] = max_secondary;

        out_response_duration[row] = response_duration;
        out_days_to_first_leakage[row] =
          reached &&
          finite_num(first_leakage_time_base) &&
          finite_num(response_time_base) ?
          first_leakage_time_base - response_time_base :
          NA_REAL;

        row++;
      }
    }
  }

  return DataFrame::create(
    Named("run_id") = out_run_id,
    Named("scenario") = out_scenario,
    Named("q_system") = out_q,
    Named("contacts_missed") = out_contacts_missed,
    Named("reached_current_generation") = out_reached,
    Named("response_time") = out_response_time,
    Named("calendar_shift") = out_calendar_shift,
    Named("first_reported_onset") =
      out_first_reported_onset,

    Named("total_potential_infections") =
      out_total_potential,
    Named("total_retained_infections") =
      out_total_retained,
    Named("max_generation_retained") =
      out_max_generation_retained,

    Named("reached_case_limit") =
      out_reached_case_limit,
    Named("reached_time_limit") =
      out_reached_time_limit,
    Named("proposed_infections_after_time_limit") =
      out_proposed_after_time_limit,

    Named("known_infected_at_response") =
      out_known_infected_at_response,
    Named("quarantined_known_infected_at_response") =
      out_quarantined_known_infected,
    Named("missed_known_infected_at_response") =
      out_missed_known_infected,
    Named("uninfected_initial_cohort_quarantined") =
      out_uninfected_initial_quarantined,
    Named("exposed_cohort_size_too_small") =
      out_exposed_cohort_too_small,

    Named("detected_cases_triggering_tracing") =
      out_detected_cases_traced,
    Named("infected_contacts_traced") =
      out_infected_contacts_traced,
    Named("downstream_uninfected_contacts_quarantined") =
      out_uninfected_contacts_traced,

    Named("post_response_transmissions") =
      out_post_response_transmissions,
    Named("post_response_transmissions_outside_quarantine") =
      out_post_response_transmissions_outside,
    Named("any_post_response_transmission_outside_quarantine") =
      out_any_leakage,
    Named("post_response_infections") =
      out_post_response_infections,
    Named("additional_infections_after_current_generation") =
      out_additional_after_current_generation,
    Named("additional_infections_after_response") =
      out_additional_after_response,

    Named("prevented_by_quarantine_total") =
      out_prevented_by_quarantine,
    Named("prevented_by_isolation_total") =
      out_prevented_by_isolation,
    Named("max_secondary_cases") =
      out_max_secondary,

    Named("response_duration_after_response_days") =
      out_response_duration,
    Named("days_to_first_leakage") =
      out_days_to_first_leakage
  );
}
')

# Transmission regimes
transmission_regimes <- tibble::tribble(
  ~transmission_regime, ~transmission_label, ~R0, ~k,
  "supercritical", "Supercritical", 2.12, 0.3157164,
  "subcritical", "Subcritical", 0.2416027, 0.3998142
)

fixed_inputs <- bind_rows(
  tibble(
    parameter = c(
      "input",
      "delay_file_used",
      "out",
      "n_runs_per_regime",
      "chunk_size",
      "q_values",
      "seed",
      "current_generation",
      "response_day",
      "response_lag_from_first_reported_onset",
      "align_response_to_first_reported_onset",
      "exposed_cohort_size",
      "backward_tracing_window",
      "daily_contact_rate",
      "contact_count_dispersion",
      "contact_multiplier",
      "notification_delay_mean",
      "notification_delay_sd",
      "max_days",
      "max_cases",
      "initial_infected",
      "p_detect",
      "detection_delay_mean",
      "detection_delay_sd",
      "isolation_delay_mean",
      "isolation_delay_sd",
      "isolation_effectiveness",
      "quarantine_duration",
      "quarantine_effectiveness",
      "quarantine_adherence"
    ),
    value = as.character(c(
      args$input,
      delay_file,
      args$out,
      args$n_runs,
      args$chunk_size,
      paste(args$q_values, collapse = ","),
      args$seed,
      args$current_generation,
      args$response_day,
      args$response_lag_from_first_reported_onset,
      args$align_response_to_first_reported_onset,
      args$exposed_cohort_size,
      args$backward_tracing_window,
      args$daily_contact_rate,
      args$contact_count_dispersion,
      args$contact_multiplier,
      args$notification_delay_mean,
      args$notification_delay_sd,
      args$max_days,
      args$max_cases,
      args$initial_infected,
      args$p_detect,
      args$detection_delay_mean,
      args$detection_delay_sd,
      args$isolation_delay_mean,
      args$isolation_delay_sd,
      args$isolation_effectiveness,
      args$quarantine_duration,
      args$quarantine_effectiveness,
      args$quarantine_adherence
    ))
  ),
  transmission_regimes %>%
    transmute(
      parameter = paste0(transmission_regime, "_R0_k"),
      value = paste0("R0=", R0, "; k=", k)
    )
)

readr::write_csv(
  fixed_inputs,
  file.path(args$out, "tables", "fixed_inputs.csv")
)

# Run simulations in checkpointed chunks
message("Running optimized two-pass dynamic tracing analysis...")
t0 <- Sys.time()

runs_list <- vector("list", nrow(transmission_regimes))

for (ii in seq_len(nrow(transmission_regimes))) {
  reg <- transmission_regimes[ii, ]
  regime_seed <- as.integer(args$seed + (ii - 1L) * 1000003L)
  
  message(
    "  ", reg$transmission_label,
    " (R0=", reg$R0,
    ", k=", reg$k, ")"
  )
  
  chunk_starts <- seq.int(
    from = 0L,
    to = args$n_runs - 1L,
    by = args$chunk_size
  )
  
  regime_chunks <- vector("list", length(chunk_starts))
  
  for (cc in seq_along(chunk_starts)) {
    run_offset <- chunk_starts[[cc]]
    chunk_n <- min(
      args$chunk_size,
      args$n_runs - run_offset
    )
    
    chunk_seed <- as.integer(
      regime_seed + run_offset * 1009L + cc * 7919L
    )
    
    chunk_file <- file.path(
      args$out,
      "chunks",
      sprintf(
        "%s_runs_%07d_%07d.csv",
        reg$transmission_regime,
        run_offset + 1L,
        run_offset + chunk_n
      )
    )
    
    message(
      "    chunk ", cc, "/", length(chunk_starts),
      ": runs ", run_offset + 1L, "-",
      run_offset + chunk_n
    )
    
    chunk_result <- run_dynamic_tracing_cpp(
      q_values = args$q_values,
      n_runs = chunk_n,
      run_offset = run_offset,
      generation_draws = generation_draws,
      incubation_draws = incubation_draws,
      current_generation = args$current_generation,
      response_day = args$response_day,
      response_lag_from_first_reported_onset =
        args$response_lag_from_first_reported_onset,
      align_response_to_first_reported_onset =
        args$align_response_to_first_reported_onset,
      exposed_cohort_size = args$exposed_cohort_size,
      backward_tracing_window = args$backward_tracing_window,
      daily_contact_rate = args$daily_contact_rate,
      contact_count_dispersion = args$contact_count_dispersion,
      contact_multiplier = args$contact_multiplier,
      notification_delay_mean = args$notification_delay_mean,
      notification_delay_sd = args$notification_delay_sd,
      max_cases = args$max_cases,
      max_days = args$max_days,
      initial_infected = args$initial_infected,
      R0 = reg$R0,
      k = reg$k,
      p_detect = args$p_detect,
      detection_delay_mean = args$detection_delay_mean,
      detection_delay_sd = args$detection_delay_sd,
      isolation_delay_mean = args$isolation_delay_mean,
      isolation_delay_sd = args$isolation_delay_sd,
      isolation_effectiveness = args$isolation_effectiveness,
      quarantine_duration = args$quarantine_duration,
      quarantine_effectiveness = args$quarantine_effectiveness,
      quarantine_adherence = args$quarantine_adherence,
      seed = chunk_seed
    ) %>%
      as_tibble() %>%
      mutate(
        transmission_regime = reg$transmission_regime,
        transmission_label = reg$transmission_label,
        R0_regime = reg$R0,
        k_regime = reg$k,
        regime_seed = regime_seed,
        chunk_id = cc,
        chunk_seed = chunk_seed
      )
    
    readr::write_csv(chunk_result, chunk_file)
    regime_chunks[[cc]] <- chunk_result
    
    rm(chunk_result)
    invisible(gc())
  }
  
  runs_list[[ii]] <- bind_rows(regime_chunks)
  rm(regime_chunks)
  invisible(gc())
}

runs <- bind_rows(runs_list) %>%
  mutate(
    global_run_id = row_number(),
    intervention_label = recode(
      scenario,
      perfect_intervention = "Perfect intervention",
      imperfect_intervention = "Imperfect intervention"
    ),
    comparison_group = factor(
      paste(transmission_label, intervention_label, sep = " — "),
      levels = c(
        "Supercritical — Perfect intervention",
        "Supercritical — Imperfect intervention",
        "Subcritical — Perfect intervention",
        "Subcritical — Imperfect intervention"
      )
    ),
    q_label = factor(
      make_q_label(q_system),
      levels = make_q_label(
        sort(unique(q_system), decreasing = TRUE)
      )
    ),
    missed_label = factor(
      make_missed_label(q_system),
      levels = make_missed_label(
        sort(unique(q_system), decreasing = TRUE)
      )
    ),
    response_duration_weeks =
      response_duration_after_response_days / 7,
    additional_response_duration_days =
      pmax(
        0,
        response_duration_after_response_days -
          args$quarantine_duration
      ),
    additional_response_duration_months =
      additional_response_duration_days / 30.4375
  )

message(
  "Simulation time: ",
  round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 2),
  " seconds"
)

readr::write_csv(
  runs,
  file.path(args$out, "tables", "run_level_results.csv")
)

runs_conditioned <- runs %>%
  filter(reached_current_generation == 1)

if (nrow(runs_conditioned) == 0) {
  stop(
    "No simulations reached the conditioning event. ",
    "Increase n_runs or reconsider the conditioning assumptions.",
    call. = FALSE
  )
}

# Computational-limit diagnostics

# These summaries do not exclude simulations that reached a limit.
# They quantify whether the limits were binding so that capped or
# time-truncated trajectories can be reported transparently.
#
summarise_limits <- function(dat) {
  dat %>%
    group_by(
      transmission_regime,
      transmission_label,
      R0_regime,
      k_regime,
      scenario,
      intervention_label,
      q_system,
      contacts_missed
    ) %>%
    summarise(
      n_runs = n(),
      n_reached_case_limit =
        sum(reached_case_limit == 1, na.rm = TRUE),
      proportion_reached_case_limit =
        mean(reached_case_limit == 1, na.rm = TRUE),
      n_reached_time_limit =
        sum(reached_time_limit == 1, na.rm = TRUE),
      proportion_reached_time_limit =
        mean(reached_time_limit == 1, na.rm = TRUE),
      n_reached_either_limit =
        sum(
          reached_case_limit == 1 |
            reached_time_limit == 1,
          na.rm = TRUE
        ),
      proportion_reached_either_limit =
        mean(
          reached_case_limit == 1 |
            reached_time_limit == 1,
          na.rm = TRUE
        ),
      total_proposed_infections_after_time_limit =
        sum(
          proposed_infections_after_time_limit,
          na.rm = TRUE
        ),
      median_proposed_infections_after_time_limit =
        median(
          proposed_infections_after_time_limit,
          na.rm = TRUE
        ),
      maximum_proposed_infections_after_time_limit =
        max(
          proposed_infections_after_time_limit,
          na.rm = TRUE
        ),
      .groups = "drop"
    )
}

limit_summary_all <- summarise_limits(runs)
limit_summary_conditioned <- summarise_limits(runs_conditioned)

readr::write_csv(
  limit_summary_all,
  file.path(
    args$out,
    "tables",
    "simulation_limit_summary_all.csv"
  )
)

readr::write_csv(
  limit_summary_conditioned,
  file.path(
    args$out,
    "tables",
    "simulation_limit_summary_conditioned.csv"
  )
)

limit_overall_all <- runs %>%
  summarise(
    analysis_set = "all_generated_runs",
    n_runs = n(),
    n_reached_case_limit =
      sum(reached_case_limit == 1, na.rm = TRUE),
    proportion_reached_case_limit =
      mean(reached_case_limit == 1, na.rm = TRUE),
    n_reached_time_limit =
      sum(reached_time_limit == 1, na.rm = TRUE),
    proportion_reached_time_limit =
      mean(reached_time_limit == 1, na.rm = TRUE),
    n_reached_either_limit =
      sum(
        reached_case_limit == 1 |
          reached_time_limit == 1,
        na.rm = TRUE
      ),
    proportion_reached_either_limit =
      mean(
        reached_case_limit == 1 |
          reached_time_limit == 1,
        na.rm = TRUE
      ),
    total_proposed_infections_after_time_limit =
      sum(
        proposed_infections_after_time_limit,
        na.rm = TRUE
      )
  )

limit_overall_conditioned <- runs_conditioned %>%
  summarise(
    analysis_set = "generation_conditioned_runs",
    n_runs = n(),
    n_reached_case_limit =
      sum(reached_case_limit == 1, na.rm = TRUE),
    proportion_reached_case_limit =
      mean(reached_case_limit == 1, na.rm = TRUE),
    n_reached_time_limit =
      sum(reached_time_limit == 1, na.rm = TRUE),
    proportion_reached_time_limit =
      mean(reached_time_limit == 1, na.rm = TRUE),
    n_reached_either_limit =
      sum(
        reached_case_limit == 1 |
          reached_time_limit == 1,
        na.rm = TRUE
      ),
    proportion_reached_either_limit =
      mean(
        reached_case_limit == 1 |
          reached_time_limit == 1,
        na.rm = TRUE
      ),
    total_proposed_infections_after_time_limit =
      sum(
        proposed_infections_after_time_limit,
        na.rm = TRUE
      )
  )

limit_overall <- bind_rows(
  limit_overall_all,
  limit_overall_conditioned
)

readr::write_csv(
  limit_overall,
  file.path(
    args$out,
    "tables",
    "simulation_limit_summary_overall.csv"
  )
)

if (any(runs_conditioned$reached_case_limit == 1, na.rm = TRUE)) {
  warning(
    "At least one conditioned realized-outbreak simulation reached max_cases = ",
    args$max_cases,
    ". Outcomes for those runs are truncated.",
    call. = FALSE
  )
}

if (any(runs_conditioned$reached_time_limit == 1, na.rm = TRUE)) {
  warning(
    "At least one conditioned realized-outbreak simulation reached max_days = ",
    args$max_days,
    ". Outcomes for those runs are right-truncated in time.",
    call. = FALSE
  )
}

conditioning_summary <- runs %>%
  group_by(
    transmission_regime,
    transmission_label,
    R0_regime,
    k_regime,
    scenario,
    intervention_label,
    q_system,
    contacts_missed
  ) %>%
  summarise(
    n_generated = n(),
    n_reached = sum(reached_current_generation == 1),
    conditioning_acceptance_probability =
      mean(reached_current_generation == 1),
    .groups = "drop"
  )

readr::write_csv(
  conditioning_summary,
  file.path(args$out, "tables", "conditioning_acceptance.csv")
)

scenario_summary <- runs_conditioned %>%
  group_by(
    transmission_regime, transmission_label, R0_regime, k_regime,
    scenario, intervention_label, comparison_group,
    q_system, contacts_missed, q_label, missed_label
  ) %>%
  summarise(
    n_runs = n(),
    leakage_events = sum(
      any_post_response_transmission_outside_quarantine,
      na.rm = TRUE
    ),
    leakage_prob = leakage_events / n_runs,
    leakage_lower = binom_ci(leakage_events, n_runs)["lower"],
    leakage_upper = binom_ci(leakage_events, n_runs)["upper"],
    
    median_additional_infections =
      median(additional_infections_after_response, na.rm = TRUE),
    q25_additional_infections =
      quant(additional_infections_after_response, 0.25),
    q75_additional_infections =
      quant(additional_infections_after_response, 0.75),
    q95_additional_infections =
      quant(additional_infections_after_response, 0.95),
    mean_additional_infections =
      mean(additional_infections_after_response, na.rm = TRUE),
    
    median_downstream_uninfected_quarantined =
      median(downstream_uninfected_contacts_quarantined, na.rm = TRUE),
    q25_downstream_uninfected_quarantined =
      quant(downstream_uninfected_contacts_quarantined, 0.25),
    q75_downstream_uninfected_quarantined =
      quant(downstream_uninfected_contacts_quarantined, 0.75),
    q95_downstream_uninfected_quarantined =
      quant(downstream_uninfected_contacts_quarantined, 0.95),
    
    median_additional_response_months =
      median(additional_response_duration_months, na.rm = TRUE),
    q25_additional_response_months =
      quant(additional_response_duration_months, 0.25),
    q75_additional_response_months =
      quant(additional_response_duration_months, 0.75),
    q95_additional_response_months =
      quant(additional_response_duration_months, 0.95),
    .groups = "drop"
  )

readr::write_csv(
  scenario_summary,
  file.path(args$out, "tables", "scenario_summary.csv")
)

if (any(
  runs_conditioned$exposed_cohort_size_too_small == 1,
  na.rm = TRUE
)) {
  warning(
    "At least one conditioned simulation had more infected people ",
    "at response than exposed_cohort_size.",
    call. = FALSE
  )
}

# Main Figure 2 and supplementary Figure S2 plotting data
main_figure_group_levels <- c(
  "Supercritical — Imperfect intervention",
  "Subcritical — Imperfect intervention"
)

si_figure_group_levels <- c(
  "Supercritical — Perfect intervention",
  "Subcritical — Perfect intervention"
)

main_figure_runs <- runs_conditioned %>%
  filter(scenario == "imperfect_intervention") %>%
  mutate(
    comparison_group = factor(
      comparison_group,
      levels = main_figure_group_levels
    )
  ) %>%
  filter(!is.na(comparison_group))

si_figure_runs <- runs_conditioned %>%
  filter(scenario == "perfect_intervention") %>%
  mutate(
    comparison_group = factor(
      comparison_group,
      levels = si_figure_group_levels
    )
  ) %>%
  filter(!is.na(comparison_group))

make_figure_summaries <- function(dat) {
  probability <- dat %>%
    group_by(
      transmission_regime,
      transmission_label,
      R0_regime,
      k_regime,
      scenario,
      intervention_label,
      comparison_group,
      q_system,
      contacts_missed,
      q_label,
      missed_label
    ) %>%
    summarise(
      n_runs = n(),
      continued_transmission_events = sum(
        any_post_response_transmission_outside_quarantine == 1,
        na.rm = TRUE
      ),
      leakage_prob = continued_transmission_events / n_runs,
      leakage_lower = binom_ci(
        continued_transmission_events,
        n_runs
      )[["lower"]],
      leakage_upper = binom_ci(
        continued_transmission_events,
        n_runs
      )[["upper"]],
      .groups = "drop"
    )
  
  severity <- dat %>%
    filter(
      is.finite(additional_infections_after_response),
      additional_infections_after_response > 0
    ) %>%
    group_by(
      transmission_regime,
      transmission_label,
      R0_regime,
      k_regime,
      scenario,
      intervention_label,
      comparison_group,
      q_system,
      contacts_missed,
      q_label,
      missed_label
    ) %>%
    summarise(
      n_runs_with_additional_infection = n(),
      median_additional_infections =
        median(additional_infections_after_response, na.rm = TRUE),
      q25_additional_infections =
        quant(additional_infections_after_response, 0.25),
      q75_additional_infections =
        quant(additional_infections_after_response, 0.75),
      q90_additional_infections =
        quant(additional_infections_after_response, 0.90),
      q95_additional_infections =
        quant(additional_infections_after_response, 0.95),
      q99_additional_infections =
        quant(additional_infections_after_response, 0.99),
      mean_additional_infections =
        mean(additional_infections_after_response, na.rm = TRUE),
      maximum_additional_infections =
        max(additional_infections_after_response, na.rm = TRUE),
      .groups = "drop"
    )
  
  operational <- dat %>%
    group_by(
      transmission_regime,
      transmission_label,
      R0_regime,
      k_regime,
      scenario,
      intervention_label,
      comparison_group,
      q_system,
      contacts_missed,
      q_label,
      missed_label
    ) %>%
    summarise(
      n_runs = n(),
      
      median_downstream_uninfected_quarantined =
        median(
          downstream_uninfected_contacts_quarantined,
          na.rm = TRUE
        ),
      q25_downstream_uninfected_quarantined =
        quant(
          downstream_uninfected_contacts_quarantined,
          0.25
        ),
      q75_downstream_uninfected_quarantined =
        quant(
          downstream_uninfected_contacts_quarantined,
          0.75
        ),
      q90_downstream_uninfected_quarantined =
        quant(
          downstream_uninfected_contacts_quarantined,
          0.90
        ),
      q95_downstream_uninfected_quarantined =
        quant(
          downstream_uninfected_contacts_quarantined,
          0.95
        ),
      q99_downstream_uninfected_quarantined =
        quant(
          downstream_uninfected_contacts_quarantined,
          0.99
        ),
      mean_downstream_uninfected_quarantined =
        mean(
          downstream_uninfected_contacts_quarantined,
          na.rm = TRUE
        ),
      probability_quarantine_ge_100 =
        mean(
          downstream_uninfected_contacts_quarantined >= 100,
          na.rm = TRUE
        ),
      probability_quarantine_ge_1000 =
        mean(
          downstream_uninfected_contacts_quarantined >= 1000,
          na.rm = TRUE
        ),
      probability_quarantine_ge_10000 =
        mean(
          downstream_uninfected_contacts_quarantined >= 10000,
          na.rm = TRUE
        ),
      probability_quarantine_ge_100000 =
        mean(
          downstream_uninfected_contacts_quarantined >= 100000,
          na.rm = TRUE
        ),
      
      median_additional_response_weeks =
        median(additional_response_duration_days / 7, na.rm = TRUE),
      q25_additional_response_weeks =
        quant(additional_response_duration_days / 7, 0.25),
      q75_additional_response_weeks =
        quant(additional_response_duration_days / 7, 0.75),
      q90_additional_response_weeks =
        quant(additional_response_duration_days / 7, 0.90),
      q95_additional_response_weeks =
        quant(additional_response_duration_days / 7, 0.95),
      q99_additional_response_weeks =
        quant(additional_response_duration_days / 7, 0.99),
      mean_additional_response_weeks =
        mean(additional_response_duration_days / 7, na.rm = TRUE),
      
      # Retain month-scale quantities for backward-compatible tables.
      median_additional_response_months =
        median(additional_response_duration_months, na.rm = TRUE),
      q25_additional_response_months =
        quant(additional_response_duration_months, 0.25),
      q75_additional_response_months =
        quant(additional_response_duration_months, 0.75),
      q90_additional_response_months =
        quant(additional_response_duration_months, 0.90),
      q95_additional_response_months =
        quant(additional_response_duration_months, 0.95),
      q99_additional_response_months =
        quant(additional_response_duration_months, 0.99),
      mean_additional_response_months =
        mean(additional_response_duration_months, na.rm = TRUE),
      probability_additional_response_ge_1_month =
        mean(additional_response_duration_months >= 1, na.rm = TRUE),
      probability_additional_response_ge_2_months =
        mean(additional_response_duration_months >= 2, na.rm = TRUE),
      probability_additional_response_ge_3_months =
        mean(additional_response_duration_months >= 3, na.rm = TRUE),
      probability_additional_response_ge_4_months =
        mean(additional_response_duration_months >= 4, na.rm = TRUE),
      probability_additional_response_ge_6_months =
        mean(additional_response_duration_months >= 6, na.rm = TRUE),
      
      .groups = "drop"
    )
  
  
  duration_positive <- dat %>%
    filter(
      is.finite(additional_response_duration_days),
      additional_response_duration_days > 0
    ) %>%
    group_by(
      transmission_regime,
      transmission_label,
      R0_regime,
      k_regime,
      scenario,
      intervention_label,
      comparison_group,
      q_system,
      contacts_missed,
      q_label,
      missed_label
    ) %>%
    summarise(
      n_runs_with_extended_response = n(),
      median_additional_response_weeks =
        median(additional_response_duration_days / 7, na.rm = TRUE),
      q25_additional_response_weeks =
        quant(additional_response_duration_days / 7, 0.25),
      q75_additional_response_weeks =
        quant(additional_response_duration_days / 7, 0.75),
      q95_additional_response_weeks =
        quant(additional_response_duration_days / 7, 0.95),
      .groups = "drop"
    )
  
  duration_probability <- dat %>%
    group_by(
      transmission_regime,
      transmission_label,
      R0_regime,
      k_regime,
      scenario,
      intervention_label,
      comparison_group,
      q_system,
      contacts_missed,
      q_label,
      missed_label
    ) %>%
    summarise(
      n_runs = n(),
      n_extended_response = sum(
        additional_response_duration_days > 0,
        na.rm = TRUE
      ),
      probability_extended_response =
        n_extended_response / n_runs,
      extended_response_lower = binom_ci(
        n_extended_response,
        n_runs
      )[["lower"]],
      extended_response_upper = binom_ci(
        n_extended_response,
        n_runs
      )[["upper"]],
      .groups = "drop"
    )
  
  merged <- probability %>%
    left_join(
      severity,
      by = c(
        "transmission_regime",
        "transmission_label",
        "R0_regime",
        "k_regime",
        "scenario",
        "intervention_label",
        "comparison_group",
        "q_system",
        "contacts_missed",
        "q_label",
        "missed_label"
      ),
      suffix = c("", "_panelB")
    ) %>%
    left_join(
      operational,
      by = c(
        "transmission_regime",
        "transmission_label",
        "R0_regime",
        "k_regime",
        "scenario",
        "intervention_label",
        "comparison_group",
        "q_system",
        "contacts_missed",
        "q_label",
        "missed_label"
      ),
      suffix = c("", "_panelsCD")
    ) %>%
    arrange(comparison_group, contacts_missed)
  
  list(
    probability = probability,
    severity = severity,
    operational = operational,
    duration_positive = duration_positive,
    duration_probability = duration_probability,
    merged = merged
  )
}

main_fig2 <- make_figure_summaries(main_figure_runs)
si_fig2 <- make_figure_summaries(si_figure_runs)

# Preserve the original object names for downstream main-figure code.
figure2_runs <- main_figure_runs
figure2_group_levels <- main_figure_group_levels
fig2_probability_data <- main_fig2$probability
severity_summary <- main_fig2$severity
operational_summary <- main_fig2$operational
figure2_manuscript_summary <- main_fig2$merged

readr::write_csv(
  main_fig2$probability,
  file.path(
    args$out,
    "tables",
    "figure2_main_panelA_probability_imperfect.csv"
  )
)

if (any(!is.finite(main_fig2$probability$leakage_lower)) ||
    any(!is.finite(main_fig2$probability$leakage_upper))) {
  stop(
    "Panel A confidence intervals contain non-finite values.",
    call. = FALSE
  )
}

readr::write_csv(
  main_fig2$probability %>%
    transmute(
      transmission_regime,
      scenario,
      q_system,
      n_runs,
      continued_transmission_events,
      leakage_prob,
      leakage_lower,
      leakage_upper,
      ci_height = leakage_upper - leakage_lower
    ),
  file.path(
    args$out,
    "tables",
    "figure2_main_panelA_CI_diagnostic.csv"
  )
)

readr::write_csv(
  main_fig2$severity,
  file.path(
    args$out,
    "tables",
    "figure2_main_panelB_conditional_infections_imperfect.csv"
  )
)

readr::write_csv(
  main_fig2$operational,
  file.path(
    args$out,
    "tables",
    "figure2_main_panelsCD_unconditional_imperfect.csv"
  )
)

readr::write_csv(
  main_fig2$duration_positive,
  file.path(
    args$out,
    "tables",
    "figure2_main_panelD_positive_duration_weeks_imperfect.csv"
  )
)

readr::write_csv(
  main_fig2$duration_probability,
  file.path(
    args$out,
    "tables",
    "figure2_main_panelD_probability_extended_response_imperfect.csv"
  )
)

readr::write_csv(
  main_fig2$merged,
  file.path(
    args$out,
    "tables",
    "figure2_main_all_manuscript_quantities_imperfect.csv"
  )
)

readr::write_csv(
  si_fig2$probability,
  file.path(
    args$out,
    "tables",
    "figureS2_panelA_probability_perfect.csv"
  )
)

readr::write_csv(
  si_fig2$severity,
  file.path(
    args$out,
    "tables",
    "figureS2_panelB_conditional_infections_perfect.csv"
  )
)

readr::write_csv(
  si_fig2$operational,
  file.path(
    args$out,
    "tables",
    "figureS2_panelsCD_unconditional_perfect.csv"
  )
)

readr::write_csv(
  si_fig2$duration_positive,
  file.path(
    args$out,
    "tables",
    "figureS2_panelD_positive_duration_weeks_perfect.csv"
  )
)

readr::write_csv(
  si_fig2$duration_probability,
  file.path(
    args$out,
    "tables",
    "figureS2_panelD_probability_extended_response_perfect.csv"
  )
)

readr::write_csv(
  si_fig2$merged,
  file.path(
    args$out,
    "tables",
    "figureS2_all_manuscript_quantities_perfect.csv"
  )
)

# Dose-response diagnostics
dose_response_diagnostics <- runs_conditioned %>%
  group_by(
    transmission_regime, transmission_label,
    scenario, intervention_label, comparison_group,
    q_system, contacts_missed, missed_label
  ) %>%
  summarise(
    n_runs = n(),
    probability_continued_transmission =
      mean(any_post_response_transmission_outside_quarantine, na.rm = TRUE),
    mean_additional_infections =
      mean(additional_infections_after_response, na.rm = TRUE),
    median_additional_infections =
      median(additional_infections_after_response, na.rm = TRUE),
    q75_additional_infections =
      quant(additional_infections_after_response, 0.75),
    q90_additional_infections =
      quant(additional_infections_after_response, 0.90),
    q95_additional_infections =
      quant(additional_infections_after_response, 0.95),
    q99_additional_infections =
      quant(additional_infections_after_response, 0.99),
    median_uninfected_quarantined =
      median(downstream_uninfected_contacts_quarantined, na.rm = TRUE),
    median_additional_response_months =
      median(additional_response_duration_months, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(transmission_regime, scenario, contacts_missed)

readr::write_csv(
  dose_response_diagnostics,
  file.path(args$out, "tables", "dose_response_diagnostics.csv")
)

# Main Figure 2 and supplementary Figure S2
if (isTRUE(args$make_plots)) {
  
  required_panelA <- c(
    "leakage_prob", "leakage_lower", "leakage_upper"
  )
  if (!all(required_panelA %in% names(main_fig2$probability))) {
    stop("Panel A CI columns are missing.", call. = FALSE)
  }
  
  if (any(main_fig2$probability$leakage_lower >
          main_fig2$probability$leakage_prob, na.rm = TRUE) ||
      any(main_fig2$probability$leakage_upper <
          main_fig2$probability$leakage_prob, na.rm = TRUE)) {
    stop(
      "Panel A confidence intervals do not contain the point estimate.",
      call. = FALSE
    )
  }
  
  
  theme_pub <- theme_classic(base_size = 17) +
    theme(
      text = element_text(size = 17),
      axis.text = element_text(colour = "black", size = 15),
      axis.title = element_text(size = 18),
      axis.title.x = element_text(margin = margin(t = 10)),
      axis.title.y = element_text(margin = margin(r = 10)),
      legend.text = element_text(size = 15),
      legend.key.size = grid::unit(0.8, "cm"),
      legend.spacing.x = grid::unit(0.35, "cm"),
      panel.grid.major.y = element_line(
        colour = "grey86",
        linewidth = 0.45
      ),
      panel.grid.minor = element_blank(),
      plot.tag = element_text(face = "bold", size = 20),
      plot.margin = margin(14, 16, 14, 16)
    )
  
  dodge_width <- 0.62
  point_size <- 4.2
  
  make_four_panel_figure <- function(
    probability_data,
    severity_data,
    operational_data,
    duration_data,
    group_values,
    group_labels
  ) {
    dodge <- position_dodge(width = dodge_width)
    
    pA <- ggplot(
      probability_data,
      aes(
        x = missed_label,
        y = leakage_prob,
        fill = comparison_group,
        group = comparison_group
      )
    ) +
      # Points first, then the exact binomial 95% CI on top.
      geom_point(
        position = dodge,
        shape = 21,
        size = point_size,
        stroke = 0.55,
        colour = "black"
      ) +
      geom_linerange(
        aes(
          ymin = leakage_lower,
          ymax = leakage_upper
        ),
        position = dodge,
        linewidth = 1.05,
        colour = "black",
        na.rm = FALSE
      ) +
      geom_errorbar(
        aes(
          ymin = leakage_lower,
          ymax = leakage_upper
        ),
        position = dodge,
        width = 0.18,
        linewidth = 1.05,
        colour = "black",
        na.rm = FALSE
      ) +
      scale_fill_manual(
        values = group_values,
        labels = group_labels,
        drop = FALSE,
        guide = guide_legend(
          title = NULL,
          nrow = 1,
          byrow = TRUE,
          override.aes = list(
            shape = 21,
            size = 4.2,
            colour = "black"
          )
        )
      ) +
      scale_x_discrete(drop = FALSE) +
      scale_y_continuous(
        breaks = seq(0, 1, by = 0.25),
        expand = expansion(mult = c(0.015, 0.025))
      ) +
      coord_cartesian(
        ylim = c(-0.015, 1.015),
        clip = "off"
      ) +
      labs(
        x = "Missed contacts",
        y = "Probability of continued transmission",
        tag = "A"
      ) +
      theme_pub +
      theme(legend.position = "bottom")
    
    pB <- ggplot(
      severity_data,
      aes(
        x = missed_label,
        y = median_additional_infections,
        fill = comparison_group,
        group = comparison_group
      )
    ) +
      geom_errorbar(
        aes(
          ymin = q25_additional_infections,
          ymax = q75_additional_infections
        ),
        position = dodge,
        width = 0.10,
        linewidth = 0.45
      ) +
      geom_point(
        position = dodge,
        shape = 21,
        size = point_size,
        stroke = 0.40,
        colour = "black",
        show.legend = FALSE
      ) +
      scale_fill_manual(
        values = group_values,
        labels = group_labels,
        drop = FALSE,
        guide = "none"
      ) +
      scale_x_discrete(drop = FALSE) +
      scale_y_continuous(
        trans = scales::log1p_trans(),
        breaks = c(0, 1, 10, 100, 1000, 5000),
        labels = scales::label_number(big.mark = ",")
      ) +
      labs(
        x = "Missed contacts",
        y = "Additional infections",
        tag = "B"
      ) +
      theme_pub
    
    pC <- ggplot(
      operational_data,
      aes(
        x = missed_label,
        y = median_downstream_uninfected_quarantined,
        fill = comparison_group,
        group = comparison_group
      )
    ) +
      geom_errorbar(
        aes(
          ymin = q25_downstream_uninfected_quarantined,
          ymax = q75_downstream_uninfected_quarantined
        ),
        position = dodge,
        width = 0.10,
        linewidth = 0.45
      ) +
      geom_point(
        position = dodge,
        shape = 21,
        size = point_size,
        stroke = 0.40,
        colour = "black",
        show.legend = FALSE
      ) +
      scale_fill_manual(
        values = group_values,
        labels = group_labels,
        drop = FALSE,
        guide = "none"
      ) +
      scale_x_discrete(drop = FALSE) +
      scale_y_continuous(
        trans = scales::log1p_trans(),
        breaks = c(0, 10, 100, 1000, 10000, 100000),
        labels = scales::label_number(
          big.mark = ",",
          accuracy = 1
        )
      ) +
      labs(
        x = "Missed contacts",
        y = "Uninfected individuals quarantined",
        tag = "C"
      ) +
      theme_pub
    
    pD <- ggplot(
      duration_data,
      aes(
        x = missed_label,
        y = median_additional_response_weeks,
        fill = comparison_group,
        group = comparison_group
      )
    ) +
      geom_errorbar(
        aes(
          ymin = q25_additional_response_weeks,
          ymax = q75_additional_response_weeks
        ),
        position = dodge,
        width = 0.10,
        linewidth = 0.60,
        colour = "black",
        na.rm = TRUE
      ) +
      geom_point(
        position = dodge,
        shape = 21,
        size = point_size,
        stroke = 0.45,
        colour = "black",
        show.legend = FALSE,
        na.rm = TRUE
      ) +
      scale_fill_manual(
        values = group_values,
        labels = group_labels,
        drop = FALSE,
        guide = "none"
      ) +
      scale_x_discrete(drop = FALSE) +
      scale_y_continuous(
        trans = scales::log1p_trans(),
        breaks = c(0, 1, 2, 4, 8, 16, 32, 64, 128),
        labels = scales::label_number(big.mark = ","),
        expand = expansion(mult = c(0.02, 0.06))
      ) +
      labs(
        x = "Missed contacts",
        y = "Additional response duration (weeks)",
        tag = "D"
      ) +
      theme_pub
    
    ((pA | pB) / (pC | pD)) +
      patchwork::plot_layout(guides = "collect") &
      theme(
        legend.position = "bottom",
        legend.box = "horizontal"
      )
  }
  
  main_group_values <- c(
    "Supercritical — Imperfect intervention" = "grey55",
    "Subcritical — Imperfect intervention" = "grey15"
  )
  
  main_group_labels <- c(
    "Supercritical — Imperfect intervention" = "High transmission",
    "Subcritical — Imperfect intervention" = "Low transmission"
  )
  
  perfect_group_values <- c(
    "Supercritical — Perfect intervention" = "white",
    "Subcritical — Perfect intervention" = "grey75"
  )
  
  perfect_group_labels <- c(
    "Supercritical — Perfect intervention" = "High transmission",
    "Subcritical — Perfect intervention" = "Low transmission"
  )
  
  fig2 <- make_four_panel_figure(
    probability_data = main_fig2$probability,
    severity_data = main_fig2$severity,
    operational_data = main_fig2$operational,
    duration_data = main_fig2$duration_positive,
    group_values = main_group_values,
    group_labels = main_group_labels
  )
  
  figS2 <- make_four_panel_figure(
    probability_data = si_fig2$probability,
    severity_data = si_fig2$severity,
    operational_data = si_fig2$operational,
    duration_data = si_fig2$duration_positive,
    group_values = perfect_group_values,
    group_labels = perfect_group_labels
  )
  
  ggsave(
    file.path(
      args$out,
      "figures",
      "Figure2_main_imperfect_scenarios.png"
    ),
    fig2,
    width = 14.5,
    height = 11.5,
    dpi = 800,
    bg = "white"
  )
  
  ggsave(
    file.path(
      args$out,
      "figures",
      "Figure2_main_imperfect_scenarios.pdf"
    ),
    fig2,
    width = 14.5,
    height = 11.5,
    bg = "white"
  )
  
  ggsave(
    file.path(
      args$out,
      "figures",
      "FigureS2_perfect_scenarios.png"
    ),
    figS2,
    width = 14.5,
    height = 11.5,
    dpi = 800,
    bg = "white"
  )
  
  ggsave(
    file.path(
      args$out,
      "figures",
      "FigureS2_perfect_scenarios.pdf"
    ),
    figS2,
    width = 14.5,
    height = 11.5,
    bg = "white"
  )
}

# Delay-distribution parameter table
gamma_parameters_from_mean_sd <- function(mean, sd) {
  tibble::tibble(
    mean = mean,
    sd = sd,
    shape = (mean / sd)^2,
    rate = mean / sd^2,
    scale = sd^2 / mean
  )
}

delay_distribution_parameters <- dplyr::bind_rows(
  gamma_parameters_from_mean_sd(
    args$detection_delay_mean,
    args$detection_delay_sd
  ) %>%
    mutate(delay = "Symptom onset to detection"),
  gamma_parameters_from_mean_sd(
    args$isolation_delay_mean,
    args$isolation_delay_sd
  ) %>%
    mutate(delay = "Detection to isolation"),
  gamma_parameters_from_mean_sd(
    args$notification_delay_mean,
    args$notification_delay_sd
  ) %>%
    mutate(delay = "Contact notification")
) %>%
  select(delay, mean, sd, shape, rate, scale)

readr::write_csv(
  delay_distribution_parameters,
  file.path(
    args$out,
    "tables",
    "gamma_delay_distribution_parameters.csv"
  )
)

# Complete paper-output report
conditioning_for_paper <- runs %>%
  group_by(
    transmission_regime,
    transmission_label,
    scenario,
    intervention_label,
    q_system,
    contacts_missed,
    q_label,
    missed_label
  ) %>%
  summarise(
    n_generated = n(),
    n_reached_generation_2 =
      sum(reached_current_generation == 1, na.rm = TRUE),
    conditioning_acceptance_probability =
      mean(reached_current_generation == 1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(
    (transmission_regime == "supercritical" &
       scenario %in% c(
         "perfect_intervention",
         "imperfect_intervention"
       )) |
      (transmission_regime == "subcritical" &
         scenario == "imperfect_intervention")
  ) %>%
  mutate(
    comparison_group = factor(
      paste(transmission_label, intervention_label, sep = " — "),
      levels = figure2_group_levels
    )
  ) %>%
  arrange(comparison_group, contacts_missed)

readr::write_csv(
  conditioning_for_paper,
  file.path(
    args$out,
    "tables",
    "conditioning_acceptance_three_groups.csv"
  )
)

fmt_pct <- function(x, digits = 1) {
  ifelse(
    is.na(x),
    "NA",
    paste0(
      formatC(100 * x, format = "f", digits = digits),
      "%"
    )
  )
}

fmt_num <- function(x, digits = 1) {
  ifelse(
    is.na(x),
    "NA",
    formatC(
      x,
      format = "f",
      digits = digits,
      big.mark = ","
    )
  )
}

fmt_int <- function(x) {
  ifelse(
    is.na(x),
    "NA",
    formatC(
      round(x),
      format = "f",
      digits = 0,
      big.mark = ","
    )
  )
}

paper_report_file <- file.path(
  args$out,
  "tables",
  "MANUSCRIPT_READY_FIGURE2_NUMBERS.txt"
)

con_paper <- file(paper_report_file, open = "wt")

writeLines(
  c(
    "MANUSCRIPT-READY FIGURE 2 RESULTS",
    "=================================",
    "",
    "Plotted groups:",
    "1. Supercritical — Perfect intervention",
    "2. Supercritical — Imperfect intervention",
    "3. Subcritical — Imperfect intervention",
    "",
    "Outcome definitions:",
    "Panel A: all generation-2-conditioned simulations.",
    "Panel B: conditional on at least one additional post-response infection.",
    "Panel C: unconditional across all generation-2-conditioned simulations.",
    "Panel D: unconditional additional duration beyond the initial 42-day quarantine.",
    ""
  ),
  con_paper
)

for (grp in figure2_group_levels) {
  writeLines(grp, con_paper)
  writeLines(
    paste(rep("-", nchar(grp)), collapse = ""),
    con_paper
  )
  
  cond_rows <- conditioning_for_paper %>%
    filter(as.character(comparison_group) == grp) %>%
    arrange(contacts_missed)
  
  if (nrow(cond_rows) > 0) {
    cond_min <- min(
      cond_rows$conditioning_acceptance_probability,
      na.rm = TRUE
    )
    cond_max <- max(
      cond_rows$conditioning_acceptance_probability,
      na.rm = TRUE
    )
    
    writeLines(
      paste0(
        "Generation-2 conditioning acceptance across tracing levels: ",
        fmt_pct(cond_min, 2),
        "–",
        fmt_pct(cond_max, 2),
        "."
      ),
      con_paper
    )
  }
  
  rows <- figure2_manuscript_summary %>%
    filter(as.character(comparison_group) == grp) %>%
    arrange(contacts_missed)
  
  for (ii in seq_len(nrow(rows))) {
    z <- rows[ii, ]
    
    writeLines(
      paste0(
        "\n",
        round(100 * z$contacts_missed),
        "% missed contacts"
      ),
      con_paper
    )
    
    writeLines(
      paste0(
        "  Panel A: probability of continued transmission = ",
        fmt_pct(z$leakage_prob, 1),
        " (95% CI ",
        fmt_pct(z$leakage_lower, 1),
        "–",
        fmt_pct(z$leakage_upper, 1),
        "; ",
        z$continued_transmission_events,
        "/",
        z$n_runs,
        ")."
      ),
      con_paper
    )
    
    if (
      is.na(z$n_runs_with_additional_infection) ||
      z$n_runs_with_additional_infection == 0
    ) {
      writeLines(
        "  Panel B: no simulations had an additional post-response infection.",
        con_paper
      )
    } else {
      writeLines(
        paste0(
          "  Panel B: among simulations with >=1 additional infection (n=",
          z$n_runs_with_additional_infection,
          "), median = ",
          fmt_int(z$median_additional_infections),
          " (IQR ",
          fmt_int(z$q25_additional_infections),
          "–",
          fmt_int(z$q75_additional_infections),
          "; 90th percentile ",
          fmt_int(z$q90_additional_infections),
          "; 95th percentile ",
          fmt_int(z$q95_additional_infections),
          "; 99th percentile ",
          fmt_int(z$q99_additional_infections),
          ")."
        ),
        con_paper
      )
    }
    
    writeLines(
      paste0(
        "  Panel C: unconditional median uninfected contacts quarantined = ",
        fmt_int(z$median_downstream_uninfected_quarantined),
        " (IQR ",
        fmt_int(z$q25_downstream_uninfected_quarantined),
        "–",
        fmt_int(z$q75_downstream_uninfected_quarantined),
        "; 95th percentile ",
        fmt_int(z$q95_downstream_uninfected_quarantined),
        "; P[>=1,000] ",
        fmt_pct(z$probability_quarantine_ge_1000, 1),
        "; P[>=10,000] ",
        fmt_pct(z$probability_quarantine_ge_10000, 1),
        "; P[>=100,000] ",
        fmt_pct(z$probability_quarantine_ge_100000, 1),
        ")."
      ),
      con_paper
    )
    
    writeLines(
      paste0(
        "  Panel D: unconditional median additional response duration = ",
        fmt_num(z$median_additional_response_months, 2),
        " months (IQR ",
        fmt_num(z$q25_additional_response_months, 2),
        "–",
        fmt_num(z$q75_additional_response_months, 2),
        "; 95th percentile ",
        fmt_num(z$q95_additional_response_months, 2),
        " months; P[>=1 month] ",
        fmt_pct(z$probability_additional_response_ge_1_month, 1),
        "; P[>=4 months] ",
        fmt_pct(z$probability_additional_response_ge_4_months, 1),
        "; P[>=6 months] ",
        fmt_pct(z$probability_additional_response_ge_6_months, 1),
        ")."
      ),
      con_paper
    )
  }
  
  writeLines("\n", con_paper)
}

writeLines(
  c(
    "FILES FOR THE PAPER",
    "-------------------",
    "figure2_panelA_probability_three_groups.csv",
    "figure2_panelB_conditional_infections_three_groups.csv",
    "figure2_panelsCD_unconditional_three_groups.csv",
    "figure2_all_manuscript_quantities_three_groups.csv",
    "conditioning_acceptance_three_groups.csv",
    "gamma_delay_distribution_parameters.csv",
    "MANUSCRIPT_READY_FIGURE2_NUMBERS.txt",
    "Figure2_three_groups_final.png",
    "Figure2_three_groups_final.pdf"
  ),
  con_paper
)

close(con_paper)

# Text-ready summary
summary_file <- file.path(
  args$out,
  "tables",
  "text_ready_summary.txt"
)

con <- file(summary_file, open = "wt")

writeLines(
  "DYNAMIC DOWNSTREAM TRACING ANALYSIS",
  con
)
writeLines(
  "===================================",
  con
)
writeLines(
  paste0(
    "Conditioned on the uncontrolled potential tree reaching generation ",
    args$current_generation,
    " by response day ",
    args$response_day,
    "."
  ),
  con
)
writeLines(
  paste0(
    "The response day was aligned to occur ",
    args$response_lag_from_first_reported_onset,
    " days after symptom onset of the first detected case: ",
    args$align_response_to_first_reported_onset,
    "."
  ),
  con
)
writeLines(
  paste0(
    "Tracing completeness q was applied to the initial known exposed ",
    "infected cohort and to downstream infected and uninfected contacts."
  ),
  con
)
writeLines("", con)

for (grp in levels(runs_conditioned$comparison_group)) {
  writeLines(grp, con)
  writeLines(paste(rep("-", nchar(grp)), collapse = ""), con)
  
  rows <- scenario_summary %>%
    filter(comparison_group == grp) %>%
    arrange(desc(q_system))
  
  for (jj in seq_len(nrow(rows))) {
    z <- rows[jj, ]
    
    writeLines(
      paste0(
        as.character(z$missed_label),
        " missed: P(leakage)=",
        sprintf("%.3f", z$leakage_prob),
        " (95% CI ",
        sprintf("%.3f", z$leakage_lower),
        "–",
        sprintf("%.3f", z$leakage_upper),
        "); median additional infections=",
        round(z$median_additional_infections),
        "; median downstream uninfected contacts quarantined=",
        round(z$median_downstream_uninfected_quarantined),
        " (IQR ",
        round(z$q25_downstream_uninfected_quarantined),
        "–",
        round(z$q75_downstream_uninfected_quarantined),
        "); median additional response duration=",
        round(z$median_additional_response_months, 1),
        " months."
      ),
      con
    )
  }
  
  writeLines("", con)
}

close(con)

message(
  "Done. Outputs written to: ",
  normalizePath(args$out, winslash = "/", mustWork = FALSE)
)

library(readr)
library(dplyr)

out_dir <- "~/Hantavirus-paper/BP results/Fig 2"

overall <- read_csv(
  file.path(out_dir, "/tables/simulation_limit_summary_overall.csv")
)

conditioned <- read_csv(
  file.path(out_dir, "/tables/simulation_limit_summary_conditioned.csv")
)

overall
conditioned %>%
  select(
    transmission_regime,
    scenario,
    q_system,
    n_runs,
    n_reached_case_limit,
    proportion_reached_case_limit,
    n_reached_time_limit,
    proportion_reached_time_limit,
    n_reached_either_limit,
    proportion_reached_either_limit
  )

library(dplyr)

runs_conditioned %>%
  filter(
    scenario == "imperfect_intervention",
    abs(q_system - 0.90) < 1e-8   # 10% of contacts missed
  ) %>%
  group_by(transmission_regime, transmission_label) %>%
  summarise(
    n_simulations = n(),
    n_quarantined_at_least_5000 = sum(
      downstream_uninfected_contacts_quarantined >= 5000,
      na.rm = TRUE
    ),
    probability_quarantined_at_least_5000 =
      n_quarantined_at_least_5000 / n_simulations,
    ci_lower = binom_ci(
      n_quarantined_at_least_5000,
      n_simulations
    )[["lower"]],
    ci_upper = binom_ci(
      n_quarantined_at_least_5000,
      n_simulations
    )[["upper"]],
    .groups = "drop"
  )

runs_conditioned %>%
  filter(
    scenario == "imperfect_intervention",
    abs(q_system - 0.90) < 1e-8   # 10% of contacts missed
  ) %>%
  group_by(transmission_regime, transmission_label) %>%
  summarise(
    n_simulations = n(),
    n_continued_transmission = sum(
      any_post_response_transmission_outside_quarantine == 1,
      na.rm = TRUE
    ),
    probability_continued_transmission =
      n_continued_transmission / n_simulations,
    ci_lower = binom_ci(
      n_continued_transmission,
      n_simulations
    )[["lower"]],
    ci_upper = binom_ci(
      n_continued_transmission,
      n_simulations
    )[["upper"]],
    .groups = "drop"
  )

library(dplyr)
library(tidyr)

panelC_summary <- runs_conditioned %>%
  filter(scenario == "imperfect_intervention") %>%
  group_by(
    transmission_regime,
    transmission_label,
    contacts_missed
  ) %>%
  summarise(
    n_simulations = n(),
    median_quarantined = median(
      downstream_uninfected_contacts_quarantined,
      na.rm = TRUE
    ),
    q25_quarantined = quantile(
      downstream_uninfected_contacts_quarantined,
      0.25,
      na.rm = TRUE,
      names = FALSE
    ),
    q75_quarantined = quantile(
      downstream_uninfected_contacts_quarantined,
      0.75,
      na.rm = TRUE,
      names = FALSE
    ),
    .groups = "drop"
  ) %>%
  mutate(
    missed_contacts = paste0(round(100 * contacts_missed), "%")
  )

panelC_summary %>%
  arrange(transmission_regime, contacts_missed) %>%
  print(n = Inf)

low_trend <- panelC_summary %>%
  filter(
    transmission_regime == "subcritical",
    abs(contacts_missed - 0.00) < 1e-8 |
      abs(contacts_missed - 0.50) < 1e-8
  ) %>%
  select(
    contacts_missed,
    median_quarantined,
    q25_quarantined,
    q75_quarantined
  ) %>%
  pivot_wider(
    names_from = contacts_missed,
    values_from = c(
      median_quarantined,
      q25_quarantined,
      q75_quarantined
    ),
    names_glue = "{.value}_{round(100 * contacts_missed)}"
  ) %>%
  mutate(
    absolute_change =
      median_quarantined_50 - median_quarantined_0,
    percent_change =
      100 * (
        median_quarantined_50 - median_quarantined_0
      ) / median_quarantined_0
  )

low_trend

