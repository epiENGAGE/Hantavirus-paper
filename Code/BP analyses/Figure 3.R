# Figure 3
parse_args <- function() {
  x <- list(
    input = getwd(),
    out = "~/Hantavirus-paper/BP results/Fig 3",
    n_runs = 10000L,
    seed = 20260529L,
    
    # Reviewer-style trigger values.
    reported_case_triggers = c(1L, 2L, 3L, 5L, 10L, 15L, 20L),
    
    # Use 95% as default, as suggested in feedback.
    q_values = c(0.95),
    
    max_days = 260,
    max_cases = 5000L,
    initial_infected = 1L,
    
    R0 = 2.12,
    k = 0.1982053,
    
    p_detect = 0.95,
    detection_delay_mean = 2.0,
    detection_delay_sd = 1.0,
    isolation_delay_mean = 0.5,
    isolation_delay_sd = 0.5,
    isolation_effectiveness = 0.99,
    
    quarantine_duration = 42.0,
    quarantine_effectiveness = 1.0,
    quarantine_adherence = 1.0,
    
    # Same duration-scaled contact-pool model as cruise script.
    backward_tracing_window = 14.0,
    daily_contact_rate = 13.50,
    contact_count_dispersion = 2.11,
    contact_multiplier = 1.0,
    
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
    else if (key %in% c("--q_values", "--q-values")) {
      x$q_values <- as.numeric(strsplit(val, ",", fixed = TRUE)[[1]])
    }
    else if (key %in% c("--max_days", "--max-days")) x$max_days <- as.numeric(val)
    else if (key %in% c("--max_cases", "--max-cases")) x$max_cases <- as.integer(val)
    else if (key %in% c("--initial_infected", "--initial-infected")) x$initial_infected <- as.integer(val)
    else if (key == "--R0") x$R0 <- as.numeric(val)
    else if (key == "--k") x$k <- as.numeric(val)
    else if (key == "--p_detect") x$p_detect <- as.numeric(val)
    else if (key %in% c("--detection_delay_mean", "--detection-delay-mean")) x$detection_delay_mean <- as.numeric(val)
    else if (key %in% c("--detection_delay_sd", "--detection-delay-sd")) x$detection_delay_sd <- as.numeric(val)
    else if (key %in% c("--isolation_delay_mean", "--isolation-delay-mean")) x$isolation_delay_mean <- as.numeric(val)
    else if (key %in% c("--isolation_delay_sd", "--isolation-delay-sd")) x$isolation_delay_sd <- as.numeric(val)
    else if (key %in% c("--isolation_effectiveness", "--isolation-effectiveness")) x$isolation_effectiveness <- as.numeric(val)
    else if (key %in% c("--quarantine_duration", "--quarantine-duration")) x$quarantine_duration <- as.numeric(val)
    else if (key %in% c("--quarantine_effectiveness", "--quarantine-effectiveness")) x$quarantine_effectiveness <- as.numeric(val)
    else if (key %in% c("--quarantine_adherence", "--quarantine-adherence")) x$quarantine_adherence <- as.numeric(val)
    else if (key %in% c("--backward_tracing_window", "--backward-tracing-window")) x$backward_tracing_window <- as.numeric(val)
    else if (key %in% c("--daily_contact_rate", "--daily-contact-rate")) x$daily_contact_rate <- as.numeric(val)
    else if (key %in% c("--contact_count_dispersion", "--contact-count-dispersion")) x$contact_count_dispersion <- as.numeric(val)
    else if (key %in% c("--contact_multiplier", "--contact-multiplier")) x$contact_multiplier <- as.numeric(val)
    else stop("Unknown argument: ", key, call. = FALSE)
    
    i <- i + 2L
  }
  
  x$reported_case_triggers <- sort(unique(x$reported_case_triggers[is.finite(x$reported_case_triggers) & x$reported_case_triggers >= 1]))
  x$q_values <- x$q_values[is.finite(x$q_values) & x$q_values >= 0 & x$q_values <= 1]
  
  if (length(x$reported_case_triggers) == 0) stop("No valid reported_case_triggers supplied.", call. = FALSE)
  if (length(x$q_values) == 0) stop("No valid q_values supplied.", call. = FALSE)
  
  x
}

required <- c("Rcpp", "dplyr", "tidyr", "ggplot2", "readr", "tibble", "scales", "patchwork")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  stop(
    "Missing packages: ", paste(missing, collapse = ", "),
    "\nInstall with install.packages(c(",
    paste(sprintf('"%s"', missing), collapse = ", "), "))",
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

quant <- function(x, p) as.numeric(stats::quantile(x, p, na.rm = TRUE, names = FALSE))
make_q_label <- function(q) paste0(round(100 * q), "%")

# Read posterior delay draws
input_dir <- normalizePath(args$input, winslash = "/", mustWork = FALSE)
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
    list.files(
      input_dir,
      pattern = "^(joint_delay_draws_for_branching|posterior_pair_level_draws)\\.csv$",
      recursive = TRUE,
      full.names = TRUE
    ),
    list.files(
      parent_dir,
      pattern = "^(joint_delay_draws_for_branching|posterior_pair_level_draws)\\.csv$",
      recursive = TRUE,
      full.names = TRUE
    )
  ))
  
  posterior_hits <- recursive_hits[basename(recursive_hits) == "posterior_pair_level_draws.csv"]
  branching_hits <- recursive_hits[basename(recursive_hits) == "joint_delay_draws_for_branching.csv"]
  
  if (length(posterior_hits) > 0) {
    delay_file <- posterior_hits[which.max(file.info(posterior_hits)$mtime)]
  } else if (length(branching_hits) > 0) {
    delay_file <- branching_hits[which.max(file.info(branching_hits)$mtime)]
  }
}

if (is.na(delay_file) || length(delay_file) == 0) {
  stop(
    "Could not find delay draws under --input or its parent: ", args$input,
    "\nLooked directly for:\n  ", paste(delay_candidates, collapse = "\n  "),
    "\nAlso tried recursive search for posterior_pair_level_draws.csv and joint_delay_draws_for_branching.csv.",
    call. = FALSE
  )
}

message("Reading delay draws from: ", normalizePath(delay_file, winslash = "/", mustWork = FALSE))
pair_draws_raw <- readr::read_csv(delay_file, show_col_types = FALSE)

if ("generation_time" %in% names(pair_draws_raw)) {
  pair_draws_raw <- pair_draws_raw %>% mutate(generation_time = as.numeric(generation_time))
} else if ("generation_interval" %in% names(pair_draws_raw)) {
  pair_draws_raw <- pair_draws_raw %>% mutate(generation_time = as.numeric(generation_interval))
} else {
  stop("Delay file does not contain generation_time or generation_interval. Columns are: ", paste(names(pair_draws_raw), collapse = ", "), call. = FALSE)
}

if ("recipient_incubation" %in% names(pair_draws_raw)) {
  pair_draws_raw <- pair_draws_raw %>% mutate(recipient_incubation = as.numeric(recipient_incubation))
} else if ("incubation" %in% names(pair_draws_raw)) {
  pair_draws_raw <- pair_draws_raw %>% mutate(recipient_incubation = as.numeric(incubation))
} else {
  stop("Delay file does not contain recipient_incubation or incubation. Columns are: ", paste(names(pair_draws_raw), collapse = ", "), call. = FALSE)
}

if ("delta" %in% names(pair_draws_raw)) {
  pair_draws_raw <- pair_draws_raw %>% mutate(delta = as.numeric(delta))
} else {
  stop("Delay file must contain delta for joint timing sampling. Columns are: ", paste(names(pair_draws_raw), collapse = ", "), call. = FALSE)
}

if ("source_incubation" %in% names(pair_draws_raw)) {
  pair_draws_raw <- pair_draws_raw %>% mutate(source_incubation = as.numeric(source_incubation))
} else {
  pair_draws_raw <- pair_draws_raw %>% mutate(source_incubation = generation_time - delta)
}

if ("serial_interval" %in% names(pair_draws_raw)) {
  pair_draws_raw <- pair_draws_raw %>% mutate(serial_interval = as.numeric(serial_interval))
} else {
  pair_draws_raw <- pair_draws_raw %>% mutate(serial_interval = delta + recipient_incubation)
}

pair_draws_joint <- pair_draws_raw %>%
  filter(
    is.finite(generation_time), generation_time > 0,
    is.finite(recipient_incubation), recipient_incubation > 0,
    is.finite(source_incubation), source_incubation > 0,
    is.finite(delta),
    is.finite(serial_interval)
  )

if (nrow(pair_draws_joint) == 0) stop("No usable joint delay draws found.", call. = FALSE)

readr::write_csv(
  pair_draws_joint,
  file.path(args$out, "tables", "joint_delay_draws_for_reported_case_trigger_used.csv")
)

generation_draws <- pair_draws_joint$generation_time
recipient_incubation_draws <- pair_draws_joint$recipient_incubation
source_incubation_draws <- pair_draws_joint$source_incubation
delta_draws <- pair_draws_joint$delta
serial_interval_draws <- pair_draws_joint$serial_interval

# C++ simulation
Rcpp::sourceCpp(code = '
#include <Rcpp.h>
#include <algorithm>
#include <numeric>
#include <cmath>
using namespace Rcpp;

inline double runif01() {
  return R::runif(0.0, 1.0);
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

inline int sample_index(int n) {
  int idx = (int)std::floor(R::runif(0.0, (double)n));
  if (idx < 0) idx = 0;
  if (idx >= n) idx = n - 1;
  return idx;
}

inline double sample_from_vec(const NumericVector& x) {
  int n = x.size();
  return x[sample_index(n)];
}

inline bool finite_num(double x) {
  return R_finite(x);
}

struct JointDelayDraw {
  double generation_time;
  double recipient_incubation;
  double delta;
  double serial_interval;
  int row_index;
};

inline int sample_joint_index_conditioned_on_source_ip(const NumericVector& source_incubation_draws,
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

inline JointDelayDraw sample_joint_delay(const NumericVector& generation_draws,
                                         const NumericVector& recipient_incubation_draws,
                                         const NumericVector& source_incubation_draws,
                                         const NumericVector& delta_draws,
                                         const NumericVector& serial_interval_draws,
                                         double source_incubation,
                                         double max_days) {
  int n = generation_draws.size();

  for (int guard = 0; guard < 200; ++guard) {
    int idx = sample_joint_index_conditioned_on_source_ip(source_incubation_draws, source_incubation);
    double inc = recipient_incubation_draws[idx];
    double delta = delta_draws[idx];
    double serial = serial_interval_draws[idx];
    double G = source_incubation + delta;

    if (R_finite(G) && G > 0.0 && G < max_days && R_finite(inc) && inc > 0.0) {
      JointDelayDraw out;
      out.generation_time = G;
      out.recipient_incubation = inc;
      out.delta = delta;
      out.serial_interval = serial;
      out.row_index = idx;
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
      out.delta = G - source_incubation;
      out.serial_interval = out.delta + inc;
      out.row_index = idx;
      return out;
    }
  }

  JointDelayDraw out;
  out.generation_time = 1.0;
  out.recipient_incubation = 21.0;
  out.delta = out.generation_time - source_incubation;
  out.serial_interval = out.delta + out.recipient_incubation;
  out.row_index = -1;
  return out;
}

struct PotentialTx {
  int source;
  int child;
  double time;
};

// [[Rcpp::export]]
DataFrame run_reported_case_trigger_cpp(NumericVector q_values,
                                        IntegerVector reported_case_triggers,
                                        int n_runs,
                                        NumericVector generation_draws,
                                        NumericVector recipient_incubation_draws,
                                        NumericVector source_incubation_draws,
                                        NumericVector delta_draws,
                                        NumericVector serial_interval_draws,
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
                                        double backward_tracing_window,
                                        double daily_contact_rate,
                                        double contact_count_dispersion,
                                        double contact_multiplier,
                                        int seed) {
  RNGScope scope;
  Environment base_env("package:base");
  Function set_seed = base_env["set.seed"];
  set_seed(seed);

  int n_q = q_values.size();
  int n_trig = reported_case_triggers.size();
  int total_rows = n_q * n_trig * n_runs;

  IntegerVector out_run_id(total_rows);
  NumericVector out_q(total_rows), out_contacts_missed(total_rows), out_response_time(total_rows);
  IntegerVector out_reported_case_trigger(total_rows), out_reached_trigger(total_rows);
  IntegerVector out_total_potential_infections(total_rows), out_total_retained_infections(total_rows);
  IntegerVector out_reported_cases_at_trigger(total_rows), out_detected_total_potential(total_rows);
  IntegerVector out_infected_quarantined_at_response(total_rows), out_future_infected_quarantined(total_rows);
  IntegerVector out_future_uninfected_contacts_quarantined(total_rows);
  IntegerVector out_total_quarantined(total_rows), out_additional_quarantine(total_rows);
  IntegerVector out_post_response_infections(total_rows), out_post_response_transmissions(total_rows);
  IntegerVector out_prevented_by_quarantine_total(total_rows), out_prevented_by_isolation_total(total_rows);
  NumericVector out_response_duration_after_response(total_rows);

  int row = 0;

  for (int qidx = 0; qidx < n_q; ++qidx) {
    double q_system = q_values[qidx];

    for (int tidx = 0; tidx < n_trig; ++tidx) {
      int trigger_count = reported_case_triggers[tidx];

      for (int rr = 0; rr < n_runs; ++rr) {
        std::vector<int> infector; infector.reserve(max_cases);
        std::vector<int> generation; generation.reserve(max_cases);
        std::vector<double> infection_time; infection_time.reserve(max_cases);
        std::vector<double> incubation_period; incubation_period.reserve(max_cases);
        std::vector<double> symptom_onset; symptom_onset.reserve(max_cases);
        std::vector<int> detected; detected.reserve(max_cases);
        std::vector<double> detection_time; detection_time.reserve(max_cases);
        std::vector<double> isolation_time; isolation_time.reserve(max_cases);
        std::vector<int> processed; processed.reserve(max_cases);
        std::vector<PotentialTx> potential_tx; potential_tx.reserve(max_cases * 2);

        auto add_case = [&](int inf, int gen, double itime, double inc) {
          double onset = itime + inc;
          int det = (runif01() < p_detect) ? 1 : 0;
          double dtime = R_PosInf;
          double isotime = R_PosInf;
          if (det == 1) {
            dtime = onset + gamma_delay_one(detection_delay_mean, detection_delay_sd, 0.0);
            isotime = dtime + gamma_delay_one(isolation_delay_mean, isolation_delay_sd, 0.0);
          }

          infector.push_back(inf);
          generation.push_back(gen);
          infection_time.push_back(itime);
          incubation_period.push_back(inc);
          symptom_onset.push_back(onset);
          detected.push_back(det);
          detection_time.push_back(dtime);
          isolation_time.push_back(isotime);
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
              serial_interval_draws,
              incubation_period[todo],
              max_days
            );
            double G = dd.generation_time;
            double child_time = infection_time[todo] + G;

            if (!finite_num(child_time) || child_time > max_days) continue;

            int child_idx = add_case(todo, generation[todo] + 1, child_time, dd.recipient_incubation);

            PotentialTx tx;
            tx.source = todo;
            tx.child = child_idx;
            tx.time = child_time;
            potential_tx.push_back(tx);
          }
        }

        int n_potential = (int)infection_time.size();

        std::vector<double> det_times;
        det_times.reserve(n_potential);

        for (int ii = 0; ii < n_potential; ++ii) {
          if (detected[ii] == 1 && finite_num(detection_time[ii])) {
            det_times.push_back(detection_time[ii]);
          }
        }

        std::sort(det_times.begin(), det_times.end());

        bool reached_trigger = ((int)det_times.size() >= trigger_count);
        double response_time = reached_trigger ? det_times[trigger_count - 1] : R_PosInf;

        std::vector<int> retained(n_potential, 0);
        std::vector<int> quarantined(n_potential, 0);
        std::vector<double> q_start(n_potential, R_PosInf);
        std::vector<double> q_end(n_potential, R_PosInf);
        std::vector<int> retained_offspring(n_potential, 0);

        for (int ii = 0; ii < initial_infected && ii < n_potential; ++ii) {
          retained[ii] = 1;
        }

        int infected_quarantined_at_response = 0;
        int future_infected_quarantined = 0;
        int future_uninfected_contacts_quarantined = 0;

        if (reached_trigger) {
          for (int ii = 0; ii < n_potential; ++ii) {
            if (infection_time[ii] <= response_time && runif01() < q_system) {
              quarantined[ii] = 1;
              q_start[ii] = response_time;
              q_end[ii] = response_time + quarantine_duration;
              infected_quarantined_at_response++;
            }
          }
        }

        std::vector<int> tx_order(potential_tx.size());
        std::iota(tx_order.begin(), tx_order.end(), 0);

        std::sort(tx_order.begin(), tx_order.end(), [&](int a, int b) {
          return potential_tx[a].time < potential_tx[b].time;
        });

        int post_response_transmissions = 0;
        int prevented_by_quarantine_total = 0;
        int prevented_by_isolation_total = 0;

        for (int oo = 0; oo < (int)tx_order.size(); ++oo) {
          PotentialTx tx = potential_tx[tx_order[oo]];

          int s = tx.source;
          int c = tx.child;
          double t = tx.time;

          if (s < 0 || s >= n_potential || c < 0 || c >= n_potential) continue;
          if (retained[s] == 0) continue;

          bool after_response = reached_trigger && t >= response_time;
          bool source_quarantined = after_response &&
                                    quarantined[s] == 1 &&
                                    finite_num(q_start[s]) &&
                                    t >= q_start[s] &&
                                    t <= q_end[s];

          bool source_isolated = after_response &&
                                 detected[s] == 1 &&
                                 finite_num(isolation_time[s]) &&
                                 t >= isolation_time[s];

          double retain_prob = 1.0;

          if (after_response) {
            if (source_quarantined && runif01() < quarantine_adherence) {
              retain_prob *= (1.0 - quarantine_effectiveness);
            }
            if (source_isolated) {
              retain_prob *= (1.0 - isolation_effectiveness);
            }
          }

          if (runif01() > retain_prob) {
            if (source_quarantined) prevented_by_quarantine_total++;
            else if (source_isolated) prevented_by_isolation_total++;
            continue;
          }

          retained[c] = 1;
          retained_offspring[s] += 1;

          if (after_response) {
            post_response_transmissions++;

            if (runif01() < q_system) {
              double trace_time = detection_time[c];
              if (!finite_num(trace_time) || trace_time < t) trace_time = t;

              if (quarantined[c] == 0) {
                future_infected_quarantined++;
              }

              quarantined[c] = 1;
              q_start[c] = trace_time;
              q_end[c] = trace_time + quarantine_duration;
            }
          }
        }

        int total_retained = 0;
        int post_response_infections = 0;
        double max_end = NA_REAL;

        for (int ii = 0; ii < n_potential; ++ii) {
          if (retained[ii] == 0) continue;

          total_retained++;

          if (reached_trigger && infection_time[ii] >= response_time) {
            post_response_infections++;
          }

          if (reached_trigger) {
            if (infection_time[ii] >= response_time) {
              if (!finite_num(max_end) || infection_time[ii] > max_end) max_end = infection_time[ii];
              if (finite_num(symptom_onset[ii]) && (!finite_num(max_end) || symptom_onset[ii] > max_end)) max_end = symptom_onset[ii];
              if (finite_num(detection_time[ii]) && (!finite_num(max_end) || detection_time[ii] > max_end)) max_end = detection_time[ii];
            }
            if (finite_num(q_end[ii]) && (!finite_num(max_end) || q_end[ii] > max_end)) max_end = q_end[ii];
          }
        }

        if (reached_trigger && daily_contact_rate > 0.0 && contact_multiplier > 0.0) {
          for (int jj = 0; jj < n_potential; ++jj) {
            if (retained[jj] == 0) continue;
            if (detected[jj] != 1 || !finite_num(detection_time[jj]) || detection_time[jj] < response_time) continue;

            bool already_quarantined_at_detection =
              quarantined[jj] == 1 &&
              finite_num(q_start[jj]) &&
              detection_time[jj] >= q_start[jj] &&
              detection_time[jj] <= q_end[jj];

            if (already_quarantined_at_detection) continue;

            double trace_start = std::max(symptom_onset[jj] - backward_tracing_window, infection_time[jj]);
            double trace_end = detection_time[jj];
            double trace_duration = std::max(0.0, trace_end - trace_start);

            int infected_contacts_in_trace_pool = 0;

            for (int cc = 0; cc < n_potential; ++cc) {
              if (retained[cc] == 0) continue;
              if (infector[cc] != jj) continue;
              if (infection_time[cc] >= trace_start && infection_time[cc] <= trace_end) {
                infected_contacts_in_trace_pool++;
              }
            }

            double mean_contacts = daily_contact_rate * contact_multiplier * trace_duration;
            double size_contacts = contact_count_dispersion * trace_duration;

            int total_traceable_contacts = infected_contacts_in_trace_pool;

            if (mean_contacts > 0.0) {
              if (size_contacts > 0.0) {
                double prob_contacts = size_contacts / (size_contacts + mean_contacts);
                total_traceable_contacts = (int)R::rnbinom(size_contacts, prob_contacts);
              } else {
                total_traceable_contacts = (int)R::rpois(mean_contacts);
              }

              if (total_traceable_contacts < infected_contacts_in_trace_pool) {
                total_traceable_contacts = infected_contacts_in_trace_pool;
              }
            }

            int uninfected_contacts = total_traceable_contacts - infected_contacts_in_trace_pool;

            int traced_contacts = 0;
            if (total_traceable_contacts > 0) {
              traced_contacts = (int)R::rbinom((double)total_traceable_contacts, q_system);
            }

            int traced_uninfected_contacts = 0;

            if (traced_contacts > 0 && uninfected_contacts > 0) {
              traced_uninfected_contacts = (int)R::rhyper(
                (double)uninfected_contacts,
                (double)infected_contacts_in_trace_pool,
                (double)traced_contacts
              );
            }

            future_uninfected_contacts_quarantined += traced_uninfected_contacts;
          }
        }

        int total_quarantined =
          infected_quarantined_at_response +
          future_infected_quarantined +
          future_uninfected_contacts_quarantined;

        double response_duration = NA_REAL;
        if (reached_trigger) {
          if (finite_num(max_end)) response_duration = std::max(0.0, max_end - response_time);
          else response_duration = 0.0;
        }

        out_run_id[row] = row + 1;
        out_q[row] = q_system;
        out_contacts_missed[row] = 1.0 - q_system;
        out_reported_case_trigger[row] = trigger_count;
        out_reached_trigger[row] = reached_trigger ? 1 : 0;
        out_response_time[row] = reached_trigger ? response_time : NA_REAL;
        out_total_potential_infections[row] = n_potential;
        out_total_retained_infections[row] = total_retained;
        out_reported_cases_at_trigger[row] = reached_trigger ? trigger_count : (int)det_times.size();
        out_detected_total_potential[row] = (int)det_times.size();
        out_infected_quarantined_at_response[row] = infected_quarantined_at_response;
        out_future_infected_quarantined[row] = future_infected_quarantined;
        out_future_uninfected_contacts_quarantined[row] = future_uninfected_contacts_quarantined;
        out_total_quarantined[row] = total_quarantined;
        out_additional_quarantine[row] = total_quarantined;
        out_post_response_infections[row] = post_response_infections;
        out_post_response_transmissions[row] = post_response_transmissions;
        out_prevented_by_quarantine_total[row] = prevented_by_quarantine_total;
        out_prevented_by_isolation_total[row] = prevented_by_isolation_total;
        out_response_duration_after_response[row] = response_duration;

        row++;
      }
    }
  }

  return DataFrame::create(
    Named("run_id") = out_run_id,
    Named("q_system") = out_q,
    Named("contacts_missed") = out_contacts_missed,
    Named("reported_case_trigger") = out_reported_case_trigger,
    Named("reached_trigger") = out_reached_trigger,
    Named("response_time") = out_response_time,
    Named("total_potential_infections") = out_total_potential_infections,
    Named("total_retained_infections") = out_total_retained_infections,
    Named("reported_cases_at_trigger") = out_reported_cases_at_trigger,
    Named("detected_total_potential") = out_detected_total_potential,
    Named("infected_quarantined_at_response") = out_infected_quarantined_at_response,
    Named("future_infected_quarantined") = out_future_infected_quarantined,
    Named("future_uninfected_contacts_quarantined") = out_future_uninfected_contacts_quarantined,
    Named("total_quarantined") = out_total_quarantined,
    Named("additional_quarantine") = out_additional_quarantine,
    Named("post_response_infections") = out_post_response_infections,
    Named("post_response_transmissions") = out_post_response_transmissions,
    Named("prevented_by_quarantine_total") = out_prevented_by_quarantine_total,
    Named("prevented_by_isolation_total") = out_prevented_by_isolation_total,
    Named("response_duration_after_response") = out_response_duration_after_response
  );
}
')

# Run simulations
transmission_regimes <- tibble::tribble(
  ~transmission_regime, ~transmission_label, ~R0,       ~k,
  "supercritical",      "High transmission", 2.12,      0.3157164,
  "subcritical",        "Low transmission",  0.2416027, 0.3998142
)

binom_ci <- function(events, n) {
  if (!is.finite(events) || !is.finite(n) || n <= 0) {
    return(c(lower = NA_real_, upper = NA_real_))
  }
  z <- stats::binom.test(events, n)$conf.int
  c(lower = unname(z[1]), upper = unname(z[2]))
}

fixed_inputs <- bind_rows(
  tibble(
    parameter = c(
      "input", "delay_file_used", "out", "n_runs_per_regime_trigger",
      "reported_case_triggers", "q_values", "seed",
      "max_days", "max_cases", "initial_infected",
      "timing_sampling", "detection_delay_distribution",
      "isolation_delay_distribution", "p_detect",
      "detection_delay_mean", "detection_delay_sd",
      "isolation_delay_mean", "isolation_delay_sd",
      "isolation_effectiveness", "quarantine_duration",
      "quarantine_effectiveness", "quarantine_adherence",
      "backward_tracing_window", "daily_contact_rate",
      "contact_count_dispersion", "contact_multiplier"
    ),
    value = as.character(c(
      args$input, delay_file, args$out, args$n_runs,
      paste(args$reported_case_triggers, collapse = ","),
      paste(args$q_values, collapse = ","),
      args$seed, args$max_days, args$max_cases, args$initial_infected,
      "joint_delta_conditioned_on_source_incubation",
      "gamma_mean_sd", "gamma_mean_sd", args$p_detect,
      args$detection_delay_mean, args$detection_delay_sd,
      args$isolation_delay_mean, args$isolation_delay_sd,
      args$isolation_effectiveness, args$quarantine_duration,
      args$quarantine_effectiveness, args$quarantine_adherence,
      args$backward_tracing_window, args$daily_contact_rate,
      args$contact_count_dispersion, args$contact_multiplier
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
  file.path(args$out, "tables", "reported_case_trigger_fixed_inputs.csv")
)

message("Running reported-case-triggered analysis for both transmission regimes...")
t0 <- Sys.time()

runs_list <- vector("list", nrow(transmission_regimes))

for (ii in seq_len(nrow(transmission_regimes))) {
  reg <- transmission_regimes[ii, ]
  regime_seed <- as.integer(args$seed + (ii - 1L) * 1000003L)
  
  message(
    "  ", reg$transmission_label,
    " (R0=", reg$R0, ", k=", reg$k, ")"
  )
  
  runs_list[[ii]] <- run_reported_case_trigger_cpp(
    q_values = args$q_values,
    reported_case_triggers = args$reported_case_triggers,
    n_runs = args$n_runs,
    generation_draws = generation_draws,
    recipient_incubation_draws = recipient_incubation_draws,
    source_incubation_draws = source_incubation_draws,
    delta_draws = delta_draws,
    serial_interval_draws = serial_interval_draws,
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
    backward_tracing_window = args$backward_tracing_window,
    daily_contact_rate = args$daily_contact_rate,
    contact_count_dispersion = args$contact_count_dispersion,
    contact_multiplier = args$contact_multiplier,
    seed = regime_seed
  ) %>%
    as_tibble() %>%
    mutate(
      transmission_regime = reg$transmission_regime,
      transmission_label = reg$transmission_label,
      R0_regime = reg$R0,
      k_regime = reg$k,
      regime_seed = regime_seed
    )
}

runs <- bind_rows(runs_list) %>%
  mutate(
    q_label = factor(
      make_q_label(q_system),
      levels = make_q_label(sort(unique(q_system), decreasing = TRUE))
    ),
    response_duration_weeks = response_duration_after_response / 7
  )

message(
  "Simulation time: ",
  round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 2),
  " seconds"
)

readr::write_csv(
  runs,
  file.path(args$out, "tables", "reported_case_trigger_run_level_results.csv")
)

trigger_attainment <- runs %>%
  group_by(
    transmission_regime, transmission_label, R0_regime, k_regime,
    reported_case_trigger, q_system, q_label
  ) %>%
  summarise(
    n_generated = n(),
    n_reached_trigger = sum(reached_trigger == 1, na.rm = TRUE),
    probability_reached_trigger = mean(reached_trigger == 1, na.rm = TRUE),
    .groups = "drop"
  )

readr::write_csv(
  trigger_attainment,
  file.path(args$out, "tables", "reported_case_trigger_attainment.csv")
)

runs_conditioned <- runs %>%
  filter(reached_trigger == 1)

if (nrow(runs_conditioned) == 0) {
  stop("No simulations reached any reported-case trigger.", call. = FALSE)
}

summary <- runs_conditioned %>%
  group_by(
    transmission_regime, transmission_label, R0_regime, k_regime,
    reported_case_trigger, q_system, q_label
  ) %>%
  summarise(
    n_runs = n(),
    
    median_total_infections = median(total_retained_infections, na.rm = TRUE),
    q25_total_infections = quant(total_retained_infections, 0.25),
    q75_total_infections = quant(total_retained_infections, 0.75),
    q90_total_infections = quant(total_retained_infections, 0.90),
    q95_total_infections = quant(total_retained_infections, 0.95),
    q99_total_infections = quant(total_retained_infections, 0.99),
    
    median_total_quarantined = median(total_quarantined, na.rm = TRUE),
    q25_total_quarantined = quant(total_quarantined, 0.25),
    q75_total_quarantined = quant(total_quarantined, 0.75),
    q90_total_quarantined = quant(total_quarantined, 0.90),
    q95_total_quarantined = quant(total_quarantined, 0.95),
    q99_total_quarantined = quant(total_quarantined, 0.99),
    
    p_infections_ge_100 = mean(total_retained_infections >= 100, na.rm = TRUE),
    p_infections_ge_500 = mean(total_retained_infections >= 500, na.rm = TRUE),
    p_infections_ge_1000 = mean(total_retained_infections >= 1000, na.rm = TRUE),
    
    p_quarantine_ge_1000 = mean(total_quarantined >= 1000, na.rm = TRUE),
    p_quarantine_ge_5000 = mean(total_quarantined >= 5000, na.rm = TRUE),
    p_quarantine_ge_10000 = mean(total_quarantined >= 10000, na.rm = TRUE),
    
    median_response_duration_weeks =
      median(response_duration_weeks, na.rm = TRUE),
    q95_response_duration_weeks =
      quant(response_duration_weeks, 0.95),
    
    .groups = "drop"
  ) %>%
  arrange(transmission_regime, q_system, reported_case_trigger)

readr::write_csv(
  summary,
  file.path(args$out, "tables", "reported_case_trigger_summary.csv")
)

# Exact binomial CIs for the main threshold probabilities.
probability_summary <- runs_conditioned %>%
  group_by(
    transmission_regime, transmission_label,
    reported_case_trigger, q_system, q_label
  ) %>%
  summarise(
    n_runs = n(),
    infections_ge_100 = sum(total_retained_infections >= 100, na.rm = TRUE),
    infections_ge_500 = sum(total_retained_infections >= 500, na.rm = TRUE),
    quarantine_ge_1000 = sum(total_quarantined >= 1000, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    p_infections_ge_100 = infections_ge_100 / n_runs,
    infections_ge_100_lower = binom_ci(infections_ge_100, n_runs)[["lower"]],
    infections_ge_100_upper = binom_ci(infections_ge_100, n_runs)[["upper"]],
    
    p_infections_ge_500 = infections_ge_500 / n_runs,
    infections_ge_500_lower = binom_ci(infections_ge_500, n_runs)[["lower"]],
    infections_ge_500_upper = binom_ci(infections_ge_500, n_runs)[["upper"]],
    
    p_quarantine_ge_1000 = quarantine_ge_1000 / n_runs,
    quarantine_ge_1000_lower = binom_ci(quarantine_ge_1000, n_runs)[["lower"]],
    quarantine_ge_1000_upper = binom_ci(quarantine_ge_1000, n_runs)[["upper"]]
  ) %>%
  ungroup()

readr::write_csv(
  probability_summary,
  file.path(
    args$out,
    "tables",
    "reported_case_trigger_threshold_probabilities_with_CI.csv"
  )
)

# Key results
key_contrasts <- summary %>%
  filter(
    abs(q_system - 0.95) < 1e-8,
    reported_case_trigger %in% c(1L, 20L)
  ) %>%
  select(
    transmission_regime, transmission_label, reported_case_trigger,
    n_runs,
    median_total_infections, q25_total_infections, q75_total_infections,
    q95_total_infections,
    median_total_quarantined, q25_total_quarantined,
    q75_total_quarantined, q95_total_quarantined,
    p_infections_ge_100, p_infections_ge_500,
    p_quarantine_ge_1000
  )

readr::write_csv(
  key_contrasts,
  file.path(args$out, "tables", "reported_case_trigger_key_contrasts.csv")
)

fmt_num <- function(x, digits = 0) {
  format(round(x, digits), big.mark = ",", scientific = FALSE, trim = TRUE)
}
fmt_pct <- function(x, digits = 1) {
  paste0(sprintf(paste0("%.", digits, "f"), 100 * x), "%")
}

text_lines <- c(
  "RESULTS-SECTION QUANTITIES",
  "All outcome summaries below condition on the requested reported-case trigger being reached.",
  ""
)

for (regime in c("supercritical", "subcritical")) {
  z <- summary %>%
    filter(
      transmission_regime == regime,
      abs(q_system - 0.95) < 1e-8,
      reported_case_trigger %in% c(1L, 20L)
    ) %>%
    arrange(reported_case_trigger)
  
  if (nrow(z) == 2) {
    a <- z[1, ]
    b <- z[2, ]
    
    text_lines <- c(
      text_lines,
      paste0(a$transmission_label, ":"),
      paste0(
        "  Trigger 1: median infections ", fmt_num(a$median_total_infections),
        " (IQR ", fmt_num(a$q25_total_infections), "--",
        fmt_num(a$q75_total_infections), "); median quarantined ",
        fmt_num(a$median_total_quarantined), " (IQR ",
        fmt_num(a$q25_total_quarantined), "--",
        fmt_num(a$q75_total_quarantined), ")."
      ),
      paste0(
        "  Trigger 20: median infections ", fmt_num(b$median_total_infections),
        " (IQR ", fmt_num(b$q25_total_infections), "--",
        fmt_num(b$q75_total_infections), "); median quarantined ",
        fmt_num(b$median_total_quarantined), " (IQR ",
        fmt_num(b$q25_total_quarantined), "--",
        fmt_num(b$q75_total_quarantined), ")."
      ),
      paste0(
        "  Probability of at least 100 infections: ",
        fmt_pct(a$p_infections_ge_100), " at trigger 1 versus ",
        fmt_pct(b$p_infections_ge_100), " at trigger 20."
      ),
      paste0(
        "  Probability of at least 500 infections: ",
        fmt_pct(a$p_infections_ge_500), " at trigger 1 versus ",
        fmt_pct(b$p_infections_ge_500), " at trigger 20."
      ),
      paste0(
        "  Probability of quarantining at least 1,000 people: ",
        fmt_pct(a$p_quarantine_ge_1000), " at trigger 1 versus ",
        fmt_pct(b$p_quarantine_ge_1000), " at trigger 20."
      ),
      ""
    )
  }
}

writeLines(
  text_lines,
  file.path(args$out, "tables", "reported_case_trigger_results_text.txt")
)
cat(paste(text_lines, collapse = "\n"), "\n")

# Figure 3: high- and low-transmission scenarios
if (isTRUE(args$make_plots)) {
  fig_dir <- file.path(args$out, "figures")
  
  q_focus <- if (any(abs(runs_conditioned$q_system - 0.95) < 1e-8)) {
    0.95
  } else {
    max(runs_conditioned$q_system, na.rm = TRUE)
  }
  
  plot_runs <- runs_conditioned %>%
    filter(abs(q_system - q_focus) < 1e-8) %>%
    mutate(
      trigger_f = factor(
        reported_case_trigger,
        levels = sort(unique(reported_case_trigger))
      ),
      transmission_label = factor(
        transmission_label,
        levels = c("High transmission", "Low transmission")
      ),
      total_retained_infections_plot = pmax(total_retained_infections, 1),
      total_quarantined_plot = pmax(total_quarantined, 1)
    )
  
  theme_pub <- theme_classic(base_size = 16) +
    theme(
      axis.text = element_text(colour = "black", size = 13),
      axis.title = element_text(size = 16),
      axis.title.x = element_text(margin = margin(t = 9)),
      axis.title.y = element_text(margin = margin(r = 9)),
      panel.grid.major.y = element_line(colour = "grey87", linewidth = 0.4),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 13),
      plot.tag = element_text(face = "bold", size = 18),
      plot.margin = margin(12, 14, 12, 14)
    )
  
  fill_values <- c(
    "High transmission" = "grey60",
    "Low transmission" = "grey15"
  )
  
  pA <- ggplot(
    plot_runs,
    aes(
      x = trigger_f,
      y = total_retained_infections_plot,
      fill = transmission_label
    )
  ) +
    geom_boxplot(
      position = position_dodge(width = 0.78),
      width = 0.66,
      linewidth = 0.42,
      colour = "black",
      outlier.shape = 16,
      outlier.size = 0.55,
      outlier.alpha = 0.18
    ) +
    scale_fill_manual(values = fill_values, drop = FALSE) +
    scale_y_log10(
      breaks = c(1, 2, 5, 10, 25, 50, 100, 250, 500,
                 1000, 2500, 5000, 10000, 25000, 50000, 100000),
      labels = scales::label_number(big.mark = ",")
    ) +
    labs(
      x = "Reported cases before intervention",
      y = "Total infections",
      tag = "A"
    ) +
    theme_pub
  
  pB <- ggplot(
    plot_runs,
    aes(
      x = trigger_f,
      y = total_quarantined_plot,
      fill = transmission_label
    )
  ) +
    geom_boxplot(
      position = position_dodge(width = 0.78),
      width = 0.66,
      linewidth = 0.42,
      colour = "black",
      outlier.shape = 16,
      outlier.size = 0.55,
      outlier.alpha = 0.18
    ) +
    scale_fill_manual(values = fill_values, drop = FALSE) +
    scale_y_log10(
      breaks = c(1, 2, 5, 10, 25, 50, 100, 250, 500,
                 1000, 2500, 5000, 10000, 25000, 50000, 100000),
      labels = scales::label_number(big.mark = ",")
    ) +
    labs(
      x = "Reported cases before intervention",
      y = "Individuals requiring quarantine",
      tag = "B"
    ) +
    theme_pub
  
  fig3 <- (pA | pB) +
    patchwork::plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  
  ggsave(
    file.path(fig_dir, "Figure3_reported_case_trigger_high_low.png"),
    fig3,
    width = 13.5,
    height = 6.3,
    dpi = 800,
    bg = "white"
  )
  
  ggsave(
    file.path(fig_dir, "Figure3_reported_case_trigger_high_low.pdf"),
    fig3,
    width = 13.5,
    height = 6.3,
    bg = "white"
  )
}

message("Done.")
message("Outputs written to: ", normalizePath(args$out, winslash = "/", mustWork = FALSE))
message(
  "Figure 3: ",
  normalizePath(
    file.path(args$out, "figures", "Figure3_reported_case_trigger_high_low.png"),
    winslash = "/", mustWork = FALSE
  )
)
