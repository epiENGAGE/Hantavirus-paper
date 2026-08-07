# Analysis for figure 6
parse_args <- function() {
  x <- list(
    input = "~/Hantavirus-paper/TD results",
    out = "~/Hantavirus-paper/BP results/Fig 6",
    n_runs = 10000L,
    seed = 20260528L,
    q_values = c(1.00, 0.95, 0.90, 0.75, 0.50),
    start_generations = 0:4,
    current_generation = 2L,
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
    isolation_effectiveness = 0.95,
    quarantine_duration = 42.0,
    quarantine_effectiveness = 1.0,
    quarantine_adherence = 1.0,
    backward_tracing_window = 14.0,
    daily_contact_rate = 13.50,
    contact_count_dispersion = 2.11,
    contact_multiplier = 1.0,
    make_plots = TRUE,
    mean_uninfected_contacts_per_detected_case = 13.50,
    uninfected_contact_dispersion = 2.11
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
    else if (key %in% c("--q_values", "--q-values")) x$q_values <- as.numeric(strsplit(val, ",", fixed = TRUE)[[1]])
    else if (key %in% c("--start_generations", "--start-generations")) x$start_generations <- as.integer(strsplit(val, ",", fixed = TRUE)[[1]])
    else if (key %in% c("--current_generation", "--current-generation")) x$current_generation <- as.integer(val)
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
  
  x$q_values <- x$q_values[is.finite(x$q_values) & x$q_values >= 0 & x$q_values <= 1]
  x$start_generations <- sort(unique(x$start_generations[is.finite(x$start_generations) & x$start_generations >= 0]))
  if (length(x$q_values) == 0) stop("No valid q_values supplied.", call. = FALSE)
  if (length(x$start_generations) == 0) stop("No valid start_generations supplied.", call. = FALSE)
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

# -----------------------------
# Read posterior delay draws
# -----------------------------
# Try the supplied directory first, then its parent. This is deliberately
# tolerant because some runs store posterior_pair_level_draws.csv in the main
# delay-output directory, while others store branching inputs in a
# bp_from_mechanistic_delays/ subdirectory.
input_dir <- normalizePath(args$input, winslash = "/", mustWork = FALSE)
parent_dir <- normalizePath(dirname(input_dir), winslash = "/", mustWork = FALSE)

delay_candidates <- unique(c(
  file.path(input_dir, "tables", "joint_delay_draws_for_branching.csv"),
  file.path(input_dir, "joint_delay_draws_for_branching.csv"),
  file.path(input_dir, "posterior_pair_level_draws.csv"),
  file.path(input_dir, "tables", "posterior_pair_level_draws.csv"),
  file.path(parent_dir, "tables", "joint_delay_draws_for_branching.csv"),
  file.path(parent_dir, "joint_delay_draws_for_branching.csv"),
  file.path(parent_dir, "posterior_pair_level_draws.csv"),
  file.path(parent_dir, "tables", "posterior_pair_level_draws.csv")
))

delay_file <- delay_candidates[file.exists(delay_candidates)][1]

# Last-resort recursive search, limited to the input directory and its parent.
# This catches layouts like */some_run_name/tables/posterior_pair_level_draws.csv.
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
  
  # Prefer the direct observed-pair posterior if both are available.
  posterior_hits <- recursive_hits[basename(recursive_hits) == "posterior_pair_level_draws.csv"]
  branching_hits <- recursive_hits[basename(recursive_hits) == "joint_delay_draws_for_branching.csv"]
  
  if (length(posterior_hits) > 0) {
    delay_file <- posterior_hits[1]
  } else if (length(branching_hits) > 0) {
    delay_file <- branching_hits[1]
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
  stop("Delay file does not contain generation_time or generation_interval.", call. = FALSE)
}

if ("recipient_incubation" %in% names(pair_draws_raw)) {
  pair_draws_raw <- pair_draws_raw %>% mutate(recipient_incubation = as.numeric(recipient_incubation))
} else if ("incubation" %in% names(pair_draws_raw)) {
  pair_draws_raw <- pair_draws_raw %>% mutate(recipient_incubation = as.numeric(incubation))
} else {
  stop("Delay file does not contain recipient_incubation or incubation.", call. = FALSE)
}

pair_draws_joint <- pair_draws_raw %>%
  filter(
    is.finite(generation_time), generation_time > 0,
    is.finite(recipient_incubation), recipient_incubation > 0
  )

if (nrow(pair_draws_joint) == 0) stop("No usable delay draws found.", call. = FALSE)

readr::write_csv(
  pair_draws_joint,
  file.path(args$out, "tables", "joint_delay_draws_for_generation_start_used.csv")
)

generation_draws <- pair_draws_joint$generation_time
incubation_draws <- pair_draws_joint$recipient_incubation

# -----------------------------
# C++ simulation core
# -----------------------------
Rcpp::sourceCpp(code = '
#include <Rcpp.h>
#include <algorithm>
#include <numeric>
using namespace Rcpp;

inline double runif01() {
  return R::runif(0.0, 1.0);
}

inline double positive_normal_one(double mean, double sd, double min_value = 0.0) {
  if (!R_finite(sd) || sd <= 0.0) return std::max(mean, min_value);
  double y = R::rnorm(mean, sd);
  int guard = 0;
  while ((!R_finite(y) || y < min_value) && guard < 100) {
    y = R::rnorm(mean, sd);
    guard++;
  }
  if (!R_finite(y) || y < min_value) y = min_value;
  return y;
}

inline double sample_from_vec(const NumericVector& x) {
  int n = x.size();
  int idx = (int)std::floor(R::runif(0.0, (double)n));
  if (idx < 0) idx = 0;
  if (idx >= n) idx = n - 1;
  return x[idx];
}

inline bool finite_num(double x) {
  return R_finite(x);
}

struct PotentialTx {
  int source;
  int child;
  double time;
};

// [[Rcpp::export]]
DataFrame run_generation_start_sensitivity_cpp(NumericVector q_values,
                                               IntegerVector start_generations,
                                               int n_runs,
                                               NumericVector generation_draws,
                                               NumericVector incubation_draws,
                                               int current_generation,
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
                                               double mean_uninfected_contacts_per_detected_case,
                                               double uninfected_contact_dispersion,
                                               int seed) {
  RNGScope scope;
  Environment base_env("package:base");
  Function set_seed = base_env["set.seed"];
  set_seed(seed);

  int n_q = q_values.size();
  int n_start = start_generations.size();
  int total_rows = n_q * n_start * n_runs;

  IntegerVector out_run_id(total_rows);
  NumericVector out_q(total_rows), out_contacts_missed(total_rows), out_response_time(total_rows);
  IntegerVector out_start_generation(total_rows), out_current_generation(total_rows);
  IntegerVector out_reached_start_generation(total_rows), out_reached_current_generation(total_rows);
  IntegerVector out_total_potential_infections(total_rows), out_total_retained_infections(total_rows);
  IntegerVector out_additional_cases_after_current_gen(total_rows), out_any_additional_cases_after_current_gen(total_rows);
  IntegerVector out_post_response_infections(total_rows), out_post_response_transmissions(total_rows);
  IntegerVector out_additional_quarantine(total_rows), out_prevented_by_quarantine_total(total_rows), out_prevented_by_isolation_total(total_rows);
  NumericVector out_response_duration_after_response(total_rows);

  int row = 0;

  for (int qidx = 0; qidx < n_q; ++qidx) {
    double q_system = q_values[qidx];

    for (int sgidx = 0; sgidx < n_start; ++sgidx) {
      int start_generation = start_generations[sgidx];

      for (int rr = 0; rr < n_runs; ++rr) {
        std::vector<int> infector; infector.reserve(max_cases);
        std::vector<int> generation; generation.reserve(max_cases);
        std::vector<double> infection_time; infection_time.reserve(max_cases);
        std::vector<double> symptom_onset; symptom_onset.reserve(max_cases);
        std::vector<int> detected; detected.reserve(max_cases);
        std::vector<double> detection_time; detection_time.reserve(max_cases);
        std::vector<double> isolation_time; isolation_time.reserve(max_cases);
        std::vector<int> processed; processed.reserve(max_cases);
        std::vector<PotentialTx> potential_tx; potential_tx.reserve(max_cases * 2);

        auto add_case = [&](int inf, int gen, double itime) {
          double inc = sample_from_vec(incubation_draws);
          double onset = itime + inc;
          int det = (runif01() < p_detect) ? 1 : 0;
          double dtime = R_PosInf;
          double isotime = R_PosInf;
          if (det == 1) {
            dtime = onset + positive_normal_one(detection_delay_mean, detection_delay_sd, 0.0);
            isotime = dtime + positive_normal_one(isolation_delay_mean, isolation_delay_sd, 0.0);
          }
          infector.push_back(inf);
          generation.push_back(gen);
          infection_time.push_back(itime);
          symptom_onset.push_back(onset);
          detected.push_back(det);
          detection_time.push_back(dtime);
          isolation_time.push_back(isotime);
          processed.push_back(0);
          return (int)infection_time.size() - 1;
        };

        for (int ii = 0; ii < initial_infected; ++ii) add_case(NA_INTEGER, 0, 0.0);

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
            double G = sample_from_vec(generation_draws);
            double child_time = infection_time[todo] + G;
            if (!finite_num(child_time) || child_time > max_days) continue;
            int child_idx = add_case(todo, generation[todo] + 1, child_time);
            PotentialTx tx;
            tx.source = todo;
            tx.child = child_idx;
            tx.time = child_time;
            potential_tx.push_back(tx);
          }
        }

        int n_potential = (int)infection_time.size();

        bool reached_start = false;
        bool reached_current = false;
        double response_time = R_PosInf;

        for (int ii = 0; ii < n_potential; ++ii) {
          if (generation[ii] == start_generation) {
            reached_start = true;
            if (infection_time[ii] < response_time) response_time = infection_time[ii];
          }
          if (generation[ii] == current_generation) reached_current = true;
        }

        std::vector<int> retained(n_potential, 0);
        std::vector<int> quarantined(n_potential, 0);
        std::vector<double> q_start(n_potential, R_PosInf);
        std::vector<double> q_end(n_potential, R_PosInf);
        std::vector<int> retained_offspring(n_potential, 0);

        for (int ii = 0; ii < initial_infected && ii < n_potential; ++ii) retained[ii] = 1;

        int additional_quarantine = 0;

        // At response start, attempt to trace/quarantine infected people who already exist.
        // q_system is tracing completeness: each already-infected person is captured with probability q.
        if (reached_start) {
          for (int ii = 0; ii < n_potential; ++ii) {
            if (infection_time[ii] <= response_time && runif01() < q_system) {
              quarantined[ii] = 1;
              q_start[ii] = response_time;
              q_end[ii] = response_time + quarantine_duration;
              additional_quarantine++;
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

          bool after_response = reached_start && t >= response_time;
          bool source_quarantined = after_response && quarantined[s] == 1 && t >= q_start[s] && t <= q_end[s];
          bool source_isolated = after_response && detected[s] == 1 && finite_num(isolation_time[s]) && t >= isolation_time[s];

          double retain_prob = 1.0;
          if (after_response) {
            if (source_quarantined && runif01() < quarantine_adherence) retain_prob *= (1.0 - quarantine_effectiveness);
            if (source_isolated) retain_prob *= (1.0 - isolation_effectiveness);
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
            // Simple forward tracing: newly infected retained cases are detected/traced with probability q
            // and then enter quarantine from their detection time.
            if (runif01() < q_system) {
              double trace_time = detection_time[c];
              if (!finite_num(trace_time) || trace_time < t) trace_time = t;
              quarantined[c] = 1;
              q_start[c] = trace_time;
              q_end[c] = trace_time + quarantine_duration;
            }
          }
        }

        int total_retained = 0;
        int additional_after_current = 0;
        int post_response_infections = 0;
        double max_end = NA_REAL;

        for (int ii = 0; ii < n_potential; ++ii) {
          if (retained[ii] == 0) continue;
          total_retained++;
          if (generation[ii] > current_generation) additional_after_current++;
          if (reached_start && infection_time[ii] >= response_time) post_response_infections++;

          if (reached_start && (infection_time[ii] >= response_time || generation[ii] > current_generation)) {
            if (!finite_num(max_end) || infection_time[ii] > max_end) max_end = infection_time[ii];
            if (finite_num(symptom_onset[ii]) && (!finite_num(max_end) || symptom_onset[ii] > max_end)) max_end = symptom_onset[ii];
            if (finite_num(detection_time[ii]) && (!finite_num(max_end) || detection_time[ii] > max_end)) max_end = detection_time[ii];
          }
          if (finite_num(q_end[ii]) && reached_start && (!finite_num(max_end) || q_end[ii] > max_end)) max_end = q_end[ii];
        }

        // Additional quarantine burden from uninfected contacts around detected post-response cases.
        // This is a stochastic proxy for contact-tracing burden and is what drives Panel C.
        if (reached_start && mean_uninfected_contacts_per_detected_case > 0.0) {
          for (int jj = 0; jj < n_potential; ++jj) {
            if (retained[jj] == 0) continue;
            if (detected[jj] != 1 || !finite_num(detection_time[jj]) || detection_time[jj] < response_time) continue;

            int uninfected_contacts = 0;
            if (uninfected_contact_dispersion > 0.0) {
              double prob_contacts = uninfected_contact_dispersion / (uninfected_contact_dispersion + mean_uninfected_contacts_per_detected_case);
              uninfected_contacts = (int)R::rnbinom(uninfected_contact_dispersion, prob_contacts);
            } else {
              uninfected_contacts = (int)R::rpois(mean_uninfected_contacts_per_detected_case);
            }

            if (uninfected_contacts > 0) {
              additional_quarantine += (int)R::rbinom((double)uninfected_contacts, q_system);
            }
          }
        }

        double response_duration = NA_REAL;
        if (reached_start) {
          if (finite_num(max_end)) response_duration = std::max(0.0, max_end - response_time);
          else response_duration = 0.0;
        }

        out_run_id[row] = row + 1;
        out_q[row] = q_system;
        out_contacts_missed[row] = 1.0 - q_system;
        out_start_generation[row] = start_generation;
        out_current_generation[row] = current_generation;
        out_reached_start_generation[row] = reached_start ? 1 : 0;
        out_reached_current_generation[row] = reached_current ? 1 : 0;
        out_response_time[row] = reached_start ? response_time : NA_REAL;
        out_total_potential_infections[row] = n_potential;
        out_total_retained_infections[row] = total_retained;
        out_additional_cases_after_current_gen[row] = additional_after_current;
        out_any_additional_cases_after_current_gen[row] = (additional_after_current > 0) ? 1 : 0;
        out_post_response_infections[row] = post_response_infections;
        out_post_response_transmissions[row] = post_response_transmissions;
        out_additional_quarantine[row] = additional_quarantine;
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
    Named("start_generation") = out_start_generation,
    Named("current_generation") = out_current_generation,
    Named("reached_start_generation") = out_reached_start_generation,
    Named("reached_current_generation") = out_reached_current_generation,
    Named("response_time") = out_response_time,
    Named("total_potential_infections") = out_total_potential_infections,
    Named("total_retained_infections") = out_total_retained_infections,
    Named("additional_cases_after_current_gen") = out_additional_cases_after_current_gen,
    Named("any_additional_cases_after_current_gen") = out_any_additional_cases_after_current_gen,
    Named("post_response_infections") = out_post_response_infections,
    Named("post_response_transmissions") = out_post_response_transmissions,
    Named("additional_quarantine") = out_additional_quarantine,
    Named("prevented_by_quarantine_total") = out_prevented_by_quarantine_total,
    Named("prevented_by_isolation_total") = out_prevented_by_isolation_total,
    Named("response_duration_after_response") = out_response_duration_after_response
  );
}
')

# -----------------------------
# Run simulations
# -----------------------------
fixed_inputs <- tibble(
  parameter = c(
    "input", "delay_file_used", "out", "n_runs", "q_values", "start_generations", "seed",
    "current_generation", "max_days", "max_cases", "initial_infected",
    "R0", "k", "p_detect", "detection_delay_mean", "detection_delay_sd",
    "isolation_delay_mean", "isolation_delay_sd", "isolation_effectiveness",
    "quarantine_duration", "quarantine_effectiveness", "quarantine_adherence",
    "mean_uninfected_contacts_per_detected_case", "uninfected_contact_dispersion"
  ),
  value = c(
    args$input, delay_file, args$out, args$n_runs,
    paste(args$q_values, collapse = ","),
    paste(args$start_generations, collapse = ","),
    args$seed,
    args$current_generation, args$max_days, args$max_cases, args$initial_infected,
    args$R0, args$k, args$p_detect, args$detection_delay_mean, args$detection_delay_sd,
    args$isolation_delay_mean, args$isolation_delay_sd, args$isolation_effectiveness,
    args$quarantine_duration, args$quarantine_effectiveness, args$quarantine_adherence,
    args$mean_uninfected_contacts_per_detected_case, args$uninfected_contact_dispersion
  )
)
readr::write_csv(fixed_inputs, file.path(args$out, "tables", "generation_start_fixed_inputs.csv"))

message("Running generation-start tracing/quarantine sensitivity analysis...")
t0 <- Sys.time()
runs <- run_generation_start_sensitivity_cpp(
  q_values = args$q_values,
  start_generations = args$start_generations,
  n_runs = args$n_runs,
  generation_draws = generation_draws,
  incubation_draws = incubation_draws,
  current_generation = args$current_generation,
  max_cases = args$max_cases,
  max_days = args$max_days,
  initial_infected = args$initial_infected,
  R0 = args$R0,
  k = args$k,
  p_detect = args$p_detect,
  detection_delay_mean = args$detection_delay_mean,
  detection_delay_sd = args$detection_delay_sd,
  isolation_delay_mean = args$isolation_delay_mean,
  isolation_delay_sd = args$isolation_delay_sd,
  isolation_effectiveness = args$isolation_effectiveness,
  quarantine_duration = args$quarantine_duration,
  quarantine_effectiveness = args$quarantine_effectiveness,
  quarantine_adherence = args$quarantine_adherence,
  mean_uninfected_contacts_per_detected_case = args$mean_uninfected_contacts_per_detected_case,
  uninfected_contact_dispersion = args$uninfected_contact_dispersion,
  seed = args$seed
) %>%
  as_tibble() %>%
  mutate(
    q_label = factor(make_q_label(q_system), levels = make_q_label(sort(unique(q_system), decreasing = TRUE))),
    response_duration_weeks = response_duration_after_response / 7
  )

message("Simulation time: ", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 2), " seconds")

readr::write_csv(runs, file.path(args$out, "tables", "generation_start_run_level_results.csv"))

runs_conditioned <- runs %>%
  filter(reached_start_generation == 1, reached_current_generation == 1)

if (nrow(runs_conditioned) == 0) {
  stop("No simulations reached both the start generation and current_generation. Increase n_runs or check parameters.", call. = FALSE)
}

summary <- runs_conditioned %>%
  group_by(q_system, q_label, start_generation, current_generation) %>%
  summarise(
    n_runs = n(),
    probability_any_additional_cases_after_current_generation = mean(any_additional_cases_after_current_gen),
    q95_total_infections = quant(total_retained_infections, 0.95),
    q95_additional_quarantine = quant(additional_quarantine, 0.95),
    q95_response_duration_weeks = quant(response_duration_weeks, 0.95),
    median_total_infections = median(total_retained_infections),
    median_additional_quarantine = median(additional_quarantine),
    median_response_duration_weeks = median(response_duration_weeks, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(start_generation, desc(q_system))

readr::write_csv(summary, file.path(args$out, "tables", "generation_start_summary.csv"))

# -----------------------------
# Figure
# -----------------------------
if (isTRUE(args$make_plots)) {
  fig_dir <- file.path(args$out, "figures")
  
  theme_pub <- theme_classic(base_size = 12) +
    theme(
      axis.text = element_text(colour = "black"),
      axis.title.x = element_text(margin = margin(t = 8)),
      axis.title.y = element_text(margin = margin(r = 8)),
      panel.grid.major.y = element_line(colour = "grey88", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      plot.tag = element_text(face = "bold", size = 13),
      plot.margin = margin(10, 12, 10, 12)
    )
  
  pA <- ggplot(summary, aes(x = start_generation, y = probability_any_additional_cases_after_current_generation,
                            colour = q_label, linetype = q_label, group = q_label)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2.3) +
    scale_x_continuous(breaks = sort(unique(summary$start_generation))) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    labs(
      x = "Generation when tracing/quarantine begins",
      y = paste0("Probability of any additional cases\nafter generation ", args$current_generation),
      colour = "Tracing completeness",
      linetype = "Tracing completeness",
      tag = "A"
    ) +
    theme_pub
  
  plot_runs <- runs_conditioned %>%
    mutate(
      start_generation_f = factor(start_generation, levels = sort(unique(start_generation))),
      q_label = factor(q_label, levels = make_q_label(sort(unique(q_system), decreasing = TRUE)))
    )
  
  pB <- ggplot(
    plot_runs,
    aes(x = start_generation_f, y = total_retained_infections, fill = q_label)
  ) +
    geom_boxplot(
      position = position_dodge(width = 0.75),
      width = 0.65,
      outlier.shape = NA,
      colour = "grey25"
    ) +
    scale_y_continuous(
      limits = c(0, 1000),
      labels = label_number(big.mark = ",")) +
    labs(
      x = "Generation when tracing/quarantine begins",
      y = "Total infections",
      fill = "Tracing completeness",
      tag = "B"
    ) +
    theme_pub
  
  pC <- ggplot(
    plot_runs,
    aes(x = start_generation_f, y = additional_quarantine, fill = q_label)
  ) +
    geom_boxplot(
      position = position_dodge(width = 0.75),
      width = 0.65,
      outlier.shape = NA,
      colour = "grey25"
    ) +
    scale_y_continuous(
      limits = c(0, 15000),
      labels = label_number(big.mark = ",")) +
    labs(
      x = "Generation when tracing/quarantine begins",
      y = "Additional quarantine",
      fill = "Tracing completeness",
      tag = "C"
    ) +
    theme_pub
  
  plot_runs$response_duration_weeks2 <- max(0, plot_runs$response_duration_weeks)
  plot_runs$response_duration_weeks <- pmax(plot_runs$response_duration_weeks - 6, 0)
  
  pD <- ggplot(
    plot_runs,
    aes(x = start_generation_f, y = response_duration_weeks, fill = q_label)
  ) +
    geom_boxplot(
      position = position_dodge(width = 0.75),
      width = 0.65,
      outlier.shape = NA,
      colour = "grey25"
    ) +
    scale_y_continuous(labels = label_number(big.mark = ",")) +
    labs(
      x = "Generation when tracing/quarantine begins",
      y = "Additional response duration (weeks)",
      fill = "Tracing completeness",
      tag = "D"
    ) +
    theme_pub
  
  fig <- (pA | pB) / (pC | pD) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  
  ggsave(file.path(fig_dir, "figure_generation_start_tracing_sensitivity.png"), fig, width = 11, height = 8.5, dpi = 600, bg = "white")
  ggsave(file.path(fig_dir, "figure_generation_start_tracing_sensitivity.pdf"), fig, width = 11, height = 8.5, bg = "white")
}

message("Done.")
message("Outputs written to: ", normalizePath(args$out, winslash = "/", mustWork = FALSE))
message("Summary: ", normalizePath(file.path(args$out, "tables", "generation_start_summary.csv"), winslash = "/", mustWork = FALSE))

# Appendix figure 2
# May need to change this wd
setwd("~/Hantavirus-paper/BP results/Fig 6")

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(scales)
})

cmd <- commandArgs(trailingOnly = TRUE)

out <- getwd()

run_file <- file.path(out, "tables", "generation_start_run_level_results.csv")
fig_dir <- file.path(out, "figures")

if (!file.exists(run_file)) {
  stop(
    "Could not find run-level results file:\n  ",
    normalizePath(run_file, winslash = "/", mustWork = FALSE),
    "\nRun the generation-start simulation first, or pass the correct output directory.",
    call. = FALSE
  )
}

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

make_q_label <- function(q) paste0(round(100 * q), "%")

runs <- readr::read_csv(run_file, show_col_types = FALSE)

required_cols <- c(
  "q_system",
  "start_generation",
  "reached_start_generation",
  "reached_current_generation",
  "total_retained_infections"
)

missing_cols <- setdiff(required_cols, names(runs))
if (length(missing_cols) > 0) {
  stop(
    "The run-level file is missing required columns: ",
    paste(missing_cols, collapse = ", "),
    call. = FALSE
  )
}

plot_runs <- runs %>%
  filter(
    reached_start_generation == 1,
    reached_current_generation == 1,
    is.finite(total_retained_infections),
    total_retained_infections > 0
  ) %>%
  mutate(
    q_label = make_q_label(q_system),
    start_generation = as.integer(start_generation)
  )

if (nrow(plot_runs) == 0) {
  stop(
    "No usable rows after filtering to reached_start_generation == 1, ",
    "reached_current_generation == 1, and total_retained_infections > 0.",
    call. = FALSE
  )
}

q_levels <- make_q_label(sort(unique(plot_runs$q_system), decreasing = TRUE))
start_levels <- sort(unique(plot_runs$start_generation))

plot_runs <- plot_runs %>%
  mutate(
    q_label = factor(q_label, levels = q_levels),
    panel_label = factor(
      paste0("Tracing starts at gen ", start_generation, "\n", q_label, " completeness"),
      levels = as.vector(t(outer(
        paste0("Tracing starts at gen ", start_levels),
        paste0(q_levels, " completeness"),
        paste,
        sep = "\n"
      )))
    )
  )

max_x <- max(plot_runs$total_retained_infections, na.rm = TRUE)

x_breaks <- c(2, 5, 10, 25, 50, 100, 250, 500, 1000, 2500, 5000)
x_breaks <- x_breaks[x_breaks <= max(max_x, 2)]

log_breaks <- exp(seq(log(1), log(max_x), length.out = 90))

p_density <- ggplot(plot_runs, aes(x = total_retained_infections)) +
  geom_histogram(
    aes(y = after_stat(density)),
    breaks = log_breaks,
    fill = "grey82",
    colour = "grey82",
    linewidth = 0.15
  ) +
  geom_density(
    linewidth = 0.9,
    colour = "black",
    adjust = 1.15
  ) +
  facet_wrap(
    ~ panel_label,
    ncol = length(q_levels),
    scales = "free_y"
  ) +
  scale_x_log10(
    breaks = x_breaks,
    labels = label_number(big.mark = ",")
  ) +
  labs(
    x = "Total infections",
    y = "Density"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text = element_text(colour = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_text(size = 15, margin = margin(t = 10)),
    axis.title.y = element_text(size = 15, margin = margin(r = 10)),
    strip.background = element_rect(fill = "grey88", colour = NA),
    strip.text = element_text(colour = "grey20", size = 9.5, face = "bold"),
    panel.spacing.x = unit(0.55, "lines"),
    panel.spacing.y = unit(0.75, "lines"),
    plot.margin = margin(10, 14, 10, 14)
  )

png_file <- file.path(fig_dir, "figure_total_infections_density_by_start_generation_and_q_CLEAN.png")
pdf_file <- file.path(fig_dir, "figure_total_infections_density_by_start_generation_and_q_CLEAN.pdf")

ggsave(
  png_file,
  p_density,
  width = 13,
  height = 10,
  dpi = 600,
  bg = "white"
)

ggsave(
  pdf_file,
  p_density,
  width = 13,
  height = 10,
  bg = "white"
)

message("Saved clean density figure:")
message("  ", normalizePath(png_file, winslash = "/", mustWork = FALSE))
message("  ", normalizePath(pdf_file, winslash = "/", mustWork = FALSE))

# Print key quantitative summaries for outbreak-size density text
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
})

cmd <- commandArgs(trailingOnly = TRUE)

run_file <- file.path(out, "tables", "generation_start_run_level_results.csv")

if (!file.exists(run_file)) {
  stop(
    "Could not find run-level results file:\n  ",
    normalizePath(run_file, winslash = "/", mustWork = FALSE),
    "\nRun the generation-start simulation first, or pass the correct output directory.",
    call. = FALSE
  )
}

runs <- readr::read_csv(run_file, show_col_types = FALSE)

required_cols <- c(
  "q_system",
  "start_generation",
  "reached_start_generation",
  "reached_current_generation",
  "total_retained_infections",
  "additional_cases_after_current_gen",
  "any_additional_cases_after_current_gen",
  "additional_quarantine",
  "response_duration_after_response"
)

missing_cols <- setdiff(required_cols, names(runs))
if (length(missing_cols) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_cols, collapse = ", "),
    call. = FALSE
  )
}

fmt_pct <- function(x, digits = 1) {
  paste0(sprintf(paste0("%.", digits, "f"), 100 * x), "\\%")
}

fmt_num <- function(x, digits = 0) {
  format(round(x, digits), big.mark = ",", scientific = FALSE, trim = TRUE)
}

fmt_dec <- function(x, digits = 1) {
  sprintf(paste0("%.", digits, "f"), x)
}

qfun <- function(x, p) {
  as.numeric(stats::quantile(x, p, na.rm = TRUE, names = FALSE))
}

make_q_label <- function(q) paste0(round(100 * q), "\\%")

dat <- runs %>%
  filter(
    reached_start_generation == 1,
    reached_current_generation == 1,
    is.finite(total_retained_infections),
    is.finite(additional_cases_after_current_gen)
  ) %>%
  mutate(
    q_label = make_q_label(q_system),
    response_duration_weeks = response_duration_after_response / 7
  )

if (nrow(dat) == 0) {
  stop(
    "No usable rows after filtering to reached_start_generation == 1 and reached_current_generation == 1.",
    call. = FALSE
  )
}

summary_by_scenario <- dat %>%
  group_by(start_generation, q_system, q_label) %>%
  summarise(
    n = n(),
    
    total_median = median(total_retained_infections, na.rm = TRUE),
    total_q75 = qfun(total_retained_infections, 0.75),
    total_q90 = qfun(total_retained_infections, 0.90),
    total_q95 = qfun(total_retained_infections, 0.95),
    total_q99 = qfun(total_retained_infections, 0.99),
    total_max = max(total_retained_infections, na.rm = TRUE),
    
    add_median = median(additional_cases_after_current_gen, na.rm = TRUE),
    add_q75 = qfun(additional_cases_after_current_gen, 0.75),
    add_q90 = qfun(additional_cases_after_current_gen, 0.90),
    add_q95 = qfun(additional_cases_after_current_gen, 0.95),
    add_q99 = qfun(additional_cases_after_current_gen, 0.99),
    add_max = max(additional_cases_after_current_gen, na.rm = TRUE),
    
    p_no_additional = mean(additional_cases_after_current_gen == 0, na.rm = TRUE),
    p_any_additional = mean(additional_cases_after_current_gen > 0, na.rm = TRUE),
    p_add_ge_5 = mean(additional_cases_after_current_gen >= 5, na.rm = TRUE),
    p_add_ge_10 = mean(additional_cases_after_current_gen >= 10, na.rm = TRUE),
    p_add_ge_25 = mean(additional_cases_after_current_gen >= 25, na.rm = TRUE),
    p_add_ge_50 = mean(additional_cases_after_current_gen >= 50, na.rm = TRUE),
    p_add_ge_100 = mean(additional_cases_after_current_gen >= 100, na.rm = TRUE),
    
    p_total_ge_25 = mean(total_retained_infections >= 25, na.rm = TRUE),
    p_total_ge_50 = mean(total_retained_infections >= 50, na.rm = TRUE),
    p_total_ge_100 = mean(total_retained_infections >= 100, na.rm = TRUE),
    p_total_ge_250 = mean(total_retained_infections >= 250, na.rm = TRUE),
    p_total_ge_500 = mean(total_retained_infections >= 500, na.rm = TRUE),
    p_total_ge_1000 = mean(total_retained_infections >= 1000, na.rm = TRUE),
    
    quarantine_median = median(additional_quarantine, na.rm = TRUE),
    quarantine_q95 = qfun(additional_quarantine, 0.95),
    
    duration_median_weeks = median(response_duration_weeks, na.rm = TRUE),
    duration_q95_weeks = qfun(response_duration_weeks, 0.95),
    
    .groups = "drop"
  ) %>%
  arrange(start_generation, desc(q_system))

cat("\n")
cat("============================================================\n")
cat("OUTBREAK-SIZE DENSITY TEXT QUANTITIES\n")
cat("============================================================\n")
cat("Run file: ", normalizePath(run_file, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("Conditioning: reached_start_generation == 1 and reached_current_generation == 1\n")
cat("============================================================\n\n")

cat("SCENARIO SUMMARY TABLE\n")
cat("------------------------------------------------------------\n")
print(
  summary_by_scenario %>%
    transmute(
      start_generation,
      q_label,
      n,
      total_median = fmt_num(total_median),
      total_q95 = fmt_num(total_q95),
      total_q99 = fmt_num(total_q99),
      add_median = fmt_num(add_median),
      add_q95 = fmt_num(add_q95),
      p_no_additional = fmt_pct(p_no_additional),
      p_add_ge_10 = fmt_pct(p_add_ge_10),
      p_add_ge_50 = fmt_pct(p_add_ge_50),
      p_total_ge_100 = fmt_pct(p_total_ge_100),
      p_total_ge_500 = fmt_pct(p_total_ge_500),
      p_total_ge_1000 = fmt_pct(p_total_ge_1000)
    ),
  n = Inf
)

cat("\n\n")
cat("THRESHOLD PROBABILITIES FOR ADDITIONAL INFECTIONS AFTER GENERATION 2\n")
cat("------------------------------------------------------------\n")
print(
  summary_by_scenario %>%
    transmute(
      start_generation,
      q_label,
      no_additional = fmt_pct(p_no_additional),
      ge_5 = fmt_pct(p_add_ge_5),
      ge_10 = fmt_pct(p_add_ge_10),
      ge_25 = fmt_pct(p_add_ge_25),
      ge_50 = fmt_pct(p_add_ge_50),
      ge_100 = fmt_pct(p_add_ge_100)
    ),
  n = Inf
)

cat("\n\n")
cat("TOTAL-INFECTION UPPER-TAIL PROBABILITIES\n")
cat("------------------------------------------------------------\n")
print(
  summary_by_scenario %>%
    transmute(
      start_generation,
      q_label,
      ge_25 = fmt_pct(p_total_ge_25),
      ge_50 = fmt_pct(p_total_ge_50),
      ge_100 = fmt_pct(p_total_ge_100),
      ge_250 = fmt_pct(p_total_ge_250),
      ge_500 = fmt_pct(p_total_ge_500),
      ge_1000 = fmt_pct(p_total_ge_1000)
    ),
  n = Inf
)

# Strong text contrasts
pick_row <- function(g, q) {
  summary_by_scenario %>%
    filter(start_generation == g, abs(q_system - q) < 1e-8) %>%
    slice(1)
}

available_q <- sort(unique(summary_by_scenario$q_system), decreasing = TRUE)
q_ref <- if (any(abs(available_q - 0.95) < 1e-8)) 0.95 else available_q[1]

early0 <- pick_row(0, q_ref)
early1 <- pick_row(1, q_ref)
delayed2 <- pick_row(2, q_ref)
delayed3 <- pick_row(3, q_ref)
delayed4 <- pick_row(4, q_ref)

cat("\n\n")
cat("READY-TO-COPY QUANTITATIVE SENTENCES\n")
cat("------------------------------------------------------------\n\n")

cat(
  "Using ",
  make_q_label(q_ref),
  " tracing completeness as the reference scenario, the probability of at least 50 additional infections after generation 2 was ",
  if (nrow(early0) == 1) fmt_pct(early0$p_add_ge_50) else "NA",
  " when tracing began at generation 0",
  if (nrow(early1) == 1) paste0(", ", fmt_pct(early1$p_add_ge_50), " when tracing began at generation 1") else "",
  if (nrow(delayed2) == 1) paste0(", ", fmt_pct(delayed2$p_add_ge_50), " when tracing began at generation 2") else "",
  if (nrow(delayed3) == 1) paste0(", ", fmt_pct(delayed3$p_add_ge_50), " when tracing began at generation 3") else "",
  if (nrow(delayed4) == 1) paste0(", and ", fmt_pct(delayed4$p_add_ge_50), " when tracing began at generation 4") else "",
  ".\n\n",
  sep = ""
)

cat(
  "At ",
  make_q_label(q_ref),
  " tracing completeness, the 95th percentile of total infections was ",
  if (nrow(early0) == 1) fmt_num(early0$total_q95) else "NA",
  " when tracing began at generation 0",
  if (nrow(early1) == 1) paste0(", ", fmt_num(early1$total_q95), " when tracing began at generation 1") else "",
  if (nrow(delayed2) == 1) paste0(", ", fmt_num(delayed2$total_q95), " when tracing began at generation 2") else "",
  if (nrow(delayed3) == 1) paste0(", ", fmt_num(delayed3$total_q95), " when tracing began at generation 3") else "",
  if (nrow(delayed4) == 1) paste0(", and ", fmt_num(delayed4$total_q95), " when tracing began at generation 4") else "",
  ".\n\n",
  sep = ""
)

# Largest increase from generation 0 to generation 4 across q values
wide_g0_g4 <- summary_by_scenario %>%
  filter(start_generation %in% c(0, 4)) %>%
  select(q_system, q_label, start_generation, p_add_ge_50, total_q95, p_total_ge_500, p_total_ge_1000) %>%
  pivot_wider(
    names_from = start_generation,
    values_from = c(p_add_ge_50, total_q95, p_total_ge_500, p_total_ge_1000),
    names_prefix = "g"
  ) %>%
  mutate(
    delta_p_add_ge_50 = p_add_ge_50_g4 - p_add_ge_50_g0,
    ratio_total_q95 = total_q95_g4 / total_q95_g0
  )

if (nrow(wide_g0_g4) > 0) {
  strongest <- wide_g0_g4 %>%
    arrange(desc(delta_p_add_ge_50)) %>%
    slice(1)
  
  cat(
    "The largest increase in the probability of at least 50 additional infections between generation-0 and generation-4 response occurred at ",
    strongest$q_label,
    " tracing completeness, rising from ",
    fmt_pct(strongest$p_add_ge_50_g0),
    " to ",
    fmt_pct(strongest$p_add_ge_50_g4),
    ". The corresponding 95th percentile of total infections changed from ",
    fmt_num(strongest$total_q95_g0),
    " to ",
    fmt_num(strongest$total_q95_g4),
    ".\n\n",
    sep = ""
  )
}

# Compare q deterioration at generation 4
g4_best <- summary_by_scenario %>%
  filter(start_generation == max(start_generation, na.rm = TRUE)) %>%
  arrange(desc(q_system)) %>%
  slice(1)

g4_worst <- summary_by_scenario %>%
  filter(start_generation == max(start_generation, na.rm = TRUE)) %>%
  arrange(q_system) %>%
  slice(1)

if (nrow(g4_best) == 1 && nrow(g4_worst) == 1) {
  cat(
    "For the latest-start scenario, generation ",
    g4_best$start_generation,
    ", lowering tracing completeness from ",
    g4_best$q_label,
    " to ",
    g4_worst$q_label,
    " increased the probability of at least 50 additional infections from ",
    fmt_pct(g4_best$p_add_ge_50),
    " to ",
    fmt_pct(g4_worst$p_add_ge_50),
    ", and increased the 95th percentile of total infections from ",
    fmt_num(g4_best$total_q95),
    " to ",
    fmt_num(g4_worst$total_q95),
    ".\n\n",
    sep = ""
  )
}

cat("DONE\n")