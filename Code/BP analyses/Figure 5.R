#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(scales)
  library(grid)
})

# ============================================================
# Figure 5:
# Threshold probabilities for the HIGH-TRANSMISSION scenario only
#
# Assumptions:
# - run_level_results.csv is in the current working directory
# - outputs are saved into ./figures and ./tables
# ============================================================

out <- "~/Hantavirus-paper/BP results/Fig 2/tables"

run_file <- file.path(out, "run_level_results.csv")
fig_dir <- file.path("~/Hantavirus-paper/BP results/Fig 5", "figures")
table_dir <- file.path("~/Hantavirus-paper/BP results/Fig 5", "tables")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

message(
  "Reading run-level results from: ",
  normalizePath(run_file, winslash = "/", mustWork = FALSE)
)

if (!file.exists(run_file)) {
  stop(
    "Run-level results file not found: ",
    normalizePath(run_file, winslash = "/", mustWork = FALSE),
    call. = FALSE
  )
}

runs <- readr::read_csv(
  run_file,
  show_col_types = FALSE
)

# ------------------------------------------------------------
# Validate required columns
# ------------------------------------------------------------

required_columns <- c(
  "transmission_regime",
  "scenario",
  "reached_current_generation",
  "q_system",
  "additional_infections_after_current_generation"
)

missing_columns <- setdiff(required_columns, names(runs))

if (length(missing_columns) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_columns, collapse = ", "),
    "\n\nAvailable columns are:\n",
    paste(names(runs), collapse = ", "),
    call. = FALSE
  )
}

make_q_label <- function(q) {
  paste0(round(100 * q), "%")
}

# ------------------------------------------------------------
# Restrict Figure 5 to:
# - high-transmission scenario
# - imperfect intervention
# - simulations reaching generation 2
# ------------------------------------------------------------

plot_runs <- runs %>%
  filter(
    transmission_regime == "supercritical",
    scenario == "imperfect_intervention",
    reached_current_generation == 1
  ) %>%
  mutate(
    q_system = as.numeric(q_system),
    q_label = factor(
      make_q_label(q_system),
      levels = make_q_label(
        sort(unique(q_system), decreasing = TRUE)
      )
    )
  )

if (nrow(plot_runs) == 0) {
  stop(
    paste0(
      "No rows remained after filtering.\n\n",
      "Expected:\n",
      '  transmission_regime == "supercritical"\n',
      '  scenario == "imperfect_intervention"\n',
      "  reached_current_generation == 1\n\n",
      "Observed transmission_regime values:\n  ",
      paste(sort(unique(runs$transmission_regime)), collapse = ", "),
      "\n\nObserved scenario values:\n  ",
      paste(sort(unique(runs$scenario)), collapse = ", ")
    ),
    call. = FALSE
  )
}

if (n_distinct(plot_runs$transmission_regime) != 1 ||
    unique(plot_runs$transmission_regime) != "supercritical") {
  stop(
    "Figure 5 contains an unintended transmission regime.",
    call. = FALSE
  )
}

if (n_distinct(plot_runs$scenario) != 1 ||
    unique(plot_runs$scenario) != "imperfect_intervention") {
  stop(
    "Figure 5 contains an unintended intervention scenario.",
    call. = FALSE
  )
}

message("Figure 5 input counts:")

print(
  plot_runs %>%
    count(
      transmission_regime,
      scenario,
      q_system,
      name = "n_conditioned_runs"
    ) %>%
    arrange(desc(q_system))
)

# ------------------------------------------------------------
# Calculate threshold probabilities
# ------------------------------------------------------------

threshold_plot_data <- plot_runs %>%
  mutate(
    `No additional cases` =
      additional_infections_after_current_generation == 0,
    
    `≥1 additional case` =
      additional_infections_after_current_generation >= 1,
    
    `≥5 additional cases` =
      additional_infections_after_current_generation >= 5,
    
    `≥10 additional cases` =
      additional_infections_after_current_generation >= 10,
    
    `≥25 additional cases` =
      additional_infections_after_current_generation >= 25,
    
    `≥50 additional cases` =
      additional_infections_after_current_generation >= 50
  ) %>%
  group_by(q_system, q_label) %>%
  summarise(
    `No additional cases` =
      mean(`No additional cases`, na.rm = TRUE),
    
    `≥1 additional case` =
      mean(`≥1 additional case`, na.rm = TRUE),
    
    `≥5 additional cases` =
      mean(`≥5 additional cases`, na.rm = TRUE),
    
    `≥10 additional cases` =
      mean(`≥10 additional cases`, na.rm = TRUE),
    
    `≥25 additional cases` =
      mean(`≥25 additional cases`, na.rm = TRUE),
    
    `≥50 additional cases` =
      mean(`≥50 additional cases`, na.rm = TRUE),
    
    n_conditioned_runs = n(),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = c(
      `No additional cases`,
      `≥1 additional case`,
      `≥5 additional cases`,
      `≥10 additional cases`,
      `≥25 additional cases`,
      `≥50 additional cases`
    ),
    names_to = "threshold",
    values_to = "probability"
  ) %>%
  mutate(
    threshold = factor(
      threshold,
      levels = c(
        "No additional cases",
        "≥1 additional case",
        "≥5 additional cases",
        "≥10 additional cases",
        "≥25 additional cases",
        "≥50 additional cases"
      )
    ),
    label = ifelse(
      is.na(probability),
      "NA",
      sub("\\.?0+$", "", sprintf("%.2f", probability))
    )
  ) %>%
  arrange(threshold, desc(q_system))

# ------------------------------------------------------------
# Save the calculated probabilities
# ------------------------------------------------------------

table_file <- file.path(
  table_dir,
  "threshold_probabilities_high_transmission.csv"
)

readr::write_csv(
  threshold_plot_data,
  table_file
)

message(
  "Saved threshold table to: ",
  normalizePath(table_file, winslash = "/", mustWork = FALSE)
)

# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------

p_threshold <- ggplot(
  threshold_plot_data,
  aes(x = q_label, y = probability)
) +
  geom_col(
    width = 0.68,
    fill = "grey75",
    colour = "grey25",
    linewidth = 0.35
  ) +
  geom_text(
    aes(label = label),
    vjust = -0.25,
    size = 3.2,
    na.rm = TRUE
  ) +
  facet_wrap(
    ~threshold,
    ncol = 3
  ) +
  scale_y_continuous(
    limits = c(0, 1.08),
    breaks = seq(0, 1, by = 0.2),
    labels = scales::number_format(accuracy = 0.1),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    x = "Tracing completeness",
    y = "Probability"
  ) +
  theme_classic(base_size = 12) +
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
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold",
      size = 11
    ),
    panel.grid.major.y = element_line(
      colour = "grey88",
      linewidth = 0.3
    ),
    panel.grid.minor = element_blank(),
    panel.spacing = grid::unit(0.6, "lines"),
    plot.margin = margin(10, 12, 10, 12)
  )

# ------------------------------------------------------------
# Save PNG and PDF
# ------------------------------------------------------------

png_file <- file.path(
  fig_dir,
  "Figure5_threshold_probabilities_high_transmission.png"
)

pdf_file <- file.path(
  fig_dir,
  "Figure5_threshold_probabilities_high_transmission.pdf"
)

ggsave(
  filename = png_file,
  plot = p_threshold,
  width = 10.5,
  height = 6.4,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = pdf_file,
  plot = p_threshold,
  width = 10.5,
  height = 6.4,
  units = "in",
  bg = "white"
)

message(
  "Saved PNG to: ",
  normalizePath(png_file, winslash = "/", mustWork = FALSE)
)

message(
  "Saved PDF to: ",
  normalizePath(pdf_file, winslash = "/", mustWork = FALSE)
)

message("Figure 5 completed successfully.")