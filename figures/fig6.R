library(tibble)
library(dplyr)
library(forcats)
library(funkyheatmap)
library(ggplot2)
library(patchwork)
library(ggh4x)

methods_order <- c("CELINA", "spVC", "C-SIDE", "MMM", "CTSV", "STANCE")

## ---------------------------------------------------------------------------
## Raw scores (1–6 rankings)
## ---------------------------------------------------------------------------

detect_raw <- tribble(
  ~method,  ~metric,            ~value,
  "CELINA", "Accuracy",              3,
  "CELINA", "FDP",       2,
  "CELINA", "KS",                    2,
  "spVC",   "Accuracy",              1,
  "spVC",   "FDP",       4,
  "spVC",   "KS",                    1,
  "C-SIDE", "Accuracy",              4,
  "C-SIDE", "FDP",       3,
  "C-SIDE", "KS",                    4,
  "MMM",    "Accuracy",              6,
  "MMM",    "FDP",       1,
  "MMM",    "KS",                    3,
  "CTSV",   "Accuracy",              2,
  "CTSV",   "FDP",       5,
  "CTSV",   "KS",                    5,
  "STANCE", "Accuracy",              5,
  "STANCE", "FDP",       6,
  "STANCE", "KS",                    4
)

robust_raw <- tribble(
  ~method,  ~metric,        ~value,
  "CELINA", "Consistency",       1,
  "spVC",   "Consistency",       5,
  "C-SIDE", "Consistency",       4,
  "MMM",    "Consistency",       2,
  "CTSV",   "Consistency",       6,
  "STANCE", "Consistency",       3
)

scale_raw <- tribble(
  ~method,  ~metric,           ~value,
  "CELINA", "Time (Gene)",       1,
  "CELINA", "Mem (Gene)",        4,
  "CELINA", "Time (Spot)",       3,
  "CELINA", "Mem (Spot)",        6,
  "CELINA", "Real Data",         2,
  "spVC",   "Time (Gene)",       2,
  "spVC",   "Mem (Gene)",        1,
  "spVC",   "Time (Spot)",       1,
  "spVC",   "Mem (Spot)",        1,
  "spVC",   "Real Data",         4,
  "C-SIDE", "Time (Gene)",       3,
  "C-SIDE", "Mem (Gene)",        6,
  "C-SIDE", "Time (Spot)",       2,
  "C-SIDE", "Mem (Spot)",        3,
  "C-SIDE", "Real Data",         1,
  "MMM",    "Time (Gene)",       4,
  "MMM",    "Mem (Gene)",        5,
  "MMM",    "Time (Spot)",       5,
  "MMM",    "Mem (Spot)",        4,
  "MMM",    "Real Data",         3,
  "CTSV",   "Time (Gene)",       6,
  "CTSV",   "Mem (Gene)",        2,
  "CTSV",   "Time (Spot)",       4,
  "CTSV",   "Mem (Spot)",        2,
  "CTSV",   "Real Data",         5,
  "STANCE", "Time (Gene)",       5,
  "STANCE", "Mem (Gene)",        3,
  "STANCE", "Time (Spot)",       6,
  "STANCE", "Mem (Spot)",        5,
  "STANCE", "Real Data",         6
)

use_raw <- tribble(
  ~method,  ~metric,        ~value,
  "CELINA", "Pkg & Doc.",        2,
  "CELINA", "Error",             1,
  "CELINA", "Output",            3,
  "spVC",   "Pkg & Doc.",        5,
  "spVC",   "Error",             4,
  "spVC",   "Output",            6,
  "C-SIDE", "Pkg & Doc.",        2,
  "C-SIDE", "Error",             5,
  "C-SIDE", "Output",            1,
  "MMM",    "Pkg & Doc.",        2,
  "MMM",    "Error",             3,
  "MMM",    "Output",            3,
  "CTSV",   "Pkg & Doc.",        6,
  "CTSV",   "Error",             1,
  "CTSV",   "Output",            5,
  "STANCE", "Pkg & Doc.",        4,
  "STANCE", "Error",             2,
  "STANCE", "Output",            3
)

## ---------------------------------------------------------------------------
## Combine domains and compute Overall as mean of all raw metrics
## ---------------------------------------------------------------------------

detect_df  <- detect_raw %>% mutate(domain = "Detection")
robust_df  <- robust_raw %>% mutate(domain = "Robustness")
scale_df   <- scale_raw  %>% mutate(domain = "Scalability")
use_df     <- use_raw    %>% mutate(domain = "Usability")

overall_df <- bind_rows(detect_raw, robust_raw, scale_raw, use_raw) %>%
  group_by(method) %>%
  summarise(metric = "Overall", value = mean(value), .groups = "drop") %>%
  mutate(domain = "Overall")

scores_long <- bind_rows(detect_df, robust_df, scale_df, use_df, overall_df)

## ---------------------------------------------------------------------------
## Colors — 1–6 gradient per domain; Robustness flat #90D5FF; Overall white
## ---------------------------------------------------------------------------

col_detect <- c("#FDE0DD", "#FCC5C0", "#FA9FB5", "#F768A1", "#C51B8A", "#7A0177")
col_robust <- c("#EFF8FF", "#BDE5FF", "#90D5FF", "#4DAFEF", "#1A7EC4", "#0A5189")
col_scale  <- c("#E5F5E0", "#C7E9C0", "#A1D99B", "#74C476", "#31A354", "#006D2C")
col_use    <- c("#FFF7BC", "#FEE391", "#FEC44F", "#FE9929", "#EC7014", "#CC4C02")

df_all <- scores_long %>%
  mutate(
    bucket = pmax(1L, pmin(6L, as.integer(round(value)))),
    fill = case_when(
      domain == "Detection"   ~ col_detect[bucket],
      domain == "Robustness"  ~ col_robust[bucket],
      domain == "Scalability" ~ col_scale[bucket],
      domain == "Usability"   ~ col_use[bucket],
      TRUE                    ~ "#FFFFFF"
    )
  ) %>%
  mutate(method = factor(method, levels = methods_order))

## ---------------------------------------------------------------------------
## Metric column order per domain
## ---------------------------------------------------------------------------

metric_orders <- list(
  Detection   = c("Accuracy", "FDP", "KS"),
  Robustness  = c("Consistency"),
  Scalability = c("Time (Gene)", "Mem (Gene)", "Time (Spot)", "Mem (Spot)", "Real Data"),
  Usability   = c("Pkg & Doc.", "Error", "Output"),
  Overall     = c("Overall")
)

## ---------------------------------------------------------------------------
## Panel builder
## ---------------------------------------------------------------------------

make_panel <- function(dom, show_y = FALSE) {
  df_dom <- df_all %>%
    filter(domain == dom) %>%
    mutate(
      method    = fct_relevel(method, rev(methods_order)),
      metric_id = match(metric, metric_orders[[dom]])
    ) %>%
    filter(!is.na(metric_id)) %>%
    droplevels()

  n_methods <- nlevels(df_dom$method)
  n_cols    <- length(metric_orders[[dom]])
  labs_df   <- df_dom %>% distinct(metric, metric_id)

  ggplot(df_dom, aes(y = method)) +
    geom_rounded_rect(
      aes(
        xmin   = metric_id - 0.48,
        xmax   = metric_id + 0.48,
        ymin   = as.numeric(method) - 0.48,
        ymax   = as.numeric(method) + 0.48,
        fill   = fill,
        radius = 0.15
      ),
      color     = "grey30",
      linewidth = 0.5
    ) +
    geom_text(
      aes(
        x     = metric_id,
        y     = as.numeric(method),
        label = ifelse(
          domain == "Overall",
          sprintf("%.2f", value),
          sprintf("%.0f", value)
        )
      ),
      size = 12
    ) +
    geom_text(
      data        = labs_df,
      aes(x = metric_id, y = n_methods + 0.5, label = metric),
      inherit.aes = FALSE,
      angle = 45, hjust = 0, vjust = 0, size = 9.5
    ) +
    ggh4x::facet_nested(cols = vars(domain)) +
    coord_cartesian(
      xlim = c(0.5, n_cols + 0.5),
      ylim = c(1, n_methods + 1.5),
      clip = "off"
    ) +
    scale_fill_identity() +
    theme(
      strip.text       = element_text(size = 30, face = "bold"),
      panel.background = element_rect(fill = "white", color = NA),
      strip.background = element_rect(fill = "grey80", color = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x      = element_blank(),
      axis.ticks       = element_blank(),
      axis.title       = element_blank(),
      axis.text.y      = if (show_y)
        element_text(size = 30, face = "bold")
      else
        element_blank(),
      plot.title       = element_blank()
    )
}

## ---------------------------------------------------------------------------
## Assemble and save
## ---------------------------------------------------------------------------

p_det  <- make_panel("Detection",   show_y = TRUE)
p_rob  <- make_panel("Robustness")
p_scal <- make_panel("Scalability")
p_usa  <- make_panel("Usability")
p_over <- make_panel("Overall")

# widths proportional to column count; Detection gets extra for y-axis labels
final_plot <- p_det + p_rob + p_scal + p_usa + p_over +
  plot_layout(
    nrow   = 1,
    widths = c(4.5, 1.5, 7.5, 4.5, 1.5)
  )

final_plot
ggsave(
  "./figures/output/fig6/fig6.pdf",
  width  = 35,
  height = 15,
  device = "pdf"
)
