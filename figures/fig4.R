source("./util.R")
library(PRROC)
library(pROC)
library(dplyr)
library(ggplot2)
library(ggh4x)
library(tidyverse)
library(patchwork)
library(ggpubr)
library(paletteer)

## ---- Helpers ------------------------------------------------------------

read_result_file <- function(base_dir, file) {
  full_path <- file.path(base_dir, file)
  if (grepl("\\.csv$", file, ignore.case = TRUE)) {
    read.csv(full_path)
  } else if (grepl("\\.xlsx$", file, ignore.case = TRUE)) {
    readxl::read_excel(full_path)
  } else {
    stop("Unsupported file type: ", file)
  }
}

compute_metrics <- function(sim_seed) {
  tags    <- sapply(strsplit(sim_seed$gene_type, ","), function(x) "ctsvg"  %in% x)
  markers <- sapply(strsplit(sim_seed$gene_type, ","), function(x) "marker" %in% x)
  fdr_nomarker        <- tags
  fdr_nomarker[!tags] <- !markers[!tags]

  if (sum(tags) > 0) {
    sim_seed$p_adj[sim_seed$p_adj == -1] <- 1
    scores   <- -log(pmax(sim_seed$p_adj, .Machine$double.xmin))
    pr       <- pr.curve(scores.class0 = scores[tags], scores.class1 = scores[!tags], curve = TRUE)
    AUPRC    <- pr$auc.integral
    pr_curve <- pr$curve

    roc_obj   <- suppressMessages(pROC::roc(response = tags, predictor = scores))
    auc_value <- as.numeric(pROC::auc(roc_obj))

    k        <- sum(tags)
    EP_value <- sum(tags[order(scores, decreasing = TRUE)[1:k]]) / k

    rejected     <- sim_seed$p_adj <= 0.05
    FDP          <- sum(rejected[!tags])            / sum(rejected)
    Power        <- sum(rejected[tags])             / sum(tags)
    FDP_nomarker <- sum(rejected[!tags & fdr_nomarker]) / sum(rejected[fdr_nomarker])
  } else {
    auc_value <- AUPRC <- pr_curve <- EP_value <- FDP <- Power <- FDP_nomarker <- NA
  }

  p_raw <- sim_seed$p_value
  p_raw <- p_raw[!tags & p_raw >= 0 & p_raw <= 1]
  ks    <- if (length(p_raw) > 0) ks.test(p_raw, "punif", 0, 1, exact = FALSE)$statistic else NA

  list(
    AUPRC = AUPRC, ks = ks, AUROC = auc_value, EP = EP_value,
    FDP = FDP, Power = Power, FDP_nomarker = FDP_nomarker,
    pr_curve = list(pr_curve), p_value = list(p_raw)
  )
}

nested_strip <- ggh4x::strip_nested(
  background_x = list(element_rect(fill = "grey90", colour = "black"), element_blank()),
  text_x       = list(element_text(size = 30, face = "bold"), element_text(size = 26)),
  by_layer_x   = TRUE
)

make_heatmap <- function(data, metric, show_y = FALSE) {
  ggplot(data, aes(
    x     = .env$metric,
    y     = display_experiment_name,
    fill  = .data[[metric]] * 100,
    label = sprintf("%.0f", .data[[metric]] * 100)
  )) +
    geom_tile(color = NA, width = 1.2, height = 0.9) +
    geom_text(size = 8) +
    scale_fill_gradient2(
      low = "#035AA6FF", mid = "white", high = "#FF9933FF",
      midpoint = 40, limits = c(0, 100), oob = scales::squish
    ) +
    ggh4x::facet_nested(
      cols = vars(test_method, dataset_name),
      scales = "free_x", space = "free_x", nest_line = FALSE,
      strip = nested_strip
    ) +
    labs(x = NULL, y = NULL) +
    theme(
      strip.text.x     = element_text(size = 20, face = "bold", margin = margin(t = 3, b = 3)),
      strip.placement  = "outside",
      panel.background = element_rect(fill = "white", color = NA, linewidth = 0.8),
      panel.spacing.x  = unit(0.45, "lines"),
      axis.text.x      = element_blank(),
      axis.ticks.x     = element_blank(),
      axis.text.y      = if (show_y) element_text(size = 26, hjust = 0) else element_blank(),
      axis.ticks.y     = element_blank(),
      legend.position  = "none",
      plot.title       = element_blank()
    )
}

## ---- Gather per-seed metrics --------------------------------------------

dim_table <- read.csv("./figures/decomposition_dim.csv")
base_dirs  <- unique(dim_table$base_dir)

pr_rows        <- list()
aggregate_list <- list()

for (base_dir in base_dirs) {
  file_methods <- c("cside","ctsv", "spvc", "celina", "stance", "mmm")
  for (file_method in file_methods) {
    all_files <- list.files(
      base_dir,
      pattern    = paste0("^", file_method, "(_[0-9]+)?\\.(xlsx|csv)$"),
      full.names = FALSE
    )
    message(paste(all_files, collapse = ", "))

    all_results <- do.call(rbind, lapply(all_files, read_result_file, base_dir = base_dir))
    temp_dim    <- dim_table[dim_table$base_dir == base_dir & dim_table$file_method == file_method, ]

    for (i in seq_len(nrow(temp_dim))) {
      idx      <- all_results$test_method     == temp_dim$test_method[i] &
                  all_results$simulation_name == temp_dim$simulation_name[i]
      temp_sim <- all_results[idx, ]

      for (seed in unique(temp_sim$seed)) {
        m <- compute_metrics(temp_sim[temp_sim$seed == seed, ])
        pr_rows[[length(pr_rows) + 1]] <- data.frame(
          simulation_name = temp_dim$simulation_name[i],
          seed            = seed,
          test_method     = temp_dim$test_method[i],
          experiment_name = temp_dim$experiment_name[i],
          dataset         = temp_dim$dataset[i],
          AUPRC           = m$AUPRC, ks = m$ks, AUROC = m$AUROC, EP = m$EP,
          FDP             = m$FDP,   Power = m$Power, FDP_nomarker = m$FDP_nomarker,
          pr_curve        = I(m$pr_curve),
          p_value         = I(m$p_value)
        )
      }
    }

    temp_results <- analyze_result(
      input_dir      = temp_dim$base_dir[i],
      save_path      = NA,
      save_sheet     = "aggregate",
      analyze_method = temp_dim$file_method[i],
      threshold      = 0.05
    ) %>%
      filter(simulation_name %in% temp_dim$simulation_name)

    aggregate_list[[length(aggregate_list) + 1]] <- temp_results
  }
}

pr_results        <- do.call(rbind, pr_rows)
aggregate_results <- do.call(rbind, aggregate_list)

## ---- Save ---------------------------------------------------------------

aggregate_df <- aggregate_results %>%
  left_join(dim_table, by = c("simulation_name", "test_method"))

## ---- Aggregate & reshape ------------------------------------------------

exp_order <- c(
  "Idealized",
  "decompose_0",
  "decompose_1", "decompose_2", "decompose_1,2", "decompose_3", "decompose_1,2,3",
  "decompose_2,4", "decompose_1,2,4", "decompose_1,2,5",
  "decompose_1,2,3,4,5",
  "Realistic"
)

exp_labels <- c(
  "  +All"         = "decompose_1,2,3,4,5",
  "Baseline"       = "decompose_0",
  "  +(i)"         = "decompose_1",
  "  +(ii)"        = "decompose_2",
  "  +(i,ii,iii)"  = "decompose_1,2",
  "  +(iv)"        = "decompose_3",
  "  +(ii,v)"      = "decompose_2,4",
  "  +(i,ii,iii,v)"   = "decompose_1,2,4",
  "  +(i,ii,iii,vi)"  = "decompose_1,2,5",
  "  +(i,ii,iii,iv)"  = "decompose_1,2,3"
)

method_labels <- c(
  "CELINA"   = "celina",
  "STANCE"   = "stance",
  "CTSV"     = "ctsv",
  "C-SIDE"   = "cside",
  "spVC"     = "spvc-noalt",
  "MMM" = "mmm"
)

agg <- pr_results %>%
  group_by(dataset, experiment_name, test_method) %>%
  summarise(
    AUPRC = mean(AUPRC, na.rm = TRUE),
    KS    = mean(ks,    na.rm = TRUE),
    EP    = mean(EP,    na.rm = TRUE),
    FDP   = mean(FDP,   na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    experiment_name         = fct_relevel(experiment_name, rev(exp_order)),
    display_experiment_name = fct_recode(experiment_name, !!!exp_labels),
    test_method  = fct_recode(factor(test_method), !!!method_labels),
    test_method  = factor(test_method, levels = c("C-SIDE","CTSV","spVC","CELINA","STANCE","MMM")),
    dataset_name = factor(dataset, levels = c("breast","ovarian","lymph"))
  )

idealized <- agg %>%
  filter(experiment_name == "Idealized") %>%
  group_by(experiment_name, test_method, display_experiment_name) %>%
  summarise(across(c(KS, EP, FDP, AUPRC), mean, na.rm = TRUE), .groups = "drop")

agg_avg <- bind_rows(
  crossing(idealized, dataset_name = factor(c("breast","lymph","ovarian"))) %>%
    mutate(dataset = as.character(dataset_name)),
  filter(agg, experiment_name != "Idealized")
) %>%
  mutate(dataset_name = fct_relevel(dataset_name, c("breast","ovarian","lymph")))

## ---- Figures ------------------------------------------------------------

ep_plot  <- make_heatmap(agg_avg, "EP",  show_y = FALSE)
fdp_plot <- make_heatmap(agg_avg, "FDP", show_y = TRUE)

ep_plot
fdp_plot

## ---- Decompositions figure (decompositions.pdf) -------------------------

#### Cell-type power panel

aggregate_all_clean <- aggregate_results %>%
  mutate(dataset_name = case_when(
    str_detect(simulation_name, "ovarian") ~ "ovarian",
    str_detect(simulation_name, "breast")  ~ "breast",
    str_detect(simulation_name, "lymph")   ~ "lymph",
    TRUE ~ "other"
  )) %>%
  mutate(experiment_name = if_else(
    str_detect(simulation_name, "decompose"),
    str_extract(simulation_name, "decompose.*"),
    "scdesign3"
  ))

ct_power <- aggregate_all_clean %>%
  filter(experiment_name %in% c("decompose_0", "decompose_1", "decompose_exp_0", "decompose_exp_1")) %>%
  # filter(test_method != "ctsv") %>%
  mutate(
    dataset_name_alt = ifelse(
      (experiment_name == "decompose_0" | experiment_name == "decompose_exp_0"),
      "Baseline",
      "Baseline+(i)"
    )
  ) %>%
  group_by(experiment_name, dataset_name_alt, cell_type, test_method) %>%
  summarise(
    FDP        = (sum(rejected_adj) - sum(tp)) / sum(rejected_adj),
    Power      = sum(tp) / sum(total_true),
    cell_prop  = mean(cell_proportion),
    total_true = sum(total_true),
    .groups    = "drop"
  ) %>%
  mutate(
    test_method = fct_relevel(factor(test_method), c("cside", "ctsv", "spvc-noalt", "celina", "stance", "mmm"))
  ) %>%
  mutate(
    test_method = fct_recode(
      test_method,
      "CELINA" = "celina",
      "STANCE" = "stance",
      "CTSV"   = "ctsv",
      "C-SIDE" = "cside",
      "spVC"   = "spvc-noalt",
      "MMM"    = "mmm"
    ),
    total_true = as.numeric(total_true)
  )

power_by_ct <- ggplot(ct_power, aes(x = cell_prop, y = Power, size = total_true, fill = dataset_name_alt)) +
  geom_point(alpha = 0.8, shape = 21, color = "black") +
  ggh4x::facet_nested(cols = vars(test_method), scales = "free_x", space = "free_x") +
  coord_cartesian(ylim = c(0, 1)) +
  scale_fill_paletteer_d(
    "tvthemes::AirNomads",
    name  = "Dataset",
    guide = guide_legend(override.aes = list(aes(fill = dataset_name_alt), shape = 21, size = 8))
  ) +
  stat_cor(
    aes(x = cell_prop, y = Power, label = ..r.label..),
    method      = "pearson",
    label.x.npc = 1, label.y.npc = 0.15,
    hjust       = 1,
    size = 8, parse = TRUE, na.rm = TRUE, inherit.aes = FALSE
  ) +
  scale_size_continuous(
    name   = "Num. True ct-SVGs",
    limits = range(ct_power$total_true, na.rm = TRUE),
    range  = c(5, 12),
    breaks = pretty(range(ct_power$total_true, na.rm = TRUE), n = 3)
  ) +
  scale_x_continuous(
    breaks = c(0.2, 0.4, 0.6),
    labels = c("20%", "40%", "60%")
  ) +
  labs(x = "Cell proportion", y = "Power") +
  theme(
    panel.grid.major        = element_blank(),
    panel.grid.minor        = element_blank(),
    legend.position         = "bottom",
    legend.direction        = "horizontal",
    legend.box              = "horizontal",
    legend.title            = element_text(size = 30, face = "bold"),
    legend.text             = element_text(size = 28),
    legend.key.size         = unit(3, "lines"),
    strip.text              = element_text(size = 32, face = "bold"),
    panel.background        = element_rect(fill = "white", color = "black", linewidth = 0.8),
    strip.background        = element_rect(fill = "grey80", color = "black"),
    axis.title              = element_text(size = 30),
    axis.text.y             = element_text(size = 26),
    axis.ticks.y            = element_line(size = 1.2),
    axis.text.x             = element_text(size = 26),
    axis.ticks.x            = element_line(size = 1.2),
    legend.background       = element_blank(),
    legend.box.background   = element_blank(),
    legend.key              = element_blank(),
    plot.margin             = margin(l = 70, t = 5, r = 5, b = 5)
  )

#### Tested-proportion panel

ord <- c(
  "Idealized",
  "decompose_0",
  "decompose_1,2",
  "decompose_1,2,3",
  "decompose_1,2,5"
)

tested_ctsvg <- aggregate_df %>%
  filter(experiment_name %in% ord) %>%
  filter(test_method %in% c("celina", "mmm")) %>%
  mutate(experiment_name = fct_relevel(experiment_name, rev(ord))) %>%
  mutate(
    display_experience_name = fct_recode(
      experiment_name,
      "Baseline"        = "decompose_0",
      "  +(i,ii,iii)"   = "decompose_1,2",
      "  +(i,ii,iii,iv)" = "decompose_1,2,3",
      "  +(i,ii,iii,vi)" = "decompose_1,2,5"
    ),
    test_method = fct_recode(
      factor(test_method, levels = c("cside","celina", "mmm")),
      "C-SIDE" = "cside",
      "CELINA" = "celina",
      "MMM" = "mmm"
    )
  ) %>%
  group_by(experiment_name, display_experience_name, test_method) %>%
  summarise(
    tested_prop = sum(true_tested) / sum(total_true),
    .groups     = "drop"
  )

tested_proportions <- ggplot(tested_ctsvg, aes(
  x     = tested_prop,
  y     = display_experience_name,
  group = test_method
)) +
  ggh4x::facet_nested(cols = vars(test_method), scales = "free_x", space = "free_x", nest_line = TRUE) +
  geom_col(position = position_dodge(width = 0.6), width = 0.7, fill = "grey40") +
  labs(x = "Proportion of Ground Truth ct-SVGs Tested", y = NULL) +
  coord_cartesian(xlim = c(0, 1.2)) +
  geom_text(
    aes(label = sprintf("%.0f%%", tested_prop * 100)),
    position = position_dodge(width = 0.8),
    hjust    = -0.15,
    size     = 9,
    color    = "black"
  ) +
  theme_bw(base_size = 18) +
  theme(
    panel.grid.major      = element_blank(),
    panel.grid.minor      = element_blank(),
    legend.position       = "none",
    legend.direction      = "horizontal",
    legend.box            = "horizontal",
    legend.title          = element_text(size = 30, face = "bold"),
    legend.text           = element_text(size = 28),
    legend.key.size       = unit(3, "lines"),
    strip.text            = element_text(size = 32, face = "bold"),
    panel.background      = element_rect(fill = "white", color = "black", linewidth = 0.8),
    strip.background      = element_rect(fill = "grey80", color = "black"),
    axis.title            = element_text(size = 30),
    axis.text             = element_blank(),
    axis.ticks            = element_blank(),
    legend.background     = element_blank(),
    legend.box.background = element_blank(),
    legend.key            = element_blank()
  )

#### Assemble and save

final_plot <-
  (
    (ep_plot           + theme(plot.margin = margin(t = 20, l = 100, b = 60))) /
    (power_by_ct       + theme(plot.margin = margin(t = 20, l = 100))) /
    (tested_proportions + theme(plot.margin = margin(t = 20, l = 100)))
  ) +
  plot_layout(heights = c(2.5, 1, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 28, face = "bold"))

ggsave("./figures/output/fig4/decompositions.pdf", final_plot, width = 28, height = 31, device = "pdf")
