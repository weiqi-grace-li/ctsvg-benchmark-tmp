source("./util.R")

library(PRROC); library(pROC)
library(ggplot2); library(tidyr); library(dplyr); library(tidyverse)
library(paletteer); library(patchwork)

# ── Constants ──────────────────────────────────────────────────────────────────

save_base_dirs <- c(
  "./data/detection_results/realistic/ovarian/",
  "./data/detection_results/realistic/lymph/",
  "./data/detection_results/realistic/breast/",
  "./data/detection_results/idealized/"
)
analyze_method <- c("stance", "celina", "spvc", "ctsv", "cside","mmm")

selected_sims <- c(
  "breast_rotate_0", "ovarian_rotate_0", "lymph_rotate_0",
  "stance_simulator_1alt_0.7",  "stance_simulator_1alt_1.5",
  "stance_simulator_alt_0.7",   "stance_simulator_alt_1.5",
  "celina_simulator_alt_I_hotspot",    "celina_simulator_alt_II_hotspot",
  "celina_simulator_alt_III_hotspot",  "celina_simulator_alt_IV_hotspot",
  "celina_simulator_alt_I_streak",     "celina_simulator_alt_II_streak",
  "celina_simulator_alt_III_streak",   "celina_simulator_alt_IV_streak",
  "celina_simulator_alt_I_gradient",   "celina_simulator_alt_II_gradient",
  "celina_simulator_alt_III_gradient", "celina_simulator_alt_IV_gradient",
  "celina_simulator_null_I", "celina_simulator_null_II",
  "celina_simulator_null_III", "celina_simulator_null_IV"
)
real_names      <- c("breast_rotate_0", "ovarian_rotate_0", "lymph_rotate_0")
unbalanced_sims <- c("stance_simulator_alt_0.7", "stance_simulator_alt_1.5")
test_methods    <- c("stance", "celina", "spvc-noalt", "ctsv", "cside","mmm")

# ── Load and compute metrics ───────────────────────────────────────────────────

plot2_data <- data.frame(
  simulation_name = character(), seed = integer(), test_method = character(),
  AUPRC = numeric(), ks = numeric(), AUROC = numeric(), EP = numeric(),
  FDP = numeric(), Power = numeric(), FDP_nomarker = numeric(),
  Type_1 = numeric(), Type_1_marker = numeric(), Type_1_null = numeric(),
  Type_1_other_ctsvg = numeric(), Type_1_other_marker = numeric(),
  Type_1_alt = numeric(), Type_1_marker_alt = numeric(), Type_1_null_alt = numeric(),
  Type_1_other_ctsvg_alt = numeric(), Type_1_other_marker_alt = numeric(),
  pr_curve = I(list()), p_value = I(list()),
  stringsAsFactors = FALSE
)

has_type <- function(gene_type_col, type) {
  sapply(strsplit(gene_type_col, ","), function(x) type %in% x)
}

for (save_base_dir in save_base_dirs) {
  for (method in analyze_method) {
    all_files <- list.files(
      save_base_dir,
      pattern    = paste0("^", method, "(_[0-9]+)?\\.(xlsx|csv)$"),
      full.names = FALSE
    )
    message(paste0(all_files, collapse = ", "))

    all_results <- do.call(rbind, lapply(all_files, function(file) {
      full_path <- file.path(save_base_dir, file)
      if      (grepl("\\.csv$",  file, ignore.case = TRUE)) read.csv(full_path)
      else if (grepl("\\.xlsx$", file, ignore.case = TRUE)) readxl::read_excel(full_path)
      else stop("Unsupported file type: ", file)
    }))

    for (file_method in unique(all_results$test_method)) {
      method_result <- subset(all_results[all_results$test_method == file_method, ],
                              simulation_name %in% selected_sims)

      for (sim in unique(method_result$simulation_name)) {
        temp_sim <- subset(method_result, simulation_name == sim)

        for (seed in unique(temp_sim$seed)) {
          temp <- temp_sim[temp_sim$seed == seed, ]

          tags         <- has_type(temp$gene_type, "ctsvg")
          markers      <- has_type(temp$gene_type, "marker")
          other_marker <- has_type(temp$gene_type, "other_marker") & !tags & !markers
          other_ctsvg  <- has_type(temp$gene_type, "other_ctsvg")  & !tags & !markers & !other_marker
          null_genes   <- !tags & !markers & !other_marker & !other_ctsvg

          fdr_nomarker        <- tags
          fdr_nomarker[!tags] <- !markers[!tags]

          if (sum(tags) > 0) {
            temp$p_adj[temp$p_adj == -1] <- 1
            scores   <- -log(pmax(temp$p_adj, .Machine$double.xmin))
            pr       <- pr.curve(scores.class0 = scores[tags], scores.class1 = scores[!tags], curve = TRUE)
            AUPRC    <- pr$auc.integral
            pr_curve <- pr$curve

            auc_value <- pROC::auc(suppressMessages(pROC::roc(response = tags, predictor = scores)))

            k        <- sum(tags)
            EP_value <- sum(tags[order(scores, decreasing = TRUE)[1:k]]) / k

            rejected_adj <- temp$p_adj <= 0.05
            FDP          <- sum(rejected_adj[!tags])                / sum(rejected_adj)
            Power        <- sum(rejected_adj[tags])                 / sum(tags)
            FDP_nomarker <- sum(rejected_adj[!tags & fdr_nomarker]) / sum(rejected_adj[fdr_nomarker])
          } else {
            auc_value <- NA; AUPRC <- NA; pr_curve <- NA
            EP_value  <- NA; FDP   <- NA; Power    <- NA; FDP_nomarker <- NA
          }

          p_raw  <- temp$p_value
          tested <- p_raw >= 0
          sig    <- p_raw < 0.05 & tested
          t1  <- function(m) sum(sig & m) / sum(m & tested)
          t1a <- function(m) sum(sig & m) / sum(m)

          p_ks <- p_raw[!tags & p_raw >= 0 & p_raw <= 1]
          ks   <- if (length(p_ks) > 0) ks.test(p_ks, "punif", 0, 1, exact = FALSE)$statistic else NA

          plot2_data <- rbind(plot2_data, data.frame(
            simulation_name = sim, seed = seed, test_method = file_method,
            AUPRC = AUPRC, ks = ks, AUROC = auc_value, EP = EP_value,
            FDP = FDP, Power = Power, FDP_nomarker = FDP_nomarker,
            Type_1 = t1(!tags), Type_1_marker = t1(markers),
            Type_1_null = t1(null_genes), Type_1_other_ctsvg = t1(other_ctsvg),
            Type_1_other_marker = t1(other_marker),
            Type_1_alt = t1a(!tags), Type_1_marker_alt = t1a(markers),
            Type_1_null_alt = t1a(null_genes), Type_1_other_ctsvg_alt = t1a(other_ctsvg),
            Type_1_other_marker_alt = t1a(other_marker),
            pr_curve = I(list(pr_curve)), p_value = I(list(p_raw))
          ))
        }
      }
    }
  }
}

# ── Shared plot helpers ────────────────────────────────────────────────────────

pal <- paletteer_d("tvthemes::AirNomads")
pal[6] = "#226060FF"

base_theme <- theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  legend.position  = "top",
  legend.direction = "horizontal",
  legend.title     = element_text(size = 20, face = "bold"),
  legend.text      = element_text(size = 20),
  legend.key.size  = unit(3, "lines"),
  strip.text       = element_text(size = 24, face = "bold"),
  panel.background = element_rect(fill = "white", color = "black", linewidth = 0.8),
  strip.background = element_rect(fill = "grey80", color = "black"),
  axis.title       = element_blank(),
  axis.text.y      = element_text(size = 22),
  axis.ticks.y     = element_line(size = 1.2),
  axis.text.x      = element_blank(),
  axis.ticks.x     = element_blank()
)

no_legend_theme   <- base_theme + theme(legend.position = "none")
titled_plot_theme <- theme(plot.title = element_text(size = 24, hjust = 0.5, face = "bold"))

shape_scale <- scale_shape_manual(
  name   = "",
  values = c("Unbalanced (Scenario 1)" = 8, "Realistic data" = 24),
  labels = c("Unbalanced (Scenario 1)" = "Unbalanced (Scenario 1)", "Realistic data" = "Realistic")
)

make_colour_scale <- function(method_levels, labels = method_levels) {
  scale_colour_manual(
    values     = pal,
    aesthetics = c("colour", "fill"),
    name       = "",
    limits     = method_levels,
    labels     = labels
  )
}

# ── Build df ───────────────────────────────────────────────────────────────────

lev_metric <- c(
  "EP", "AUPRC", "FDP", "KS",
  "Type_1", "Type_1_marker", "Type_1_other_marker", "Type_1_other_ctsvg", "Type_1_null",
  "Type_1_alt", "Type_1_marker_alt", "Type_1_other_marker_alt", "Type_1_other_ctsvg_alt", "Type_1_null_alt"
)

df <- plot2_data %>%
  filter(simulation_name %in% selected_sims, test_method %in% test_methods) %>%
  mutate(
    experiment_name = if_else(simulation_name %in% real_names, "real", "pure"),
    test_method     = factor(test_method, levels = c("cside", "ctsv", "spvc-noalt","celina", "stance", "mmm"))
  ) %>%
  group_by(simulation_name, test_method, experiment_name) %>%
  summarise(
    AUPRC                   = mean(AUPRC, na.rm = TRUE),
    KS                      = mean(ks, na.rm = TRUE),
    EP                      = mean(EP, na.rm = TRUE),
    FDP                     = mean(FDP, na.rm = TRUE),
    Type_1                  = mean(Type_1, na.rm = TRUE),
    Type_1_marker           = mean(Type_1_marker, na.rm = TRUE),
    Type_1_null             = mean(Type_1_null, na.rm = TRUE),
    Type_1_other_ctsvg      = mean(Type_1_other_ctsvg, na.rm = TRUE),
    Type_1_other_marker     = mean(Type_1_other_marker, na.rm = TRUE),
    Type_1_alt              = mean(Type_1_alt, na.rm = TRUE),
    Type_1_marker_alt       = mean(Type_1_marker_alt, na.rm = TRUE),
    Type_1_null_alt         = mean(Type_1_null_alt, na.rm = TRUE),
    Type_1_other_ctsvg_alt  = mean(Type_1_other_ctsvg_alt, na.rm = TRUE),
    Type_1_other_marker_alt = mean(Type_1_other_marker_alt, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols      = c(EP, AUPRC, KS, FDP,
                  Type_1, Type_1_marker, Type_1_other_marker, Type_1_other_ctsvg, Type_1_null,
                  Type_1_alt, Type_1_marker_alt, Type_1_other_marker_alt, Type_1_other_ctsvg_alt, Type_1_null_alt),
    names_to  = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = fct_relevel(metric, lev_metric),
    metric_display = fct_recode(
      metric,
      "Type I Error Only Tested"                = "Type_1",
      "Type I Error Only Tested (marker)"       = "Type_1_marker",
      "Type I Error Only Tested (null)"         = "Type_1_null",
      "Type I Error Only Tested (other ctsvg)"  = "Type_1_other_ctsvg",
      "Type I Error Only Tested (other marker)" = "Type_1_other_marker",
      "Type I Error"                            = "Type_1_alt",
      "Type I Error (marker)"                   = "Type_1_marker_alt",
      "Type I Error (null)"                     = "Type_1_null_alt",
      "Type I Error (other ctsvg)"              = "Type_1_other_ctsvg_alt",
      "Type I Error (other marker)"             = "Type_1_other_marker_alt"
    ),
    test_method = fct_recode(
      test_method,
      "spVC"   = "spvc-noalt",
      "STANCE" = "stance",
      "CELINA" = "celina",
      "C-SIDE" = "cside",
      "CTSV"   = "ctsv",
      "MMM"    = "mmm"
    )
  )

# ── Overall plot ───────────────────────────────────────────────────────────────

overall_metric <- c("EP", "AUPRC", "FDP", "KS", "Type_1_alt")

pure_overall <- subset(df, experiment_name == "pure") %>% filter(!is.na(value), metric %in% overall_metric)
real_overall <- subset(df, experiment_name == "real") %>% filter(!is.na(value), metric %in% overall_metric)

pure_overall$dataset_type <- ifelse(pure_overall$simulation_name %in% unbalanced_sims,
                                    "Unbalanced (Scenario 1)", "Balanced (Scenario 2-6)")
real_overall$dataset_type <- "Realistic data"

method_levels <- unique(pure_overall$test_method)

ref_lines <- expand.grid(metric_display = c("FDP", "Type I Error"), test_method = method_levels) %>%
  mutate(metric_display = factor(metric_display), yint = 0.05, ylab = 0.025, label = "target 0.05")

overall <- ggplot(pure_overall, aes(x = test_method, y = value, group = test_method)) +
  geom_boxplot(aes(fill = test_method), width = 0.55, alpha = 0.7,
               outlier.colour = NA, outlier.fill = NA, outlier.shape = 21, outlier.size = 2.5) +
  facet_wrap(~ metric_display, nrow = 1, scales = "fixed") +
  geom_point(data = real_overall,
             aes(x = test_method, y = value, fill = test_method, shape = dataset_type),
             color = "black", size = 5, stroke = 1.2) +
  geom_hline(data = ref_lines, aes(yintercept = yint), linetype = "dashed", color = "grey20", size = 0.6) +
  geom_text(data = ref_lines, aes(y = ylab, label = label),
            x = 3.5, y = -0.04, hjust = 0.5, color = "grey20", size = 7) +
  coord_cartesian(ylim = c(-0.05, 1.05)) +
  shape_scale +
  make_colour_scale(method_levels) +
  guides(colour = "none", fill = "none", shape = "none") +
  base_theme

overall

# ── EP and FDP subfigures ──────────────────────────────────────────────────────

sim_dim <- read.csv("./simulators/idealized/simulator_setup.csv") %>%
  mutate(simulation_name = paste0(
    sim_name,
    ifelse(!is.na(phi)     & phi     != "", paste0("_", phi),     ""),
    ifelse(!is.na(scene)   & scene   != "", paste0("_", scene),   ""),
    ifelse(!is.na(pattern) & pattern != "", paste0("_", pattern), "")
  ))

ep <- pure_overall %>%
  filter(metric == "EP") %>%
  left_join(sim_dim[, c("simulation_name", "pattern_plot", "scenario_plot")], by = "simulation_name") %>%
  mutate(pattern_plot = factor(pattern_plot, levels = c("hotspot", "streak", "gradient")))
ep$scenario_display <- paste0("Scenario ", ep$scenario_plot)

EP <- ggplot(subset(ep, dataset_type != "Unbalanced (Scenario 1)"),
             aes(x = test_method, y = value, group = test_method)) +
  geom_boxplot(aes(fill = test_method), width = 0.55, alpha = 0.8,
               outlier.colour = NA, outlier.fill = NA, outlier.shape = 21, outlier.size = 2.5) +
  facet_wrap(~ pattern_plot, nrow = 1, scales = "fixed") +
  geom_point(data = subset(ep, dataset_type == "Unbalanced (Scenario 1)"),
             aes(x = test_method, y = value, fill = test_method, color = test_method, shape = dataset_type),
             size = 5, stroke = 1.2) +
  coord_cartesian(ylim = c(-0.05, 1.05)) +
  ggtitle("EP") +
  shape_scale +
  make_colour_scale(method_levels) +
  guides(colour = guide_legend(order = 1), fill = "none",
         shape  = guide_legend(order = 2, override.aes = list(size = 6))) +
  no_legend_theme + titled_plot_theme

EP
ggsave("./figures/output/fig2/ep.pdf", EP, width = 9, height = 6.5, device = "pdf")


fdp <- pure_overall %>%
  filter(metric == "FDP") %>%
  left_join(sim_dim[, c("simulation_name", "pattern_plot", "scenario_plot")], by = "simulation_name")
fdp$scenario_display <- paste0("Scenario ", fdp$scenario_plot)
fdp <- fdp %>%
  mutate(scenario_display = fct_recode(
    factor(scenario_display),
    "Scenario 1"   = "Scenario 1",
    "Scenario 3-4" = "Scenario 3",
    "Scenario 3-4" = "Scenario 4",
    "Scenario 5-6" = "Scenario 5",
    "Scenario 5-6" = "Scenario 6"
  ))

FDP <- ggplot(subset(fdp, scenario_display != "Scenario 1"),
              aes(x = test_method, y = value, group = test_method)) +
  geom_boxplot(aes(fill = test_method), width = 0.55, alpha = 0.8,
               outlier.colour = NA, outlier.fill = NA, outlier.shape = 21, outlier.size = 2.5) +
  geom_point(data = subset(fdp, dataset_type == "Unbalanced (Scenario 1)"),
             aes(x = test_method, y = value, fill = test_method, color = test_method, shape = dataset_type),
             size = 5, stroke = 1.2) +
  geom_hline(data = ref_lines, aes(yintercept = yint), linetype = "dashed", color = "grey20", size = 0.6) +
  geom_text(data = ref_lines, aes(y = ylab, label = label),
            x = 3.5, y = -0.04, hjust = 0.5, color = "grey20", size = 7) +
  facet_wrap(~ scenario_display, nrow = 1, scales = "fixed") +
  coord_cartesian(ylim = c(-0.05, 1.05)) +
  ggtitle("FDP") +
  shape_scale +
  make_colour_scale(method_levels) +
  guides(colour = guide_legend(order = 1), fill = "none",
         shape  = guide_legend(order = 2, override.aes = list(size = 6))) +
  no_legend_theme + titled_plot_theme

FDP
ggsave("./figures/output/fig2/fdp.pdf", FDP, width = 9, height = 6.5, device = "pdf")


panel <- (overall / (EP + FDP)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 28, face = "bold"))
panel
ggsave("./figures/output/fig2/overall.pdf", panel, width = 22, height = 15, device = "pdf")

# ── Type I error ───────────────────────────────────────────────────────────────

recode_type1 <- function(df) {
  df %>% mutate(metric_display = fct_recode(
    metric,
    "Overall"      = "Type_1",
    "Marker"       = "Type_1_marker",
    "Other ctSVGs" = "Type_1_other_ctsvg",
    "Other Marker" = "Type_1_other_marker",
    "Null Genes"   = "Type_1_null"
  ))
}

Type_1_metric <- c("Type_1", "Type_1_marker", "Type_1_other_ctsvg", "Type_1_other_marker", "Type_1_null")

pure_t1 <- subset(df, experiment_name == "pure") %>% filter(!is.na(value), metric %in% Type_1_metric) %>% recode_type1()
real_t1 <- subset(df, experiment_name == "real") %>% filter(!is.na(value), metric %in% Type_1_metric) %>% recode_type1()

pure_t1$dataset_type <- ifelse(pure_t1$simulation_name %in% unbalanced_sims,
                               "Unbalanced (Scenario 1)", "Balanced (Scenario 2-6)")
real_t1$dataset_type <- "Realistic data"
method_levels        <- unique(pure_t1$test_method)

ref_lines_t1 <- expand.grid(
  metric_display = c("Overall", "Marker", "Other ctSVGs", "Other Marker", "Null Genes"),
  test_method    = method_levels
) %>% mutate(metric_display = factor(metric_display), yint = 0.05, ylab = 0.025, label = "target 0.05")

type_1 <- ggplot(subset(pure_t1, dataset_type == "Balanced (Scenario 2-6)"),
                 aes(x = test_method, y = value, group = test_method)) +
  geom_boxplot(aes(fill = test_method), width = 0.55, alpha = 0.7,
               outlier.colour = NA, outlier.fill = NA, outlier.shape = 21, outlier.size = 2.5) +
  facet_wrap(~ metric_display, nrow = 1, scales = "fixed") +
  geom_point(data = subset(pure_t1, dataset_type == "Unbalanced (Scenario 1)"),
             aes(x = test_method, y = value, fill = test_method, color = test_method, shape = dataset_type),
             size = 5, stroke = 1.2) +
  geom_point(data = real_t1,
             aes(x = test_method, y = value, fill = test_method, shape = dataset_type),
             color = "black", size = 5, stroke = 1.2) +
  geom_hline(data = ref_lines_t1, aes(yintercept = yint), linetype = "dashed", color = "grey20", size = 0.6) +
  geom_text(data = ref_lines_t1, aes(y = ylab, label = label),
            x = 3.5, hjust = 0, color = "grey20", size = 7) +
  coord_cartesian(ylim = c(-0.05, 1.05)) +
  ggtitle("Type I Error (Only Tested Genes)") +
  shape_scale +
  make_colour_scale(method_levels) +
  guides(colour = "none",
         fill   = guide_legend(order = 1),
         shape  = guide_legend(order = 2, override.aes = list(size = 6))) +
  base_theme + titled_plot_theme

type_1

# ── Type I error (alt denominator) ────────────────────────────────────────────

recode_type1_alt <- function(df) {
  df %>% mutate(metric_display = fct_recode(
    metric,
    "Overall"      = "Type_1_alt",
    "Marker"       = "Type_1_marker_alt",
    "Other ctSVG"  = "Type_1_other_ctsvg_alt",
    "Other Marker" = "Type_1_other_marker_alt",
    "Null Gene"    = "Type_1_null_alt"
  ))
}

Type_1_metric_alt <- c("Type_1_alt", "Type_1_marker_alt", "Type_1_null_alt")

pure_t1a <- subset(df, experiment_name == "pure") %>% filter(!is.na(value), metric %in% Type_1_metric_alt) %>% recode_type1_alt()
real_t1a <- subset(df, experiment_name == "real") %>% filter(!is.na(value), metric %in% Type_1_metric_alt) %>% recode_type1_alt()

pure_t1a$dataset_type <- ifelse(pure_t1a$simulation_name %in% unbalanced_sims,
                                "Unbalanced (Scenario 1)", "Balanced (Scenario 2-6)")
real_t1a$dataset_type <- "Realistic data"
method_levels         <- unique(pure_t1a$test_method)

ref_lines_t1a <- expand.grid(
  metric_display = c("Overall", "Marker", "Null Gene"),
  test_method    = method_levels
) %>% mutate(metric_display = factor(metric_display), yint = 0.05, ylab = 0.025, label = "target 0.05")

type_1_alt <- ggplot(pure_t1a, aes(x = test_method, y = value, group = test_method)) +
  geom_boxplot(aes(fill = test_method), width = 0.55, alpha = 0.7,
               outlier.colour = NA, outlier.fill = NA, outlier.shape = 21, outlier.size = 2.5) +
  facet_wrap(~ metric_display, nrow = 1, scales = "fixed") +
  geom_point(data = subset(pure_t1a, dataset_type == "Unbalanced (Scenario 1)"),
             aes(x = test_method, y = value, fill = test_method, color = test_method, shape = dataset_type),
             size = 5, stroke = 1.2) +
  geom_point(data = real_t1a,
             aes(x = test_method, y = value, fill = test_method, shape = dataset_type),
             color = "black", size = 5, stroke = 1.2) +
  geom_hline(data = ref_lines_t1a, aes(yintercept = yint), linetype = "dashed", color = "grey20", size = 0.6) +
  geom_text(data = ref_lines_t1a, aes(y = ylab, label = label),
            x = 3.5, y = -0.04, hjust = 0.5, color = "grey20", size = 7) +
  coord_cartesian(ylim = c(-0.05, 1.05)) +
  ggtitle("Type I Error") +
  shape_scale +
  make_colour_scale(method_levels) +
  guides(colour = "none",
         fill   = guide_legend(order = 1),
         shape  = guide_legend(order = 2, override.aes = list(size = 6))) +
  no_legend_theme + titled_plot_theme

type_1_alt
# ggsave("./figures/output/fig2/type1.pdf", type_1_alt, width = 20.1, height = 6.6, device = "pdf")
ggsave("./figures/output/fig2/type1_exp.pdf", type_1_alt, width = 20.1, height = 6.6, device = "pdf")

type_1_panel <- (type_1 / type_1_alt) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 28, face = "bold"))
type_1_panel
# ggsave("./figures/output/fig2/type1.pdf", type_1_panel, width = 20, height = 15, device = "pdf")
ggsave("./figures/output/fig2/type1_exp.pdf", type_1_panel, width = 20, height = 15, device = "pdf")

# ── Supplementary: C-SIDE p-value distribution ────────────────────────────────

selected_scenarios <- c(
  "stance_simulator_alt_1.5", "stance_simulator_1alt_1.5",
  "celina_simulator_alt_I_hotspot", "celina_simulator_alt_II_streak",
  "celina_simulator_alt_II_gradient", "celina_simulator_null_IV",
  "ovarian_scdesign3_remove300"
)

bw <- 0.05

p_dist_cside <- subset(plot2_data, test_method == "cside") %>%
  filter(simulation_name %in% selected_scenarios, seed == 42) %>%
  mutate(scenario_display = fct_recode(
    factor(simulation_name),
    "Idealized Scenario 1" = "stance_simulator_alt_1.5",
    "Idealized Scenario 2" = "stance_simulator_1alt_1.5",
    "Idealized Scenario 3" = "celina_simulator_alt_I_hotspot",
    "Idealized Scenario 4" = "celina_simulator_alt_II_streak",
    "Idealized Scenario 5" = "celina_simulator_alt_II_gradient",
    "Idealized Scenario 6" = "celina_simulator_null_IV",
    "Realistic"            = "ovarian_scdesign3_remove300"
  )) %>%
  select(scenario_display, p_value) %>%
  unnest(p_value)

p_binned <- p_dist_cside %>%
  mutate(bin = cut(p_value, breaks = seq(0, 1, by = bw), include.lowest = TRUE, right = TRUE)) %>%
  count(scenario_display, bin, name = "n") %>%
  group_by(scenario_display) %>%
  mutate(prop = n / sum(n), xmid = (as.numeric(bin) - 1) * bw + bw / 2) %>%
  ungroup() %>%
  mutate(scenario_display = fct_relevel(scenario_display, c(
    "Idealized Scenario 1", "Idealized Scenario 2", "Idealized Scenario 3",
    "Idealized Scenario 4", "Idealized Scenario 5", "Idealized Scenario 6",
    "Realistic"
  )))

cside_p <- ggplot(p_binned, aes(x = xmid, y = prop)) +
  geom_col(width = bw, fill = "#FF9933FF", color = "black") +
  facet_wrap(~ scenario_display, nrow = 3, scales = "fixed") +
  scale_x_continuous(breaks = seq(0, 1, by = 0.05)) +
  scale_y_continuous(breaks = seq(0, 0.5, by = 0.05)) +
  labs(x = "p-value", y = "Proportion") +
  theme_minimal(base_size = 14) +
  no_legend_theme + titled_plot_theme

ggsave("./figures/output/fig2/cside_p_dist.pdf", cside_p, width = 18, height = 20, device = "pdf")
