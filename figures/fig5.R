library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(forcats)
library(scales)
library(patchwork)
library(paletteer)
library(ggh4x)
library(tibble)

source("./util.R")  # provides collect_details

## ---- Configuration -------------------------------------------------------

base_dir      <- "./data/detection_results/realistic"
mod_base_dir  <- "./data/detection_results/realistic"
methods       <- c("cside", "ctsv", "spvc", "celina", "stance", "mmm")
top_number    <- seq(10, 200, 10)

method_labels <- c(
  "C-SIDE" = "cside",
  "CTSV"   = "ctsv",
  "spVC"   = "spvc",
  "CELINA" = "celina",
  "STANCE" = "stance",
  "MMM"    = "mmm"
)

# Method color palette — shared by panels a/b; STANCE-mod shown in gray
fig2_pal      <- as.character(paletteer::paletteer_d("tvthemes::AirNomads"))
fig2_pal[3]   <- "#FFD700"
fig2_pal[6]   <- "#226060FF"
method_colors <- c(
  setNames(fig2_pal[1:5], c("C-SIDE", "CTSV", "spVC", "CELINA", "STANCE")),
  "STANCE-mod" = "gray15",
  "MMM"        = fig2_pal[6]
)


## ---- Helpers ---------------------------------------------------------------

read_method_files <- function(dataset, method, dir = base_dir) {
  ext       <- if (method %in% c("celina", "stance", "mmm")) "xlsx" else "csv"
  dir_path  <- file.path(dir, dataset)
  all_files <- list.files(dir_path, full.names = TRUE)
  matched   <- all_files[
    grepl(paste0("^", method, "(_[0-9]+)?\\.", ext, "$"), basename(all_files)) &
    !grepl("runtime", basename(all_files))
  ]
  if (length(matched) == 0) return(NULL)
  rows <- lapply(matched, function(f) {
    if (ext == "csv") read.csv(f) else as.data.frame(readxl::read_excel(f))
  })
  do.call(rbind, rows)
}

read_stance_mod_files <- function(dataset) {
  df <- read_method_files(dataset, "stance", dir = mod_base_dir)
  if (is.null(df)) return(NULL)
  df <- df[df$test_method == "stance-mod", ]
  df[grepl(paste0("^", dataset, "_rotation_[0-9]+$"), df$simulation_name), ]
}

# fill_zero: when TRUE, seeds with 0 tested genes contribute 0 instead of NA
# (used for C-SIDE which reports 0 tested pairs rather than skipping)
compute_cat_curve <- function(data, k_values, fill_zero = FALSE) {
  data      <- data[data$p_value != -1, ]
  data$pair <- paste(data$cell_type, data$gene_name, sep = "-")
  rotations <- unique(data$simulation_name)
  seeds     <- unique(data$seed)

  seed_cats <- lapply(seeds, function(s) {
    sd        <- data[data$seed == s, ]
    rot_lists <- lapply(rotations, function(rot) {
      rd <- sd[sd$simulation_name == rot, ]
      rd$pair[order(rd$p_value)]
    })
    min_tested <- min(sapply(rot_lists, length))
    if (min_tested == 0) return(rep(if (fill_zero) 0 else NA_real_, length(k_values)))
    sapply(k_values, function(k) {
      eff_k <- min(k, min_tested)
      length(Reduce(intersect, lapply(rot_lists, head, eff_k))) / eff_k
    })
  })
  rowMeans(do.call(cbind, seed_cats), na.rm = TRUE)
}

compute_ep_df <- function(data) {
  data %>%
    mutate(
      p_adj    = ifelse(p_adj == -1, 1, p_adj),
      is_ctsvg = sapply(strsplit(gene_type, ","), function(x) "ctsvg" %in% x)
    ) %>%
    group_by(simulation_name, seed) %>%
    group_modify(function(sub, ...) {
      scores <- -log(pmax(sub$p_adj, .Machine$double.xmin))
      k      <- sum(sub$is_ctsvg)
      if (k == 0) return(data.frame(ep = NA_real_))
      data.frame(ep = sum(sub$is_ctsvg[order(scores, decreasing = TRUE)[1:k]]) / k)
    }) %>%
    ungroup()
}

## ---- Compute CAT: ovarian only -------------------------------------------

cat_rows <- list()

for (m in methods) {
  data <- read_method_files("ovarian", m)
  if (is.null(data) || nrow(data) == 0) {
    message(sprintf("No data: ovarian / %s", m))
    next
  }
  message(sprintf("CAT: ovarian / %s", m))
  cat_rows[[length(cat_rows) + 1]] <- data.frame(
    method   = m,
    k        = top_number,
    cat      = compute_cat_curve(data, top_number, fill_zero = (m == "cside")),
    group_id = m
  )
}

data <- tryCatch(read_stance_mod_files("ovarian"), error = function(e) NULL)
if (!is.null(data) && nrow(data) > 0) {
  message("CAT: ovarian / stance-mod")
  cat_rows[[length(cat_rows) + 1]] <- data.frame(
    method   = "stance",
    k        = top_number,
    cat      = compute_cat_curve(data, top_number),
    group_id = "STANCE-mod"
  )
}

cat_df <- do.call(rbind, cat_rows) %>%
  mutate(
    method   = fct_recode(factor(method, levels = methods), !!!method_labels),
    group_id = if_else(
      group_id %in% methods,
      as.character(fct_recode(factor(group_id, levels = methods), !!!method_labels)),
      group_id
    ),
    group_id = factor(group_id, levels = names(method_colors))
  )

## ---- Compute EP: ovarian only --------------------------------------------

ep_rows <- list()

for (m in methods) {
  data <- read_method_files("ovarian", m)
  if (is.null(data) || nrow(data) == 0) next
  message(sprintf("EP: ovarian / %s", m))
  ep_data <- compute_ep_df(as.data.frame(data))
  ep_rows[[length(ep_rows) + 1]] <- data.frame(
    method_raw      = m,
    simulation_name = ep_data$simulation_name,
    seed            = ep_data$seed,
    ep              = ep_data$ep
  )
}

data <- tryCatch(read_stance_mod_files("ovarian"), error = function(e) NULL)
if (!is.null(data) && nrow(data) > 0) {
  message("EP: ovarian / stance-mod")
  ep_data <- compute_ep_df(as.data.frame(data))
  ep_rows[[length(ep_rows) + 1]] <- data.frame(
    method_raw      = "stance-mod",
    simulation_name = ep_data$simulation_name,
    seed            = ep_data$seed,
    ep              = ep_data$ep
  )
}

ep_df <- do.call(rbind, ep_rows) %>%
  filter(!is.na(ep)) %>%
  mutate(
    method = case_when(
      method_raw == "stance-mod" ~ "STANCE-mod",
      TRUE ~ as.character(
        fct_recode(factor(method_raw, levels = methods), !!!method_labels)
      )
    ),
    method        = factor(method, levels = names(method_colors)),
    rotation_type = if_else(grepl("_0$", simulation_name), "Original", "Rotated")
  )

ep_df_avg <- ep_df %>%
  group_by(method, simulation_name, rotation_type) %>%
  summarise(ep = mean(ep, na.rm = TRUE), .groups = "drop")

## ---- Runtime data (panels c/d) -------------------------------------------

runtime_files <- c(
  "./data/detection_results/scalability/runtime.csv",
  "./data/detection_results/scalability/runtime_stance.csv",
  "./data/detection_results/scalability/runtime_celina.csv",
  "./data/detection_results/scalability/runtime_ctsv.csv",
  "./data/detection_results/scalability/runtime_all_cside.csv",
  "./data/detection_results/scalability/runtime_all_spvc.csv",
  "./data/detection_results/scalability/runtime_all_mmm.csv"
)

runtime_results <- lapply(runtime_files, function(file) {
  message(paste0("Reading file- ", file))
  tmp <- read.csv(file) %>% filter(status == "ok")
  if ("parameter_setting" %in% names(tmp)) tmp <- select(tmp, -parameter_setting)
  tmp
}) %>% bind_rows()

runtime_df <- runtime_results %>%
  mutate(
    dataset = sub("-\\d+$", "", dataset),
    n_side  = as.numeric(sub(".*n_side:\\s*(\\d+).*", "\\1", dataset)),
    n_genes = ifelse(
      grepl("n_gene", dataset),
      as.numeric(sub(".*n_gene:\\s*(\\d+).*", "\\1", dataset)),
      100
    ),
    elapsed = elapsed / 60,
    peak_mb = peak_mb / 1024,
    method  = factor(method, levels = c("cside", "ctsv", "spvc-noalt", "celina", "stance", "mmm"))
  ) %>%
  mutate(
    type = ifelse(n_side == 30, "Varying Genes (900 Spots)", "Varying Spots (100 Genes)")
  ) %>%
  select(method, elapsed, peak_mb, n_side, n_genes, type)

runtime_long <- runtime_df %>%
  pivot_longer(cols = c(elapsed, peak_mb), names_to = "metric", values_to = "value") %>%
  mutate(
    label  = ifelse(type == "Varying Genes (900 Spots)", n_genes, n_side^2),
    method = factor(method, c("cside", "ctsv", "spvc-noalt", "celina", "stance", "mmm"))
  ) %>%
  mutate(
    method = fct_recode(method,
      "CELINA" = "celina", "STANCE" = "stance", "CTSV" = "ctsv",
      "C-SIDE" = "cside",  "spVC"   = "spvc-noalt", "MMM" = "mmm"
    ),
    metric = fct_recode(factor(metric),
      "Runtime"      = "elapsed",
      "Peak Memory"  = "peak_mb"
    )
  )

cutoff <- 15
hybrid_trans <- trans_new(
  "hybrid",
  transform = function(y) ifelse(y <= cutoff, y, cutoff + log(y - cutoff + 1)),
  inverse   = function(z) ifelse(z <= cutoff, z, cutoff + (exp(z - cutoff) - 1))
)

ann_runtime <- runtime_long %>% filter(metric == "Runtime")     %>% distinct(metric, type)
ann_peak    <- runtime_long %>% filter(metric == "Peak Memory") %>% distinct(metric, type)

pane_c_base <- ggplot(runtime_long, aes(x = label, y = value, color = method, group = method)) +
  geom_line(linewidth = 1.5) +
  # geom_point(size = 4) +
  facet_grid(
    rows   = vars(metric),
    cols   = vars(type),
    scales = "free"
  ) +
  ggh4x::facetted_pos_scales(
    y = list(
      metric == "Runtime" ~ scale_y_continuous(
        trans  = hybrid_trans,
        breaks = c(0, 5, 10, 60, 1440),
        labels = c("0 Min", "5 Min", "10 Min", "1 Hour", "1 Day"),
        limits = c(0, 2880)
      ),
      metric == "Peak Memory" ~ scale_y_continuous(
        breaks = c(0, 5, 10, 15, 20, 25),
        labels = c("0 GB", "5 GB", "10 GB", "15 GB", "20 GB", "25 GB"),
        limits = c(0, 25)
      )
    )
  ) +
  scale_color_manual(values = fig2_pal, name = "Methods") +
  theme(
    panel.grid.major     = element_blank(),
    panel.grid.minor     = element_blank(),
    legend.position      = "none",
    legend.direction     = "horizontal",
    legend.box           = "horizontal",
    # legend.title         = element_text(size = 30, face = "bold"),
    # legend.text          = element_text(size = 28),
    # legend.key.size      = unit(3, "lines"),
    strip.text           = element_text(size = 32, face = "bold"),
    panel.background     = element_rect(fill = "white", color = "black", linewidth = 0.8),
    strip.background     = element_rect(fill = "grey90", color = "black"),
    axis.text.y          = element_text(size = 26),
    axis.ticks.y         = element_line(linewidth = 1.2),
    axis.text.x          = element_text(size = 26),
    axis.ticks.x         = element_line(linewidth = 1.2),
    legend.background    = element_blank(),
    legend.box.background = element_blank(),
    legend.key           = element_blank(),
    axis.title           = element_blank(),
    plot.title           = element_blank()
  )

pane_c <- pane_c_base +
  geom_hline(
    data        = ann_peak,
    aes(yintercept = 16),
    linetype    = "dashed", color = "black", linewidth = 1.2,
    inherit.aes = FALSE
  ) +
  geom_text(
    data        = ann_peak,
    aes(x = 0, y = 16, label = "16 GB Limit"),
    hjust = -0.01, vjust = -0.5, size = 10,
    inherit.aes = FALSE
  ) +
  geom_hline(
    data        = ann_runtime,
    aes(yintercept = 120),
    linetype    = "dashed", color = "black", linewidth = 1.2,
    inherit.aes = FALSE
  ) +
  geom_text(
    data        = ann_runtime,
    aes(x = 0, y = 160, label = "2 Hours"),
    hjust = -0.01, vjust = 0, size = 10,
    inherit.aes = FALSE
  )

## ---- Decomposition runtime data (panel d) --------------------------------

decomposition_files <- c(
  "./data/detection_results/decomposed/breast/runtime_all_celina.csv", 
  "./data/detection_results/decomposed/breast/runtime_all_cside.csv", 
  "./data/detection_results/decomposed/breast/runtime_all_ctsv.csv", 
  "./data/detection_results/decomposed/breast/runtime_all_mmm.csv", 
  "./data/detection_results/decomposed/breast/runtime_all_spvc.csv", 
  "./data/detection_results/decomposed/breast/runtime_all_stance.csv", 
  "./data/detection_results/decomposed/lymph/runtime_all_celina.csv", 
  "./data/detection_results/decomposed/lymph/runtime_all_cside.csv", 
  "./data/detection_results/decomposed/lymph/runtime_all_ctsv.csv", 
  "./data/detection_results/decomposed/lymph/runtime_all_mmm.csv", 
  "./data/detection_results/decomposed/lymph/runtime_all_spvc.csv", 
  "./data/detection_results/decomposed/lymph/runtime_all_stance.csv",
  "./data/detection_results/decomposed/ovarian/runtime_all_celina.csv", 
  "./data/detection_results/decomposed/ovarian/runtime_all_cside.csv", 
  "./data/detection_results/decomposed/ovarian/runtime_all_ctsv.csv", 
  "./data/detection_results/decomposed/ovarian/runtime_all_mmm.csv", 
  "./data/detection_results/decomposed/ovarian/runtime_all_spvc.csv", 
  "./data/detection_results/decomposed/ovarian/runtime_all_stance.csv"
)

selected_exp <- c(
  "0", "exp_0",
  "1,2", "exp_1,2",
  "1,2,3", "exp_1,2,3",
  "1,2,5", "exp_1,2,5",
  "1,2,3,4,5", "exp_1,2,3,4,5"
)

decompose_results <- lapply(decomposition_files, read.csv) %>% bind_rows()

decompose_df <- decompose_results %>%
  mutate(
    dataset_name    = str_extract(dataset, "^[^_]+"),
    experiment_name = str_extract(dataset, "(?<=decompose(_exp)?_).*")
  ) %>%
  mutate(
    experiment_name = ifelse(experiment_name %in% c("", "-", NA), NA, experiment_name),
    experiment_name = sub("-\\d+$", "", experiment_name)
  ) %>%
  filter(experiment_name %in% selected_exp) %>%
  filter(status == "ok") %>%
  filter(method != "spvc-gam") %>%
  mutate(
    experiment_name = factor(experiment_name, rev(selected_exp))
  ) %>%
  mutate(
    experiment_name_display = fct_recode(factor(experiment_name),
      "Baseline"       = "0",       "Baseline"       = "exp_0",
      "  +(i,ii,iii)"  = "1,2",     "  +(i,ii,iii)"  = "exp_1,2",
      "  +(i,ii,iii,iv)" = "1,2,3", "  +(i,ii,iii,iv)" = "exp_1,2,3",
      "  +(i,ii,iii,vi)" = "1,2,5", "  +(i,ii,iii,vi)" = "exp_1,2,5",
      "  +All"         = "1,2,3,4,5", "  +All"        = "exp_1,2,3,4,5"
    )
  ) %>%
  mutate(
    test_method = fct_recode(factor(method),
      "CELINA" = "celina", "STANCE" = "stance", "CTSV" = "ctsv",
      "C-SIDE" = "cside",  "spVC"   = "spvc-noalt", "MMM"    = "mmm"
    ),
    test_method = factor(test_method, levels = c("C-SIDE", "CTSV", "spVC", "CELINA", "STANCE", "MMM"))
  ) %>%
  group_by(dataset_name, experiment_name_display, test_method) %>%
  summarise(cpu_time = mean(elapsed * ncores, na.rm = TRUE) / 60, .groups = "drop") %>%
  mutate(dataset_name = factor(dataset_name, c("breast", "ovarian", "lymph")))

ratio_df <- decompose_df %>%
  group_by(dataset_name, test_method) %>%
  mutate(baseline_cpu = cpu_time[experiment_name_display == "Baseline"][1]) %>%
  ungroup() %>%
  mutate(cpu_ratio = cpu_time / baseline_cpu)

pane_d <- ggplot(ratio_df,
  aes(
    x     = "CPU Runtime",
    y     = experiment_name_display,
    fill  = log1p(cpu_ratio),
    label = sprintf("%.0f", ceiling(cpu_ratio))
  )
) +
  geom_tile(color = NA, width = 1.2, height = 0.9) +
  geom_text(size = 8) +
  scale_fill_gradient2(
    low      = "#035AA6FF",
    mid      = "white",
    high     = "#FF9933FF",
    midpoint = log1p(2),
    limits   = range(log1p(c(0, 450)))
  ) +
  ggh4x::facet_nested(
    cols      = vars(test_method, dataset_name),
    scales    = "free_x",
    space     = "free_x",
    nest_line = FALSE,
    strip     = ggh4x::strip_nested(
      background_x = list(
        element_rect(fill = "grey90", colour = "black"),
        element_blank()
      ),
      text_x = list(
        element_text(size = 30, face = "bold"),
        element_text(size = 26)
      ),
      by_layer_x = TRUE
    )
  ) +
  labs(x = NULL, y = NULL, fill = "CPU Time (Ratio to Baseline)") +
  theme(
    strip.placement   = "outside",
    panel.background  = element_rect(fill = "white", color = NA, linewidth = 0.8),
    panel.spacing.x   = unit(0.45, "lines"),
    axis.text.x       = element_blank(),
    axis.ticks.x      = element_blank(),
    axis.text.y       = element_blank(),
    axis.ticks.y      = element_blank(),
    legend.position   = "bottom",
    legend.text       = element_blank(),
    plot.title        = element_blank(),
    legend.key.width  = unit(3, "cm"),
    legend.title      = element_text(size = 26, face = "bold")
  )

## ---- QC runtime table (panel e) — QC1 only, no facet --------------------

DATASET_ORDER    <- c("breast", "kidney_inter", "kidney_core", "brain", "lung")
METHOD_ORDER     <- c("C-SIDE", "CTSV", "spVC", "CELINA", "STANCE", "MMM")
dataset_labels_7 <- c(
  breast       = "Breast",
  kidney_inter = "Kidney Inter",
  kidney_core  = "Kidney Core",
  brain        = "DLPFC",
  lung         = "Lung"
)
VAR_LEVELS <- c("# genes", "# spots", METHOD_ORDER)

qc_meta <- tibble(
  dataset  = DATASET_ORDER,
  f1_genes = c(8384, 15426, 14850, 15124, 17553),
  f1_spots = c(345L, 2048L, 2917L, 3639L, 3848L)
)

cpu_raw <- tribble(
  ~qc,  ~dataset,        ~method,    ~value,
  "f1", "breast",        "C-SIDE",    3,
  "f1", "breast",        "CTSV",      219,
  "f1", "breast",        "spVC",      3,
  "f1", "breast",        "CELINA",    0,
  "f1", "breast",        "STANCE",    0,
  "f1", "breast",        "MMM",       1,
  "f1", "kidney_inter",  "C-SIDE",    9,
  "f1", "kidney_inter",  "CTSV",      1241,
  "f1", "kidney_inter",  "spVC",      576,
  "f1", "kidney_inter",  "CELINA",    6,
  "f1", "kidney_inter",  "STANCE",    1513,
  "f1", "kidney_inter",  "MMM",       38,
  "f1", "kidney_core",   "C-SIDE",    10,
  "f1", "kidney_core",   "CTSV",      990,
  "f1", "kidney_core",   "spVC",      371,
  "f1", "kidney_core",   "CELINA",    49,
  "f1", "kidney_core",   "STANCE",    NA,
  "f1", "kidney_core",   "MMM",       122,
  "f1", "brain",         "C-SIDE",    10,
  "f1", "brain",         "CTSV",      NA,
  "f1", "brain",         "spVC",      1054,
  "f1", "brain",         "CELINA",    11,
  "f1", "brain",         "STANCE",    NA,
  "f1", "brain",         "MMM",       209,
  "f1", "lung",          "C-SIDE",    12,
  "f1", "lung",          "CTSV",      NA,
  "f1", "lung",          "spVC",      NA,
  "f1", "lung",          "CELINA",    120,
  "f1", "lung",          "STANCE",    NA,
  "f1", "lung",          "MMM",       107
)

COL_ERROR <- "grey75"
COL_META  <- "white"

# Match pane_d: log1p scale, midpoint = log1p(2), max = log1p(400)
lmid_qc <- log1p(2)
lmax_qc <- log1p(400)
cpu_ramp    <- colorRamp(c("#035AA6FF", "white", "#FF9933FF"))
cpu_col_fn  <- function(x) {
  lx   <- log1p(x)
  t    <- ifelse(lx <= lmid_qc,
                 0.5 * lx / lmid_qc,
                 0.5 + 0.5 * (lx - lmid_qc) / (lmax_qc - lmid_qc))
  cols <- cpu_ramp(pmax(0, pmin(1, t)))
  rgb(cols[, 1], cols[, 2], cols[, 3], maxColorValue = 255)
}

qc1 <- bind_rows(
  qc_meta |> transmute(
    dataset, variable = "# genes",
    label = formatC(f1_genes, format = "d", big.mark = ","),
    fill  = COL_META, bold = FALSE
  ),
  qc_meta |> transmute(
    dataset, variable = "# spots",
    label = formatC(f1_spots, format = "d", big.mark = ","),
    fill  = COL_META, bold = FALSE
  ),
  cpu_raw |> filter(qc == "f1") |> mutate(
    variable = method,
    is_error = is.na(value),
    label    = if_else(is_error, ">1,920", as.character(as.integer(value))),
    fill     = if_else(is_error, COL_ERROR, cpu_col_fn(coalesce(value, 0))),
    bold     = TRUE
  ) |> select(dataset, variable, label, fill, bold)
) |> mutate(
  variable = factor(variable, levels = VAR_LEVELS),
  dataset  = factor(dataset, levels = rev(DATASET_ORDER),
                    labels = rev(dataset_labels_7[DATASET_ORDER]))
)

pane_e <- ggplot(qc1, aes(x = variable, y = dataset, fill = fill)) +
  geom_tile(color = "white", linewidth = 0) +
  geom_text(data = filter(qc1, !bold), aes(label = label), size = 8) +
  geom_text(data = filter(qc1,  bold), aes(label = label), size = 8) +
  facet_wrap(~ "CPU Time") +
  scale_fill_identity() +
  scale_x_discrete(position = "top") +
  scale_y_discrete(expand = expansion(add = c(0.5, 1.5))) +
  labs(x = NULL, y = NULL) +
  theme(
    panel.grid       = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    panel.border     = element_rect(fill = NA, color = "black", linewidth = 0.8),
    strip.text       = element_text(size = 32, face = "bold"),
    strip.background = element_rect(fill = "grey90", color = "black"),
    axis.text      = element_blank(),
    axis.text.y      = element_blank(),
    axis.ticks       = element_blank()
  )

## ---- Real-data stacked bar (panel f) — fig7 pane_b ----------------------

ALL_METHODS    <- c("cside", "ctsv", "spvc", "celina", "stance", "mmm")
METHOD_LABELS  <- c("cside" = "C-SIDE", "ctsv" = "CTSV", "spvc-noalt" = "spVC",
                    "celina" = "CELINA", "stance" = "STANCE", "mmm" = "MMM")
DATASETS       <- c("breast", "brain_filtered_0.1", "lung_filtered_0.1",
                    "kidney_inter_filtered_0.1", "kidney_core_filtered_0.1")

real_data_results <- collect_details(
  input_dir      = c(
    "./data/detection_results/real_data/breast/",
    "./data/detection_results/real_data/dlpfc/",
    "./data/detection_results/real_data/kidney/",
    "./data/detection_results/real_data/lung/"
  ),
  analyze_method = ALL_METHODS
)

real_data_df <- dplyr::bind_rows(
  lapply(Filter(Negate(is.null), real_data_results), function(df) {
    dplyr::mutate(df, tested = (p_value >= 0))
  })
)

filter_2 <- real_data_df %>%
  filter(simulation_name %in% DATASETS) %>%
  mutate(
    p_value = replace(p_value, is.na(p_value), -1),
    p_adj   = replace(p_adj,   is.na(p_adj),   -1),
    tested   = p_value >= 0,
    rejected = tested & p_adj < 0.05
  )

pair_wide <- filter_2 %>%
  select(simulation_name, gene_name, cell_type, test_method, rejected) %>%
  pivot_wider(names_from = test_method, values_from = rejected, values_fill = FALSE)

tested_wide <- filter_2 %>%
  select(simulation_name, gene_name, cell_type, test_method, tested) %>%
  pivot_wider(names_from = test_method, values_from = tested, values_fill = FALSE)

for (m in ALL_METHODS) {
  if (!m %in% colnames(pair_wide))   pair_wide[[m]]   <- FALSE
  if (!m %in% colnames(tested_wide)) tested_wide[[m]] <- FALSE
}

NON_STANCE <- c("cside", "ctsv", "spvc-noalt", "celina", "mmm")

ds_all5 <- filter_2 %>%
  filter(test_method %in% NON_STANCE) %>%
  distinct(simulation_name, test_method) %>%
  group_by(simulation_name) %>%
  summarise(n_methods = n_distinct(test_method), .groups = "drop") %>%
  filter(n_methods == 5) %>%
  pull(simulation_name)

stacked_colors <- c(
  "Not tested"       = "grey85",
  "Tested, not sig." = "grey65",
  "Significant"      = "grey35"
)
cat_levels <- c("Not tested", "Tested, not sig.", "Significant")

f7_pb_rows <- list()
for (ds in DATASETS) {
  pw      <- pair_wide   %>% filter(simulation_name == ds)
  tw      <- tested_wide %>% filter(simulation_name == ds)
  n_total <- nrow(pw)
  if (n_total == 0) next
  for (m in NON_STANCE) {
    if (!m %in% colnames(pw)) pw[[m]] <- FALSE
    if (!m %in% colnames(tw)) tw[[m]] <- FALSE
  }
  n_rej <- rowSums(pw[, NON_STANCE, drop = FALSE])
  n_tst <- rowSums(tw[, NON_STANCE, drop = FALSE])
  for (thresh in 1:5) {
    lbl <- if (thresh == 5) "All" else as.character(thresh)
    if (thresh == 5 && !ds %in% ds_all5) next
    f7_pb_rows[[length(f7_pb_rows) + 1]] <- data.frame(
      dataset        = ds,
      threshold      = lbl,
      not_tested     = sum(n_tst <  thresh) / n_total,
      tested_not_sig = sum(n_tst >= thresh & n_rej <  thresh) / n_total,
      significant    = sum(n_tst >= thresh & n_rej >= thresh) / n_total
    )
  }
}

f7_pb_df <- do.call(rbind, f7_pb_rows) %>%
  group_by(threshold) %>%
  summarise(across(c(not_tested, tested_not_sig, significant), mean), .groups = "drop") %>%
  pivot_longer(c(not_tested, tested_not_sig, significant),
               names_to = "category", values_to = "prop") %>%
  mutate(
    category  = factor(case_when(
      category == "not_tested"     ~ "Not tested",
      category == "tested_not_sig" ~ "Tested, not sig.",
      TRUE                         ~ "Significant"
    ), levels = cat_levels),
    threshold = factor(threshold, levels = c("1", "2", "3", "4", "All"))
  )

fig7_theme <- theme(
  panel.grid.major  = element_blank(),
  panel.grid.minor  = element_blank(),
  strip.text        = element_text(size = 32, face = "bold"),
  panel.background  = element_rect(fill = "white", color = "black", linewidth = 0.8),
  strip.background  = element_rect(fill = "grey90", color = "black"),
  legend.position   = "bottom",
  legend.direction  = "horizontal",
  legend.title      = element_blank(),
  legend.text       = element_text(size = 26),
  legend.key.size   = unit(2, "lines"),
  axis.title.x      = element_text(size = 32, margin = margin(b = 15)),
  axis.title.y      = element_text(size = 32, margin = margin(t = 15)),
  axis.text.x       = element_text(size = 26, hjust = 0.5),
  axis.ticks.x      = element_blank()
)

pane_f <- ggplot(f7_pb_df, aes(x = threshold, y = prop, fill = category)) +
  geom_col(position = "stack", width = 0.65, color = NA, linewidth = 0) +
  coord_flip() +
  facet_wrap(~ "Proportion of gene-cell-type pairs") +
  scale_fill_manual(name = NULL, values = stacked_colors,
                    guide = guide_legend(reverse = TRUE)) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(x = "Min. # methods agreeing", y = NULL) +
  fig7_theme +
  theme(
    axis.text.y  = element_text(size = 26),
    axis.title.y = element_text(size = 32, margin = margin(r = 15)),
    axis.title.x = element_text(size = 32, margin = margin(t = 15))
  )

## ---- Shared theme (panels a/b) --------------------------------------------

shared_theme <- theme(
  panel.grid.major  = element_blank(),
  panel.grid.minor  = element_blank(),
  legend.position   = "none",
  legend.direction  = "horizontal",
  legend.title      = element_text(size = 26, face = "bold"),
  legend.text       = element_text(size = 26),
  legend.key.size   = unit(1.5, "lines"),
  strip.text        = element_text(size = 32, face = "bold"),
  strip.background  = element_rect(fill = "grey90", color = "black"),
  panel.background  = element_rect(fill = "white", color = "black", linewidth = 0.8),
  axis.title.x      = element_text(size = 32),
  axis.title.y      = element_text(size = 32),
  axis.text.x       = element_text(size = 26),
  axis.text.y       = element_text(size = 26),
  # axis.line         = element_line(color = "black", linewidth = 0.8)
)

## ---- Panel a: CAT — one solid line per method/variant, ovarian only ------

pane_a <- ggplot(cat_df, aes(
  x     = k,
  y     = cat,
  color = group_id,
  group = group_id
)) +
  geom_line(linewidth = 1.5, na.rm = TRUE) +
  scale_color_manual(name = "Method", values = method_colors) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    labels = c("0", "0.25", "0.50", "0.75", "1.00")
  ) +
  facet_wrap(~ "CAT") +
  labs(x = "Top M gene-cell-type pairs", y = NULL) +
  shared_theme

## ---- Panel b: EP — one point per (method x rotation), ovarian only ------

pane_b <- ggplot(ep_df_avg, aes(x = method, y = ep, fill = method, shape = rotation_type)) +
  geom_point(color = "black", size = 9, stroke = 1.2) +
  scale_fill_manual(name = "Method", values = method_colors, guide = "none") +
  scale_shape_manual(
    name   = "",
    values = c("Original" = 24, "Rotated" = 21),
    labels = c("Original" = "Original (rotate 0°)", "Rotated" = "Rotated (30°/60°/90°)")
  ) +
  scale_y_continuous(
    # breaks = c(0, 0.05, 0.10, 0.15, 0.2, 0.25)
    limits = c(0, 1),
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    labels = c("0", "0.25", "0.50", "0.75", "1.00")
  ) +
  facet_wrap(~ "EP") +
  labs(x = NULL, y = NULL) +
  shared_theme +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank()
  )

## ---- Assemble and save -----------------------------------------------------

row1 <- (
  pane_a + theme(plot.margin = margin(t = 5, b = 10, l = 10, r = 10))
) | (
  pane_b + theme(plot.margin = margin(t = 5, b = 10, l = 20,  r = 10))
)

row_ef <- (
  (pane_e + theme(plot.margin = margin(t = 10, b = 10, l = 25, r = 25))) |
  (pane_f + theme(plot.margin = margin(t = 10, b = 10, l =  25, r = 10)))
) + plot_layout(widths = c(1.2, 1))

combined <- row1 /
  (pane_c + theme(plot.margin = margin(t = 100, l = 100, b = 15))) /
  (pane_d + theme(plot.margin = margin(t = 20, l = 100, b = 15))) /
  row_ef +
  plot_layout(heights = c(1, 2, 1, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 28, face = "bold"))

ggsave("./figures/output/fig5/fig5.pdf", combined,
       width = 28, height = 34, device = "pdf")
