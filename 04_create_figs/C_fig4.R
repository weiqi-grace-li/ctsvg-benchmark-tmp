## --- gather data 
source("./scDesign_simulators/simulation_util.R")
library(PRROC)
library(pROC)


dim_table = read.csv("./figures/decomposition_dim.csv")

pr_results = data.frame(
  simulation_name = character(),
  seed            = integer(),
  test_method     = character(),
  experiment_name = character(),
  dataset = character(),
  AUPRC           = numeric(),
  ks = numeric(),
  AUROC = numeric(), 
  EP = numeric(),
  FDP  = numeric(), 
  Power = numeric(),
  FDP_nomarker = numeric(),
  pr_curve = I(list()),
  p_value = I(list()),
  stringsAsFactors = FALSE
)

aggregate_results = data.frame()
base_dirs = unique(dim_table$base_dir)

for (base_dir in base_dirs){
    # file_methods = unique(dim_table$file_method[which(dim_table$base_dir == base_dir)])
    file_methods = c("spvc")
    for (file_method in file_methods){
      all_files <- list.files(
        base_dir,
        pattern = paste0(".*", file_method, "_?([0-9]*)\\.(xlsx|csv)$"),
        full.names = FALSE
      )
      message(paste0(all_files, collapse = ", "))
      
      all_results <- do.call(
        rbind,
        lapply(all_files, function(file) {
          full_path <- file.path(base_dir, file)
          
          if (grepl("\\.csv$", file, ignore.case = TRUE)) {
            read.csv(full_path)
          } else if (grepl("\\.xlsx$", file, ignore.case = TRUE)) {
            readxl::read_excel(full_path)
          } else {
            stop("Unsupported file type: ", file)
          }
        })
      )
      
      temp_dim = dim_table[which(dim_table$base_dir == base_dir & dim_table$file_method == file_method),]
      
      for (i in 1:nrow(temp_dim)){
        idx = ((all_results$test_method == temp_dim[i, "test_method"]) & (all_results$simulation_name == temp_dim[i, "simulation_name"]))
        temp_sim = all_results[idx, ]
        seeds = unique(temp_sim$seed)
        for (seed in seeds){
          temp_sim_seed = temp_sim[which(temp_sim$seed == seed),]
          tags = sapply(strsplit(temp_sim_seed$gene_type, ","), function(x) "ctsvg" %in% x)
          markers = sapply(strsplit(temp_sim_seed$gene_type, ","), function(x) "marker" %in% x)
          fdr_nomarker = tags
          fdr_nomarker[which(!tags)] = !markers[which(!tags)]
          
          if (sum(tags) > 0){
            # calculate AUPRC
            temp_sim_seed$p_adj[which(temp_sim_seed$p_adj == -1)] = 1 # change untested to p value 1 
            # scores = 1-temp_sim_seed$p_adj # large score more significance
            p_adj_clean = pmax(temp_sim_seed$p_adj, .Machine$double.xmin)
            scores = -log(p_adj_clean)
            pr = pr.curve(scores.class0 = scores[tags], scores.class1 = scores[!tags], curve = TRUE)
            AUPRC = pr$auc.integral
            pr_curve = pr$curve
            
            # calculate AUROC 
            roc_obj = suppressMessages(pROC::roc(response=tags, predictor = scores))
            auc_value = pROC::auc(roc_obj) 
            # calculate EP 
            k = sum(tags)
            top_K = order(scores, decreasing = TRUE)[1:k]
            EP_value = sum(tags[top_K])/k
            
            # calculate FDP 
            rejected_adj = (temp_sim_seed$p_adj <=0.05)
            FDP =  sum(rejected_adj[!tags])/sum(rejected_adj)
            Power = sum(rejected_adj[tags])/sum(tags)
            
            # calculate FDP no marker 
            FDP_nomarker =  sum(rejected_adj[(!tags & fdr_nomarker)])/sum(rejected_adj[fdr_nomarker])
            
          } else {
            auc_value = NA
            AUPRC = NA
            pr_curve = NA 
            EP_value = NA
            FDP = NA
            Power = NA
            FDP_nomarker = NA
          }
          
          # calculate K-S distance, only do this if any genes are tested 
          p_raw = temp_sim_seed$p_value
          p_raw = p_raw[!tags & p_raw >= 0 & p_raw <= 1]
          
          if(length(p_raw) > 0){
            ks = ks.test(p_raw, "punif", 0, 1, exact = FALSE)
            ks = ks$statistic
          } else {
            ks = NA
          }
          
          temp_result = data.frame(
            simulation_name = temp_dim[i, "simulation_name"], 
            seed = seed, 
            test_method = temp_dim[i, "test_method"], 
            experiment_name = temp_dim[i, "experiment_name"], 
            dataset = temp_dim[i, "dataset"], 
            AUPRC = AUPRC, 
            ks = ks, 
            AUROC = auc_value, 
            EP = EP_value,
            FDP = FDP, 
            Power = Power, 
            FDP_nomarker = FDP_nomarker,
            pr_curve = I(list(pr_curve)), 
            p_value = I(list(p_raw))
          )
          
          pr_results = rbind(pr_results, temp_result)
          
        } 
      }
      # for each method + base_dir combo 
      temp_results = analyze_result(
        input_dir = temp_dim[i, "base_dir"],
        save_path = NA,
        save_sheet = "aggregate",
        analyze_method = temp_dim[i, "file_method"], 
        threshold = c(0.05)
      )
      
      # filter simulation names 
      temp_results = temp_results %>%
        filter(simulation_name %in% dim_table$simulation_name[which(dim_table$base_dir == base_dir & dim_table$file_method == file_method)])
      
      if (nrow(aggregate_results) == 0){
        aggregate_results = temp_results
      } else {
        aggregate_results = rbind(aggregate_results, temp_results)
      }
    }
}


library(dplyr)
library(stringr)
library(ggplot2)
library(tidyverse)

# inputs ---------------------------------------------------------------
# pr_results: data.frame with columns simulation_name, seed, test_method, AUPRC
# simulation_names: character vector of simulation names to keep

# 1) filter ------------------------------------------------------------
df <- pr_results 
aggregate_df <- aggregate_results %>%
  left_join(dim_table, by = c("simulation_name", "test_method"))

source("./utility/analyze_result.R")
write_sheet(
  df[,c("simulation_name", "dataset", "experiment_name", "seed", "test_method", "AUPRC", 
        "ks", "AUROC", "EP", "FDP", "Power", "FDP_nomarker")], 
  save_path = "./Results/AUPRC_KS.xlsx", 
  sheet_name = "aggregate"
)

write_sheet(
  aggregate_df, 
  save_path = "./Results/overall_decomposition_results.xlsx", 
  sheet_name = "aggregate"
)

# 4) aggregate mean AUPRC ---------------------------------------------
agg <- df %>%
  group_by(dataset, experiment_name, test_method) %>%
  summarise(
    AUPRC = mean(AUPRC, na.rm = TRUE), 
    KS = mean(ks, na.rm = TRUE),
    EP = mean(EP, na.rm = TRUE),
    FDP = mean(FDP, na.rm = TRUE),
    .groups = "drop")

ord = c(
          "Idealized",
          # "decompose_exp_-1",
          "decompose_0",
          "decompose_1", "decompose_2","decompose_1,2", "decompose_3", "decompose_1,2,3",
          "decompose_2,4", "decompose_1,2,4", "decompose_1,2,5",
          "decompose_1,2,3,4,5",
          "Realistic"
)
agg <- agg %>%
  mutate(experiment_name = fct_relevel(experiment_name, rev(ord))) %>%
  mutate(
    display_experience_name = fct_recode(experiment_name, 
       "  +All" = "decompose_1,2,3,4,5",
       "Baseline" = "decompose_0",
       "  +(i)" = "decompose_1", 
       "  +(ii)" = "decompose_2",
       "  +(i,ii,iii)" = "decompose_1,2", 
       "  +(iv)" = "decompose_3", 
       "  +(ii,v)" = "decompose_2,4", 
       "  +(i,ii,iii,v)" = "decompose_1,2,4",
       "  +(i,ii,iii,vi)" = "decompose_1,2,5", 
       "  +(i,ii,iii,iv)" = "decompose_1,2,3"                                          
    )
  )

agg <- agg %>%
  mutate(
    test_method = fct_recode(
      factor(test_method),
      "CELINA" = "celina", 
      "STANCE" = "stance", 
      "CTSV" = "ctsv",
      "C-SIDE" = "cside",
      "spVC"  = "spvc-original",
      "spVC-gam" = "spvc-gam"
    )
  ) %>%
  mutate(
    test_method = factor(
      test_method,
      levels = c("CELINA","STANCE","C-SIDE","CTSV","spVC", "spVC-gam")
    ),
    dataset_name = factor(dataset, levels = c("breast","ovarian","lymph"))
  )
paletteer::paletteer_d("tvthemes::AirNomads")
paletteer::paletteer_d("fishualize::Aluterus_scriptus")
# -------- plot (rows = experiment, columns grouped by test -> dataset) --------
# two-level column headers: first = test_name, second = dataset_name
## ---- let's see 
library(ggplot2)
library(ggh4x)
library(patchwork)

methods <- levels(agg$test_method)

## 1. Build one plot per method -------------------------------------------

# plot_list <- lapply(methods, function(m) {
#   ggplot(
#     subset(agg, test_method == m),
#     aes(
#       x = "EP",
#       y = display_experience_name,
#       fill = EP * 100,
#       label = sprintf("%.0f", EP * 100)
#     )
#   ) +
#     geom_tile(color = NA, width = 1.2, height = 0.9) +
#     geom_text(size = 8) +
#     scale_fill_gradient2(
#       low  = "#035AA6FF",
#       mid = "white",
#       high = "#FF9933FF",
#       midpoint = 40,
#       limits   = c(0, 100),
#       oob      = scales::squish
#     ) +
#     # # within a single method, you only need dataset_name on the x facet
#     # ggh4x::facet_nested(cols = vars(test_method, dataset_name),  # <-- key: nested facets
#     #                     scales = "free_x", space = "free_x", nest_line = TRUE) +
#     # --- NESTED FACETS: first layer = test_method, second = dataset_name ---
#     ggh4x::facet_nested(
#       cols      = vars(test_method, dataset_name),
#       scales    = "free_x",
#       space     = "free_x",
#       nest_line = FALSE,
#       strip     = ggh4x::strip_nested(
#         background_x = list(
#           element_rect(fill = "grey90", colour = "black"),  # layer 1 boxed
#           element_blank()                                   # layer 2 no box
#         ),
#         text_x = list(
#           element_text(size = 30, face = "bold"),           # layer 1 text
#           element_text(size = 26)                           # layer 2 text
#         ),
#         by_layer_x = TRUE
#       )
#     ) +
#     labs(
#       x = NULL, y = NULL,
#       fill = "Early Precision"
#       ,title = m
#     ) +
#     # theme_minimal(base_size = 12) +
#     theme(
#       strip.text.x     = element_text(size = 20, face = "bold",
#                                       margin = margin(t = 3, b = 3)),
#       strip.placement  = "outside",
#       panel.background = element_rect(fill = "white", color = NA, linewidth = 0.8),
#       panel.spacing.x  = unit(0.45, "lines"),
#       axis.text.x      = element_blank(),
#       axis.ticks.x     = element_blank(),
#       axis.text.y      = element_blank(),
#       axis.ticks.y     = element_blank(),
#       legend.position = "none",
#       plot.title = element_blank()
#       
#     )
# })
# 
# ## 2. Put y-axis labels back only on the first plot -----------------------
# # 
# # plot_list[[1]] <- plot_list[[1]] +
# #   theme(
# #     axis.text.y  = element_text(size = 24, face = "bold", hjust = 0, color = "black", margin = margin(r = 0)),
# #     axis.ticks.y = element_blank()
# #   )
# 
# ## 3. Insert spacers between plots and combine with shared legend ---------
# 
# spacer <- plot_spacer()
# 
# plots_with_gaps <- list()
# for (i in seq_along(plot_list)) {
#   plots_with_gaps[[length(plots_with_gaps) + 1]] <- plot_list[[i]]
#   if (i < length(plot_list)) {
#     plots_with_gaps[[length(plots_with_gaps) + 1]] <- spacer
#   }
# }
# 
# combined <- wrap_plots(
#   plots_with_gaps,
#   nrow   = 1,
#   widths = rep(c(1, 0), length.out = length(plots_with_gaps))  # 0.25 controls gap width
# ) +
#   plot_layout(guides = "collect") &                                # <-- single shared legend
#   theme(
#     legend.position = "none",
#     plot.margin = margin(l = 1, t = 5, r = 5, b = 5)
#   )

### make the idealized simulations all have the same results 
idealized = agg %>%
  filter(experiment_name == "Idealized")%>% 
  group_by(experiment_name, test_method, display_experience_name)%>%
  summarise(
    KS = mean(KS, na.rm = TRUE), 
    EP = mean(EP, na.rm = TRUE), 
    FDP = mean(FDP, na.rm = TRUE),
    AUPRC = mean(AUPRC, na.rm = TRUE),
    .group = "drop"
  )
idealized_avg = crossing(idealized, dataset_name = c("breast", "lymph", "ovarian"))
idealized_avg$dataset_name = factor(idealized_avg$dataset_name)
idealized_avg$dataset = as.character(idealized_avg$dataset_name)
others = agg %>% 
  filter(experiment_name != "Idealized") 
agg_avg = rbind(idealized_avg[, colnames(others)], others) %>%
  mutate(dataset_name = fct_relevel(dataset_name, c("breast", "ovarian", "lymph")))


combined = ggplot(agg_avg,
  aes(
    x = "EP",
    y = display_experience_name,
    fill = EP * 100,
    label = sprintf("%.0f", EP * 100)
  )
) +
  geom_tile(color = NA, width = 1.2, height = 0.9) +
  geom_text(size = 8) +
  scale_fill_gradient2(
    low  = "#035AA6FF",
    mid = "white",
    high = "#FF9933FF",
    midpoint = 40,
    limits   = c(0, 100),
    oob      = scales::squish
  ) +
  # # within a single method, you only need dataset_name on the x facet
  # ggh4x::facet_nested(cols = vars(test_method, dataset_name),  # <-- key: nested facets
  #                     scales = "free_x", space = "free_x", nest_line = TRUE) +
  # --- NESTED FACETS: first layer = test_method, second = dataset_name ---
  ggh4x::facet_nested(
    cols      = vars(test_method, dataset_name),
    scales    = "free_x",
    space     = "free_x",
    nest_line = FALSE,
    strip     = ggh4x::strip_nested(
      background_x = list(
        element_rect(fill = "grey90", colour = "black"),  # layer 1 boxed
        element_blank()                                   # layer 2 no box
      ),
      text_x = list(
        element_text(size = 30, face = "bold"),           # layer 1 text
        element_text(size = 26)                           # layer 2 text
      ),
      by_layer_x = TRUE
    )
  ) +
  labs(
    x = NULL, y = NULL
  ) +
  # theme_minimal(base_size = 12) +
  theme(
    strip.text.x     = element_text(size = 20, face = "bold",
                                    margin = margin(t = 3, b = 3)),
    strip.placement  = "outside",
    panel.background = element_rect(fill = "white", color = NA, linewidth = 0.8),
    panel.spacing.x  = unit(0.45, "lines"),
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    axis.text.y      = element_blank(),
    axis.ticks.y     = element_blank(),
    legend.position = "none",
    plot.title = element_blank()

  )
combined
ggsave("../Paper Writing/plots/figure4/spvc_only.pdf", combined, width = 10, height = 12.85, device = "pdf")

FDP = ggplot(agg_avg,
                  aes(
                    x = "FDP",
                    y = display_experience_name,
                    fill = FDP * 100,
                    label = sprintf("%.0f", FDP * 100)
                  )
) +
  geom_tile(color = NA, width = 1.2, height = 0.9) +
  geom_text(size = 8) +
  scale_fill_gradient2(
    low  = "#035AA6FF",
    mid = "white",
    high = "#FF9933FF",
    midpoint = 40,
    limits   = c(0, 100),
    oob      = scales::squish
  ) +
  # # within a single method, you only need dataset_name on the x facet
  # ggh4x::facet_nested(cols = vars(test_method, dataset_name),  # <-- key: nested facets
  #                     scales = "free_x", space = "free_x", nest_line = TRUE) +
  # --- NESTED FACETS: first layer = test_method, second = dataset_name ---
  ggh4x::facet_nested(
    cols      = vars(test_method, dataset_name),
    scales    = "free_x",
    space     = "free_x",
    nest_line = FALSE,
    strip     = ggh4x::strip_nested(
      background_x = list(
        element_rect(fill = "grey90", colour = "black"),  # layer 1 boxed
        element_blank()                                   # layer 2 no box
      ),
      text_x = list(
        element_text(size = 30, face = "bold"),           # layer 1 text
        element_text(size = 26)                           # layer 2 text
      ),
      by_layer_x = TRUE
    )
  ) +
  labs(
    x = NULL, y = NULL
  ) +
  # theme_minimal(base_size = 12) +
  theme(
    strip.text.x     = element_text(size = 20, face = "bold",
                                    margin = margin(t = 3, b = 3)),
    strip.placement  = "outside",
    panel.background = element_rect(fill = "white", color = NA, linewidth = 0.8),
    panel.spacing.x  = unit(0.45, "lines"),
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    axis.text.y      = element_text(size = 26, hjust = 0),
    axis.ticks.y     = element_blank(),
    legend.position = "none",
    plot.title = element_blank()
    
  )

FDP

ggsave("../Paper Writing/plots/supfigs/FDP.pdf", FDP, width = 28, height = 12.85, device = "pdf")

# 
# library(ggplot2)
# library(ggh4x)
# library(patchwork)
# 
# methods <- levels(agg$test_method)
# 
# ## 1. Build one plot per method -------------------------------------------
# 
# plot_list <- lapply(methods, function(m) {
#   ggplot(
#     subset(agg, test_method == m),
#     aes(
#       x = "FDP",
#       y = display_experience_name,
#       fill = FDP * 100,
#       label = sprintf("%.2f", FDP * 100)
#     )
#   ) +
#     geom_tile(color = NA, width = 1.2, height = 0.9) +
#     geom_text(size = 4) +
#     scale_fill_gradient2(
#       low  = "#035AA6FF",
#       high = "#FF9933FF",
#       midpoint = 50,
#       limits   = c(0, 100),
#       oob      = scales::squish
#     ) +
#     # within a single method, you only need dataset_name on the x facet
#     ggh4x::facet_nested(
#       cols      = vars(dataset_name),
#       scales    = "free_x",
#       space     = "free_x",
#       nest_line = TRUE
#     ) +
#     labs(
#       x = NULL, y = NULL,
#       fill = "Early Precision"
#       ,title = m
#     ) +
#     theme_minimal(base_size = 12) +
#     theme(
#       panel.grid       = element_blank(),
#       strip.text.x     = element_text(size = 12, face = "bold",
#                                       margin = margin(t = 3, b = 3)),
#       strip.placement  = "outside",
#       strip.background = element_rect(fill = NA, color = NA),
#       panel.spacing.x  = unit(0.45, "lines"),
#       axis.text.x      = element_blank(),
#       axis.ticks.x     = element_blank(),
#       ggh4x.facet.nestline = element_line(linewidth = 1.2, colour = "black"),
#       axis.text.y      = element_blank(),   # remove y for now; add on first plot
#       axis.ticks.y     = element_blank(),
#       plot.title = element_text(
#         hjust = 0.5, face = "bold", size = 16
#       )
#       
#     )
# })
# 
# ## 2. Put y-axis labels back only on the first plot -----------------------
# 
# plot_list[[1]] <- plot_list[[1]] +
#   theme(
#     axis.text.y  = element_text(size = 12, face = "bold"),
#     axis.ticks.y = element_blank()
#   )
# 
# ## 3. Insert spacers between plots and combine with shared legend ---------
# 
# spacer <- plot_spacer()
# 
# plots_with_gaps <- list()
# for (i in seq_along(plot_list)) {
#   plots_with_gaps[[length(plots_with_gaps) + 1]] <- plot_list[[i]]
#   if (i < length(plot_list)) {
#     plots_with_gaps[[length(plots_with_gaps) + 1]] <- spacer
#   }
# }
# 
# combined <- wrap_plots(
#   plots_with_gaps,
#   nrow   = 1,
#   widths = rep(c(1, 0.02), length.out = length(plots_with_gaps))  # 0.25 controls gap width
# ) +
#   plot_layout(guides = "collect") &                                # <-- single shared legend
#   theme(
#     legend.position = "none"
#   )
# combined
# 
# ggsave("../Paper Writing/plots/supfigs/FDP.pdf", width = 14.98, height = 7.52, device = "pdf")
# 
# 
# 
# # --- check some curve plot 
# stance_breast = df[which(df$dataset == "breast" & df$experiment_name == "decompose_1,2" &
#          df$seed == 42 & df$test_method == "stance")
#        ,"pr_curve"][[1]]
# stance_lymph = df[which(df$dataset == "lymph" & df$experiment_name == "decompose_1,2" &
#                           df$seed == 42 & df$test_method == "stance")
#                   ,"pr_curve"][[1]]
# stance_ovarian = df[which(df$dataset == "ovarian" & df$experiment_name == "decompose_1,2" &
#                           df$seed == 42 & df$test_method == "stance")
#                   ,"pr_curve"][[1]]
# 
# plot(stance_breast[,1], stance_breast[,2], type = "l", col = "red")
# lines(stance_lymph[,1], stance_lymph[,2], type = "l", col = "green")
# lines(stance_ovarian[,1], stance_ovarian[,2], type = "l", col = "orange")
# 
# ctsv_lymph = df[which(df$dataset == "lymph" & df$experiment_name == "scdesign3" &
#                            df$seed == 42 & df$test_method == "ctsv")
#                    ,"p_value"][[1]]
# hist(ctsv_lymph)
# 
# 
# celina = df[which(df$dataset == "breast" & df$experiment_name == "decompose_0" &
#                            df$seed == 42 & df$test_method == "celina")
#                    ,"p_value"][[1]]
