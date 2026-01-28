source("./scDesign_simulators/simulation_util.R")
source("./utility/analyze_result.R")

## --- gather data 
library(PRROC)
library(pROC)
library(ggplot2)
library(tidyr)
library(dplyr)
library(tidyverse)

# save_base_dirs = c(
#   "./scDesign_simulators/ovarian_small/ctsvg_result/scdesign3_seeds/",
#   "./scDesign_simulators/lymph_small/ctsvg_result/scdesign3_seeds/",
#   "./scDesign_simulators/breast_small/ctsvg_result/scdesign3_seeds/",
#   "./Results/pure_simulation/"
# )

save_base_dirs = c(
  "./Results/rotation/realistic/ovarian/",
  "./Results/rotation/realistic/lymph/",
  "./Results/rotation/realistic/breast/",
  "./Results/pure_simulation/"
)
analyze_method = c("stance", "celina", "spvc", "ctsv", "cside", "spvc-gam")
pr_results = data.frame(
  simulation_name = character(),
  seed            = integer(),
  test_method     = character(),
  AUPRC           = numeric(),
  ks = numeric(),
  AUROC = numeric(), 
  EP = numeric(),
  FDP  = numeric(), 
  Power = numeric(),
  FDP_nomarker = numeric(),
  # --- set of type 1 
  Type_1 = numeric(),
  Type_1_marker = numeric(),
  Type_1_null = numeric(),
  Type_1_other_ctsvg = numeric(),
  Type_1_other_marker = numeric(),
  Type_1_alt = numeric(),
  Type_1_marker_alt = numeric(),
  Type_1_null_alt = numeric(),
  Type_1_other_ctsvg_alt = numeric(),
  Type_1_other_marker_alt = numeric(),
  pr_curve = I(list()),
  p_value = I(list()),
  stringsAsFactors = FALSE
)

## for spvc, we only collect spvc-original and the experiment name with _spvcdebug.  

selected_sims = c(
  "breast_rotation_0", "ovarian_rotation_0", "lymph_rotation_0", # realistic results 
  "stance_simulator_1alt_0.7",
  "stance_simulator_1alt_1.5",
  "stance_simulator_alt_0.7",
  "stance_simulator_alt_1.5",
  "celina_simulator_alt_I_hotspot",
  "celina_simulator_alt_II_hotspot",
  "celina_simulator_alt_III_hotspot",
  "celina_simulator_alt_IV_hotspot",
  "celina_simulator_alt_I_streak",
  "celina_simulator_alt_II_streak",
  "celina_simulator_alt_III_streak",
  "celina_simulator_alt_IV_streak",
  "celina_simulator_alt_I_gradient",
  "celina_simulator_alt_II_gradient",
  "celina_simulator_alt_III_gradient",
  "celina_simulator_alt_IV_gradient",
  "celina_simulator_null_I",
  "celina_simulator_null_II",
  "celina_simulator_null_III",
  "celina_simulator_null_IV"
)

for (save_base_dir in save_base_dirs){
  for (method in analyze_method){
    all_files <- list.files(
      save_base_dir,
      pattern = paste0("^[^$]*", method, "_?([0-9]*)\\.(xlsx|csv)$"),
      full.names = FALSE
    )
    message(paste0(all_files, collapse = ", "))
    
    all_results <- do.call(
      rbind,
      lapply(all_files, function(file) {
        full_path <- file.path(save_base_dir, file)
        
        if (grepl("\\.csv$", file, ignore.case = TRUE)) {
          read.csv(full_path)
        } else if (grepl("\\.xlsx$", file, ignore.case = TRUE)) {
          readxl::read_excel(full_path)
        } else {
          stop("Unsupported file type: ", file)
        }
      })
    )
    
    file_methods = unique(all_results$test_method)
    for (file_method in file_methods){
      # if (file_method == "spvc-gam") next 
      method_result = all_results[which(all_results$test_method == file_method),]
      method_result = subset(method_result, method_result$simulation_name %in% selected_sims)
      simulations = unique(method_result$simulation_name)
      for (sim in simulations){
        temp_sim = subset(method_result, simulation_name == sim)
        seeds = unique(temp_sim$seed)
        for (seed in seeds){
          # get the table and the tag 
          temp_sim_seed = temp_sim[which(temp_sim$seed == seed),]
          tags = sapply(strsplit(temp_sim_seed$gene_type, ","), function(x) "ctsvg" %in% x)
          markers = sapply(strsplit(temp_sim_seed$gene_type, ","), function(x) "marker" %in% x)
          other_marker = (sapply(strsplit(temp_sim_seed$gene_type, ","), function(x) "other_marker" %in% x) & !tags & !markers)
          other_ctsvg = (sapply(strsplit(temp_sim_seed$gene_type, ","), function(x) "other_ctsvg" %in% x) & !tags & !markers & !other_marker)
          null_genes = (!tags & !markers & !other_marker & !other_ctsvg)
          
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
          
         
          
          # Type 1 
          p_raw = temp_sim_seed$p_value
          
          Type_1 = sum(((p_raw < 0.05) & (p_raw >= 0) & (!tags)))/sum(!tags & (p_raw >= 0))
          Type_1_marker = sum(((p_raw < 0.05) & (p_raw >= 0) & (markers)))/sum(markers & (p_raw >= 0))
          Type_1_null = sum(((p_raw < 0.05) & (p_raw >= 0) & (null_genes)))/sum(null_genes & (p_raw >= 0))
          Type_1_other_ctsvg = sum(((p_raw < 0.05) & (p_raw >= 0) & (other_ctsvg)))/sum(other_ctsvg & (p_raw >= 0))
          Type_1_other_marker = sum(((p_raw < 0.05) & (p_raw >= 0) & (other_marker)))/sum(other_marker & (p_raw >= 0))
          
          Type_1_alt = sum(((p_raw < 0.05) & (p_raw >= 0) & (!tags)))/sum(!tags)
          Type_1_marker_alt = sum(((p_raw < 0.05) & (p_raw >= 0) & (markers)))/sum(markers)
          Type_1_null_alt = sum(((p_raw < 0.05) & (p_raw >= 0) & (null_genes)))/sum(null_genes)
          Type_1_other_ctsvg_alt = sum(((p_raw < 0.05) & (p_raw >= 0) & (other_ctsvg)))/sum(other_ctsvg)
          Type_1_other_marker_alt = sum(((p_raw < 0.05) & (p_raw >= 0) & (other_marker)))/sum(other_marker)
          
          # KS distance 
          p_raw = p_raw[!tags & p_raw >= 0 & p_raw <= 1]
          if(length(p_raw) > 0){
            ks = ks.test(p_raw, "punif", 0, 1, exact = FALSE)
            ks = ks$statistic
          } else {
            ks = NA
          }
        
          temp_result = data.frame(
            simulation_name = sim, 
            seed = seed, 
            test_method = file_method, 
            AUPRC = AUPRC, 
            ks = ks, 
            AUROC = auc_value, 
            EP = EP_value,
            FDP = FDP, 
            Power = Power, 
            FDP_nomarker = FDP_nomarker,
            Type_1 = Type_1,
            Type_1_marker = Type_1_marker, 
            Type_1_null = Type_1_null, 
            Type_1_other_ctsvg = Type_1_other_ctsvg, 
            Type_1_other_marker = Type_1_other_marker, 
            Type_1_alt = Type_1_alt,
            Type_1_marker_alt = Type_1_marker_alt,
            Type_1_null_alt = Type_1_null_alt,
            Type_1_other_ctsvg_alt = Type_1_other_ctsvg_alt,
            Type_1_other_marker_alt = Type_1_other_marker_alt,
            pr_curve = I(list(pr_curve)), 
            p_value = I(list(p_raw))
          )
          
          pr_results = rbind(pr_results, temp_result)
        } # seed
      } # simulation 
    } # method
  } # method_file 
}

write_sheet(
  pr_results[,c("simulation_name", "seed", "test_method", "AUPRC", "ks", "AUROC", "EP", "FDP", "FDP_nomarker", "Power", "Type_1", "Type_1_marker", 
                "Type_1_alt", "Type_1_marker_alt")],
  save_path = "./Results/overall_plot.xlsx",
  sheet_name = "aggregate"
)

str(pr_results)
plot_sim_names = selected_sims

test_methods = c(
  "spvc-original", 
  "spvc-gam",
  "stance", 
  "celina", 
  "cside",
  "ctsv"
)

real_names = c(
  "breast_rotation_0", "ovarian_rotation_0", "lymph_rotation_0"
)

lev_metric = c(
  "EP","AUPRC","FDP","KS",
  "Type_1", "Type_1_marker", "Type_1_other_marker", "Type_1_other_ctsvg", "Type_1_null", 
  "Type_1_alt", "Type_1_marker_alt", "Type_1_other_marker_alt", "Type_1_other_ctsvg_alt", "Type_1_null_alt"
)

df <- pr_results %>%
  filter(simulation_name %in% plot_sim_names) %>%
  filter(test_method %in% test_methods) %>%
  mutate(
    experiment_name = if_else(simulation_name %in% real_names, "real", "pure"),
    test_method  = factor(test_method,  levels = c("celina","stance", "cside","ctsv", "spvc-gam", "spvc-original"))
  )  %>%
  group_by(simulation_name, test_method, experiment_name) %>%
  summarise(
    AUPRC = mean(AUPRC, na.rm = TRUE), 
    KS = mean(ks, na.rm = TRUE),
    EP = mean(EP, na.rm = TRUE),
    FDP = mean(FDP, na.rm = TRUE), 
    Type_1 = mean(Type_1, na.rm = TRUE),
    Type_1_marker = mean(Type_1_marker, na.rm = TRUE),
    Type_1_null = mean(Type_1_null, na.rm = TRUE),
    Type_1_other_ctsvg = mean(Type_1_other_ctsvg, na.rm = TRUE),
    Type_1_other_marker = mean(Type_1_other_marker, na.rm = TRUE),
    Type_1_alt = mean(Type_1_alt, na.rm = TRUE),
    Type_1_marker_alt = mean(Type_1_marker_alt, na.rm = TRUE),
    Type_1_null_alt = mean(Type_1_null_alt, na.rm = TRUE),
    Type_1_other_ctsvg_alt = mean(Type_1_other_ctsvg_alt, na.rm = TRUE),
    Type_1_other_marker_alt = mean(Type_1_other_marker_alt, na.rm = TRUE),
    .groups = "drop")%>%
  pivot_longer(cols = c(
    EP, AUPRC, KS, FDP,
    Type_1, Type_1_marker, Type_1_other_marker, Type_1_other_ctsvg, Type_1_null,
    Type_1_alt, Type_1_marker_alt, Type_1_other_marker_alt, Type_1_other_ctsvg_alt, Type_1_null_alt
  ),
 names_to = "metric", values_to = "value") %>%
  mutate(
    metric = fct_relevel(metric, lev_metric)
  ) %>% 
  mutate(
    metric_display = fct_recode(
      metric, 
      "Type I Error Only Tested" = "Type_1", 
      "Type I Error Only Tested (marker)" = "Type_1_marker", 
      "Type I Error Only Tested (null)" = "Type_1_null", 
      "Type I Error Only Tested (other ctsvg)" = "Type_1_other_ctsvg", 
      "Type I Error Only Tested (other marker)" = "Type_1_other_marker",
      "Type I Error" = "Type_1_alt", 
      "Type I Error (marker)" = "Type_1_marker_alt", 
      "Type I Error (null)" = "Type_1_null_alt", 
      "Type I Error (other ctsvg)" = "Type_1_other_ctsvg_alt", 
      "Type I Error (other marker)" = "Type_1_other_marker_alt"
    ),
    test_method = fct_recode(test_method, "spVC" = "spvc-original", "STANCE" = "stance",  "CELINA" = "celina", "C-SIDE" = "cside",  "CTSV" = "ctsv")
  )

real = subset(df, experiment_name == "real")
pure = subset(df, experiment_name == "pure")

lev_metric = c(
  "EP","AUPRC","FDP","KS","Type I Error", "Type I Error (marker)", 
  "Type I Error (null)", "Type I Error (other_ctsvg)", 
  "Type I Error (other_marker)"
)

overall_metric = c(
  "EP","AUPRC","FDP","KS","Type_1_alt"
)

# (Optional) lock x order too, if you have a preferred method order:
# lev_method <- c("CSIDE","CELINA","CTSV","spVC","Method5")
# pure_long$test_method <- factor(pure_long$test_method, levels = lev_method)
# real_long$test_method <- factor(real_long$test_method, levels = lev_method)

methods <- unique(pure$test_method)

ref_lines <- expand.grid(
  metric_display = c("FDP", "Type I Error"),   # the panels that need the line
  test_method = methods           # one line per method
) %>%
  mutate(
    metric_display = factor(metric_display, levels = lev_metric),
    yint   = 0.05,
    ylab   = 0.025,
    label  = "target 0.05"
  )

pure = pure %>%
  filter(!is.na(value)) %>%
  filter(metric %in% overall_metric)
real = real %>% 
  filter(!is.na(value))%>%
  filter(metric %in% overall_metric)

pure$dataset_type <- ifelse(pure$simulation_name %in% c("stance_simulator_alt_0.7","stance_simulator_alt_1.5"),
                            "Unbalanced (Scenario 1)", "Balanced (Scenario 2-6)")
real$dataset_type <- "Realistic data"   # ensure the column exists in `real`

method_levels <- unique(pure$test_method)  

library(paletteer)

pal <- paletteer_d("tvthemes::AirNomads")

overall = ggplot(pure, aes(x = test_method, y = value, group = test_method)) +
  geom_boxplot(
    aes(fill = test_method),
    width = 0.55, 
    alpha = 0.7, 
    outlier.colour = NA, 
    outlier.fill = NA, 
    outlier.shape = 21, 
    outlier.size = 2.5
  ) +
  facet_wrap(~ metric_display, nrow = 1, scales = "fixed") +
  # geom_point(
  #   data = subset(pure, dataset_type == "Unbalanced (Scenario 1)"),
  #   aes(x = test_method, y = value, fill = test_method, color = test_method, shape = dataset_type),
  #   # color = "black",         # outline
  #   size = 5,stroke = 1.2 
  # )+
  geom_point(
    data = real,
    aes(x = test_method, y = value, fill = test_method, shape = dataset_type),
    color = "black",         # outline
    size = 5,stroke = 1.2 
  )+
  geom_hline(
    data = ref_lines, aes(yintercept = yint),
    linetype = "dashed", color = "grey20", size = 0.6
  ) +
  geom_text(
    data = ref_lines, aes(y = ylab, label = label),
    x = 3.5, y = -0.04, hjust = 0.5, color = "grey20", size = 7
  ) +
  coord_cartesian(ylim = c(-0.05, 1.05))+
  # now shape has mapped values; specify keys
  scale_shape_manual(
    name   = "",
    values = c(
      "Unbalanced (Scenario 1)" = 8,
      "Realistic data" = 24
    ), 
    label = c(
      "Unbalanced (Scenario 1)" = "Unbalanced (Scenario 1)",
      "Realistic data" = "Realistic"
    )
  ) +
  scale_colour_manual(
    values     = pal,
    aesthetics = c("colour", "fill"),
    name       = "",
    limits     = method_levels,
    labels     = method_labels
  )+
  guides(
    colour = "none",
    fill   = guide_legend(order = 1),  # no separate fill legend
    shape  = guide_legend(order = 2, override.aes = list(size = 6))
  )+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    legend.direction = "horizontal",
    legend.title = element_text(size = 20, face = "bold"),
    legend.text  = element_text(size = 20),
    legend.key.size = unit(3, "lines"),
    strip.text = element_text(size = 24, face = "bold"), # for fecet heading 
    panel.background = element_rect(fill = "white", color = "black", linewidth = 0.8),
    strip.background = element_rect(fill = "grey80", color = "black"),
    axis.title = element_blank(),
    axis.text.y = element_text(size = 22),
    axis.ticks.y = element_line(size = 1.2), 
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank()
  )

overall
# ggsave("../Paper Writing/plots/figure2/all_metrics_with_spvcgam.pdf", overall, width = 20.1, height = 6.6, device = "pdf")
# ggsave("../Paper Writing/plots/figure2/all_metrics.pdf", overall, width = 20.1, height = 6.6, device = "pdf")


## trying to get to subfigures 
library(dplyr)
library(ggplot2)
library(patchwork)

# Balanced scenarios
sim_dim = read.csv("./utility/simulator_setup.csv")
sim_dim <- sim_dim %>%
  mutate(
    simulation_name = paste0(
      sim_name,
      ifelse(!is.na(phi)   & phi   != "", paste0("_", phi), ""),
      ifelse(!is.na(scene) & scene != "", paste0("_", scene), ""),
      ifelse(!is.na(pattern) & pattern != "", paste0("_", pattern), "")
    )
  )

ep <- pure %>%
  filter(metric == "EP") %>%
  left_join(sim_dim[, c("simulation_name", "pattern_plot", "scenario_plot")], by ="simulation_name")


ep$scenario_display = paste0("Scenario ", ep$scenario_plot)
ep = ep %>% 
  mutate(
    pattern_plot = factor(pattern_plot, levels = c("hotspot", "streak", "gradient"))
  )

EP = ggplot(subset(ep, dataset_type != "Unbalanced (Scenario 1)"), aes(x = test_method, y = value, group = test_method)) +
  geom_boxplot(
    aes(fill = test_method),
    width = 0.55,
    alpha = 0.8,
    outlier.colour = NA,
    outlier.fill = NA,
    outlier.shape = 21,
    outlier.size = 2.5
  ) +
  facet_wrap(~ pattern_plot, nrow = 1, scales = "fixed") +
  geom_point(
    data = subset(ep, dataset_type == "Unbalanced (Scenario 1)"),
    aes(x = test_method, y = value, fill = test_method, color = test_method, shape = dataset_type),
    # color = "black",         # outline
    size = 5,stroke = 1.2 
  )+
  coord_cartesian(ylim = c(-0.05, 1.05))+
  ggtitle("EP")+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    legend.direction = "horizontal",
    legend.title = element_text(size = 20, face = "bold"),
    legend.text  = element_text(size = 20),
    legend.key.size = unit(3, "lines"),
    strip.text = element_text(size = 24, face = "bold"), # for fecet heading 
    panel.background = element_rect(fill = "white", color = "black", linewidth = 0.8),
    strip.background = element_rect(fill = "grey80", color = "black"),
    axis.title = element_blank(),
    axis.text.y = element_text(size = 22),
    axis.ticks.y = element_line(size = 1.2), 
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(), 
    plot.title = element_text(size = 24, hjust = 0.5, face = "bold")
  )+
  scale_shape_manual(
    name   = "",
    values = c(
      "Unbalanced (Scenario 1)" = 8,
      "Realistic data" = 24
    ), 
    label = c(
      "Unbalanced (Scenario 1)" = "Unbalanced (Scenario 1)",
      "Realistic data" = "Realistic"
    )
  ) +
  scale_colour_manual(
    values     = pal,
    aesthetics = c("colour", "fill"),
    name       = "",
    limits     = method_levels,
    labels     = method_labels
  )+
  guides(
    colour =guide_legend(order = 1) ,
    fill   = "none",  # no separate fill legend
    shape  = guide_legend(order = 2, override.aes = list(size = 6))
  )


EP

ggsave("../Paper Writing/plots/figure2/ep.pdf", EP, width = 9, height = 6.5, device = "pdf")

fdp <- pure %>%
  filter(metric == "FDP") %>%
  left_join(sim_dim[, c("simulation_name", "pattern_plot", "scenario_plot")], by ="simulation_name")


fdp$scenario_display = paste0("Scenario ", fdp$scenario_plot)

fdp = fdp %>% 
  mutate(
    scenario_display = factor(scenario_display)
  ) %>% 
  mutate(
    scenario_display = fct_recode(
      scenario_display, 
      "Scenario 1" = "Scenario 1",
      "Scenario 3-4" = "Scenario 3",
      "Scenario 3-4" = "Scenario 4",
      "Scenario 5-6" = "Scenario 5",
      "Scenario 5-6" = "Scenario 6",
    )
  )

FDP = ggplot(subset(fdp, scenario_display != "Scenario 1"), aes(x = test_method, y = value, group = test_method)) +
  geom_boxplot(
    aes(fill = test_method),
    width = 0.55,
    alpha = 0.8,
    outlier.colour = NA,
    outlier.fill = NA,
    outlier.shape = 21,
    outlier.size = 2.5
  ) +
  geom_point(
    data = subset(fdp, dataset_type == "Unbalanced (Scenario 1)"),
    aes(x = test_method, y = value, fill = test_method, color = test_method, shape = dataset_type),
    # color = "black",         # outline
    size = 5,stroke = 1.2 
  )+
  geom_hline(
    data = ref_lines, aes(yintercept = yint),
    linetype = "dashed", color = "grey20", size = 0.6
  ) +
  geom_text(
    data = ref_lines, aes(y = ylab, label = label),
    x = 3.5, y = -0.04,hjust = 0.5, color = "grey20", size = 7
  ) +
  facet_wrap(~ scenario_display, nrow = 1, scales = "fixed") +
  coord_cartesian(ylim = c(-0.05, 1.05))+
  ggtitle("FDP")+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    legend.direction = "horizontal",
    legend.title = element_text(size = 20, face = "bold"),
    legend.text  = element_text(size = 20),
    legend.key.size = unit(3, "lines"),
    strip.text = element_text(size = 24, face = "bold"), # for fecet heading 
    panel.background = element_rect(fill = "white", color = "black", linewidth = 0.8),
    strip.background = element_rect(fill = "grey80", color = "black"),
    axis.title = element_blank(),
    axis.text.y = element_text(size = 22),
    axis.ticks.y = element_line(size = 1.2), 
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(), 
    plot.title = element_text(size = 24, hjust = 0.5, face = "bold")
  )+
  scale_shape_manual(
    name   = "",
    values = c(
      "Unbalanced (Scenario 1)" = 8,
      "Realistic data" = 24
    ), 
    label = c(
      "Unbalanced (Scenario 1)" = "Unbalanced (Scenario 1)",
      "Realistic data" = "Realistic"
    )
  ) +
  scale_colour_manual(
    values     = pal,
    aesthetics = c("colour", "fill"),
    name       = "",
    limits     = method_levels,
    labels     = method_labels
  )+
  guides(
    colour =guide_legend(order = 1) ,
    fill   = "none",  # no separate fill legend
    shape  = guide_legend(order = 2, override.aes = list(size = 6))
  )

FDP

ggsave("../Paper Writing/plots/figure2/fdp.pdf", FDP, width = 9, height = 6.5, device = "pdf")

p_all <- (overall / (EP + FDP)) +
  plot_annotation(tag_levels = "a")

panel = p_all & theme(
  plot.tag = element_text(size = 28, face = "bold"),
)

panel 

ggsave("../Paper Writing/plots/figure2/overall.pdf", panel, width = 20, height = 15, device = "pdf")



### type 1 error 
real = subset(df, experiment_name == "real") 
pure = subset(df, experiment_name == "pure")
Type_1_metric = c(
  "Type_1", "Type_1_marker", "Type_1_other_ctsvg", "Type_1_other_marker", "Type_1_null"
)
pure = pure %>%
  filter(!is.na(value)) %>%
  filter(metric %in% Type_1_metric) %>%
  mutate(
    metric_display = fct_recode(
      metric, 
      "Overall" ="Type_1" ,
      "Marker" ="Type_1_marker" ,
      "Other ctSVGs" ="Type_1_other_ctsvg" ,
      "Other Marker" ="Type_1_other_marker",
      "Null Genes" ="Type_1_null" 
    )
  )
real = real %>% 
  filter(!is.na(value))%>%
  filter(metric %in% Type_1_metric) %>%
  mutate(
    metric_display = fct_recode(
      metric, 
      "Overall" ="Type_1" ,
      "Marker" ="Type_1_marker" ,
      "Other ctSVGs" ="Type_1_other_ctsvg" ,
      "Other Marker" ="Type_1_other_marker",
      "Null Genes" ="Type_1_null" 
    )
  )

pure$dataset_type <- ifelse(pure$simulation_name %in% c("stance_simulator_alt_0.7","stance_simulator_alt_1.5"),
                            "Unbalanced (Scenario 1)", "Balanced (Scenario 2-6)")
real$dataset_type <- "Realistic data"   # ensure the column exists in `real`

method_levels <- unique(pure$test_method)  

method_labels <- c("spVC" = "spVC Balanced (Scenario 2–6)")
library(paletteer)

pal <- paletteer_d("tvthemes::AirNomads")

type_1 = ggplot(subset(pure, dataset_type == "Balanced (Scenario 2-6)"), aes(x = test_method, y = value, group = test_method)) +
  geom_boxplot(
    aes(fill = test_method),
    width = 0.55, 
    alpha = 0.7, 
    outlier.colour = NA, 
    outlier.fill = NA, 
    outlier.shape = 21, 
    outlier.size = 2.5
  ) +
  facet_wrap(~ metric_display, nrow = 1, scales = "fixed") +
  geom_point(
    data = subset(pure, dataset_type == "Unbalanced (Scenario 1)"),
    aes(x = test_method, y = value, fill = test_method, color = test_method, shape = dataset_type),
    # color = "black",         # outline
    size = 5,stroke = 1.2 
  )+
  geom_point(
    data = real,
    aes(x = test_method, y = value, fill = test_method, shape = dataset_type),
    color = "black",         # outline
    size = 5,stroke = 1.2 
  )+
  geom_hline(
    data = ref_lines, aes(yintercept = yint),
    linetype = "dashed", color = "grey20", size = 0.6
  ) +
  geom_text(
    data = ref_lines, aes(y = ylab, label = label),
    x = 3.5, hjust = 0, color = "grey20", size = 7
  ) +
  coord_cartesian(ylim = c(-0.05, 1.05))+
  # now shape has mapped values; specify keys
  scale_shape_manual(
    name   = "",
    values = c(
      "Unbalanced (Scenario 1)" = 8,
      "Realistic data" = 24
    ), 
    label = c(
      "Unbalanced (Scenario 1)" = "Unbalanced (Scenario 1)",
      "Realistic data" = "Realistic"
    )
  ) +
  scale_colour_manual(
    values     = pal,
    aesthetics = c("colour", "fill"),
    name       = "",
    limits     = method_levels,
    labels     = method_labels
  )+
  guides(
    colour = "none",
    fill   = guide_legend(order = 1),  # no separate fill legend
    shape  = guide_legend(order = 2, override.aes = list(size = 6))
  )+
  ggtitle("Type I Error (Only Tested Genes)")+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.title = element_text(size = 20, face = "bold"),
    legend.text  = element_text(size = 20),
    legend.key.size = unit(3, "lines"),
    strip.text = element_text(size = 24, face = "bold"), # for fecet heading 
    panel.background = element_rect(fill = "white", color = "black", linewidth = 0.8),
    strip.background = element_rect(fill = "grey80", color = "black"),
    axis.title = element_blank(),
    axis.text.y = element_text(size = 22),
    axis.ticks.y = element_line(size = 1.2), 
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(), 
    plot.title = element_text(size = 24, face = "bold", hjust = 0.5)
  )

type_1


## type 1 error alt
real = subset(df, experiment_name == "real") 
pure = subset(df, experiment_name == "pure")
Type_1_metric_alt = c(
  "Type_1_alt", "Type_1_marker_alt", "Type_1_other_ctsvg_alt", "Type_1_other_marker_alt", "Type_1_null_alt"
)
pure = pure %>%
  filter(!is.na(value)) %>%
  filter(metric %in% Type_1_metric_alt) %>%
  mutate(
    metric_display = fct_recode(
      metric, 
      "Overall" ="Type_1_alt" ,
      "Marker" ="Type_1_marker_alt" ,
      "Other ctSVG" ="Type_1_other_ctsvg_alt" ,
      "Other Marker" ="Type_1_other_marker_alt",
      "Null Gene" ="Type_1_null_alt" 
    )
  )
real = real %>% 
  filter(!is.na(value))%>%
  filter(metric %in% Type_1_metric_alt) %>%
  mutate(
    metric_display = fct_recode(
      metric, 
      "Overall" ="Type_1_alt" ,
      "Marker" ="Type_1_marker_alt" ,
      "Other ctSVG" ="Type_1_other_ctsvg_alt" ,
      "Other Marker" ="Type_1_other_marker_alt",
      "Null Gene" ="Type_1_null_alt" 
    )
  )

pure$dataset_type <- ifelse(pure$simulation_name %in% c("stance_simulator_alt_0.7","stance_simulator_alt_1.5"),
                            "Unbalanced (Scenario 1)", "Balanced (Scenario 2-6)")
real$dataset_type <- "Realistic data"   # ensure the column exists in `real`

method_levels <- unique(pure$test_method)  

method_labels <- c("spVC" = "spVC Balanced (Scenario 2–6)")
library(paletteer)

pal <- paletteer_d("tvthemes::AirNomads")

ref_lines <- expand.grid(
  metric_display = c("Overall", "Marker", "Other Marker", "Other ctSVG", "Null Gene"),   # the panels that need the line
  test_method = method_levels           # one line per method
) %>%
  mutate(
    metric_display = factor(metric_display),
    yint   = 0.05,
    ylab   = 0.025,
    label  = "target 0.05"
  )


type_1_alt = ggplot(pure, aes(x = test_method, y = value, group = test_method)) +
  geom_boxplot(
    aes(fill = test_method),
    width = 0.55, 
    alpha = 0.7, 
    outlier.colour = NA, 
    outlier.fill = NA, 
    outlier.shape = 21, 
    outlier.size = 2.5
  ) +
  facet_wrap(~ metric_display, nrow = 1, scales = "fixed") +
  geom_point(
    data = subset(pure, dataset_type == "Unbalanced (Scenario 1)"),
    aes(x = test_method, y = value, fill = test_method, color = test_method, shape = dataset_type),
    # color = "black",         # outline
    size = 5,stroke = 1.2 
  )+
  geom_point(
    data = real,
    aes(x = test_method, y = value, fill = test_method, shape = dataset_type),
    color = "black",         # outline
    size = 5,stroke = 1.2 
  )+
  geom_hline(
    data = ref_lines, aes(yintercept = yint),
    linetype = "dashed", color = "grey20", size = 0.6
  ) +
  geom_text(
    data = ref_lines, aes(y = ylab, label = label),
    x = 3.5, y = -0.04, hjust = 0.5, color = "grey20", size = 7
  ) +
  coord_cartesian(ylim = c(-0.05, 1.05))+
  # now shape has mapped values; specify keys
  scale_shape_manual(
    name   = "",
    values = c(
      "Unbalanced (Scenario 1)" = 8,
      "Realistic data" = 24
    ), 
    label = c(
      "Unbalanced (Scenario 1)" = "Unbalanced (Scenario 1)",
      "Realistic data" = "Realistic"
    )
  ) +
  scale_colour_manual(
    values     = pal,
    aesthetics = c("colour", "fill"),
    name       = "",
    limits     = method_levels,
    labels     = method_labels
  )+
  guides(
    colour = "none",
    fill   = guide_legend(order = 1),  # no separate fill legend
    shape  = guide_legend(order = 2, override.aes = list(size = 6))
  )+
  ggtitle("Type I Error")+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    legend.direction = "horizontal",
    legend.title = element_text(size = 20, face = "bold"),
    legend.text  = element_text(size = 20),
    legend.key.size = unit(3, "lines"),
    strip.text = element_text(size = 24, face = "bold"), # for fecet heading 
    panel.background = element_rect(fill = "white", color = "black", linewidth = 0.8),
    strip.background = element_rect(fill = "grey80", color = "black"),
    axis.title = element_blank(),
    axis.text.y = element_text(size = 22),
    axis.ticks.y = element_line(size = 1.2), 
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(),
    plot.title = element_text(size = 24, face = "bold", hjust = 0.5)
  )

type_1_alt

ggsave("../Paper Writing/plots/figure2/type1.pdf", type_1_alt, width = 20.1, height = 6.6, device = "pdf")

type_1_all <- (type_1 / type_1_alt) +
  plot_annotation(tag_levels = "a")

type_1_panel = type_1_all & theme(
  plot.tag = element_text(size = 28, face = "bold"),
)

type_1_panel
  

ggsave("../Paper Writing/plots/figure2/type1.pdf", type_1_panel, width = 20, height = 15, device = "pdf")


### supplementary CSIDE p value distribution 

cside = subset(pr_results, test_method == "cside")
selected_scenarios = c(
  "stance_simulator_alt_1.5", "stance_simulator_1alt_1.5", 
  "celina_simulator_alt_I_hotspot", "celina_simulator_alt_II_streak", "celina_simulator_alt_II_gradient", "celina_simulator_null_IV", 
  "ovarian_scdesign3_remove300"
)
p_dist_cside <- cside %>%
  filter(simulation_name %in% selected_scenarios, seed == 42) %>%
  mutate(scenario_display = fct_recode(
    factor(simulation_name), 
     "Idealized Scenario 1" = "stance_simulator_alt_1.5", 
     "Idealized Scenario 2" = "stance_simulator_1alt_1.5", 
     "Idealized Scenario 3" = "celina_simulator_alt_I_hotspot", 
     "Idealized Scenario 4" = "celina_simulator_alt_II_streak", 
     "Idealized Scenario 5" = "celina_simulator_alt_II_gradient", 
     "Idealized Scenario 6" = "celina_simulator_null_IV", 
     "Realistic" = "ovarian_scdesign3_remove300"
  )) %>%   # rename or label scenario
  select(scenario_display, p_value) %>%
  unnest(p_value)


bw <- 0.05
breaks <- seq(0, 1, by = bw)

p_binned <- p_dist_cside %>%
  mutate(bin = cut(p_value, breaks = breaks, include.lowest = TRUE, right = TRUE)) %>%
  count(scenario_display, bin, name = "n") %>%
  group_by(scenario_display) %>%
  mutate(prop = n / sum(n),
         xmid = (as.numeric(bin) - 1) * bw + bw/2) %>%
  ungroup() %>%
  mutate(
    scenario_display = fct_relevel(
      scenario_display, c(
        "Idealized Scenario 1", "Idealized Scenario 2", "Idealized Scenario 3", 
        "Idealized Scenario 4", "Idealized Scenario 5", "Idealized Scenario 6", 
        "Realistic"
      )  
    )
  )


cside_p = ggplot(p_binned, aes(x = xmid, y = prop)) +
  geom_col(width = bw, fill = "#FF9933FF", color = "black") +
  facet_wrap(~ scenario_display, nrow = 3, scales = "fixed") +
  scale_x_continuous(breaks = seq(0, 1, by = 0.05)) +
  scale_y_continuous(breaks = seq(0, 0.5, by = 0.05)) +
  labs(x = "p-value", y = "Proportion") +
  theme_minimal(base_size = 14) + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    legend.direction = "horizontal",
    legend.title = element_text(size = 20, face = "bold"),
    legend.text  = element_text(size = 20),
    legend.key.size = unit(3, "lines"),
    strip.text = element_text(size = 24, face = "bold"), # for fecet heading 
    panel.background = element_rect(fill = "white", color = "black", linewidth = 0.8),
    strip.background = element_rect(fill = "grey80", color = "black"),
    axis.title = element_blank(),
    axis.text.y = element_text(size = 22),
    axis.ticks.y = element_line(size = 1.2), 
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(),
    plot.title = element_text(size = 24, face = "bold", hjust = 0.5)
  )

ggsave("../Paper Writing/plots/figure2/cside_p_dist.pdf", cside_p, width = 18, height = 20, device = "pdf")
