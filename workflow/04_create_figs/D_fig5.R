runtime_files = c(
  "./Results/scalability/runtime.csv"
  ,"./Results/scalability/runtime_stance.csv"
  ,"./Results/scalability/runtime_celina.csv"
  ,"./Results/scalability/runtime_cside.csv"
  ,"./Results/scalability/runtime_spvc.csv"
  ,"./Results/scalability/runtime_ctsv.csv"
)

pal <- paletteer::paletteer_d("tvthemes::AirNomads")

# Example: change the 3rd color to something else
pal[3] <- "#FFD700"

runtime_results = data.frame()

for (file in runtime_files){
  if (nrow(runtime_results) == 0){
    runtime_results = read.csv(file)
  } else {
    runtime_results = rbind(runtime_results, read.csv(file))
  }
}
str(runtime_results)

runtime_df <- runtime_results %>%
  mutate(
    n_side  = as.numeric(sub(".*n_side:\\s*(\\d+).*", "\\1", dataset)),
    n_genes = ifelse(
      grepl("n_gene", dataset),
      as.numeric(sub(".*n_gene:\\s*(\\d+).*", "\\1", dataset)),
      100
    ),
    elapsed = elapsed/60, 
    peak_mb = peak_mb/1024,
    method  = factor(method,  levels = c("stance","celina", "cside","ctsv", "spvc-gam", "spvc-original"))
  ) %>%
  mutate(
    type = ifelse(
      n_side == 30, "Varying Genes (900 Spots)", "Varying Spots (100 Genes)"
    )
  )%>%
  filter(status == "ok") %>%
  filter(method != "spvc-gam") %>%
  select(method, elapsed, peak_mb, n_side, n_genes, type)

runtime_long <- runtime_df %>%
  pivot_longer(
    cols = c(elapsed, peak_mb),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    label = ifelse(type == "Varying Genes (900 Spots)", n_genes, n_side^2), 
    method = factor(method, c("celina", "stance", "cside", "ctsv", "spvc-original"))
  )  %>%
  mutate(
    method = fct_recode(method, 
        "CELINA" = "celina",
        "STANCE" = "stance",
        "CTSV" = "ctsv",
        "C-SIDE" = "cside",
        "spVC"  = "spvc-original"
    ), 
    metric = fct_recode(factor(metric), 
        "Runtime" = "elapsed", 
        "Peak Memory" = "peak_mb"
    )
  )

## define scale for elapsed 
library(scales)

cutoff <- 15  # example: linear up to 15 min, then log

hybrid_trans <- trans_new(
  "hybrid",
  transform = function(y) ifelse(y <= cutoff, 
                                 y, 
                                 cutoff + log(y - cutoff + 1)),
  inverse   = function(z) ifelse(z <= cutoff, 
                                 z, 
                                 cutoff + (exp(z - cutoff) - 1))
)


## --- test variable gene at spot 900
top = ggplot(runtime_long, aes(x = label, y = value, color = method, group = method)) +
  geom_line(size = 1.5) +
  geom_point(size = 4) +
  facet_grid(
    rows = vars(metric),      # <-- two rows: Genes / Spots
    cols = vars(type),    # <-- two columns: elapsed / peak_mb
    scales = "free" 
    # ,switch = "y"
  )+
  ggh4x::facetted_pos_scales(
    y = list(
      metric == "Runtime" ~ scale_y_continuous(
        trans = hybrid_trans,
        breaks = c(0, 5, 10, 60, 1440),
        labels = c("0 Min",  "5 Min", "10 Min", "1 Hour", "1 Day"),
        limits = c(0, 2880)
      ),
      metric == "Peak Memory" ~ scale_y_continuous(
        breaks = c(0, 5, 10, 15),
        labels = c("0 GB", "5 GB", "10 GB", "15 GB"),
        limits = c(0, 20)
      )
    )
  ) +
  scale_color_manual(
    # "tvthemes::AirNomads",
    values = pal,
    name = "Methods"
  ) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.title = element_text(size = 30, face = "bold"),
    legend.text  = element_text(size = 28),
    legend.key.size = unit(3, "lines"),
    strip.text = element_text(size = 32, face = "bold"), # for fecet heading 
    panel.background = element_rect(fill = "white", color = "black", linewidth = 0.8),
    strip.background = element_rect(fill = "grey80", color = "black"),
    axis.text.y = element_text(size = 26),
    axis.ticks.y = element_line(size = 1.2), 
    axis.text.x = element_text(size = 26), 
    axis.ticks.x = element_line(size = 1.2),
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.key = element_blank(), 
    axis.title = element_blank(), 
    plot.title = element_blank()
  )

ann_runtime <- runtime_long |>
  filter(metric == "Runtime") |>
  distinct(metric, type)      # one row per Runtime facet

ann_peak <- runtime_long |>
  filter(metric == "Peak Memory") |>
  distinct(metric, type)      # one row per Peak Memory facet


final_top <- top +
  geom_hline(
    data = ann_peak,
    aes(yintercept = 16),
    linetype = "dashed",
    color = "black",
    linewidth = 1.2,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = ann_peak,
    aes(x = 0, y = 16, label = "16 GB Limit"),
    hjust = -0.01,
    vjust = -0.5,
    size = 10,
    # fontface = "bold",
    inherit.aes = FALSE
  ) +
  geom_hline(
    data = ann_runtime,
    aes(yintercept = 120),
    linetype = "dashed",
    color = "black",
    linewidth = 1.2,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = ann_runtime,
    aes(x = 0, y = 160, label = "2 Hours"),
    hjust = -0.01,
    vjust = 0,
    size = 10,
    # fontface = "bold",
    inherit.aes = FALSE
  )

#### scalability across decomposition experiments 
decomposition_files = c(
  "./scDesign_simulators/lymph_small/ctsvg_result/decompose_remove300/runtime.csv", 
  "./scDesign_simulators/ovarian_small/ctsvg_result/decompose_v3/runtime.csv", 
  "./scDesign_simulators/breast_small/ctsvg_result/decompose_v3/runtime.csv"
)

selected_exp = c(
  "0",  "exp_0",
  # "1",  "exp_1",
  # "2", "exp_2",
  "1,2", "exp_1,2", 
  # "3",  "exp_3",
  "1,2,3","exp_1,2,3",
  # "2,4", "exp_2,4",
  # "1,2,4",  "exp_1,2,4",
  "1,2,5",    "exp_1,2,5",
  "1,2,3,4,5", "exp_1,2,3,4,5"
)

decompose_results = data.frame()

for (file in decomposition_files){
  if (nrow(decompose_results) == 0){
    decompose_results = read.csv(file)
  } else {
    decompose_results = rbind(decompose_results, read.csv(file))
  }
}

decompose_df = decompose_results %>%
  mutate(
    # dataset = the prefix before first "_"
    dataset_name = str_extract(dataset, "^[^_]+"),
    
    # extract everything after "decompose_" or "decompose_exp_"
    experiment_name = str_extract(dataset, "(?<=decompose(_exp)?_).*")
  ) %>%
  mutate(
    # clean up weird empty cases
    experiment_name = ifelse(experiment_name == "", NA, experiment_name),
    experiment_name = ifelse(experiment_name == "-", NA, experiment_name)
  ) %>%
  filter(experiment_name %in% selected_exp) %>%
  filter(status == "ok") %>%
  filter(method != "spvc-gam") %>%
  mutate(
    experiment_name = factor(experiment_name, rev(selected_exp))
  ) %>%
  mutate(
    experiment_name_display = fct_recode(factor(experiment_name), 
      "Baseline" = "0", 
      "  +(i)" = "1", 
      "  +(ii)" = "2",
      "  +(i,ii,iii)" = "1,2", 
      "  +(iv)" = "3", 
      "  +(ii,v)" = "2,4", 
      "  +(i,ii,iii,v)" = "1,2,4",
      "  +(i,ii,iii,vi)" = "1,2,5", 
      "  +(i,ii,iii,iv)" = "1,2,3", 
      "Baseline" = "exp_0", 
      "  +(i)" = "exp_1", 
      "  +(ii)" = "exp_2",
      "  +(i,ii,iii)" = "exp_1,2", 
      "  +(iv)" = "exp_3", 
      "  +(ii,v)" = "exp_2,4", 
      "  +(i,ii,iii,v)" = "exp_1,2,4",
      "  +(i,ii,iii,vi)" = "exp_1,2,5", 
      "  +(i,ii,iii,iv)" = "exp_1,2,3", 
      "  +All" = "exp_1,2,3,4,5", 
      "  +All" = "1,2,3,4,5" 
      
    )
  ) %>%
  mutate(
    test_method = fct_recode(
      factor(method),
      "CELINA" = "celina", 
      "STANCE" = "stance", 
      "CTSV" = "ctsv",
      "C-SIDE" = "cside",
      "spVC"  = "spvc-original"
  )) %>% 
  mutate(
    test_method = factor(
      test_method,
      levels = c("CELINA","STANCE","C-SIDE","CTSV","spVC")
    )
  ) %>% 
  group_by(dataset_name, experiment_name_display, test_method) %>%
  summarise(
    cpu_time = mean(elapsed * ncores, na.rm = TRUE)/60,
    .group = "drop"
  ) %>%
  mutate(
    dataset_name = factor(dataset_name, c("breast", "ovarian", "lymph"))
  )

ratio_df <- decompose_df %>%
  # extract baseline cpu_time per dataset + method
  group_by(dataset_name, test_method) %>%
  mutate(
    baseline_cpu = cpu_time[experiment_name_display == "Baseline"][1]
  ) %>%
  ungroup() %>%
  mutate(
    cpu_ratio = cpu_time / baseline_cpu
  )


decomposed = ggplot(ratio_df,
       aes(
         x = "CPU Runtime",
         y = experiment_name_display,
         fill = log1p(cpu_ratio),
         label = sprintf("%.0f", ceiling(cpu_ratio))
       )
) +
  geom_tile(color = NA, width = 1.2, height = 0.9) +
  geom_text(size = 8) +
  scale_fill_gradient2(
    low  = "#035AA6FF",
    mid = "white",
    high = "#FF9933FF",
    midpoint = log1p(2),
    limits   = range(log1p(c(0, 400)))
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
    x = NULL, y = NULL,
    fill = "CPU Time (Ratio to Baseline)"
  ) +
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
    legend.position = "bottom",
    legend.text = element_blank(),
    plot.title = element_blank(),
    legend.key.width = unit(3, "cm"), 
    legend.title = element_text(size = 26, face = "bold")
    
  )


final = (final_top + theme(plot.margin = margin(t = 20, l = 100)) )/
  (decomposed    + theme(plot.margin = margin(t = 20, l = 100)))+
  plot_layout(heights = c(2.5, 1))+
  plot_annotation(tag_levels = "a") &
  theme(
    plot.tag = element_text(size = 28, face = "bold")
  )

final

ggsave("../Paper Writing/plots/scalability/scalable.pdf", final, width = 24, height = 20, device = "pdf")
