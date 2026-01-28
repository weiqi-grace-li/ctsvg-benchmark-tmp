## analyze celina result
library(dplyr)
library(tidyr)
library(forcats)
path = "./xenium-publication/formal_synthetic/breast_cancer_demon/experiment_result_v2/"

files = c(
  "experiment_celina.xlsx",
  "experiment_stance.xlsx",
  "experiment_cside.csv",
  "experiment_spvc.csv",
  "experiment_ctsv.csv"
)

top_number = seq(10, 1000, 10)
seeds = c(74, 153, 321, 561, 997)

load("./xenium-publication/formal_synthetic/breast_cancer_demon/packaged/breast_ori.RData")
cell_thres = 0.05

### ---------------- ###
### functions ###
### ---------------- ###
top_n_pairs = function(result_table, n){
  if (nrow(result_table) > n){
    return(rownames(result_table)[order(result_table$p_adj)[1:n]]) 
  } else{
    return(rownames(result_table))
  }
}

overlap = function(a, b, N){
  if (length(a) == N & length(b) == N){
    overlap = length(intersect(a, b))/N
  } else{
    overlap = length(intersect(a, b))/max(length(a), length(b))
  }
}

overlap_table = function(sig_genes, N, seeds){
  temp_overlap = rep(NA, 6)
  names(temp_overlap) = c("visium_pseudob", "visium_pseudob2", "pseudob_pseudob2", "syn_visium", "syn_pseudob", "syn")
  temp_overlap["visium_pseudob"] = overlap(sig_genes$visium, sig_genes$pseudob, N)
  temp_overlap["visium_pseudob2"] = overlap(sig_genes$visium, sig_genes$pseudob2, N)
  temp_overlap["pseudob_pseudob2"] = overlap(sig_genes$pseudob, sig_genes$pseudob2, N)
  
  overlap_syn_visium = c()
  overlap_syn_pseudob = c()
  overlap_syn = c()
  
  for (seed_select in seeds){
    sig = sig_genes[[paste0("A_seed", seed_select)]]
    overlap_syn_visium = c(overlap_syn_visium, overlap(sig_genes$visium, sig, N))  
    overlap_syn_pseudob = c(overlap_syn_pseudob, overlap(sig_genes$pseudob, sig, N))
  }
  temp_overlap["syn_visium"] = mean(overlap_syn_visium)
  temp_overlap["syn_pseudob"] = mean(overlap_syn_pseudob)
  
  compare_set = combn(seeds, 2)
  for (j in 1:ncol(compare_set)){
    sig_a = sig_genes[[paste0("A_seed", compare_set[1,j])]]
    sig_b = sig_genes[[paste0("A_seed", compare_set[2,j])]]
    overlap_syn = c(overlap_syn, overlap(sig_a, sig_b, N))
  }
  temp_overlap["syn"] = mean(overlap_syn)
  return(temp_overlap)
}

### -------------------- ###
### create final results ###
### -------------------- ###
overlap_results = data.frame(
  test_method = character(), 
  top_N = integer(), 
  cell_prop = character(), # either full, or above 0.05
  visium_pseudob = numeric(),
  visium_pseudob2 = numeric(),
  pseudob_pseudob2 = numeric(),
  syn_visium = numeric(),
  syn_pseudob = numeric(),
  syn = numeric()
)


### ---------------- ###
### main result loop ###
### ---------------- ###

for (file in files){
  
  ## --- read experiment results 
  full_path <- paste0(path, file)
  
  if (grepl("\\.csv$", file, ignore.case = TRUE)) {
    temp_result = read.csv(full_path)
  } else if (grepl("\\.xlsx$", file, ignore.case = TRUE)) {
    temp_result = readxl::read_excel(full_path)
  } else {
    stop("Unsupported file type: ", file)
  }  
  
  temp_result = subset(temp_result, p_value >= 0)
  
  file_methods = unique(temp_result$test_method)
  print(paste0("File ", file, " has method, ", paste0(file_methods, collapse = ", ")))
  
  
  ## --- separate experiment results by method (only necessary for spVC)
  for (file_method in file_methods){
    
    ## --- put results in list 
    # organize results into list with visium, xenium and synthetic 
    temp_result_method = subset(temp_result, test_method == file_method)
    temp_result_method = as.data.frame(temp_result_method)
    temp_result_list = list(
      visium = subset(temp_result_method, simulation_name == "visium"),
      pseudoa = subset(temp_result_method, simulation_name == "pseudoa"),
      pseudob = subset(temp_result_method, simulation_name == "pseudob"),
      pseudob2 = subset(temp_result_method, simulation_name == "pseudob2")
    )
    temp_result_list_0.05 = list(
      visium = subset(temp_result_method, ((simulation_name == "visium") & (cell_proportion > cell_thres))),
      pseudoa = subset(temp_result_method, ((simulation_name == "pseudoa") & (cell_proportion > cell_thres))),
      pseudob = subset(temp_result_method, ((simulation_name == "pseudob") & (cell_proportion > cell_thres))),
      pseudob2 = subset(temp_result_method, ((simulation_name == "pseudob2") & (cell_proportion > cell_thres)))
    )
    seeds = c(74, 153, 321, 561)
    for (seed_select in seeds){
      temp_result_list[[paste0("A_seed", seed_select)]] = subset(temp_result_method, simulation_name == paste0("A_seed", seed_select))
      temp_result_list[[paste0("C_seed", seed_select)]] = subset(temp_result_method, simulation_name == paste0("C_seed", seed_select))
      
      temp_result_list_0.05[[paste0("A_seed", seed_select)]] = subset(temp_result_method, simulation_name == paste0("A_seed", seed_select) & cell_proportion > cell_thres)
      temp_result_list_0.05[[paste0("C_seed", seed_select)]] = subset(temp_result_method, simulation_name == paste0("C_seed", seed_select) & cell_proportion > cell_thres)
    }
    
    # set all rownames as gene-cell-type pair 
    for (table in names(temp_result_list)){
      rownames(temp_result_list[[table]]) = paste0(temp_result_list[[table]]$cell_type, "-", temp_result_list[[table]]$gene_name)
      rownames(temp_result_list_0.05[[table]]) = paste0(temp_result_list_0.05[[table]]$cell_type, "-", temp_result_list_0.05[[table]]$gene_name)
      # if celina or stance, add p adj 
      if (file_method %in% c("celina", "stance")){
        temp_result_list[[table]]$p_adj = p.adjust(temp_result_list[[table]]$p_value, "BH")
        temp_result_list_0.05[[table]]$p_adj = p.adjust(temp_result_list_0.05[[table]]$p_value, "BH")
      }
      
    }
    
    ## ---- Fill the overlap table
    for (N in top_number){
      sig_genes = list(
        visium = top_n_pairs(temp_result_list$visium, N),
        pseudob = top_n_pairs(temp_result_list$pseudob, N), 
        pseudob2 = top_n_pairs(temp_result_list$pseudob2, N)
      )
      
      sig_genes_0.05 = list(
        visium = top_n_pairs(temp_result_list_0.05$visium, N),
        pseudob = top_n_pairs(temp_result_list_0.05$pseudob, N), 
        pseudob2 = top_n_pairs(temp_result_list_0.05$pseudob2, N)
      )
      for (seed_select in seeds){
        sig_genes[[paste0("A_seed", seed_select)]] = top_n_pairs(temp_result_list[[paste0("A_seed", seed_select)]], N)
        sig_genes_0.05[[paste0("A_seed", seed_select)]] = top_n_pairs(temp_result_list_0.05[[paste0("A_seed", seed_select)]], N)
      }
      
      overlap_results = rbind(
        overlap_results, 
        data.frame(
          t(c(test_method = file_method, top_N = N, cell_prop = "full", overlap_table(sig_genes, N, seeds)))
        )
      )

      overlap_results = rbind(
        overlap_results, 
        data.frame(
          t(c(test_method = file_method, top_N = N, cell_prop = "0.05", overlap_table(sig_genes_0.05, N, seeds)))
        )
      )
    }
  }
}

overlap_summary = overlap_results %>%
  filter(cell_prop == "full") %>%
  filter(test_method != "spvc-gam") %>%
  pivot_longer(cols = c(visium_pseudob, visium_pseudob2, pseudob_pseudob2, syn_visium, syn_pseudob, syn),
               names_to = "metric", values_to = "value") %>%
  mutate(
    value = as.numeric(value), 
    top_N = as.integer(top_N),
    test_method = factor(test_method, c("celina", "stance", "cside", "ctsv", "spvc-original"))
  )%>%
  filter(top_N <=200) %>%
  mutate(
    test_method = fct_recode(
      test_method, 
      "CELINA" = "celina", 
      "STANCE" = "stance", 
      "C-SIDE" = "cside", 
      "CTSV" = "ctsv", 
      "spVC" = "spvc-original"
    )
  )
  

overlap_set1 = overlap_summary %>%
  filter(metric %in% c("visium_pseudob", "visium_pseudob2", "pseudob_pseudob2", "syn_visium")) %>%
  mutate(display_metric = fct_relevel(metric, rev(c("pseudob_pseudob2", "visium_pseudob", "visium_pseudob2", "syn_visium")))) %>%
  mutate(
    display_metric = fct_recode(
      metric, 
      "Xenium Rep 1 vs. Rep 2: Replicate"= "pseudob_pseudob2",
      "Xenium Rep 1 vs. Visium Rep: Replicate + Platform + Pseudo-spot"= "visium_pseudob",
      "Xenium Rep 2 vs. Visium Rep: Replicate + Platform + Pseudo-spot"= "visium_pseudob2",
      "Simulated vs. Visium Rep: Replicate + Platform + Pseudo-spot + scDesign3"= "syn_visium" 
    )
  ) 
  

validdation = ggplot(overlap_set1, aes(
  x = top_N, 
  y = value,
  color = display_metric,
  linetype = display_metric,
  group = display_metric
)) +
  geom_line(size = 1) +
  facet_wrap(~test_method, nrow = 1, scales = "fixed") +
  scale_color_manual(
   name = "Comparisons",
   values = c(
      "Xenium Rep 1 vs. Rep 2: Replicate"= "black",
      "Xenium Rep 1 vs. Visium Rep: Replicate + Platform + Pseudo-spot"= "black",
      "Xenium Rep 2 vs. Visium Rep: Replicate + Platform + Pseudo-spot"   = "black", 
      "Simulated vs. Visium Rep: Replicate + Platform + Pseudo-spot + scDesign3" = "coral"
    )
  ) +
  scale_linetype_manual(
  name = "Comparisons",
  values = c(
    "Xenium Rep 1 vs. Rep 2: Replicate" = "solid",
    "Xenium Rep 1 vs. Visium Rep: Replicate + Platform + Pseudo-spot" = "dashed",
    "Xenium Rep 2 vs. Visium Rep: Replicate + Platform + Pseudo-spot"   = "dashed", 
    "Simulated vs. Visium Rep: Replicate + Platform + Pseudo-spot + scDesign3" = "solid"
  )) +
  guides(
  color = guide_legend(nrow = 3, byrow = TRUE), 
  linetype = guide_legend(nrow = 3, byrow = TRUE)
  )+
  labs(
    x = "Top N significant gene-cell-type pairs", 
    y = "Overlap (%)"
  )+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.box = "vertical",          # stack legends vertically
    legend.title.position = "top",    # title above each legend key row
    legend.direction = "horizontal",
    legend.title = element_text(size = 20, face = "bold"),
    legend.text  = element_text(size = 20),
    legend.key.size = unit(1.5, "lines"),
    strip.text = element_text(size = 24, face = "bold"), # for fecet heading 
    panel.background = element_rect(fill = "white", color = "black", linewidth = 0.8),
    strip.background = element_rect(fill = "grey80", color = "black"),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x  = element_text(size = 20),
    axis.text.y  = element_text(size = 20)  
  )

overlap_set2 = overlap_summary %>%
  filter(metric %in% c("syn_pseudob", "syn")) 
  
ggplot(overlap_set2, aes(
    x = top_N, 
    y = value,
    color = metric,
    linetype = metric,
    group = metric
  )) +
  geom_line(size = 1) +
  facet_wrap(~test_method, nrow = 1, scales = "fixed") +
  scale_color_manual(
    name = "Comparisons",
    values = c(
      "syn"= "black",
      "syn_pseudob" = "coral"
    )
  ) +
  guides(
    color = guide_legend(nrow = 3, byrow = TRUE), 
  )+
  labs(
    x = "Top N significant gene-cell-type pairs", 
    y = "Overlap (%)"
  )+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.box = "vertical",          # stack legends vertically
    legend.title.position = "top",    # title above each legend key row
    legend.direction = "horizontal",
    legend.title = element_text(size = 20, face = "bold"),
    legend.text  = element_text(size = 20),
    legend.key.size = unit(1.5, "lines"),
    strip.text = element_text(size = 24, face = "bold"), # for fecet heading 
    panel.background = element_rect(fill = "white", color = "black", linewidth = 0.8),
    strip.background = element_rect(fill = "grey80", color = "black"),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x  = element_text(size = 20),
    axis.text.y  = element_text(size = 20)  
  )


validdation

ggsave("../Paper Writing/plots/figure3/pseudo_spot_validation.pdf", width = 21.16, height = 8, device = "pdf")
