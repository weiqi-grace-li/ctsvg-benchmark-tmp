### --------------------------------------------------------------------- ###
library(ggplot2)
library(arrow)
library(Seurat)
library(Matrix)
library(readr)
source("./util.R")
source("./simulators/realistic/pseudo_spot_simulator.R")

preprocessed_path = "./data/realistic_preprocessed/breast/"

save_base_dir = "./data/realistic_packaged/breast/"

scdesign_raw_path = paste0(save_base_dir, "scdesign_raw.RData")
scdesign_dev_path = paste0(save_base_dir, "scdesign_dev.RData")

ncores = 36
seed = 42
### --------------------------------------------------------------------- ###
### STEP 1: build breast_small_ori from preprocessed data ###
### --------------------------------------------------------------------- ###
breast_small_ori = build_realistic_ori(preprocessed_path)
save(breast_small_ori, file = paste0(save_base_dir, "breast_small_ori.RData"))

### --------------------------------------------------------------------- ###
### STEP 2: get scdesign 3 running ###
### --------------------------------------------------------------------- ###
scdesign_raw = scdesign3_fit_raw(
  data_ori = breast_small_ori,
  save_path_raw = scdesign_raw_path,
  save_path_dev = scdesign_dev_path,
  seed = seed,
  ncores = ncores,
  cell_count_threshold = 100
)

### --------------------------------------------------------------------- ###
### STEP 3: select ctsvg ###
### --------------------------------------------------------------------- ###
load(scdesign_dev_path)  # -> dev_table
special_genes = select_ctsvg_by_deviance(
  data_ori = breast_small_ori,
  dev_table = dev_table,
  ct_prop_threshold = 0.03,
  ctsvg_top_n = 50,
  save_dir = save_base_dir,
  plot_dir = paste0(save_base_dir, "ctsvg_selection/plots/")
)

### --------------------------------------------------------------------- ###
### STEP 4: generate scdesign3 newdata ###
### --------------------------------------------------------------------- ###
set.seed(253)
seeds = c(42, sample(1:1000, 5, replace = FALSE)) # 42, 734, 60, 67, 798, 419
seeds = seeds[1:5]

for (seed in seeds){
  scdesign_new = generate_scdesign_newdata(
    gene_names = breast_small_ori$gene_names,
    special_genes = special_genes,
    scdesign_raw_path = scdesign_raw_path,
    save_path = paste0(save_base_dir, "scdesign_new_v2_", seed, ".RData"),
    seed = seed,
    ncores = ncores
  )

  breast_new = generate_simulator(
    data_ori = breast_small_ori,
    scdesign_new_path = paste0(save_base_dir, "scdesign_new_v2_", seed, ".RData"),
    save_path = paste0(save_base_dir, "data_new_v3_", seed, ".RData"),
    seed,
    special_genes
  )
}

### --------------------------------------------------------------------- ###
### extreme_genes_over300.RData -- needed to run test methods downstream ###
### --------------------------------------------------------------------- ###
genes_over300 = compute_genes_over300(save_base_dir, seeds, save_dir = save_base_dir)

### --------------------------------------------------------------------- ###
### Create boundary file (interactive -- cannot be run unattended) ###
### --------------------------------------------------------------------- ###
# identify() requires manual clicking -- not scriptable. Use the boundary RData
# published on Zenodo instead of hand-picking your own.
# plot(breast_small_ori$pseudo_meta[,1], breast_small_ori$pseudo_meta[,2], pch = ".")
# breast_boundary = identify(breast_small_ori$pseudo_meta, order = TRUE)
# breast_boundary = as.matrix(breast_small_ori$pseudo_meta)[breast_boundary$ind[order(breast_boundary$order)], ]
# points(breast_boundary, type = "l", col = "blue")
# save(breast_boundary, file = paste0(save_base_dir, "breast_boundary.RData"))
load(paste0(save_base_dir, "breast_boundary.RData"))
