### --------------------------------------------------------------------- ###
library(ggplot2)
library(arrow)
library(Seurat)
library(Matrix)
library(readr)
source("./util.R")
source("./simulators/realistic/pseudo_spot_simulator.R")

preprocessed_path = "./data/realistic_preprocessed/ovarian/"

save_base_dir = "./data/realistic_packaged/ovarian/"

scdesign_raw_path = paste0(save_base_dir, "scdesign_raw.RData")
scdesign_dev_path = paste0(save_base_dir, "scdesign_dev.RData")

ncores = 36
ncores_small = 30

seed = 42

### --------------------------------------------------------------------- ###
### STEP 1: build ovarian_small_ori from preprocessed data ###
### --------------------------------------------------------------------- ###
ovarian_small_ori = build_realistic_ori(preprocessed_path)
save(ovarian_small_ori, file = paste0(save_base_dir, "ovarian_small_ori.RData"))


### --------------------------------------------------------------------- ###
### STEP 2: get scdesign 3 running ###
### --------------------------------------------------------------------- ###
scdesign_raw = scdesign3_fit_raw(
  data_ori = ovarian_small_ori,
  save_path_raw = scdesign_raw_path,
  save_path_dev = scdesign_dev_path,
  seed = seed,
  ncores = ncores,
  ncores_small = ncores_small,
  cell_count_threshold=100
)


### --------------------------------------------------------------------- ###
### STEP 3: select ctsvg ###
### --------------------------------------------------------------------- ###
load(scdesign_dev_path)  # -> dev_table
special_genes = select_ctsvg_by_deviance(
  data_ori = ovarian_small_ori,
  dev_table = dev_table,
  ct_prop_threshold = 0.04,
  ctsvg_top_n = 100,
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
    gene_names = ovarian_small_ori$gene_names,
    special_genes = special_genes,
    scdesign_raw_path = scdesign_raw_path,
    save_path = paste0(save_base_dir, "scdesign_new_v2_", seed, ".RData"),
    seed = seed,
    ncores = ncores
  )

  ovarian_new = generate_simulator(
    data_ori = ovarian_small_ori,
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
# plot(ovarian_small_ori$pseudo_meta[,1], ovarian_small_ori$pseudo_meta[,2], pch = ".")
# ovarian_boundary = identify(ovarian_small_ori$pseudo_meta, order = TRUE)
# ovarian_boundary = as.matrix(ovarian_small_ori$pseudo_meta)[ovarian_boundary$ind[order(ovarian_boundary$order)], ]
# points(ovarian_boundary, type = "l", col = "blue")
# save(ovarian_boundary, file = paste0(save_base_dir, "ovarian_boundary.RData"))
load(paste0(save_base_dir, "ovarian_boundary.RData"))
