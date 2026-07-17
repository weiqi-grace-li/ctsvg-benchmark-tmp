### --------------------------------------------------------------------- ###
library(ggplot2)
library(arrow)
library(Seurat)
library(Matrix)
library(readr)
source("./util.R")
source("./simulators/realistic/pseudo_spot_simulator.R")

preprocessed_path = "./data/realistic_preprocessed/validation/"

save_base_dir = "./data/realistic_packaged/validation/"

scdesign_raw_path = paste0(save_base_dir, "scdesign_raw.RData")
scdesign_dev_path = paste0(save_base_dir, "scdesign_dev.RData")

ncores = 256

#### ---------------------------------------- #####
#### Assemble original list  and handle data  #####
#### ---------------------------------------- #####
validation_ori = build_realistic_ori(preprocessed_path)
save(validation_ori, file = paste0(save_base_dir, "validation_ori.RData"))

##### ------------------------------------------------  #####
##### Select ctSVG genes  #####
##### ------------------------------------------------  #####

# all_genes_special marks every gene as ct_svg for every cell type, purely to stop
# generate_scdesign_newdata() from flattening any gene's mean -- validation wants full
# realistic per-gene signal preserved, not a null/alt ground-truth split.
all_genes_special = setNames(
  lapply(names(validation_ori$cell_type_proportion), function(ct) list(ct_svg = validation_ori$gene_names, marker = c(""))),
  names(validation_ori$cell_type_proportion)
)

##### -------------------------------------------------------  #####
##### PART 2: fit scDesign3 (single fit, resampled per seed)   #####
##### -------------------------------------------------------  #####
scdesign_raw = scdesign3_fit_raw(
  data_ori = validation_ori,
  save_path_raw = scdesign_raw_path,
  save_path_dev = scdesign_dev_path,
  seed = 42,
  ncores = ncores,
  cell_count_threshold = 50
)

##### ------------------------------------------------------------------  #####
##### PART 3: resample per seed and package into pseudo_spot_simulator   #####
##### ------------------------------------------------------------------  #####
seeds = c(74, 153, 321, 561, 997)
for (seed in seeds){
  generate_scdesign_newdata(
    gene_names = validation_ori$gene_names,
    special_genes = all_genes_special,
    scdesign_raw_path = scdesign_raw_path,
    save_path = paste0(save_base_dir, "scdesign_new_", seed, ".RData"),
    seed = seed,
    ncores = ncores
  )

  validation_new = generate_simulator(
    data_ori = validation_ori,
    scdesign_new_path = paste0(save_base_dir, "scdesign_new_", seed, ".RData"),
    save_path = paste0(save_base_dir, "validation_new_seed", seed, ".RData"),
    seed = seed,
    special_genes = list()
  )
}

### --------------------------------------------------------------------- ###
### Create boundary file (interactive -- cannot be run unattended) ###
### --------------------------------------------------------------------- ###
# identify() requires manual clicking -- not scriptable. Use the boundary RData
# published on Zenodo instead of hand-picking your own.
# plot(validation_ori$pseudo_meta[,1], validation_ori$pseudo_meta[,2], pch = ".")
# validation_boundary = identify(validation_ori$pseudo_meta, order = TRUE)
# validation_boundary = as.matrix(validation_ori$pseudo_meta)[validation_boundary$ind[order(validation_boundary$order)], ]
# points(validation_boundary, type = "l", col = "blue")
# save(validation_boundary, file = paste0(save_base_dir, "validation_boundary.RData"))
load(paste0(save_base_dir, "validation_boundary.RData"))
