### --------------------------------------------------------------------- ###
library(ggplot2)
library(arrow)
library(Seurat)
library(Matrix)
library(readr)
library(Triangulation)
source("./scDesign_simulators/simulation_util.R")
source("./scDesign_simulators/pseudo_spot_simulator.R")

xenium_path = "./xenium-publication/formal_synthetic/breast_cancer_small/"

save_base_dir = "./scDesign_simulators/breast_small/"

scdesign_raw_path = paste0(save_base_dir, "packaged/scdesign_raw.RData")
scdesign_new_path = paste0(save_base_dir, "packaged/scdesign_new_v2.RData")
# v1: with special_genes v1 
# v2: with special_genes v2 
data_new_path = paste0(save_base_dir, "packaged/data_new_v2.RData")

ncores = 48
ncores_cside = 30
seed = 42

### --------------------------------------------------------------------- ###
### STEP 4: generate scdesign3 newdata ###
### --------------------------------------------------------------------- ###
# load(paste0(save_base_dir, "packaged/special_genes_v2.RData"))
set.seed(253)
seeds = c(42, sample(1:1000, 5, replace = FALSE)) # 42, 734, 60, 67, 798, 419
seeds = seeds[1:5]


### --------------------------------------------------------------------- ###
### STEP 5: rerun all ctsvg algorithms ###
### --------------------------------------------------------------------- ###
# --- let's figure out the collection of all genes in each seed > 300 cell counts
# for breast cancer, we have MUC6, OPRPN over 300
# genes_over300 = c()
# for (seed in seeds){
#   load(paste0(save_base_dir, "packaged/data_new_v2_", seed, ".RData"))
#   temp_genes = names(which(apply(data_new$cell_counts, MARGIN = 1, FUN = max) > 300))
#   ct_svg = unique(unlist(lapply(data_new$special_genes, FUN = function(x){x$ct_svg})))
#   if (any(temp_genes %in% ct_svg)){
#     stop("Genes in ct_svgs!")
#   }
#   genes_new = setdiff(genes_over300, temp_genes)
#   print(paste0("New genes detected in dataset, ", paste0(genes_new, collapse = ", ")))
#   genes_over300 = c(genes_over300, temp_genes)
# }
# genes_over300 = unique(genes_over300)
# save(genes_over300, file = paste0(save_base_dir, "packaged/extreme_genes_over300.RData"))
load(paste0(save_base_dir, "packaged/extreme_genes_over300.RData"))

# --- run scdesign3 simulated datasets
for (seed in seeds){
  
  load(paste0(save_base_dir, "packaged/data_new_v3_", seed, ".RData"))
  
  # empty spots or genes
  
  # genes_to_save = which(rowSums(data_new$spot_counts) > 0) # genes to save 1
  # genes to save 2, rule out genes with more than 1K count per cell.
  if (any(rownames(data_new$cell_counts) != rownames(data_new$spot_counts))){
    stop("Order of genes are different between cell count matrix and spot count matrix!")
  }
  genes_to_save = names(which( (!rownames(data_new$cell_counts) %in% genes_over300) & (rowSums(data_new$spot_counts) > 0) ))
  genes_to_remove = setdiff(data_new$gene_names, genes_to_save)
  ct_svgs = unique(unlist(lapply(data_new$special_genes, FUN = function(x) {x$ct_svg})))
  if (any(genes_to_remove %in% ct_svgs)){
    stop("ct_svg gene has usually high new count but cannot be removed!")
  } else{
    message(paste0("Removing gene due to high usually high new cell count: ", paste(genes_to_remove, collapse = ", ")))
  }
  
  # spot comp
  spot_comp = data_new$spot_composition
  save_cell_types = names(which(data_new$cell_type_proportion > 0.01))
  
  if (any(rownames(data_new$cell_metadata) != colnames(data_new$cell_counts))){
    stop("Cell order doesn't match for cell metadata and cell count data")
  }
  save_cells = which((data_new$cell_metadata$type %in% save_cell_types) & colSums(data_new$cell_counts[genes_to_save,]) > 0)
  
  spot_comp = spot_comp[, save_cell_types]
  
  if (any(colnames(data_new$spot_counts) != rownames(spot_comp))){
    stop("Spot counts and spot comp spot order don't match!")
  }
  
  spot_to_save = which(colSums(data_new$spot_counts[genes_to_save,]) > 0 & rowSums(spot_comp) > 0)
  
  # renormalize
  spot_comp = t(apply(spot_comp[spot_to_save, ], MARGIN = 1, FUN = function(x) x/sum(x)))
  if (any(round(rowSums(spot_comp), 2) != 1)){
    stop("Spot compotion doesn't add up to 1 for some rows!")
  }
  
  # create sim_obj
  sim_obj = list()
  sim_obj[["gene_names"]] = data_new$gene_names[which(data_new$gene_names %in% genes_to_save)]
  sim_obj[["cell_metadata"]] = data_new$cell_metadata[save_cells,]
  sim_obj[["special_genes"]] = data_new$special_genes
  sim_obj[["seed"]] = data_new$seed
  sim_obj[["cell_type_proportion"]] = data_new$cell_type_proportion[save_cell_types]
  sim_obj[["dispersion"]] = data_new$dispersion
  sim_obj[["spot_composition"]] = spot_comp
  str(sim_obj)
  
  # gather spvc boundary
  # Tr.cell <- TriMesh(breast_boundary, n = 2)
  # Tr.cell$Tr = Tr.cell$Tr[-8, ]
  # save(Tr.cell, file = paste0(save_base_dir, "packaged/Tr_cell.RData"))
  # load(paste0(save_base_dir, "packaged/breast_boundary.RData"))
  load(paste0(save_base_dir, "packaged/Tr_cell.RData"))
  
  run_all_tests(
    sp_count = data_new$spot_counts[genes_to_save, spot_to_save],
    sp_comp = spot_comp,
    sp_coords = data_new$spot_coords[spot_to_save,],
    
    # feed newly generated cell counts, not original cells
    sc_count = data_new$cell_counts[genes_to_save,save_cells],
    sc_metadata = data_new$cell_metadata[save_cells,],
    spvc_tri = Tr.cell,
    sim_name = "breast_realistic_v3",
    # save_dir = paste0(save_base_dir, "ctsvg_result/scdesign3_seeds_v3/"),
    # time_save_path = paste0(save_base_dir, "ctsvg_result/scdesign3_seeds_v3/runtime.csv"),
    save_dir = "./Results/rotation/realistic/",
    time_save_path = "./Results/rotation/idealized/runtime.csv",
    sim_obj = sim_obj,
    ncores = ncores,
    ncores_cside = ncores_cside,
    run_method = c("spvc", "cside", "stance", "celina", "ctsv", "spvc-gam")
  )
  
}

### --------------------------------------------------------------------- ###
### STEP 6: save final results ###
### --------------------------------------------------------------------- ###
analyze_result(
  input_dir = "./Results/rotation/realistic/breast/",
  save_path ="./Results/rotation/realistic/breast/breast_rotation_summary.xlsx",
  save_sheet = "aggregate",
  threshold = c(0.05)
)

