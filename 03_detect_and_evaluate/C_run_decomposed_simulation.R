library(ggplot2)
library(arrow)
library(Seurat)
library(Matrix)
library(Triangulation)
library(readr)
source("./scDesign_simulators/simulation_util.R")
source("./scDesign_simulators/pseudo_spot_simulator.R")
source("./scDesign_simulators/scDesign_decompose.R")

xenium_path = "./xenium-publication/formal_synthetic/breast_cancer_small/"

save_base_dir = "./scDesign_simulators/breast_small/"

scdesign_raw_path = paste0(save_base_dir, "packaged/scdesign_raw.RData")
scdesign_new_path = paste0(save_base_dir, "packaged/scdesign_new_v2.RData")
# v1: with special_genes v1 
# v2: with special_genes v2 
# data_new_path = paste0(save_base_dir, "packaged/data_new_v2.RData") # 20251118 
data_new_path = paste0(save_base_dir, "packaged/data_new_v3_42.RData") 
ncores = 96
ncores_cside = 10
seed = 42

# load data new 
load(data_new_path)
load(paste0(save_base_dir, "packaged/breast_small_ori.RData"))
### --------------------------------------------------------------------- ###
### Decomposition 1: Cell type proportions  ###
### --------------------------------------------------------------------- ###
set.seed(253)
seeds = c(42, sample(1:1000, 5, replace = FALSE)) # 42, 734, 60, 67, 798
seeds = seeds[1:5]

# -- load genes to be removed
load(paste0(save_base_dir, "packaged/extreme_genes_over300.RData"))

for (seed in seeds){
  if (seed %in% c(42, 734)){
    sim_types = list(
      c()
    )
  } else if (seed %in% c(60, 67, 798)){
    sim_types = list(
      c(0), c(1), c(2), c(3),
      c(1, 2),
      c(2, 4),
      c(1, 2, 4),
      c(1, 2, 5),
      c(1, 2, 3),
      c()
    )
  }
  
  for (sim in 1:length(sim_types)){
    
    print(paste0("Now running replicate: ", which(seeds == seed), " | sim_type: ", paste(sim_types[[sim]], collapse = ",")))
    
    decompose = decompose_simulator$new(data_new, data_ori = breast_small_ori, dispersion = 0.7, seed = seed, sim_type = sim_types[[sim]], n_markers = 25)
    
    # --- create Tr.Cell for spvc 
    coords = decompose$cell_metadata
    boundary = matrix(c(min(coords$x), min(coords$y), max(coords$x), min(coords$y), max(coords$x), max(coords$y), min(coords$x), max(coords$y)), ncol = 2, byrow = TRUE)
    Tr.cell = TriMesh(boundary, n = 2)
    
    if (2 %in% sim_types[[sim]]){
      load(paste0(save_base_dir, "packaged/Tr_cell.RData"))
    }
    
    # --- adjustments 
    # genes to save 2, rule out genes with more than 1K count per cell. 
    if (any(rownames(decompose$cell_counts) != rownames(decompose$spot_counts))){
      stop("Order of genes are different between cell count matrix and spot count matrix!")
    }
    genes_to_save = names(which( (!rownames(decompose$cell_counts) %in% genes_over300) & (rowSums(decompose$spot_counts) > 0) ))
    genes_to_remove = setdiff(decompose$gene_names, genes_to_save)
    ct_svgs = unique(unlist(lapply(decompose$special_genes, FUN = function(x) {x$ct_svg})))
    if (any(genes_to_remove %in% ct_svgs)){
      stop("ct_svg gene has usually high new count but cannot be removed!")
    } else{
      message(paste0("Removing gene due to high usually high new cell count: ", paste(genes_to_remove, collapse = ", ")))
    }
    
    # spot comp
    spot_comp = decompose$spot_composition
    save_cell_types = names(which(decompose$cell_type_proportion > 0.01))
    
    if (any(rownames(decompose$cell_metadata) != colnames(decompose$cell_counts))){
      stop("Cell order doesn't match for cell metadata and cell count data")
    }
    save_cells = which((decompose$cell_metadata$type %in% save_cell_types) & colSums(decompose$cell_counts[genes_to_save,]) > 0)
    
    spot_comp = spot_comp[, save_cell_types]
    
    if (any(colnames(decompose$spot_counts) != rownames(spot_comp))){
      stop("Spot counts and spot comp spot order don't match!")
    }
    
    spot_to_save = which(colSums(decompose$spot_counts[genes_to_save,]) > 0 & rowSums(spot_comp) > 0)
    
    # renormalize
    spot_comp = t(apply(spot_comp[spot_to_save, ], MARGIN = 1, FUN = function(x) x/sum(x)))
    if (any(round(rowSums(spot_comp), 2) != 1)){
      stop("Spot compotion doesn't add up to 1 for some rows!")
    }
    
    # create sim_obj
    sim_obj = list()
    sim_obj[["gene_names"]] = decompose$gene_names[which(decompose$gene_names %in% genes_to_save)]
    sim_obj[["cell_metadata"]] = decompose$cell_metadata[save_cells,]
    sim_obj[["special_genes"]] = decompose$special_genes
    sim_obj[["seed"]] = decompose$seed
    sim_obj[["cell_type_proportion"]] = decompose$cell_type_proportion[save_cell_types]
    sim_obj[["dispersion"]] = decompose$dispersion[genes_to_save, spot_to_save]
    sim_obj[["spot_composition"]] = spot_comp
    str(sim_obj)
    
    # --- run the tests 
    run_all_tests(
      sp_count = decompose$spot_counts[genes_to_save,spot_to_save],
      sp_comp = spot_comp,
      sp_coords = decompose$spot_coords[spot_to_save,],
      
      sc_count = decompose$cell_counts[genes_to_save,save_cells],
      sc_metadata = decompose$cell_metadata[save_cells,],
      
      spvc_tri = Tr.cell,
      sim_name = paste0("breast_decompose_exp_", paste(sim_types[[sim]], collapse = ",")),
      save_dir = paste0(save_base_dir, "ctsvg_result/decompose_v3/"),
      time_save_path = paste0(save_base_dir, "ctsvg_result/decompose_v3/runtime.csv"),
      sim_obj = sim_obj,
      ncores = ncores, 
      ncores_cside =  ncores_cside,
      # run_method = run_method
      run_method = c("spvc-gam")
    )   
  }
}

### --------------------------------------------------------------------- ###
### Decomposition 2: use realistic location  ###
### --------------------------------------------------------------------- ###
analyze_result(
  input_dir = paste0(save_base_dir, "ctsvg_result/decompose_v3/"),
  save_path = paste0(save_base_dir, "ctsvg_result/decompose_v3/summary_result.xlsx"),
  save_sheet = "aggregate",
  threshold = c(0.05)
)

