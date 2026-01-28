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
### STEP 1: First collect data_ori, for scdesign3 to run first ###
### --------------------------------------------------------------------- ###
# # --- 1. read in original cell data 
# 
# # cell meta info 
# cell_meta = as.data.frame(read_parquet(paste0(path, "cell_meta.parquet")))
# str(cell_meta)
# rownames(cell_meta) = cell_meta[, "__index_level_0__"]
# cell_meta = cell_meta[, c("x", "y", "type")]
# 
# # cell count info 
# cell_count = Read10X(paste0(path, "cell_counts_10x/"))
# 
# stopifnot(length(intersect(rownames(cell_meta), colnames(cell_count)))==length(rownames(cell_meta)))
# stopifnot(length(intersect(rownames(cell_meta), colnames(cell_count)))==length(colnames(cell_count)))
# 
# # align meta and count 
# cell_names = intersect(rownames(cell_meta), colnames(cell_count))
# cell_meta = cell_meta[cell_names, ]
# cell_count = cell_count[, cell_names]
# table(cell_meta$type)
# 
# # read in pseudo_count information 
# pseudo_count = as.data.frame(read_parquet(paste0(path, "pseudo_count.parquet")))
# rownames(pseudo_count) = pseudo_count[,"__index_level_0__"]
# pseudo_count = pseudo_count[,-which(colnames(pseudo_count)=="__index_level_0__")]
# pseudo_count = t(as.matrix(pseudo_count))
# which(rowSums(pseudo_count)==0)
# pseudo_count = pseudo_count[-which(rowSums(pseudo_count)==0),]
# which(colSums(pseudo_count)==0)
# 
# 
# # drop empty cells and empty genes 
# which(rowSums(cell_count) == 0)
# which(colSums(cell_count) == 0)
# 
# # assemble general form 
# 
# # --- 2. cell gene drop, cell spot prop, pseudo_meta
# # read cell by gene drop fraction 
# cell_gene_drop = readMM(paste0(path, "gene_dropout_frac.mtx.gz"))
# rownames(cell_gene_drop) = read_lines(paste0(path, "gene_dropout_frac.rows.tsv.gz"))
# colnames(cell_gene_drop) = read_lines(paste0(path, "gene_dropout_frac.cols.tsv.gz"))
# cell_gene_drop = as(cell_gene_drop, "CsparseMatrix")
# gene_cell_drop = t(cell_gene_drop)
# 
# # read cell by gene keep fraction
# cell_gene_keep = readMM(paste0(path, "gene_keep_frac.mtx.gz"))
# rownames(cell_gene_keep) = read_lines(paste0(path, "gene_keep_frac.rows.tsv.gz"))
# colnames(cell_gene_keep) = read_lines(paste0(path, "gene_keep_frac.cols.tsv.gz"))
# cell_gene_keep = as(cell_gene_keep, "CsparseMatrix")
# gene_cell_keep = t(cell_gene_keep)
# 
# # read cell by gene keep mask 
# keep_mask = readMM(paste0(path, "mask_keep.mtx.gz"))
# rownames(keep_mask) = read_lines(paste0(path, "mask_keep.rows.tsv.gz"))
# colnames(keep_mask) = read_lines(paste0(path, "mask_keep.cols.tsv.gz"))
# keep_mask = as(keep_mask, "CsparseMatrix")
# keep_mask = t(keep_mask)
# 
# # read cell by gene fallback mask 
# fallback_mask = readMM(paste0(path, "mask_fallback.mtx.gz"))
# rownames(fallback_mask) = read_lines(paste0(path, "mask_fallback.rows.tsv.gz"))
# colnames(fallback_mask) = read_lines(paste0(path, "mask_fallback.cols.tsv.gz"))
# fallback_mask = as(fallback_mask, "CsparseMatrix")
# fallback_mask = t(fallback_mask)
# 
# # read cell fallback vector 
# cell_fallback <- read.delim(paste0(path, "cell_fallback.tsv"), sep = "\t", stringsAsFactors = FALSE)
# cell_fallback <- setNames(cell_fallback$fallback, cell_fallback$cell_id)
# 
# # read cell by spot proportion 
# cell_spot_prop <- readMM(paste0(path, "cell_spot_prop.mtx.gz"))
# rownames(cell_spot_prop) = read_lines(paste0(path, "cell_spot_prop.rows.tsv.gz"))
# colnames(cell_spot_prop) = read_lines(paste0(path, "cell_spot_prop.cols.tsv.gz"))
# cell_spot_prop = as(cell_spot_prop, "CsparseMatrix")
# table(round(rowSums(cell_spot_prop), 2)) # vast majority is 1, handful cells not = 1 because part of it belongs to spots not in overlapping region   
# 
# # read in pseudo meta information 
# pseudo_meta <- as.data.frame(read_parquet(paste0(path, "pseudo_meta.parquet")))
# rownames(pseudo_meta) = pseudo_meta[,"__index_level_0__"]
# pseudo_meta = pseudo_meta[,c("x_transformed", "y_transformed")]
# colnames(pseudo_meta) = c("x", "y")
# 
# 
# # --- 3. align cell, gene, spot order 
# # align cell orders 
# stopifnot(nrow(cell_spot_prop)==length(cell_names))
# cell_spot_prop = cell_spot_prop[cell_names, ]
# 
# stopifnot(ncol(gene_cell_drop)==length(cell_names))
# gene_cell_drop = gene_cell_drop[,cell_names]
# 
# stopifnot(ncol(keep_mask)==length(cell_names))
# keep_mask = keep_mask[,cell_names]
# 
# stopifnot(ncol(fallback_mask)==length(cell_names))
# fallback_mask = fallback_mask[,cell_names]
# 
# stopifnot(length(cell_fallback)==length(cell_names))
# cell_fallback = cell_fallback[cell_names]
# 
# stopifnot(ncol(gene_cell_keep)==length(cell_names))
# gene_cell_keep = gene_cell_keep[,cell_names]
# 
# 
# rownames(cell_meta) = colnames(cell_count) = rownames(cell_spot_prop) =  paste0("Cell", cell_names)
# colnames(gene_cell_drop) = colnames(gene_cell_keep) = colnames(keep_mask) = colnames(fallback_mask) = names(cell_fallback) =  paste0("Cell", cell_names)
# 
# # align gene order 
# genes = intersect(intersect(intersect(rownames(cell_count), rownames(gene_cell_keep)), rownames(pseudo_count)), rownames(gene_cell_drop))
# cell_count = cell_count[genes,]
# gene_cell_keep = gene_cell_keep[genes, ]
# gene_cell_drop = gene_cell_drop[genes, ]
# pseudo_count = pseudo_count[genes, ]
# 
# # align spot order 
# spot = intersect(intersect(colnames(cell_spot_prop), rownames(pseudo_meta)), colnames(pseudo_count))
# cell_spot_prop = cell_spot_prop[,spot] 
# pseudo_meta = pseudo_meta[spot,]
# pseudo_count = pseudo_count[,spot]
# 
# # check pseudo_count
# which(rowSums(pseudo_count)==0)
# which(colSums(pseudo_count)==0)
# 
# 
# # --- 4. create breast_ori object 
# breast_small_ori = pseudo_data_ori$new(cell_count, cell_meta, pseudo_meta, pseudo_count, cell_spot_prop, gene_cell_keep = gene_cell_keep, 
#                                        gene_cell_drop = gene_cell_drop,
#                                        mask_keep = keep_mask, mask_fallback= fallback_mask, cell_fallback= cell_fallback)
# 
# save(breast_small_ori, file = paste0(path, "breast_small_ori.RData"))
load(paste0(save_base_dir, "packaged/breast_small_ori.RData"))

### --------------------------------------------------------------------- ###
### STEP 2: get scdesign 3 running ###
### --------------------------------------------------------------------- ###
# scdesign_raw = scdesign3_fit_raw(
#   data_ori = breast_small_ori, 
#   save_path_raw = scdesign_raw_path, 
#   save_path_marginal = scdesign_marginal_path, 
#   seed = seed, 
#   ncores = ncores, 
#   cell_count_threshold=100
# )

### --------------------------------------------------------------------- ###
### STEP 3: select ctsvg ###
### --------------------------------------------------------------------- ###

# we have run ctsvg selection results previously, we will just directly use it
# select_result = ctsvg_result(input_dir = paste0(save_base_dir, "ctsvg_selection/"), threshold = 0.01, top_n = 50)
# special_genes = construct_special_genes(
#   selected_list = select_result$ctsvg_table,
#   cell_types = names(breast_small_ori$cell_type_proportion), # special genes contain all cell types
#   save_dir = paste0(save_base_dir, "packaged/"),
#   data_ori = breast_small_ori,
#   agreement = 2,
#   cell_detect_threshold = 5,
#   fallback = 0,
#   plot_dir = paste0(save_base_dir, "ctsvg_selection/plots/")
# )


### --------------------------------------------------------------------- ###
### STEP 4: generate scdesign3 newdata ###
### --------------------------------------------------------------------- ###
# load(paste0(save_base_dir, "packaged/special_genes_v2.RData"))
set.seed(253)
seeds = c(42, sample(1:1000, 5, replace = FALSE)) # 42, 734, 60, 67, 798, 419
seeds = seeds[1:5]

# for (seed in seeds){
#   # scdesign_new = generate_scdesign_newdata(
#   #   gene_names = breast_small_ori$gene_names,
#   #   special_genes = special_genes,
#   #   scdesign_raw_path = scdesign_raw_path,
#   #   save_path = paste0(save_base_dir, "packaged/scdesign_new_v2_", seed, ".RData"),
#   #   seed = seed,
#   #   ncores = ncores
#   # )
# 
#   breast_new = generate_simulator(
#     data_ori = breast_small_ori,
#     scdesign_new_path = paste0(save_base_dir, "packaged/scdesign_new_v2_", seed, ".RData"),
#     # save_path = paste0(save_base_dir, "packaged/data_new_v2_", seed, ".RData"), ## 20251118 change pseudo spot 
#     save_path = paste0(save_base_dir, "packaged/data_new_v3_", seed, ".RData"),
#     seed,
#     special_genes,
#     strat = "A"
#   )
# }

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
### STEP 5.2: Rotation results ###
### save to Results/rotation/realistic/
### --------------------------------------------------------------------- ###
load(paste0(save_base_dir, "packaged/extreme_genes_over300.RData"))

# --- run scdesign3 simulated datasets 
# for (seed in seeds){
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
  for (rotate in c(0, 30, 60, 90)){
    message(paste0("Running replicate ", which(seeds == seed), " | rotate = ", rotate))
    load(paste0(save_base_dir, "packaged/breast_boundary.RData"))
    # rotation 
    # pos.rotated <- rotate_points(data_new$spot_coords[spot_to_save,], angle_degrees = rotate)
    # boundary.rotated = rotate_points(breast_boundary, angle_degrees = rotate)
    
    # rotation + translation 
    pos.rotated <- rotate_translate_points(data_new$spot_coords[spot_to_save,], angle_degrees = rotate, type = "C")
    boundary.rotated = rotate_translate_points(breast_boundary, angle_degrees = rotate, type = "C")
    
    Tr.cell.rotated = TriMesh(boundary.rotated, n = 2)
    # Tr.cell.rotated$Tr = Tr.cell.rotated$Tr[-10,] # 20260107 revision, automatic remove degenerate triangle 
    Tr.cell.rotated = remove_colinear(Tr.cell.rotated)
    
    save_dir = "./Results/rotation/realistic/breast/"
    sim_name = paste0("breast_rotation_translation_C_", rotate)
    
    run_all_tests(
      sp_count = data_new$spot_counts[genes_to_save, spot_to_save],
      sp_comp = spot_comp,
      sp_coords = pos.rotated,

      # feed newly generated cell counts, not original cells
      sc_count = data_new$cell_counts[genes_to_save,save_cells],
      sc_metadata = data_new$cell_metadata[save_cells,],
      spvc_tri = Tr.cell.rotated,
      sim_name = sim_name,
      save_dir = save_dir,
      time_save_path = "./Results/rotation/realistic/breast/runtime.csv",
      sim_obj = sim_obj,
      ncores = ncores,
      ncores_cside = ncores_cside,
      run_method = c("celina", "cside", "stance", "ctsv", "spvc", "spvc-gam")
    )
    
    library(STANCE)
    source("./utility/stance_util.R")
    # create stance object 
    gene_thres = 0.05; spot_thres = 0; normalized = F; correction = F
    stance_obj = creatSTANCEobject(counts = data_new$spot_counts[genes_to_save, spot_to_save], pos = pos.rotated, prop = spot_comp)
    stance_obj = data_preprocess(object = stance_obj, gene.threshold = gene_thres, spot.threshold = spot_thres, normalized = normalized)
    
    # swop to distance perserve normalization 
    stance_obj@location = dist_perserve_trans(pos.rotated)
    
    # run STANCE individual test 
    stance_obj = build_kernelMatrix(object = stance_obj)
    stance_obj = runTest1(object = stance_obj, correction = correction, pv.adjust = "BY")
    gene.list = rownames(stance_obj@gene_expression)
    utSVG.list = gene.list[stance_obj@Test_1$p_value_adj < 0.05] # use adjusted p value
    mySTANCE = runTest2(object = stance_obj, Genes_to_test = utSVG.list, Cell_types_to_test = NULL, # test for all cell types 
                        correction = correction, ncores = ncores)
    
    save_file = retrieve_file_start_with(save_dir, "stance")
    save_file = ifelse(is.na(save_file), "stance.xlsx", save_file)
    stance_details_save(
      sim_obj, mySTANCE, sim_name, "oracle", 
      file_path = paste0(save_dir, save_file)
      , test_method = "stance-mod"
    )  
  }
}



### --------------------------------------------------------------------- ###
### STEP 6: save final results ###
### --------------------------------------------------------------------- ###
# analyze_result(
#   input_dir = paste0(save_base_dir, "ctsvg_result/scdesign3_seeds/"),
#   save_path = paste0(save_base_dir, "ctsvg_result/scdesign3_seeds/breast_cancer_scdesign3_results.xlsx"),
#   save_sheet = "aggregate",
#   threshold = c(0.05)
# )
# 
analyze_result(
  input_dir = "./Results/rotation/realistic/breast/",
  save_path ="./Results/rotation/realistic/breast/breast_rotation_summary.xlsx",
  save_sheet = "aggregate",
  threshold = c(0.05)
)










# library(Seurat)
# library(SeuratObject)
#
# load(paste0(save_base_dir, "packaged/scdesign_dev.RData"))
# dev_table
# non_zero_prop = matrix(NA, ncol = ncol(dev_table), nrow = nrow(dev_table))
# rownames(non_zero_prop) = rownames(dev_table)
# colnames(non_zero_prop) = colnames(dev_table)
# non_zero_count =  non_zero_prop
# moran = non_zero_prop
# for (ct in colnames(non_zero_prop)){
#   idx = which(breast_small_ori$cell_metadata$type == ct)
#   for (gene in rownames(non_zero_prop)){
#     non_zero_prop[gene, ct] = length(which(breast_small_ori$cell_counts_ori[gene, idx] >0))/length(idx)
#     non_zero_count[gene, ct] = length(which(breast_small_ori$cell_counts_ori[gene, idx] >0))
#   }
#   features = FindSpatiallyVariableFeatures(
#     breast_small_ori$cell_counts_ori[,idx],
#     spatial.location = breast_small_ori$cell_metadata[idx,],
#     selection.method = "moransi",nfeatures = 200
#   )
#   moran[,ct] = features[rownames(moran), "p.value"]
# }
#
# hist(non_zero_prop, breaks = seq(0, 1, 0.01))
#
# cumsum(sort(breast_small_ori$cell_type_proportion, decreasing = TRUE))
# selected_cells = names(which(breast_small_ori$cell_type_proportion>0.03))
# dev_table_1 = dev_table[,selected_cells]
# dev_table_1[which((non_zero_prop[, selected_cells] < 0.01) | (non_zero_count[,selected_cells] < 100))] = 0
# moran_thres = sort(moran)[200]
# dev_table_1[which(moran[,selected_cells] > moran_thres)] = 0
# ord = order(dev_table_1, decreasing = TRUE)
# top_idx = ord[1:50]
# dev_table_1[top_idx]
# coords = arrayInd(top_idx, dim(dev_table_1))
# rownames(coords) = NULL
# res <- vector("list", length(breast_small_ori$cell_type_proportion))
# names(res) = names(breast_small_ori$cell_type_proportion)
# for (i in seq_len(50)) {
#   col <- colnames(dev_table_1)[coords[i, 2]]
#   row <- rownames(dev_table_1)[coords[i, 1]]
#   res[[col]] = c(res[[col]], row)
# }
#
# plot_dir = paste0(save_base_dir, "ctsvg_selection/plots/")
# height = 4
# width = 3.5
# dpi = 300
#
# special_genes = res
# save(special_genes, file = paste0(save_base_dir, "packaged/special_genes_v2.RData"))
#
# for (ct in names(res)){
#   # - double check whether the gene in the cell type has at least over cel_detect threshold non_zero count
#   idx_ct = which(breast_small_ori$cell_metadata$type == ct)
#   ctsvg = res[[ct]]
#
#   # - plot
#   if (!is.null(plot_dir)){
#     dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
#     col_map <- viridis_pal(option = "magma")(10)
#
#     for (gene in ctsvg){
#       base_df <- breast_small_ori$cell_metadata[idx_ct, c("x","y"), drop = FALSE]
#       expr <- as.numeric(breast_small_ori$cell_counts_ori[gene, idx_ct])
#       z    <- as.numeric(scale(expr))
#       df   <- cbind(base_df, z = z)
#
#       p <- ggplot(df, aes(x, y, colour = z)) +
#         geom_point(size = 2) +
#         stat_density_2d(aes(x = x, y = y), colour = "black", linewidth = 0.25, alpha = 0.4) +
#         scale_colour_gradientn(colors = col_map, name = "z") +
#         coord_equal() + theme_bw(base_size = 11) +
#         theme(panel.grid = element_blank(), legend.position = "right") +
#         labs(title = paste0(ct, " — ", gene), x = NULL, y = NULL)
#
#       safe_ct <- gsub("[^[:alnum:]_]+", "-", ct)
#       safe_g  <- gsub("[^[:alnum:]_]+", "-", gene)
#       fn <- file.path(plot_dir, sprintf("[%.2f_%.2f_%d]_%s_%s.png", dev_table[gene, ct], non_zero_prop[gene, ct], non_zero_count[gene, ct],safe_ct, safe_g))
#       ggsave(fn, p, width = width, height = height, dpi = dpi, bg = "white")
#     }
#   }
# }
#
#
#
#

# load(paste0(save_base_dir, "packaged/special_genes_v2.RData"))
# for (name in names(special_genes)){
#   special_genes[[name]] = list(ct_svg = special_genes[[name]], marker = NULL)
# }
