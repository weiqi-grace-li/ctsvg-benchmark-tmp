library(Seurat)

## --- process visium data
lung_visium <- Load10X_Spatial(
  data.dir = "./real_data/lung/raw/",  # folder that contains the .h5 and spatial/
  filename = "CytAssist_FFPE_Human_Lung_Squamous_Cell_Carcinoma_filtered_feature_bc_matrix.h5"
)

# save(lung_visium, file = "./real_data/lung/processed/lung_visium.RData")

load("./real_data/lung/processed/lung_visium.RData")
spot_counts <- GetAssayData(lung_visium, assay = "Spatial", layer = "counts")
dim(spot_counts)
spot_coords <- GetTissueCoordinates(lung_visium)

# step 1. align gene names to ensembl ID 
library(AnnotationDbi)
library(org.Hs.eg.db)
library(Matrix)

sym <- rownames(spot_counts)
ens <- mapIds(
  org.Hs.eg.db,
  keys      = sym,
  keytype   = "SYMBOL",
  column    = "ENSEMBL",
  multiVals = "first"
)
length(sym)
length(unique(sym))
length(ens)
length(unique(ens))

keep <- !is.na(ens)
spot_counts <- spot_counts[keep, ]
rownames(spot_counts) <- ens[keep]

spot_counts = as(rowsum(as.matrix(spot_counts), group = ens[keep]), "dgCMatrix") # collapse symbols with same ensembl id
dim(spot_counts)


# step 2. preprocess spot count data 
spot_counts_nonzero = spot_counts>0
gene.thres = 20 
spot.thres = 20 

# 2.1 remove any genes that has less than 20 non-zero expression spots  
genes_to_save = which(rowSums(spot_counts_nonzero) >= gene.thres)

# 2.2 remove any spots that has less than 20 non-zero expression genes 
spots_to_save = which(colSums(spot_counts_nonzero) >= spot.thres)

# 2.3 get filtered spot_counts 
spot_counts = spot_counts[genes_to_save, spots_to_save]
dim(spot_counts)

# 2.4 gather x, y coordinates 
spot_coords <- GetTissueCoordinates(lung_visium)
spot_coords <- spot_coords[spots_to_save, c("x", "y"), drop = FALSE]
all(rownames(spot_coords)==colnames(spot_counts))


## --- process single cell reference 
library(zellkonverter)
library(SpatialExperiment)
library(SingleCellExperiment)
library(Triangulation)
source("./scDesign_simulators/simulation_util.R")

lung_ref <- readH5AD("./real_data/lung/raw/c1870f1f-ca36-4d96-b03b-7dc0e96d83ee.h5ad", reader = "R")

cell_counts = assay(lung_ref, "count")
# filter 1: save only genes that appear in spatial data
overlapping_genes = intersect(rownames(spot_counts), rownames(cell_counts))
cell_counts = cell_counts[overlapping_genes, ]
dim(cell_counts)

# filter 2: use threshold 
cell_counts_nonzero = cell_counts > 0 
genes_to_save = which(rowSums(cell_counts_nonzero) >= gene.thres) # doesn't filter out anything 
cells_to_save = which(colSums(cell_counts_nonzero) >= spot.thres) # doesn't filter out anything 
cell_counts = cell_counts[genes_to_save, cells_to_save]

# 3.3 filter 2 save only cell types with more than 25 cells (doesn't filter out anything)
cell_metadata = colData(lung_ref)[cells_to_save,"cell_type"]
Sys.setlocale("LC_CTYPE", "en_US.UTF-8")
levels(cell_metadata) <- gsub("/", "_", levels(cell_metadata))
cell_types_to_save = names(table(cell_metadata))[which(table(cell_metadata)>=25)] # non is less than 25
cells_to_save2 = cell_metadata %in% cell_types_to_save
cell_metadata = droplevels(cell_metadata[cells_to_save2])
cell_counts   = cell_counts[, cells_to_save2]

cell_metadata = as.data.frame(cell_metadata)
colnames(cell_metadata) = "type"
rownames(cell_metadata) = colnames(cell_counts)

# step 4. gather everything for running RCTD 
source("./utility/RCTD_util.R")

str(spot_counts)
str(spot_coords)
str(cell_counts)
str(cell_metadata)

type = cell_metadata$type 
names(type) = rownames(cell_metadata)

# weights = run_RCTD(
#   ref_count = round(cell_counts, 0),
#   ref_type = type,
#   spa_coords = as.data.frame(spot_coords),
#   spa_counts = spot_counts,
#   max_cores = 30,
#   doublet_mode = "full"
# )

# sort(colSums(weights)/sum(weights), decreasing = TRUE)
# cell_type_weights_to_save = which(colSums(weights)/sum(weights) > 0.01)
# weight_clean = weights[,cell_type_weights_to_save]
# weight_clean = weight_clean/rowSums(weight_clean)
# save(weight_clean, file = "./real_data/lung/processed/weight_clean.RData")
load("./real_data/lung/processed/weight_clean.RData")

# filter 3 save only cell types that remain in the weight matrix 
cells_to_save3 = cell_metadata$type %in% colnames(weight_clean)
cell_metadata = cell_metadata[cells_to_save3, , drop = FALSE]
cell_counts = cell_counts[, cells_to_save3]

# assemble the entire dataset 
# remaining_spots = rownames(weight_clean)
# lung = list()
# lung$dispersion = 0 
# lung$seed = 0 
# 
# lung$cell_type_proportion = colSums(weight_clean)/sum(weight_clean)
# names(lung$cell_type_proportion) = colnames(weight_clean)
# lung$cell_counts = cell_counts
# lung$cell_metadata = cell_metadata
# lung$spot_counts = spot_counts[,remaining_spots]
# lung$spot_coords = spot_coords[remaining_spots,]
# lung$spot_composition = weight_clean
# special_genes = list()
# for (ct in names(lung$cell_type_proportion)){
#   special_genes[[ct]] = list(
#     marker = c(""), ct_svg = c("")
#   )
# }
# lung$special_genes = special_genes
# lung$gene_names = rownames(spot_counts)
# save(lung, file = "./real_data/lung/processed/lung.RData")
load("./real_data/lung/processed/lung.RData")

# spot_coords = lung$spot_coords 
# save(spot_coords, file = "./real_data/lung/processed/spot_coords.RData")
# load("./real_data/lung/processed/spot_coords.RData")
# plot(spot_coords[,"x"], spot_coords[,"y"], pch = ".")
# boundary.idx <- identify(spot_coords, order = TRUE)
# lung_boundary = spot_coords[boundary.idx$ind[order(boundary.idx$order)], ]
# points(lung_boundary, type = "l", col = "blue")
# Tr.cell <- TriMesh(lung_boundary, n = 2)
# Tr.cell = remove_colinear(Tr.cell)
# save(Tr.cell, file = "./real_data/lung/processed/Tr.cell.RData")

load("./real_data/lung/processed/Tr.cell.RData")
source("./scDesign_simulators/simulation_util.R")

# reduce the size of single cell reference, cap each cell type at 10,000
# set.seed(1)
# selected_cells = c()
# for (ct in colnames(lung$spot_composition)){
#   cell_idx = which(lung$cell_metadata$type == ct)
#   if (length(cell_idx) <= 10000){
#     selected_cells = c(selected_cells, cell_idx)
#     message(paste0("selecting ", ct, ": kept all ", length(cell_idx)))
#   } else {
#     random_idx = sample(cell_idx, 10000, replace = FALSE)
#     selected_cells = c(selected_cells, random_idx)
#     message(paste0("selecting ", ct, ": cap at 10,000"))
#   }
# }
# 
# cell_metadata_reduced = cell_metadata[selected_cells, , drop = FALSE]
# cell_counts_reduced = cell_counts[, selected_cells]
# save(cell_metadata_reduced, file = "./real_data/lung/processed/cell_metadata_reduced.RData")
# save(cell_counts_reduced, file = "./real_data/lung/processed/cell_counts_reduced.RData")

load("./real_data/lung/processed/cell_metadata_reduced.RData")
load("./real_data/lung/processed/cell_counts_reduced.RData")

run_all_tests(
  sp_count = as.matrix(lung$spot_counts),
  sp_comp = as.matrix(lung$spot_composition),
  sp_coords = as.matrix(lung$spot_coords),
  
  # feed newly generated cell counts, not original cells 
  sc_count = as.matrix(cell_counts_reduced),
  sc_metadata = cell_metadata_reduced,
  spvc_tri = Tr.cell,
  sim_name = "lung",
  save_dir = "./real_data/lung/results/",
  time_save_path = "./real_data/lung/results/runtime.csv",
  sim_obj = lung,
  ncores = 40,
  ncores_cside = 40, 
  run_method = c("spvc"), 
  timeout_sec = 48*3600
)


## filtered version 

# gene level filtering: 
expressed_prop = 0.1
total_spots = ncol(lung$spot_counts)
min_spot = ceiling(expressed_prop*total_spots)

spot_counts_filter = lung$spot_counts>0
genes_to_save = which(rowSums(spot_counts_filter)>min_spot)
length(genes_to_save)

# # spot level filtering: 先不做特别激进的spot level filtering好了 (didn't filter out any)
min_count = 100
spots_to_save = which(colSums(lung$spot_counts) > min_count)
which(colSums(lung$spot_counts) < min_count)

# cell type level
cell_types_to_save = names(lung$cell_type_proportion)[which(lung$cell_type_proportion>0.03)]
cells_to_save = rownames(cell_metadata_reduced)[which(cell_metadata_reduced$type %in% cell_types_to_save)]

lung_filtered = lung

overlapping_gene = intersect(rownames(lung$spot_counts)[genes_to_save], rownames(cell_counts_reduced))
lung_filtered$cell_counts = cell_counts_reduced[overlapping_gene, cells_to_save]
lung_filtered$cell_metadata = cell_metadata_reduced[cells_to_save,,drop = FALSE]
lung_filtered$spot_counts = lung$spot_counts[genes_to_save, ]
lung_filtered$gene_names = rownames(lung_filtered$spot_counts)
lung_filtered$special_genes = lung$special_genes[cell_types_to_save]
lung_filtered$cell_type_proportion = lung$cell_type_proportion[cell_types_to_save]
lung_filtered$spot_composition = sweep(
  lung$spot_composition[, cell_types_to_save],
  1,
  rowSums(lung$spot_composition[, cell_types_to_save]),
  "/"
)

str(lung_filtered)
# save(lung_filtered, file = "./real_data/lung/processed/lung_filtered.RData") # filtered gene threshold was 0.5 
# load("./real_data/lung/processed/lung_filtered.RData")
# save(lung_filtered, file = "./real_data/lung/processed/lung_filtered_v3.RData") # filtered gene threshold was 0.1 