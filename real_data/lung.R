library(Seurat)

## -------------------------------------------##
## Step 1: load and process spatial data ##
## -------------------------------------------##

## --- 1.1 process visium data
lung_visium <- Load10X_Spatial(
  data.dir = "./data/real_raw/lung/",  # folder that contains the .h5 and spatial/
  filename = "CytAssist_FFPE_Human_Lung_Squamous_Cell_Carcinoma_filtered_feature_bc_matrix.h5"
)

## - spot counts
spot_counts <- GetAssayData(lung_visium, assay = "Spatial", layer = "counts")
dim(spot_counts)

## - coordinates
spot_coords <- GetTissueCoordinates(lung_visium)

## --- 1.2 GENE ID Alignment (spot count)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(Matrix)

# - convert gene symbols to Ensembl IDs
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

# - drop unmapped symbols
keep <- !is.na(ens)
spot_counts <- spot_counts[keep, ]
rownames(spot_counts) <- ens[keep]

# - collapse any duplicate Ensembl IDs
spot_counts = as(rowsum(as.matrix(spot_counts), group = ens[keep]), "dgCMatrix") # collapse symbols with the same Ensembl ID
dim(spot_counts)


# --- step 1.3 basic filtering on spatial dataset
spot_counts_nonzero = spot_counts>0
gene.thres = 20
spot.thres = 20

# - gene filter: remove any genes that have fewer than 20 non-zero expression spots
genes_to_save = which(rowSums(spot_counts_nonzero) >= gene.thres)

# - spot filter: remove any spots that have fewer than 20 non-zero expression genes
spots_to_save = which(colSums(spot_counts_nonzero) >= spot.thres)

spot_counts = spot_counts[genes_to_save, spots_to_save]
dim(spot_counts)

# --- step 1.4 spatial locations for remaining spots
spot_coords <- GetTissueCoordinates(lung_visium)
spot_coords <- spot_coords[spots_to_save, c("x", "y"), drop = FALSE]
all(rownames(spot_coords)==colnames(spot_counts))


## -------------------------------------------##
## Step 2: load and process reference data ##
## -------------------------------------------##

## --- step 2.1 load single cell reference
library(zellkonverter)
library(SpatialExperiment)
library(SingleCellExperiment)
library(Triangulation)

lung_ref <- readH5AD("./data/real_raw/lung/c1870f1f-ca36-4d96-b03b-7dc0e96d83ee.h5ad", reader = "R")

# - cell counts
cell_counts = assay(lung_ref, "count")

## --- step 2.2 basic filters

# - overlap filter: save only genes that appear in spatial data
overlapping_genes = intersect(rownames(spot_counts), rownames(cell_counts))
cell_counts = cell_counts[overlapping_genes, ]
dim(cell_counts)

# - gene filter: remove any genes that have fewer than 20 cells (none fall below this threshold for this dataset)
cell_counts_nonzero = cell_counts > 0
genes_to_save = which(rowSums(cell_counts_nonzero) >= gene.thres)

# - cell filter: remove any cells that have fewer than 20 non-zero expression genes (none fall below this threshold for this dataset)
cells_to_save = which(colSums(cell_counts_nonzero) >= spot.thres)

cell_counts = cell_counts[genes_to_save, cells_to_save]

# - cell type filter: remove any cell types with fewer than 25 cells in reference (none fall below this threshold for this dataset)
cell_metadata = colData(lung_ref)[cells_to_save,"cell_type"]
Sys.setlocale("LC_CTYPE", "en_US.UTF-8")
levels(cell_metadata) <- gsub("/", "_", levels(cell_metadata))
cell_types_to_save = names(table(cell_metadata))[which(table(cell_metadata)>=25)]
cells_to_save2 = cell_metadata %in% cell_types_to_save
cell_metadata = droplevels(cell_metadata[cells_to_save2])
cell_counts   = cell_counts[, cells_to_save2]

cell_metadata = as.data.frame(cell_metadata)
colnames(cell_metadata) = "type"
rownames(cell_metadata) = colnames(cell_counts)


## -------------------------------------------##
## Step 3: running RCTD  ##
## -------------------------------------------##
# --- step 3.1 gather everything for running RCTD
source("./method_wrappers/RCTD_util.R")

str(spot_counts)
str(spot_coords)
str(cell_counts)
str(cell_metadata)

type = cell_metadata$type
names(type) = rownames(cell_metadata)

# --- step 3.2 run RCTD
weights = run_RCTD(
  ref_count = round(cell_counts, 0),
  ref_type = type,
  spa_coords = as.data.frame(spot_coords),
  spa_counts = spot_counts,
  max_cores = 30,
  doublet_mode = "full"
)

# --- step 3.3 normalize, remove cell types with less than 1% presence, renormalize
sort(colSums(weights)/sum(weights), decreasing = TRUE)
cell_type_weights_to_save = which(colSums(weights)/sum(weights) > 0.01)
weight_clean = weights[,cell_type_weights_to_save]
weight_clean = weight_clean/rowSums(weight_clean)

## -------------------------------------------##
## Step 4: Assemble spatial, reference, composition  ##
## -------------------------------------------##

# --- step 4.1: remove any filtered cell types from single cell reference
cells_to_save3 = cell_metadata$type %in% colnames(weight_clean)
cell_metadata = cell_metadata[cells_to_save3, , drop = FALSE]
cell_counts = cell_counts[, cells_to_save3]

# --- step 4.2: assemble into one object
remaining_spots = rownames(weight_clean)
lung = list()
lung$dispersion = 0
lung$seed = 0

lung$cell_type_proportion = colSums(weight_clean)/sum(weight_clean)
names(lung$cell_type_proportion) = colnames(weight_clean)
lung$cell_counts = cell_counts
lung$cell_metadata = cell_metadata
lung$spot_counts = spot_counts[,remaining_spots]
lung$spot_coords = spot_coords[remaining_spots,]
lung$spot_composition = weight_clean
special_genes = list()
for (ct in names(lung$cell_type_proportion)){
  special_genes[[ct]] = list(
    marker = c(""), ct_svg = c("")
  )
}
lung$special_genes = special_genes
lung$gene_names = rownames(spot_counts)

# - reduce cell reference numbers
# reduce the size of single cell reference, cap each cell type at 10,000
set.seed(1)
selected_cells = c()
for (ct in colnames(lung$spot_composition)){
  cell_idx = which(lung$cell_metadata$type == ct)
  if (length(cell_idx) <= 10000){
    selected_cells = c(selected_cells, cell_idx)
    message(paste0("selecting ", ct, ": kept all ", length(cell_idx)))
  } else {
    random_idx = sample(cell_idx, 10000, replace = FALSE)
    selected_cells = c(selected_cells, random_idx)
    message(paste0("selecting ", ct, ": cap at 10,000"))
  }
}

cell_metadata_reduced = cell_metadata[selected_cells, , drop = FALSE]
cell_counts_reduced = cell_counts[, selected_cells]

# - reassemble into filter1 with reduced cell reference
filter1 = lung
filter1$cell_counts = cell_counts_reduced
filter1$cell_metadata = cell_metadata_reduced

# --- final packaged output: lightly-filtered packaged dataset (reduced cell reference)
save(filter1, file = "./data/real_packaged/lung/filter1.RData")

# --- Step 4.3: create tr.cell object for spvc
spot_coords = lung$spot_coords
plot(spot_coords[,"x"], spot_coords[,"y"], pch = ".")
boundary.idx <- identify(spot_coords, order = TRUE)
lung_boundary = spot_coords[boundary.idx$ind[order(boundary.idx$order)], ]
points(lung_boundary, type = "l", col = "blue")
Tr.cell <- TriMesh(lung_boundary, n = 2)
Tr.cell = remove_colinear(Tr.cell)

# --- final packaged output: triangulation mesh for spvc
save(Tr.cell, file = "./data/real_packaged/lung/Tr.cell.RData")


## -------------------------------------------##
## Step 5: Assemble a more aggressively filtered dataset   ##
## -------------------------------------------##

# - gene level filtering:
expressed_prop = 0.1
total_spots = ncol(lung$spot_counts)
min_spot = ceiling(expressed_prop*total_spots)

spot_counts_filter = lung$spot_counts>0
genes_to_save = which(rowSums(spot_counts_filter)>min_spot)
length(genes_to_save)

# - spot-level filtering: a conservative minimum-count threshold; verified below that it removes no spots for this dataset
min_count = 100
spots_to_save = which(colSums(lung$spot_counts) > min_count)
which(colSums(lung$spot_counts) < min_count)

# - cell type level filtering:
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

# --- final packaged output: aggressively gene/spot/cell-type-filtered packaged dataset
filter2 = lung_filtered
save(filter2, file = "./data/real_packaged/lung/filter2.RData")
