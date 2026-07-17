library(spatialLIBD)

## -------------------------------------------##
## Step 1 & 2: load and process spatial data / cell reference  ##
## -------------------------------------------##

# --- Step 1.1 & 2.1: fetch and load raw data (dlpfc has no local raw folder; fetched via spatialLIBD)
zip_path <- spatialLIBD::fetch_data(type = "spatialDLPFC_snRNAseq")  # downloads a .zip
unz_dir  <- unzip(zip_path, exdir = tempdir())  # ephemeral staging, not a packaged deliverable

ref <- HDF5Array::loadHDF5SummarizedExperiment(
  file.path(tempdir(), "sce_DLPFC_annotated")
)

spe <- fetch_data(type = "spe")
brain.select <- spe[,spe$sample_id == "151673"]

library(zellkonverter)
library(SpatialExperiment)
library(SingleCellExperiment)
library(Triangulation)

spot_counts = assay(brain.select, "counts")
cell_counts = assay(ref, "counts")

# --- step 1.2: Gene ID alignment (cell reference only)
library(AnnotationDbi)
library(org.Hs.eg.db)

sym <- rownames(cell_counts)
ens <- mapIds(
  org.Hs.eg.db,
  keys      = sym,
  keytype   = "SYMBOL",
  column    = "ENSEMBL",
  multiVals = "first"
)
# - drop unmapped symbols
keep <- !is.na(ens)
cell_counts <- cell_counts[keep, ]
rownames(cell_counts) <- ens[keep]

# - collapse any duplicate Ensembl IDs
cell_counts = rowsum(cell_counts, group = ens[keep])


# --- step 1.3: Basic filtering on spatial dataset
gene.thres = 20
spot.thres = 20

# - gene filter: remove any genes expressed in fewer than 20 spots
dim(spot_counts)
spot_counts_filter = spot_counts > 0
genes_to_save = which(rowSums(spot_counts_filter) >= gene.thres)

# - spot filter: remove any spots with fewer than 20 expressed genes
spots_to_save = which(colSums(spot_counts_filter) >= spot.thres)

# - get filtered spot_counts
spot_counts = spot_counts[genes_to_save, spots_to_save]
dim(spot_counts)

# --- step 1.4: Spatial locations for remaining spots
spot_coords = as.data.frame(spatialCoords(brain.select))
spot_coords = spot_coords[spots_to_save, ]
colnames(spot_coords) = c("x", "y")

# --- step 2.2: Basic filtering on single cell reference data
cell_counts_filter = cell_counts > 0

# - gene filter: remove any genes expressed in fewer than 20 cells
genes_to_save = which(rowSums(cell_counts_filter) >= gene.thres)

# - cell filter: remove any cells with fewer than 20 expressed genes
cells_to_save = which(colSums(cell_counts_filter) >= spot.thres)

# - overlap filter: save only genes that appear in spatial data
overlapping_genes = intersect(names(genes_to_save), rownames(spot_counts))

cell_counts = cell_counts[overlapping_genes, cells_to_save]

all(rownames(colData(ref))==colnames(cell_counts))

# - cell type labels from broad annotation
cell_metadata = colData(ref)[cells_to_save,"cellType_broad_k"]
Sys.setlocale("LC_CTYPE", "en_US.UTF-8")
levels(cell_metadata) <- gsub("/", "_", levels(cell_metadata))

# - cell type filter: remove any cell types with fewer than 25 cells in reference
cell_types_to_save = names(table(cell_metadata))[which(table(cell_metadata)>=25)]
cells_to_save2 = cell_metadata %in% cell_types_to_save
cell_metadata = droplevels(cell_metadata[cells_to_save2])
cell_counts   = cell_counts[, cells_to_save2]

cell_metadata = as.data.frame(cell_metadata)
colnames(cell_metadata) = "type"
rownames(cell_metadata) = colnames(cell_counts)

## -------------------------------------------##
## Step 3: Run RCTD  ##
## -------------------------------------------##

# --- Step 3.1: gather everything for running RCTD
source("./method_wrappers/RCTD_util.R")

str(spot_counts)
str(spot_coords)
str(cell_counts)
str(cell_metadata)

type = cell_metadata$type
names(type) = rownames(cell_metadata)

# --- Step 3.2: Run RCTD in full mode
weights = run_RCTD(
  ref_count = cell_counts,
  ref_type = type,
  spa_coords = as.data.frame(spot_coords),
  spa_counts = spot_counts,
  max_cores = 30,
  doublet_mode = "full"
)

# --- Step 3.3: normalize, remove cell types with less than 1% presence in dataset, renormalize
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
# - reduce the size of single cell reference, cap each cell type at 10,000
set.seed(1)
selected_cells = c()
for (ct in colnames(weight_clean)){
  cell_idx = which(cell_metadata$type == ct)
  if (length(cell_idx) <= 10000){
    selected_cells = c(selected_cells, cell_idx)
    message(paste0("selecting ", ct, ": kept all ", length(cell_idx)))
  } else {
    random_idx = sample(cell_idx, 10000, replace = FALSE)
    selected_cells = c(selected_cells, random_idx)
    message(paste0("selecting ", ct, ": cap at 10,000"))
  }
}

remaining_spots = rownames(weight_clean)
brain = list()
brain$dispersion = 0
brain$seed = 0
brain$cell_type_proportion = colSums(weight_clean)/sum(weight_clean)
names(brain$cell_type_proportion) = colnames(weight_clean)
brain$cell_counts = cell_counts[,selected_cells]
brain$cell_metadata = cell_metadata[selected_cells, , drop = FALSE]
brain$spot_counts = spot_counts[,remaining_spots]
brain$spot_coords = spot_coords[remaining_spots,]
brain$spot_composition = weight_clean
special_genes = list()
for (ct in names(brain$cell_type_proportion)){
  special_genes[[ct]] = list(
    marker = c(""), ct_svg = c("")
  )
}
brain$special_genes = special_genes
brain$gene_names = rownames(spot_counts)
str(brain)

# --- final packaged output: lightly-filtered packaged dataset
filter1 = brain
save(filter1, file = "./data/real_packaged/dlpfc/filter1.RData")

# --- Step 4.3: create tr.cell object for spvc
# get only the coordinates and the decomposition
plot(brain$spot_coords[,"x"], brain$spot_coords[,"y"], pch = ".")
boundary.idx <- identify(brain$spot_coords, order = TRUE)
brain_bounary = brain$spot_coords[boundary.idx$ind[order(boundary.idx$order)], ]
points(brain_bounary, type = "l", col = "blue")
Tr.cell <- TriMesh(brain_bounary, n = 2)
Tr.cell = remove_colinear(Tr.cell)

# --- final packaged output: triangulation mesh for spvc
save(Tr.cell, file = "./data/real_packaged/dlpfc/Tr.cell.RData")


## -------------------------------------------##
## Step 5: Assemble a more aggressively filtered dataset  ##
## -------------------------------------------##

# - gene level filtering:
expressed_prop = 0.1
total_spots = ncol(brain$spot_counts)
min_spot = ceiling(expressed_prop*total_spots)

spot_counts_filter = brain$spot_counts>0
genes_to_save = which(rowSums(spot_counts_filter)>min_spot)
length(genes_to_save)

# - spot-level filtering: a conservative minimum-count threshold; verified below that it removes no spots for this dataset
min_count = 100
spots_to_save = which(colSums(brain$spot_counts) > min_count)
which(colSums(brain$spot_counts) < min_count)

brain_filtered = brain

overlapping_gene = intersect(rownames(brain$spot_counts)[genes_to_save], rownames(brain$cell_counts))
brain_filtered$cell_counts = brain_filtered$cell_counts[overlapping_gene, ]
brain_filtered$spot_counts = brain_filtered$spot_counts[genes_to_save, ]
brain_filtered$gene_names = rownames(brain_filtered$spot_counts)
str(brain_filtered)

# --- final packaged output: aggressively gene/spot-filtered packaged dataset
filter2 = brain_filtered
save(filter2, file = "./data/real_packaged/dlpfc/filter2.RData")
