library(zellkonverter)
library(SpatialExperiment)
library(SingleCellExperiment)
library(Triangulation)

## -------------------------------------------##
## Step 1 & 2: load and process spatial data / cell reference  ##
## -------------------------------------------##

# --- Step 1.1: select one patient's spatial data
inter.spe <- readH5AD("./data/real_raw/kidney/visium_merge_inter_upload.h5ad", reader = "R")
inter.spe
table(inter.spe$donor_id)

inter.select = inter.spe[,inter.spe$donor_id == "PD47171"]

# --- Step 2.1: load and select the same patient's single-cell reference
# (independent of kidney_core.R, so this script is self-contained from raw data)
ref.sce <- readH5AD("./data/real_raw/kidney/RCC_upload_final_raw_counts.h5ad", reader = "R")
ref.select = ref.sce[,ref.sce$patient == "PD47171"]

spot_counts = assay(inter.select, "X")
cell_counts = assay(ref.select, "X")

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
spot_coords = reducedDim(inter.select, "spatial")
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

all(rownames(colData(ref.select))==colnames(cell_counts))

# - cell type labels from broad annotation
cell_metadata = colData(ref.select)[cells_to_save,"broad_type"]
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
weight_broad = weights[,cell_type_weights_to_save]
weight_broad = weight_broad/rowSums(weight_broad)


## -------------------------------------------##
## Step 4: Assemble spatial, reference, composition  ##
## -------------------------------------------##

# --- step 4.1: remove any filtered cell types from single cell reference
cells_to_save3 = cell_metadata$type %in% colnames(weight_broad)
cell_metadata = cell_metadata[cells_to_save3, , drop = FALSE]
cell_counts = cell_counts[, cells_to_save3]

# --- step 4.2: assemble into one object
remaining_spots = rownames(weight_broad)
kidney_inter = list()
kidney_inter$dispersion = 0
kidney_inter$seed = 0
kidney_inter$cell_type_proportion = colSums(weight_broad)/sum(weight_broad)
names(kidney_inter$cell_type_proportion) = colnames(weight_broad)
kidney_inter$cell_counts = cell_counts
kidney_inter$cell_metadata = cell_metadata
kidney_inter$spot_counts = spot_counts[,remaining_spots]
kidney_inter$spot_coords = spot_coords[remaining_spots,]
kidney_inter$spot_composition = weight_broad
special_genes = list()
for (ct in names(kidney_inter$cell_type_proportion)){
  special_genes[[ct]] = list(
    marker = c(""), ct_svg = c("")
  )
}
kidney_inter$special_genes = special_genes
kidney_inter$gene_names = rownames(spot_counts)
str(kidney_inter)

# --- final packaged output: lightly-filtered packaged dataset
filter1 = kidney_inter
save(filter1, file = "./data/real_packaged/kidney_inter/filter1.RData")

# --- Step 4.3: create tr.cell object for spvc
plot(kidney_inter$spot_coords[,"x"], kidney_inter$spot_coords[,"y"], pch = ".")
boundary.idx <- identify(kidney_inter$spot_coords, order = TRUE)
kidney_inter_bounary = kidney_inter$spot_coords[boundary.idx$ind[order(boundary.idx$order)], ]
points(kidney_inter_bounary, type = "l", col = "blue")
Tr.cell <- TriMesh(kidney_inter_bounary, n = 2)
Tr.cell = remove_colinear(Tr.cell)

# --- final packaged output: triangulation mesh for spvc
save(Tr.cell, file = "./data/real_packaged/kidney_inter/Tr.cell.RData")


## -------------------------------------------##
## Step 5: Assemble a more aggressively filtered dataset  ##
## -------------------------------------------##

# - gene level filtering:
expressed_prop = 0.1
total_spots = ncol(kidney_inter$spot_counts)
min_spot = ceiling(expressed_prop*total_spots)

spot_counts_filter = kidney_inter$spot_counts>0
genes_to_save = which(rowSums(spot_counts_filter)>min_spot)
length(genes_to_save)

# - spot-level filtering: a conservative minimum-count threshold; verified below that it removes no spots for this dataset
min_count = 100
spots_to_save = which(colSums(kidney_inter$spot_counts) > min_count)
which(colSums(kidney_inter$spot_counts) < min_count)

kidney_inter_filtered = kidney_inter

overlapping_gene = intersect(rownames(kidney_inter$spot_counts)[genes_to_save], rownames(kidney_inter$cell_counts))
kidney_inter_filtered$cell_counts = kidney_inter_filtered$cell_counts[overlapping_gene, ]
kidney_inter_filtered$spot_counts = kidney_inter_filtered$spot_counts[genes_to_save, ]
kidney_inter_filtered$gene_names = rownames(kidney_inter_filtered$spot_counts)
str(kidney_inter_filtered)

# --- final packaged output: aggressively gene/spot-filtered packaged dataset
filter2 = kidney_inter_filtered
save(filter2, file = "./data/real_packaged/kidney_inter/filter2.RData")
