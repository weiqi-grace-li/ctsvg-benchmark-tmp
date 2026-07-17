library(dplyr)

## -------------------------------------------##
## Step 1: load and process spatial data  ##
## -------------------------------------------##

# --- step 1.1: load raw spatial data files
coords <- read.delim( # spot coordinates, in pixels
  "./data/real_raw/breast/A1_labeled_coordinates.tsv",
  header = TRUE,
  sep = "\t"
)
str(coords)

library(data.table)
counts <- fread( # spot counts, indexed by grid coordinates
  "./data/real_raw/breast/count-matrices/A1.tsv.gz"
  )
str(counts)

selection <- fread( # mapping between grid and pixel coordinates
  "./data/real_raw/breast/A1_selection.tsv.gz"
)
str(selection)

# - spot composition is provided by the paper (Supplementary Table 8);
#   no RCTD deconvolution needed for this dataset
comp = readxl::read_excel("./data/real_raw/breast/supp_8.xlsx", sheet = "major_A1") |>
  as.data.frame()
rownames(comp) = comp[,1]
comp = comp[,-1]

# - align pixel coordinates to grid coordinates
coords_align = coords %>%
  left_join(selection, by = c("x" = "new_x", "y" = "new_y"))

str(coords_align)
str(coords)
rownames(coords_align) = paste0(coords_align$x.y, "x", coords_align$y.y)

spots_to_save = intersect(rownames(coords_align), unique(counts$V1))

# - spot_coords: pixel locations of retained spots
spot_coords = coords_align[spots_to_save, c("pixel_x.x", "pixel_y.x")]
spot_coords = as.matrix(spot_coords)
colnames(spot_coords) = c("x", "y")
str(spot_coords)

# - spot_composition: pre-provided cell type proportions
spot_composition = as.matrix(comp[spots_to_save,])
str(spot_composition)

# - spot_counts: gene x spot count matrix (transposed from fread)
spot_counts = as.data.frame(counts)
rownames(spot_counts) = spot_counts$V1
spot_counts = spot_counts[spots_to_save,-1]
spot_counts = t(as.matrix(spot_counts))
str(spot_counts)

# --- step 1.2: Gene ID alignment
# Note: unlike lung/kidney/dlpfc, the breast spatial and SC data share the same
# gene naming convention; no symbol-to-Ensembl conversion is needed.

# --- step 1.3: Basic filtering on spatial dataset
spot_counts_nonzero = spot_counts>0
gene.thres = 20
spot.thres = 20

# - gene filter: remove any genes that have fewer than 20 non-zero expression spots
genes_to_save = which(rowSums(spot_counts_nonzero) >= gene.thres)

# - spot filter: remove any spots that have fewer than 20 non-zero expression genes
spots_to_save = which(colSums(spot_counts_nonzero) >= spot.thres)

# - get filtered spot_counts
spot_counts = spot_counts[genes_to_save, spots_to_save]
dim(spot_counts)

# --- step 1.4: Spatial locations for remaining spots
spot_coords <- spot_coords[spots_to_save, c("x", "y"), drop = FALSE]
all(rownames(spot_coords)==colnames(spot_counts))
spot_composition = spot_composition[spots_to_save,]


## -------------------------------------------##
## Step 2: load and process reference data  ##
## -------------------------------------------##

# --- step 2.1: load single cell reference (HER2+ patients only, following paper)
ref.meta <- read.csv(
  "./data/real_raw/breast/sc/Whole_miniatlas_meta.csv",
  header = TRUE
)
str(ref.meta)
her2_idx = which(ref.meta$subtype=="HER2+")
table(ref.meta[her2_idx, "celltype_major"])
her2_barcode = ref.meta$NAME[her2_idx]

library(Matrix)
counts <- readMM("./data/real_raw/breast/sc/matrix.mtx.gz")
barcodes <- readLines("./data/real_raw/breast/sc/barcodes.tsv.gz")
genes <- read.delim("./data/real_raw/breast/sc/features.tsv.gz", header = FALSE)

gene_names <- genes$V1   # features.tsv.gz's first two columns both store the gene symbol
rownames(counts) <- gene_names
colnames(counts) <- barcodes

cell_counts <- counts[, her2_barcode, drop = FALSE]

# --- step 2.2: Basic filtering on single cell reference data
cell_counts_nonzero = cell_counts > 0

# - gene filter: keep only genes expressed in >=20 cells and overlapping with spatial data
genes_to_save = rownames(cell_counts)[which(rowSums(cell_counts_nonzero) >= gene.thres)]
genes_to_save = intersect(genes_to_save, rownames(spot_counts))

# - cell filter: remove any cells with fewer than 20 expressed genes
cells_to_save = which(colSums(cell_counts_nonzero) >= spot.thres)
cell_counts = cell_counts[genes_to_save, cells_to_save]

# - cell type labels and recoding to match composition categories
cell_metadata = ref.meta$celltype_major[her2_idx]
cell_metadata <- as.character(cell_metadata)
cell_metadata[cell_metadata %in% c("Normal Epithelial", "Cancer Epithelial")] <- "Epithelial" # recode to match decomposition results
cell_metadata[cell_metadata %in% c("Plasmablasts")] <- "Plasma Cells" # recode to match decomposition results
cell_metadata <- factor(cell_metadata)
Sys.setlocale("LC_CTYPE", "en_US.UTF-8")
levels(cell_metadata) <- gsub("/", "_", levels(cell_metadata))

# - cell type filter: remove any cell types with fewer than 25 cells (none fall below this threshold for this dataset)
cell_types_to_save = names(table(cell_metadata))[which(table(cell_metadata)>=25)]
cells_to_save2 = cell_metadata %in% cell_types_to_save
cell_metadata = droplevels(cell_metadata[cells_to_save2])
cell_counts   = cell_counts[, cells_to_save2]

cell_metadata = as.data.frame(cell_metadata)
colnames(cell_metadata) = "type"
rownames(cell_metadata) = colnames(cell_counts)


## -------------------------------------------##
## Step 3: Process spot composition  ##
## -------------------------------------------##

# Note: for breast, cell type proportions are provided by the paper (supp_8.xlsx)
# rather than estimated via RCTD.  The step below mirrors the RCTD post-processing
# in other datasets: filter to cell types with >1% overall presence, then renormalize.

# --- step 3.3: filter cell types with less than 1% presence, renormalize
sort(colSums(spot_composition)/sum(spot_composition), decreasing = TRUE)
cell_type_weights_to_save = which(colSums(spot_composition)/sum(spot_composition) > 0.01)
spot_composition = spot_composition[,cell_type_weights_to_save]
spot_composition = spot_composition/rowSums(spot_composition)


## -------------------------------------------##
## Step 4: Assemble spatial, reference, composition  ##
## -------------------------------------------##

# --- step 4.1: remove any filtered cell types from single cell reference
cells_to_save3 = cell_metadata$type %in% colnames(spot_composition)
cell_metadata = cell_metadata[cells_to_save3, , drop = FALSE]
cell_counts = cell_counts[, cells_to_save3]

# --- step 4.2: assemble into one object
breast = list()
breast$dispersion = 0
breast$seed = 0

breast$cell_type_proportion = colSums(spot_composition)/sum(spot_composition)
names(breast$cell_type_proportion) = colnames(spot_composition)
breast$cell_counts = cell_counts
breast$cell_metadata = cell_metadata
breast$spot_counts = spot_counts
breast$spot_coords = spot_coords
breast$spot_composition = spot_composition
special_genes = list()
for (ct in names(breast$cell_type_proportion)){
  special_genes[[ct]] = list(
    marker = c(""), ct_svg = c("")
  )
}
breast$special_genes = special_genes
breast$gene_names = rownames(spot_counts)
str(breast)

# --- final packaged output for breast: no filter2 (no aggressive-filter stage for this dataset)
filter1 = breast
save(filter1, file = "./data/real_packaged/breast/filter1.RData")

# --- Step 4.3: create tr.cell object for spvc
# load("./data/real_packaged/breast/filter1.RData")
# plot(breast$spot_coords[,"x"], breast$spot_coords[,"y"], pch = ".")
# boundary.idx <- identify(breast$spot_coords, order = TRUE)
# breast_bounary = breast$spot_coords[boundary.idx$ind[order(boundary.idx$order)], ]
# points(breast_bounary, type = "l", col = "blue")
# Tr.cell <- TriMesh(breast_bounary, n = 2)
# Tr.cell = remove_colinear(Tr.cell)
# --- final packaged output: triangulation mesh for spvc
# save(Tr.cell, file = "./data/real_packaged/breast/Tr.cell.RData")
