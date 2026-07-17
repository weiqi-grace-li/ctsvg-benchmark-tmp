library(Matrix)

source("./util.R")

ncores     <- 50
run_method <- c("celina", "stance", "spvc", "ctsv", "cside", "mmm")
save_dir <- paste0("./data/detection_results/realistic_validation/")
seeds        <- c(74, 153, 321, 561, 997)
sim_save_dir <- "./data/realistic_packaged/validation/"
path         <- "./data/realistic_preprocessed/validation/"

# -------------------------------------------------------
# Load count data
# -------------------------------------------------------
visium    <- t(as.matrix(read.csv(paste0(path, "visium_count.csv"), row.names = 1)))
pseudob   <- t(read.csv(paste0(path, "pseudo_count.csv"), row.names = 1))
pseudob2  <- t(read.csv(paste0(path, "pseudo_count_rep2B.csv"), row.names = 1))
spot_meta <- read.csv(paste0(path, "pseudo_meta.csv"), row.names = 1)

synthetic    <- list()
for (seed in seeds) {
  load(paste0(sim_save_dir, "validation_new_seed", seed, ".RData"))  # -> data_new
  synthetic[[paste0("A_seed", seed)]] <- data_new$spot_counts
}

# -------------------------------------------------------
# Clean empty spots / genes
# -------------------------------------------------------
remove_empty <- function(spot_counts) {
  empty_genes <- which(rowSums(spot_counts) == 0)
  message(paste0(length(empty_genes), " genes out of ", nrow(spot_counts), " has no reads..."))
  empty_spots <- which(colSums(spot_counts) == 0)
  message(paste0(length(empty_spots), " spots out of ", ncol(spot_counts), " has no reads..."))
  cleaned <- spot_counts
  if (length(empty_genes) > 0) cleaned <- cleaned[-empty_genes, ]
  if (length(empty_spots) > 0) cleaned <- cleaned[, -empty_spots]
  Matrix(cleaned, sparse = TRUE)
}

visium_clean   <- remove_empty(visium)
pseudob_clean  <- remove_empty(pseudob)
pseudob2_clean <- remove_empty(pseudob2)

synthetic_clean           <- list()
genes_synthetic_clean     <- c()
spotnames_synthetic_clean <- c()
for (name in names(synthetic)) {
  synthetic_clean[[name]] <- remove_empty(synthetic[[name]])
  if (length(synthetic_clean) == 1) {
    genes_synthetic_clean     <- rownames(synthetic_clean[[name]])
    spotnames_synthetic_clean <- colnames(synthetic_clean[[name]])
  } else {
    genes_synthetic_clean     <- intersect(genes_synthetic_clean, rownames(synthetic_clean[[name]]))
    spotnames_synthetic_clean <- intersect(spotnames_synthetic_clean, colnames(synthetic_clean[[name]]))
  }
}

# -------------------------------------------------------
# Intersect common genes and spots across all datasets
# -------------------------------------------------------
common_genes <- Reduce(intersect, list(
  rownames(visium_clean),
  rownames(pseudob_clean), rownames(pseudob2_clean),
  genes_synthetic_clean
))
message(paste0("Selected ", length(common_genes), " genes..."))

common_spots <- Reduce(intersect, list(
  colnames(visium_clean),
  colnames(pseudob_clean), colnames(pseudob2_clean),
  spotnames_synthetic_clean
))
message(paste0("Selected ", length(common_spots), " spots..."))

visium_clean   <- visium_clean[common_genes, common_spots]
pseudob_clean  <- pseudob_clean[common_genes, common_spots]
pseudob2_clean <- pseudob2_clean[common_genes, common_spots]
for (name in names(synthetic_clean)) {
  synthetic_clean[[name]] <- synthetic_clean[[name]][common_genes, common_spots]
}

# -------------------------------------------------------
# Shared spatial inputs (from C_seed997 reference)
# -------------------------------------------------------
spot_comp        <- data_new$spot_composition[common_spots, ]
spot_coords      <- as.matrix(spot_meta[common_spots, 1:2])
colnames(spot_coords) <- c("x", "y")

sc_count    <- data_new$cell_counts_ori
sc_metadata <- data_new$cell_metadata

# -------------------------------------------------------
# sim_obj
# -------------------------------------------------------
sim_obj <- list(
  gene_names           = common_genes,
  cell_metadata        = data_new$cell_metadata,
  special_genes        = setNames(
    lapply(unique(data_new$cell_metadata$type), function(ct)
      list(marker = c("1"), ctsvg = c("1"))),
    unique(data_new$cell_metadata$type)
  ),
  seed                 = -1,
  cell_type_proportion = data_new$cell_type_proportion,
  dispersion           = -1,
  spot_composition     = spot_comp
)

# -------------------------------------------------------
# spVC triangulation (loaded only when needed)
# -------------------------------------------------------
spvc_tri <- NULL
if (any(c("spvc") %in% run_method)) {
  library(Triangulation)
  load(paste0(sim_save_dir, "validation_boundary.RData"))
  spvc_tri <- TriMesh(breast_boundary, n = 2)
}

# -------------------------------------------------------
# Collect all count datasets
# -------------------------------------------------------
count_data_list <- list(
  visium   = as.matrix(visium_clean),  
  pseudob  = as.matrix(pseudob_clean),  
  pseudob2 = as.matrix(pseudob2_clean)
)
for (name in names(synthetic_clean)) {
  count_data_list[[name]] <- as.matrix(synthetic_clean[[name]])
}

# -------------------------------------------------------
# Run all methods on each dataset
# -------------------------------------------------------
for (name in names(count_data_list)) {
  message(paste0("=== Running dataset: ", name, " ==="))
  run_all_tests(
    sp_count    = count_data_list[[name]],
    sp_comp     = spot_comp,
    sp_coords   = spot_coords,
    sc_count    = sc_count,
    sc_metadata = sc_metadata,
    spvc_tri    = spvc_tri,
    sim_name    = name,
    save_dir    = save_dir,
    sim_obj     = sim_obj,
    ncores      = ncores,
    run_method  = run_method
  )
}
