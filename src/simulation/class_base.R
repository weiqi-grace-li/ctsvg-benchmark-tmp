# classes/BaseDataset.R

# If you use sparse matrices anywhere:
suppressPackageStartupMessages({
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package 'Matrix' is required for BaseDataset classes.")
  }
})

# ---- BaseDataset: the shared payload across all datasets ----
BaseDataset <- setRefClass(
  "BaseDataset",
  fields = list(
    # core payload
    gene_names        = "character",
    cell_counts       = "ANY",         # matrix or dgCMatrix
    cell_metadata     = "data.frame",  # must have rownames = cells
    spot_counts       = "ANY",         # matrix or dgCMatrix (you mostly use matrix)
    spot_coords       = "matrix",      # spots x 2
    spot_composition  = "ANY",         # matrix (spots x celltypes)
    # optional identifiers
    cell_names        = "character",
    spot_names        = "character",
    dataset_id        = "character"
  ),
  methods = list(
    initialize = function(
    gene_names,
    cell_counts,
    cell_metadata,
    spot_counts,
    spot_coords,
    spot_composition,
    cell_names = character(),
    spot_names = character(),
    dataset_id = ""
    ) {
      .self$gene_names       <- gene_names
      .self$cell_counts      <- cell_counts
      .self$cell_metadata    <- cell_metadata
      .self$spot_counts      <- spot_counts
      .self$spot_coords      <- spot_coords
      .self$spot_composition <- spot_composition
      .self$cell_names       <- cell_names
      .self$spot_names       <- spot_names
      .self$dataset_id       <- dataset_id
      .self$validate_basic()
      invisible(.self)
    },
    
    validate_basic = function() {
      # gene_names
      stopifnot(is.character(.self$gene_names), length(.self$gene_names) > 0)
      
      # cell metadata must have rownames
      stopifnot(is.data.frame(.self$cell_metadata))
      if (is.null(rownames(.self$cell_metadata)) || anyNA(rownames(.self$cell_metadata))) {
        stop("BaseDataset: cell_metadata must have non-NA rownames (cell IDs).")
      }
      
      # spot_coords: matrix with 2 columns
      stopifnot(is.matrix(.self$spot_coords), ncol(.self$spot_coords) == 2)
      if (is.null(rownames(.self$spot_coords)) || anyNA(rownames(.self$spot_coords))) {
        stop("BaseDataset: spot_coords must have non-NA rownames (spot IDs).")
      }
      
      invisible(TRUE)
    },
    
    get_n_genes = function() length(.self$gene_names),
    get_n_cells = function() nrow(.self$cell_metadata),
    get_n_spots = function() nrow(.self$spot_coords)
  )
)
