# src/methods/wrappers/ctsv.R

#' Run CTSV for cell-type-specific spatial variation
#'
#' @param sp_counts numeric matrix (genes x spots)
#' @param sp_coords numeric matrix/data.frame (spots x 2)
#' @param sp_comp numeric matrix/data.frame (spots x cell types), same order as CTSV expects
#' @param ncores integer number of cores
#'
#' @return CTSV object (with pval/qval slots)
run_ctsv <- function(
    sp_counts,
    sp_coords,
    sp_comp,
    ncores = 4
) {
  if (!requireNamespace("CTSV", quietly = TRUE)) {
    stop("Package 'CTSV' is required but not installed.")
  }
  if (!requireNamespace("SpatialExperiment", quietly = TRUE)) {
    stop("Package 'SpatialExperiment' is required but not installed.")
  }
  if (!requireNamespace("BiocParallel", quietly = TRUE)) {
    stop("Package 'BiocParallel' is required but not installed.")
  }
  
  if (is.data.frame(sp_counts)) sp_counts <- as.matrix(sp_counts)
  if (is.data.frame(sp_comp)) sp_comp <- as.matrix(sp_comp)
  if (is.data.frame(sp_coords)) sp_coords <- as.matrix(sp_coords)
  
  if (nrow(sp_coords) != ncol(sp_counts)) {
    stop(sprintf(
      "sp_coords must have nrow == ncol(sp_counts). Got nrow(sp_coords)=%d, ncol(sp_counts)=%d.",
      nrow(sp_coords), ncol(sp_counts)
    ))
  }
  if (nrow(sp_comp) != ncol(sp_counts)) {
    stop(sprintf(
      "sp_comp must have nrow == ncol(sp_counts). Got nrow(sp_comp)=%d, ncol(sp_counts)=%d.",
      nrow(sp_comp), ncol(sp_counts)
    ))
  }
  
  # build SpatialExperiment
  spe <- SpatialExperiment::SpatialExperiment(
    assays = list(counts = sp_counts),
    spatialCoords = sp_coords
  )
  
  # CTSV parallel backend
  bp <- BiocParallel::MulticoreParam(
    workers = ncores,
    progressbar = TRUE
  )
  
  # run CTSV
  ctsv_obj <- CTSV::CTSV(
    spe = spe,
    prop = sp_comp,
    num_core = ncores,
    BPPARAM = bp
  )
  
  ctsv_obj
}
