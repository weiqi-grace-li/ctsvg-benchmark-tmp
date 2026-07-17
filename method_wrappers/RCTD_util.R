library(spacexr)

###------------------###
### RCTD Run Wrapper ###
###------------------###

#' Run RCTD with provided reference and spatial data
#'
#' @param ref_count [genes x n_cells] matrix of single-cell RNA-seq counts
#' @param ref_type vector of cell type labels (length = n_cells, names = cell names)
#' @param spa_coords [spots x 2] matrix of spatial coordinates (rownames = spot names, colnames = c("x", "y"))
#' @param spa_counts [genes x n_spots] matrix of spatial transcriptomics counts
#' @param max_cores number of CPU cores for parallelization
#' @param doublet_mode one of "full", "doublet", or "multi" for RCTD mode
#'
#' @return normalized cell type weights for each spot [n_spots x n_cell_types]
run_RCTD = function(ref_count, ref_type, spa_coords, spa_counts, max_cores = 4, doublet_mode = "full") {
  # Create reference and spatial RNA objects
  ref = Reference(ref_count, ref_type)
  puck = SpatialRNA(spa_coords, spa_counts)

  # Build and run RCTD object
  myRCTD = create.RCTD(puck, ref, max_cores = max_cores)
  myRCTD = run.RCTD(myRCTD, doublet_mode = doublet_mode)

  # Normalize cell type weights to sum to 1 per spot
  weights = normalize_weights(myRCTD@results$weights)  # [n_spots x n_cell_types]
  return(weights)
}
