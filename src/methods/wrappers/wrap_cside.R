# src/methods/wrappers/cside.R

# NOTE:
# - This wrapper assumes you already have spot compositions (sp_comp) to import.
# - It creates an RCTD object mainly as a container required by CSIDE.
# - It includes your quadrant-based cell-type filter + iterative weight-threshold filtering.

#' Internal: filter cell types for CSIDE based on region counts + weight threshold
#'
#' @param myRCTD spacexr RCTD object (with @spatialRNA and @results$weights populated)
#' @param cell_type_threshold numeric, minimum total cells per cell type (overall)
#' @param weight_threshold numeric in [0,1], threshold for pixel/spot weights used in filtering
#' @param doublet_mode logical, pass-through to spacexr::aggregate_cell_types
#'
#' @return character vector of cell types to keep
cside_celltype_filter <- function(myRCTD, cell_type_threshold, weight_threshold, doublet_mode = FALSE) {
  if (!requireNamespace("spacexr", quietly = TRUE)) {
    stop("Package 'spacexr' is required but not installed.")
  }
  
  # quadrant thresholding
  region_thresh <- cell_type_threshold / 4
  barcodes <- colnames(myRCTD@spatialRNA@counts)
  coords <- myRCTD@spatialRNA@coords[barcodes, , drop = FALSE]
  
  medx <- stats::median(coords$x)
  medy <- stats::median(coords$y)
  
  r1 <- barcodes[coords$x < medx & coords$y < medy]
  cell_type_filter <- spacexr::aggregate_cell_types(myRCTD, r1, doublet_mode = doublet_mode) >= region_thresh
  
  r2 <- barcodes[coords$x < medx & coords$y > medy]
  cell_type_filter <- cell_type_filter & (spacexr::aggregate_cell_types(myRCTD, r2, doublet_mode = doublet_mode) >= region_thresh)
  
  r3 <- barcodes[coords$x > medx & coords$y > medy]
  cell_type_filter <- cell_type_filter & (spacexr::aggregate_cell_types(myRCTD, r3, doublet_mode = doublet_mode) >= region_thresh)
  
  r4 <- barcodes[coords$x > medx & coords$y < medy]
  cell_type_filter <- cell_type_filter & (spacexr::aggregate_cell_types(myRCTD, r4, doublet_mode = doublet_mode) >= region_thresh)
  
  # overall threshold
  cell_type_count <- spacexr::aggregate_cell_types(myRCTD, barcodes, doublet_mode = doublet_mode)
  cell_types_default <- names(which(cell_type_count >= cell_type_threshold))
  
  # keep only those passing quadrant filter
  cell_types <- intersect(cell_types_default, names(which(cell_type_filter)))
  
  # iterative filtering by weight_threshold
  my_beta <- as.matrix(sweep(myRCTD@results$weights, 1, rowSums(myRCTD@results$weights), "/"))
  
  while (TRUE) {
    if (length(cell_types) == 0) break
    
    # spacexr:::filter_barcodes_cell_types is internal; keep as in your original code.
    res <- spacexr:::filter_barcodes_cell_types(
      barcodes = barcodes,
      cell_types = cell_types,
      my_beta = my_beta[barcodes, , drop = FALSE],
      thresh = weight_threshold
    )
    
    cell_types_remain <- names(which(colSums(res$my_beta) >= cell_type_threshold))
    diff_types <- setdiff(cell_types, cell_types_remain)
    
    if (length(diff_types) == 0) {
      break
    } else {
      warning(
        paste0(
          "CSIDE cell-type filter: removing cell types: ",
          paste(diff_types, collapse = ", "),
          " because they did not contain sufficient pixels passing total cell type weight >= ",
          weight_threshold,
          ". Consider lowering weight_threshold or removing these cell types."
        )
      )
    }
    cell_types <- cell_types_remain
  }
  
  cell_types
}

#' Run CSIDE using imported spot compositions (no deconvolution)
#'
#' @param sp_counts numeric matrix (genes x spots)
#' @param sp_coords data.frame/matrix (spots x 2) with columns named x,y preferred
#' @param sp_comp numeric matrix/data.frame (spots x cell_types) proportions or weights to import
#' @param sc_counts numeric matrix (genes x cells) scRNA reference counts (required by create.RCTD)
#' @param sc_metadata data.frame with rownames = cell IDs, and a column `type` for cell types
#' @param df integer, df for run.CSIDE.nonparam
#' @param cell_type_threshold numeric
#' @param ncores integer
#' @param weight_threshold numeric in [0,1]
#' @param doublet_mode logical for CSIDE (you used FALSE)
#' @param gene_threshold numeric gene threshold for CSIDE
#' @param fdr numeric fdr argument for CSIDE (used internally; you later recompute BH yourself)
#' @param cell_min numeric passed to create.RCTD as CELL_MIN_INSTANCE
#' @param counts_min numeric passed to create.RCTD as counts_MIN
#' @param fc_cutoff_reg numeric passed to create.RCTD
#' @param fc_cutoff numeric passed to create.RCTD
#'
#' @return RCTD object (possibly with @de_results populated if CSIDE ran)
run_cside_import <- function(
    sp_counts,
    sp_coords,
    sp_comp,
    sc_counts,
    sc_metadata,
    df = 15,
    cell_type_threshold = 0,
    ncores = 4,
    weight_threshold = 0.8,
    doublet_mode = FALSE,
    gene_threshold = 5e-5,
    fdr = 0.05,
    cell_min = 25,
    counts_min = 10,
    fc_cutoff_reg = exp(-Inf),
    fc_cutoff = exp(-Inf)
) {
  if (!requireNamespace("spacexr", quietly = TRUE)) {
    stop("Package 'spacexr' is required but not installed.")
  }
  
  if (is.data.frame(sp_counts)) sp_counts <- as.matrix(sp_counts)
  if (is.data.frame(sc_counts)) sc_counts <- as.matrix(sc_counts)
  if (is.data.frame(sp_comp)) sp_comp <- as.matrix(sp_comp)
  if (is.data.frame(sp_coords)) sp_coords <- as.data.frame(sp_coords)
  
  # minimal checks
  if (nrow(sp_coords) != ncol(sp_counts)) {
    stop(sprintf("sp_coords must have nrow == ncol(sp_counts). Got nrow(sp_coords)=%d, ncol(sp_counts)=%d.",
                 nrow(sp_coords), ncol(sp_counts)))
  }
  if (nrow(sp_comp) != ncol(sp_counts)) {
    stop(sprintf("sp_comp must have nrow == ncol(sp_counts). Got nrow(sp_comp)=%d, ncol(sp_counts)=%d.",
                 nrow(sp_comp), ncol(sp_counts)))
  }
  if (!all(c("x", "y") %in% colnames(sp_coords))) {
    # spacexr SpatialRNA expects coords with x,y columns
    colnames(sp_coords)[1:2] <- c("x", "y")
  }
  if (!("type" %in% colnames(sc_metadata))) {
    stop("sc_metadata must contain a column named 'type' with cell type labels.")
  }
  
  # create SpatialRNA puck
  puck <- spacexr::SpatialRNA(coords = sp_coords, counts = sp_counts)
  
  # reference (not used for weights since we import, but required by create.RCTD)
  sc_types <- factor(sc_metadata$type)
  names(sc_types) <- rownames(sc_metadata)
  reference <- spacexr::Reference(sc_counts, sc_types)
  
  myRCTD <- spacexr::create.RCTD(
    spatialRNA = puck,
    reference = reference,
    max_cores = ncores,
    CELL_MIN_INSTANCE = cell_min,
    counts_MIN = counts_min,
    fc_cutoff_reg = fc_cutoff_reg,
    fc_cutoff = fc_cutoff
  )
  
  # import user-provided weights/proportions
  myRCTD <- spacexr::import_weights(myRCTD, sp_comp)
  
  # ensure mode
  myRCTD@config$RCTDmode <- "full"
  
  # pre-check: cell types worth testing
  cell_types_check <- cside_celltype_filter(
    myRCTD = myRCTD,
    cell_type_threshold = cell_type_threshold,
    weight_threshold = weight_threshold,
    doublet_mode = doublet_mode
  )
  
  if (length(cell_types_check) == 0) {
    return(myRCTD)
  }
  
  # run CSIDE with error-handling similar to your original
  myRCTD_temp <- try(
    spacexr::run.CSIDE.nonparam(
      myRCTD,
      df = df,
      cell_types = NULL,               # avoid conflict with internal selection
      gene_threshold = gene_threshold,
      cell_type_threshold = cell_type_threshold,
      fdr = fdr,
      doublet_mode = doublet_mode,
      weight_threshold = weight_threshold
    ),
    silent = TRUE
  )
  
  if (inherits(myRCTD_temp, "try-error")) {
    warning("run.CSIDE.nonparam failed; retrying with more aggressive filtering (default create.RCTD args).")
    myRCTD2 <- spacexr::create.RCTD(spatialRNA = puck, reference = reference, max_cores = ncores)
    
    myRCTD_temp2 <- try(
      spacexr::run.CSIDE.nonparam(
        myRCTD2,
        df = df,
        cell_types = NULL,
        gene_threshold = gene_threshold,
        cell_type_threshold = cell_type_threshold,
        fdr = fdr,
        doublet_mode = doublet_mode,
        weight_threshold = weight_threshold
      ),
      silent = TRUE
    )
    
    if (inherits(myRCTD_temp2, "try-error")) {
      warning("run.CSIDE.nonparam failed even after retry; returning RCTD object without de_results.")
      return(myRCTD)
    } else {
      return(myRCTD_temp2)
    }
  }
  
  myRCTD_temp
}
