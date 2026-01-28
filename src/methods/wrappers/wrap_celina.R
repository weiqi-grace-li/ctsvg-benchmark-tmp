# src/methods/wrappers/celina.R

#' Run CELINA on spot-level counts + locations + cell-type proportions
#'
#' @param count numeric matrix (genes x spots) spot gene expression matrix
#' @param pos numeric matrix/data.frame (spots x 2) spatial coordinates with colnames c("x","y") preferred
#' @param comp numeric matrix/data.frame (spots x cell_types) cell-type proportions per spot
#' @param sc_count numeric matrix (genes x cells) scRNA-seq counts
#' @param sc_type vector of length ncol(sc_count) giving cell type label per cell (names optional)
#' @param approximation logical, whether to use CELINA approximation in kernel calculation
#' @param ncores integer number of cores
#' @param threshold numeric gene filtering threshold for preprocess_input()
#' @param cell_types_to_test character vector; if NULL uses colnames(comp)
#' @param project character project name stored in CELINA object
#'
#' @return CELINA object (with results in @result)
run_celina_spot <- function(
    count,
    pos,
    comp,
    sc_count,
    sc_type,
    approximation = FALSE,
    ncores = 4,
    threshold = 5e-5,
    cell_types_to_test = NULL,
    project = "Sim"
) {
  if (!requireNamespace("CELINA", quietly = TRUE)) {
    stop("Package 'CELINA' is required but not installed.")
  }
  
  # --- basic checks
  if (is.data.frame(count)) count <- as.matrix(count)
  if (is.data.frame(sc_count)) sc_count <- as.matrix(sc_count)
  if (is.data.frame(comp)) comp <- as.matrix(comp)
  if (is.data.frame(pos)) pos <- as.matrix(pos)
  
  if (is.null(colnames(comp))) stop("comp must have colnames = cell type names.")
  if (nrow(pos) != ncol(count)) {
    stop(sprintf("pos must have nrow == ncol(count). Got nrow(pos)=%d, ncol(count)=%d.",
                 nrow(pos), ncol(count)))
  }
  if (nrow(comp) != ncol(count)) {
    stop(sprintf("comp must have nrow == ncol(count). Got nrow(comp)=%d, ncol(count)=%d.",
                 nrow(comp), ncol(count)))
  }
  if (nrow(sc_count) != nrow(count)) {
    stop(sprintf("sc_count must have same genes (rows) as count. Got nrow(sc_count)=%d, nrow(count)=%d.",
                 nrow(sc_count), nrow(count)))
  }
  if (length(sc_type) != ncol(sc_count)) {
    stop(sprintf("sc_type length must equal ncol(sc_count). Got length(sc_type)=%d, ncol(sc_count)=%d.",
                 length(sc_type), ncol(sc_count)))
  }
  
  # --- preprocess sc_type names (CELINA expects labels aligned to scRNA cells)
  if (is.null(names(sc_type))) {
    names(sc_type) <- colnames(sc_count)
  } else {
    # if names exist but don't match, don't silently reorder; just warn
    if (!all(names(sc_type) %in% colnames(sc_count))) {
      warning("Names(sc_type) do not match colnames(sc_count). Proceeding with sc_type as-is.")
    }
  }
  
  # --- cell types to test
  if (is.null(cell_types_to_test)) {
    cell_types_to_test <- unique(colnames(comp))
  }
  
  # --- create and run CELINA
  celina_obj <- CELINA::Create_Celina_Object(
    celltype_mat = t(comp),                # CELINA expects celltype x spot
    gene_expression_mat = count,            # genes x spots
    location = pos,                         # spots x 2
    covariates = NULL,
    project = project
  )
  
  celina_obj <- CELINA::preprocess_input(
    celina_obj,
    cell_types_to_test = cell_types_to_test,
    scRNA_count = sc_count,
    sc_cell_type_labels = as.factor(sc_type),
    threshold = threshold
  )
  
  celina_obj = CELINA::Calculate_Kernel(celina_obj, approximation = approximation)
  
  celina_obj = CELINA::Testing_interaction_all(
    celina_obj,
    celltype_to_test = cell_types_to_test,
    num_cores = ncores
  )
  
  celina_obj
}
