# src/methods/wrappers/stance.R

#' Run STANCE (two-stage) ctSVG test
#'
#' @param count numeric matrix (genes x spots)
#' @param pos numeric matrix/data.frame (spots x 2)
#' @param comp numeric matrix/data.frame (spots x cell_types) cell-type proportions per spot
#' @param gene_thres numeric, gene filtering threshold (STANCE::data_preprocess)
#' @param spot_thres numeric/integer, spot filtering threshold (STANCE::data_preprocess)
#' @param normalized logical, whether input is already normalized (STANCE::data_preprocess)
#' @param correction logical, whether to apply correction in STANCE tests
#' @param ncores integer, number of cores for stage-2 test
#' @param utsvg_thres numeric, adjusted p-value threshold for selecting utSVGs from stage 1
#' @param topn_utsvg integer or NULL; if not NULL, take top N utSVGs instead of using utsvg_thres
#' @param pv_adjust character, adjustment method used by runTest1 (you used "BY")
#'
#' @return STANCE object with Test_1 and Test_2 populated
run_stance <- function(
    count,
    pos,
    comp,
    gene_thres = 0.05,
    spot_thres = 10,
    normalized = FALSE,
    correction = FALSE,
    ncores = 4,
    utsvg_thres = 0.05,
    topn_utsvg = NULL,
    pv_adjust = "BY"
) {
  if (!requireNamespace("STANCE", quietly = TRUE)) {
    stop("Package 'STANCE' is required but not installed.")
  }
  
  if (is.data.frame(count)) count <- as.matrix(count)
  if (is.data.frame(comp)) comp <- as.matrix(comp)
  if (is.data.frame(pos)) pos <- as.matrix(pos)
  
  if (nrow(pos) != ncol(count)) {
    stop(sprintf("pos must have nrow == ncol(count). Got nrow(pos)=%d, ncol(count)=%d.",
                 nrow(pos), ncol(count)))
  }
  if (nrow(comp) != ncol(count)) {
    stop(sprintf("comp must have nrow == ncol(count). Got nrow(comp)=%d, ncol(count)=%d.",
                 nrow(comp), ncol(count)))
  }
  
  # create stance object
  stance_obj <- STANCE::creatSTANCEobject(counts = count, pos = pos, prop = comp)
  
  # preprocess
  stance_obj <- STANCE::data_preprocess(
    object = stance_obj,
    gene.threshold = gene_thres,
    spot.threshold = spot_thres,
    normalized = normalized
  )
  
  # stage 1: global (utSVG) screening
  stance_obj <- STANCE::build_kernelMatrix(object = stance_obj)
  stance_obj <- STANCE::runTest1(object = stance_obj, correction = correction, pv.adjust = pv_adjust)
  
  gene.list <- rownames(stance_obj@gene_expression)
  
  # pick utSVGs for stage 2
  if (is.null(topn_utsvg)) {
    utSVG.list <- gene.list[stance_obj@Test_1$p_value_adj < utsvg_thres]
  } else {
    topn_utsvg <- min(as.integer(topn_utsvg), length(gene.list))
    ord <- order(stance_obj@Test_1$p_value_adj)
    utSVG.list <- gene.list[ord[seq_len(topn_utsvg)]]
    best_i <- ord[1]
    message(paste0(
      "Selecting only top ", topn_utsvg,
      " | best gene ", gene.list[best_i],
      " | best p_adj ", stance_obj@Test_1$p_value_adj[best_i]
    ))
  }
  
  message("##-----")
  message(paste0(
    "utSVG selected ", length(utSVG.list),
    " | ", length(gene.list), " genes for stage 2 test."
  ))
  message("##-----")
  
  # stage 2: ct-specific test (all cell types by default)
  stance_obj <- STANCE::runTest2(
    object = stance_obj,
    Genes_to_test = utSVG.list,
    Cell_types_to_test = NULL,
    correction = correction,
    ncores = ncores
  )
  
  stance_obj
}
