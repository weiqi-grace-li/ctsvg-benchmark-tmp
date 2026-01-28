# src/methods/wrappers/spvc.R

# This file contains:
# 1) test.spVC.alt
# 2) run_spvc_alt 

# ----------------------------
# spVC entry: test.spVC.alt
# ----------------------------
test.spVC.alt <- function(Y, X = NULL, S, V, Tr, para.cores = 1, scaleX = FALSE,
                          subset = 1:nrow(Y), p.adjust.method = "BH",
                          p.adjust.thresh = 0.05, linear.fit = FALSE,
                          reduced.only = FALSE, twostep = TRUE,
                          filter.min.nonzero = 100, filter.spot.counts = 100,
                          fix.constant = NULL, fix.varying = NULL,
                          size.factors = NULL, ori = FALSE) {
  # data prep ----
  # standardize location points and boundary
  min.x <- min(V[, 1]); max.x <- max(V[, 1])
  min.y <- min(V[, 2]); max.y <- max(V[, 2])
  V[, 1] <- (V[, 1] - min.x) / max.x
  V[, 2] <- (V[, 2] - min.y) / max.y
  S[, 1] <- (S[, 1] - min.x) / max.x
  S[, 2] <- (S[, 2] - min.y) / max.y
  
  # filter spots with low gene expression
  idx.s <- which(colSums(Y) >= filter.spot.counts)
  
  S <- S[idx.s, ]
  ind <- inVT(V, Tr, S[, 1], S[, 2])$ind.inside
  cat("spVC model will use ", length(ind) / ncol(Y) * 100, "% of the original data.\n")
  S.est <- S[ind, ]
  
  if (is.null(X)) {
    X.est <- matrix(1, nrow = length(ind), ncol = 1)
  } else {
    X <- matrix(X, nrow = ncol(Y))
    X <- X[idx.s, ]
    X <- matrix(X, nrow = length(idx.s))
    if (is.null(colnames(X))) colnames(X) <- paste0("X", 1:ncol(X))
    X.est <- cbind(1, X[ind, ])
    colnames(X.est) <- c("0", colnames(X))
    if (scaleX == TRUE) {
      X.est[, -1] <- scale(X.est[, -1])
    }
  }
  colnames(X.est)[1] <- "0"
  
  Y <- Y[, idx.s]
  Y.est <- Y[, ind]
  
  if (is.null(size.factors)) {
    size.factors <- colSums(Y.est)
    size.factors <- size.factors / median(size.factors)
  } else {
    size.factors <- size.factors[ind]
  }
  
  d <- 2
  basis.cell <- basis(V = V, Tr = Tr, d = d, r = 1, Z = as.matrix(S.est))
  
  B <- basis.cell$B
  Q2 <- basis.cell$Q2
  
  TV.all <- list()
  Mat <- BPST::build(d - 1)
  for (k in 1:nrow(Tr)) {
    TV.all[[k]] <- locTV(V[Tr[k, 1], ], V[Tr[k, 2], ], V[Tr[k, 3], ], Mat, d)
  }
  P <- bdiag(TV.all)
  BQ2 <- as.matrix(B %*% Q2)
  BQ2.center <- scale(BQ2, scale = FALSE)
  PQ2 <- as.matrix(crossprod(Q2, P) %*% Q2)
  
  dat.fit <- as.data.frame(X.est)
  pen.list <- list()
  p.X <- ncol(X.est)
  for (ii in 1:p.X) {
    dat.fit[[p.X + ii]] <- kr(matrix(X.est[, ii], ncol = 1), BQ2.center[, -1])
    pen.list[[ii]] <- list(PQ2[2:ncol(Q2), 2:ncol(Q2)])
  }
  names(dat.fit)[1:p.X] <- paste0("beta_", colnames(X.est))
  names(dat.fit)[p.X + 1:p.X] <- paste0("gamma_", colnames(X.est))
  names(pen.list) <- paste0("gamma_", colnames(X.est))
  
  if (is.null(rownames(Y.est))) rownames(Y.est) <- paste0("gene", 1:nrow(Y.est))
  idx <- names(which(apply(
    Y.est[subset, , drop = FALSE], 1,
    function(x) sum(x != 0) > filter.min.nonzero
  )))
  cat("Conducting tests for", length(idx), " genes.\n")
  print.idx <- 1:length(idx)
  names(print.idx) <- idx
  
  # full spVC and evaluate the significance of each individual component ----
  results.full <- NULL
  if (twostep == FALSE) {
    # remove gamma_0 and beta_0 unless ori=TRUE
    if (ori) {
      formula.full <- stats::as.formula(paste0("Y ~ 0 + ", paste0(names(dat.fit), collapse = " + ")))
    } else {
      remove <- which(names(dat.fit) %in% c("beta_0", "gamma_0"))
      formula.full <- stats::as.formula(paste0("Y ~ 0 + ", paste0(names(dat.fit)[-remove], collapse = " + ")))
    }
    
    results.full <- parallel::mclapply(
      idx,
      mc.cores = para.cores,
      FUN = function(x) {
        if (print.idx[x] %% 500 == 0) {
          cat("Fitting Full Model for Gene", print.idx[x], "out of", length(idx), "genes.\n")
        }
        fit.spVC.alt(
          formula.full,
          Y.iter = as.vector(Y.est[x, ]),
          dat.fit = dat.fit,
          size.factors = size.factors,
          pen.list = pen.list
        )
      }
    )
    names(results.full) <- idx
  }
  
  # prepare output ----
  if (twostep == FALSE) {
    results <- list(
      results.full = results.full,
      BQ2.center.est = colMeans(BQ2)
    )
  } else {
    results <- list(
      results.varying = NULL,
      results.reduced = NULL,
      BQ2.center.est = colMeans(BQ2)
    )
  }
  
  results
}

# ----------------------------
# Wrapper: run_spvc_alt
# ----------------------------

#' Run spVC (alt) using triangulation inputs
#'
#' @param sp_counts numeric matrix (genes x spots)
#' @param sp_coords numeric matrix/data.frame (spots x 2)
#' @param sp_comp numeric matrix/data.frame (spots x cell types) or NULL
#' @param Tr.cell list with $V and $Tr
#' @param ncores integer
#' @param twostep logical; FALSE for full model (matches your current usage)
#' @param filter_min_nonzero integer
#' @param filter_spot_counts integer
#' @param ori logical passed to test.spVC.alt
#' @param scaleX logical
#' @param subset integer vector of gene indices to test
#'
#' @return list from test.spVC.alt (e.g., results.full)
run_spvc_alt <- function(
    sp_counts,
    sp_coords,
    sp_comp,
    Tr.cell,
    ncores = 4,
    twostep = FALSE,
    filter_min_nonzero = 100,
    filter_spot_counts = 100,
    ori = FALSE,
    scaleX = FALSE,
    subset = NULL
) {
  # dependency checks for functions used inside test.spVC.alt
  # (We keep these checks here so errors are clearer at runtime.)
  pkgs <- c("BPST", "Matrix", "parallel")
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      stop("Package '", p, "' is required but not installed.")
    }
  }
  # inVT, basis, locTV, bdiag, kr, fit.spVC.alt are assumed to exist
  needed_funs <- c("inVT", "basis", "locTV", "bdiag", "kr", "fit.spVC.alt")
  missing <- needed_funs[!vapply(needed_funs, exists, logical(1), mode = "function", inherits = TRUE)]
  if (length(missing) > 0) {
    stop(
      "Missing required spVC helper functions: ",
      paste(missing, collapse = ", "),
      ". If these are in your spVC source files, source them before calling run_spvc_alt()."
    )
  }
  
  if (is.data.frame(sp_counts)) sp_counts <- as.matrix(sp_counts)
  if (is.data.frame(sp_coords)) sp_coords <- as.matrix(sp_coords)
  if (!is.null(sp_comp) && is.data.frame(sp_comp)) sp_comp <- as.matrix(sp_comp)
  
  if (nrow(sp_coords) != ncol(sp_counts)) {
    stop(sprintf("sp_coords must have nrow == ncol(sp_counts). Got %d vs %d.",
                 nrow(sp_coords), ncol(sp_counts)))
  }
  if (!is.null(sp_comp) && nrow(sp_comp) != ncol(sp_counts)) {
    stop(sprintf("sp_comp must have nrow == ncol(sp_counts). Got %d vs %d.",
                 nrow(sp_comp), ncol(sp_counts)))
  }
  if (is.null(Tr.cell$V) || is.null(Tr.cell$Tr)) {
    stop("Tr.cell must contain $V and $Tr.")
  }
  
  Tr.cell = remove_colinear(Tr.cell)
  
  # Avoid side effects: test.spVC.alt rescales V and S in-place
  V <- as.matrix(Tr.cell$V)
  Tr <- as.matrix(Tr.cell$Tr)
  S <- as.matrix(sp_coords)
  
  if (is.null(subset)) subset <- seq_len(nrow(sp_counts))
  
  test.spVC.alt(
    Y = sp_counts,
    X = sp_comp,
    S = S,
    V = V,
    Tr = Tr,
    para.cores = ncores,
    scaleX = scaleX,
    subset = subset,
    twostep = twostep,
    filter.min.nonzero = filter_min_nonzero,
    filter.spot.counts = filter_spot_counts,
    ori = ori
  )
}


are_colinear = function(A, B, C, tol = 1e-6){
  abs(det(rbind(B-A, C-A))) < tol
}

remove_colinear = function(Tr.cell){
  save = c()
  for (i in 1:nrow(Tr.cell$Tr)){
    A = Tr.cell$V[Tr.cell$Tr[i, 1], ]
    B = Tr.cell$V[Tr.cell$Tr[i, 2], ]
    C = Tr.cell$V[Tr.cell$Tr[i, 3], ]
    if (!are_colinear(A, B, C)){
      save = c(save, i)
    }
  }
  Tr.cell$Tr = Tr.cell$Tr[save, ]
  return(Tr.cell)
}






