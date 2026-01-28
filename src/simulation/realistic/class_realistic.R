pseudo_simulator <- setRefClass("pseudo_simulator",
                                fields = list(
                                  dispersion = "dgCMatrix",
                                  seed = "numeric",
                                  cell_type_proportion = "numeric",
                                  cell_counts = "dgCMatrix",
                                  cell_counts_ori = "dgCMatrix",
                                  cell_counts_predrop = "dgCMatrix", 
                                  cell_metadata = "data.frame",
                                  spot_counts = "matrix",
                                  spot_counts_predrop = "matrix", 
                                  spot_coords = "matrix",
                                  spot_composition = "matrix",
                                  gene_names = "character",
                                  spot_names = "character", 
                                  cell_names = "character",
                                  mu_matrix = "dgCMatrix",
                                  special_genes = "list"
                                  # copula_matrix = "matrix" too big give up 
                                ),
                                
                                methods = list(
                                  initialize = function(data_ori, newcount_list, para_list, seed, special_genes, strat = "C") {
                                    .self$seed = seed
                                    .self$cell_type_proportion = data_ori$cell_type_proportion 
                                    .self$cell_counts_ori = data_ori$cell_counts_ori
                                    .self$cell_metadata = data_ori$cell_metadata
                                    .self$gene_names = data_ori$gene_names
                                    .self$cell_names = data_ori$cell_names
                                    .self$spot_names = data_ori$spot_names
                                    .self$special_genes = special_genes
                                    .self$merge(data_ori, newcount_list, para_list, strat)
                                    .self$form_spot(data_ori)
                                  },
                                  
                                  #' @description Run full simulation.
                                  merge = function(data_ori, newcount_list, para_list, strat){
                                    
                                    # --- cell count matrices 
                                    # replace cell_counts_ori with new counts
                                    newcounts = do.call(cbind, lapply(newcount_list, function(n) n[.self$gene_names, , drop = FALSE]))
                                    newcounts <- as(Matrix::Matrix(newcounts, sparse = TRUE), "CsparseMatrix") # should be sparse, but guard if it's not 
                                    # .self$cell_counts_predrop = .self$cell_counts_ori 
                                    # # both cell_counts_ori and cell_counts_pre_drop are dgcMatrix, assign by column 
                                    # .self$cell_counts_predrop[,match(colnames(newcounts), colnames(.self$cell_counts_predrop))] = newcounts
                                    
                                    .self$cell_counts_predrop = replace_columns_sparse(data_ori$cell_counts_ori, newcounts)
                                    
                                    message(sprintf("Generated cell counts %d, making up %d cells...", ncol(newcounts), (ncol(.self$cell_counts_predrop)-ncol(newcounts))))
                                    
                                    
                                    # --- cell count matrix with drop 
                                    stopifnot(all(dim(data_ori$gene_cell_drop) == dim(.self$cell_counts_predrop)))
                                    if (strat == "A"){
                                      .self$cell_counts = .self$cell_counts_predrop * (1-data_ori$gene_cell_drop)
                                    } else if (strat == "B"){
                                      .self$cell_counts = .self$cell_counts_predrop * data_ori$gene_cell_keep
                                    } else{
                                      K <- (.self$cell_counts_predrop * data_ori$gene_cell_keep) * data_ori$mask_keep
                                      F <- ((.self$cell_counts_predrop %*% Matrix::Diagonal(x = data_ori$cell_fallback)) *
                                              data_ori$mask_fallback)
                                      # A <- Matrix::drop0(K + F)     # dgCMatrix (sparse)
                                      # A@x <- round(A@x)             # R rounds .5 to even
                                      # .self$cell_counts = Matrix::drop0(A)
                                      .self$cell_counts <- Matrix::drop0(K + F)
                                    }
                                    
                                    # --- mu matrix 
                                    mu_matrix = t(do.call(rbind, lapply(para_list, function(x) x$mean_mat[, .self$gene_names, drop=FALSE])))
                                    mu_matrix = as(Matrix::Matrix(mu_matrix, sparse = TRUE), "CsparseMatrix")
                                    
                                    # .self$mu_matrix = .self$cell_counts_ori
                                    # .self$mu_matrix[,match(colnames(mu_matrix), colnames(.self$mu_matrix))] = mu_matrix
                                    .self$mu_matrix = replace_columns_sparse(data_ori$cell_counts_ori, mu_matrix)
                                    
                                    # --- dispersion
                                    # dispersion = t(do.call(rbind, lapply(para_list, function(x) x$sigma_mat[, .self$gene_names, drop=FALSE]))) ## 20251118 revise 
                                    
                                    dispersion = t(
                                      do.call(
                                        rbind, 
                                        lapply(para_list, function(x){
                                          sigma_sub <- x$sigma_mat[, .self$gene_names, drop = FALSE]
                                          
                                          # start from all zeros (or all NA, depending on what you want)
                                          inv_sigma <- matrix(
                                            0,
                                            nrow = nrow(sigma_sub),
                                            ncol = ncol(sigma_sub),
                                            dimnames = dimnames(sigma_sub)
                                          )
                                          
                                          # logical mask: only invert entries that are non-NA and > 0
                                          idx <- !is.na(sigma_sub) & sigma_sub > 0
                                          
                                          inv_sigma[idx] <- 1 / sigma_sub[idx]
                                          
                                          inv_sigma
                                        }
                                        )
                                      )
                                    )
                                    
                                    dispersion = Matrix::Matrix(dispersion, sparse = TRUE)
                                    dispersion = as(dispersion, "CsparseMatrix")
                                    # initialize dispersion to -1 
                                    
                                    # .self$dispersion <- Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0),
                                    #                                          dims = dim(.self$cell_counts_ori),
                                    #                                          dimnames = dimnames(.self$cell_counts_ori))
                                    
                                    # .self$dispersion[,match(colnames(dispersion), colnames(.self$dispersion))] = dispersion  
                                    zero_matrix <- Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0),
                                                                        dims = dim(.self$cell_counts_ori),
                                                                        dimnames = dimnames(.self$cell_counts_ori))
                                    .self$dispersion = replace_columns_sparse(zero_matrix, dispersion)
                                  },
                                  
                                  
                                  form_spot = function(data_ori) {
                                    # spot counts
                                    .self$spot_counts = round(as.matrix(.self$cell_counts %*% data_ori$cell_spot_prop), 0)
                                    .self$spot_counts_predrop = as.matrix(.self$cell_counts_predrop %*% data_ori$cell_spot_prop)
                                    
                                    .self$spot_coords = as.matrix(data_ori$pseudo_meta[.self$spot_names, ])
                                    cell_spot = data_ori$cell_spot_prop
                                    rownames(cell_spot) = data_ori$cell_metadata$type[match(rownames(cell_spot), rownames(data_ori$cell_metadata))]
                                    
                                    # library(Matrix.utils)
                                    # spot_comp = aggregate.Matrix(spot_comp, groupings = rownames(spot_comp), fun = "sum")
                                    # .self$spot_composition = t(apply(spot_comp, MARGIN=2, FUN=function(x) x/sum(x)))
                                    
                                    # set groups as factors 
                                    group <- factor(rownames(cell_spot), levels = names(.self$cell_type_proportion))
                                    
                                    # creates a n_cells by n_cell_type matrix of 1s 
                                    S <- sparseMatrix(
                                      i = seq_len(nrow(cell_spot)),
                                      j = as.integer(group),
                                      x = 1,
                                      dims = c(nrow(cell_spot), nlevels(group))
                                    )
                                    # this is n_cell_type x n_cells, n_cells x n_spots 
                                    spot_comp = t(S) %*% cell_spot
                                    total = colSums(spot_comp)
                                    scale = 1/ifelse(total>0, total, 1)
                                    .self$spot_composition = as.matrix(t(spot_comp %*% Matrix::Diagonal(x = scale)))
                                    colnames(.self$spot_composition) = levels(group)
                                    rownames(.self$spot_composition) = colnames(cell_spot)
                                  }, 
                                  sim_plot = function(gene_name, desc = ""){
                                    expr = log1p(.self$cell_counts[gene_name, ])
                                    expr = scales::rescale(expr)
                                    loc = .self$cell_metadata[, c("x", "y")]
                                    cell_type = .self$cell_metadata$type 
                                    
                                    # 2. Create "Overall" dataframe
                                    df_all = tibble(
                                      Expression = expr,
                                      x = loc[, 1],
                                      y = loc[, 2],
                                      cell_type = "Overall"
                                    )
                                    
                                    # 3. Create per-cell-type dataframe
                                    df_by_type = tibble(
                                      Expression = expr,
                                      x = loc[, 1],
                                      y = loc[, 2],
                                      cell_type = cell_type
                                    )
                                    
                                    # 4. Combine and plot
                                    df_plot = bind_rows(df_all, df_by_type)
                                    
                                    ordered_levels = c("Overall", names(sort(.self$cell_type_proportion, decreasing = TRUE)))
                                    
                                    # Apply factor level ordering
                                    df_plot$cell_type = factor(df_plot$cell_type, levels = ordered_levels)
                                    
                                    overall = ggplot(df_plot, aes(x = x, y = y, color = Expression)) +
                                      geom_point(size = 0.1) +
                                      facet_wrap(~cell_type, nrow = 1) +
                                      scale_colour_gradientn(colors = viridis_pal(option = "magma")(10)) +
                                      coord_fixed(ratio = 1) +
                                      theme(axis.text.x = element_text(angle = 45)) +
                                      ggtitle(paste("Simulated Expression of", gene_name, "(Overall + Cell Types)"))
                                    
                                    
                                    expr = log1p(.self$cell_counts_predrop[gene_name, ])
                                    expr = scales::rescale(expr)
                                    loc = .self$cell_metadata[, c("x", "y")]
                                    cell_type = .self$cell_metadata$type 
                                    
                                    # 2. Create "Overall" dataframe
                                    df_all = tibble(
                                      Expression = expr,
                                      x = loc[, 1],
                                      y = loc[, 2],
                                      cell_type = "Overall"
                                    )
                                    
                                    # 3. Create per-cell-type dataframe
                                    df_by_type = tibble(
                                      Expression = expr,
                                      x = loc[, 1],
                                      y = loc[, 2],
                                      cell_type = cell_type
                                    )
                                    
                                    # 4. Combine and plot
                                    df_plot = bind_rows(df_all, df_by_type)
                                    
                                    ordered_levels = c("Overall", names(sort(.self$cell_type_proportion, decreasing = TRUE)))
                                    
                                    # Apply factor level ordering
                                    df_plot$cell_type = factor(df_plot$cell_type, levels = ordered_levels)
                                    
                                    predrop = ggplot(df_plot, aes(x = x, y = y, color = Expression)) +
                                      geom_point(size = 0.1) +
                                      facet_wrap(~cell_type, nrow = 1) +
                                      scale_colour_gradientn(colors = viridis_pal(option = "magma")(10)) +
                                      coord_fixed(ratio = 1) +
                                      theme(axis.text.x = element_text(angle = 45)) +
                                      ggtitle(paste("Pre-drop Simulated Expression of", gene_name, "(Overall + Cell Types)"))
                                    
                                    # original 
                                    expr = log1p(.self$cell_counts_ori[gene_name, ])
                                    expr =  scales::rescale(expr)
                                    loc = .self$cell_metadata[, c("x", "y")]
                                    cell_type = .self$cell_metadata$type
                                    
                                    # 2. Create "Overall" dataframe
                                    df_all = tibble(
                                      Expression = expr,
                                      x = loc[, 1],
                                      y = loc[, 2],
                                      cell_type = "Overall"
                                    )
                                    
                                    # 3. Create per-cell-type dataframe
                                    df_by_type = tibble(
                                      Expression = expr,
                                      x = loc[, 1],
                                      y = loc[, 2],
                                      cell_type = cell_type
                                    )
                                    
                                    # 4. Combine and plot
                                    df_plot = bind_rows(df_all, df_by_type)
                                    
                                    ordered_levels = c("Overall", names(sort(.self$cell_type_proportion, decreasing = TRUE)))
                                    
                                    # Apply factor level ordering
                                    df_plot$cell_type = factor(df_plot$cell_type, levels = ordered_levels)
                                    
                                    original = ggplot(df_plot, aes(x = x, y = y, color = Expression)) +
                                      geom_point(size = 0.1) +
                                      facet_wrap(~cell_type, nrow = 1) +
                                      scale_colour_gradientn(colors = viridis_pal(option = "magma")(10)) +
                                      coord_fixed(ratio = 1) +
                                      theme(axis.text.x = element_text(angle = 45)) +
                                      ggtitle(paste("Raw Expression of", gene_name, "(Overall + Cell Types)"))
                                    
                                    
                                    # mu 
                                    expr = log1p(.self$mu_matrix[gene_name, ])
                                    expr =  scales::rescale(expr)
                                    loc = .self$cell_metadata[, c("x", "y")]
                                    cell_type = .self$cell_metadata$type
                                    
                                    # 2. Create "Overall" dataframe
                                    df_all = tibble(
                                      Expression = expr,
                                      x = loc[, 1],
                                      y = loc[, 2],
                                      cell_type = "Overall"
                                    )
                                    
                                    # 3. Create per-cell-type dataframe
                                    df_by_type = tibble(
                                      Expression = expr,
                                      x = loc[, 1],
                                      y = loc[, 2],
                                      cell_type = cell_type
                                    )
                                    
                                    # 4. Combine and plot
                                    df_plot = bind_rows(df_all, df_by_type)
                                    
                                    ordered_levels = c("Overall", names(sort(.self$cell_type_proportion, decreasing = TRUE)))
                                    
                                    # Apply factor level ordering
                                    df_plot$cell_type = factor(df_plot$cell_type, levels = ordered_levels)
                                    
                                    mu = ggplot(df_plot, aes(x = x, y = y, color = Expression)) +
                                      geom_point(size = 0.1) +
                                      facet_wrap(~cell_type, nrow = 1) +
                                      scale_colour_gradientn(colors = viridis_pal(option = "magma")(10)) +
                                      coord_fixed(ratio = 1) +
                                      theme(axis.text.x = element_text(angle = 45)) +
                                      ggtitle(paste("Mu of", gene_name, "(Overall + Cell Types)"))
                                    
                                    
                                    original / mu / predrop / overall
                                    
                                    
                                  }
                                )
)


# Replace many columns of a dgCMatrix by rebuilding once (fast, sparse-safe).
# original_counts : dgCMatrix (genes × all_cells)
# replacement_counts : Csparse/dgCMatrix (genes × subset_of_cells_to_replace)
replace_columns_sparse <- function(original_counts, replacement_counts) {
  # 1) map replacement columns into original column positions
  original_cols    <- colnames(original_counts)
  replacement_cols <- colnames(replacement_counts)
  col_map          <- match(replacement_cols, original_cols)  # target positions
  present          <- !is.na(col_map)
  if (!any(present)) return(original_counts)
  if (sum(present) == ncol(original_counts)) return(replacement_counts[,colnames(original_counts)])
  
  cols_to_replace <- col_map[present]
  
  # 2) take ORIGINAL entries that are NOT in columns being replaced
  original_T <- as(original_counts, "dgTMatrix")  # triplet form
  keep_entry <- !(original_T@j + 1L %in% cols_to_replace)
  i_keep <- original_T@i[keep_entry] + 1L
  j_keep <- original_T@j[keep_entry] + 1L
  x_keep <- original_T@x[keep_entry]
  
  # 3) take REPLACEMENT entries and remap their column indices
  repl_T <- as(replacement_counts[, present, drop = FALSE], "dgTMatrix")
  i_repl <- repl_T@i + 1L
  j_repl <- cols_to_replace[repl_T@j + 1L]  # map into original positions
  x_repl <- repl_T@x
  
  # 4) build the final sparse matrix once
  combined <- Matrix::sparseMatrix(
    i = c(i_keep, i_repl),
    j = c(j_keep, j_repl),
    x = c(x_keep, x_repl),
    dims = dim(original_counts),
    dimnames = dimnames(original_counts)
  )
  as(combined, "dgCMatrix")
}
