# src/methods/formatters/cside_format.R

cside_details_to_df <- function(
    sim_obj,
    cside_obj,
    sim_name,
    deconv_method,
    cell_threshold,
    df,
    weight_threshold
) {
  all_gene <- sim_obj$gene_names
  
  cside_results <- cside_obj@de_results$all_gene_list
  cside_convmat <- cside_obj@de_results$gene_fits$con_mat
  
  cell_types <- unique(sim_obj$cell_metadata$type)
  results <- list()
  
  for (ct in cell_types) {
    # sometimes cell types are not tested
    if (ct %in% names(cside_results)) {
      result <- cside_results[[ct]]
      gene_tested <- rownames(result)
      
      # create BH adjusted p value (exactly as your original)
      p_adj <- p.adjust(result$p_val, method = "BH")
      result <- cbind(result, p_adj)
      
      # make sure we have found the same significant genes as cside (exactly as your original)
      sig_genes_user <- sort(rownames(result)[p_adj < 0.05])
      sig_genes_cside <- sort(rownames(cside_obj@de_results$sig_gene_list[[ct]]))
      if (any(sig_genes_user != sig_genes_cside)) {
        stop("BH adjusted pvalue result doesn't match CSIDE")
      }
    } else {
      result <- NULL
      gene_tested <- character(0)
    }
    
    # gene categories (same as your original)
    marker <- sim_obj$special_genes[[ct]]$marker
    ctsvg <- sim_obj$special_genes[[ct]]$ct_svg
    other_marker <- c()
    other_ctsvg <- c()
    for (oct in setdiff(cell_types, ct)) {
      other_marker <- c(other_marker, sim_obj$special_genes[[oct]]$marker)
      other_ctsvg <- c(other_ctsvg, sim_obj$special_genes[[oct]]$ct_svg)
    }
    check_list <- list(marker = marker, ctsvg = ctsvg, other_marker = other_marker, other_ctsvg = other_ctsvg)
    
    # create result table and prefill items (exactly 15 columns as your original)
    result_ct <- as.data.frame(matrix(NA, nrow = length(all_gene), ncol = 15))
    colnames(result_ct) <- c(
      "simulation_name", "test_method", "deconv_method", "cell_threshold", "df", "weight_threshold",
      "seed", "cell_type", "cell_proportion",
      "gene_name", "dispersion", "gene_type", "p_value", "p_adj", "convergence"
    )
    
    result_ct$simulation_name <- sim_name
    result_ct$test_method <- "cside"
    result_ct$deconv_method <- deconv_method
    result_ct$cell_threshold <- cell_threshold
    result_ct$df <- df
    result_ct$weight_threshold <- weight_threshold
    result_ct$seed <- sim_obj$seed
    result_ct$cell_type <- ct
    result_ct$cell_proportion <- sim_obj$cell_type_proportion[ct]
    result_ct$gene_name <- all_gene
    
    if (length(sim_obj$dispersion) == 1) {
      result_ct$dispersion <- sim_obj$dispersion
    } else {
      result_ct$dispersion <- rowMeans(sim_obj$dispersion[all_gene, ])
    }
    
    # generate gene type and fill p-values (exactly as your original)
    for (gene in all_gene) {
      gene_labels <- names(check_list)[sapply(check_list, function(lst) gene %in% lst)]
      gene_label <- if (length(gene_labels) == 0) "Null" else paste(gene_labels, collapse = ",")
      
      idx <- which(result_ct$gene_name == gene)
      result_ct[idx, "gene_type"] <- gene_label
      
      if (gene %in% gene_tested) {
        result_ct[idx, "p_value"] <- result[gene, "p_val"]
        result_ct[idx, "p_adj"] <- result[gene, "p_adj"]
      } else {
        result_ct[idx, "p_value"] <- -1
        result_ct[idx, "p_adj"] <- -1
      }
      
      if ((gene %in% rownames(cside_convmat)) & (ct %in% colnames(cside_convmat))) {
        result_ct[idx, "convergence"] <- cside_convmat[gene, ct]
      } else {
        result_ct[idx, "convergence"] <- "unknown"
      }
    }
    
    results[[ct]] <- result_ct
  }
  
  long_df <- do.call(rbind, results)
  long_df
}
