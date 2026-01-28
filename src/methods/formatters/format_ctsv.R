# src/methods/formatters/ctsv_format.R

ctsv_details_to_df <- function(
    sim_obj,
    ctsv_obj,
    sim_name,
    deconv_method
) {
  all_gene <- sim_obj$gene_names
  
  # Reconstruct per-cell-type result tables (matches your original)
  K <- ncol(sim_obj$spot_composition)
  
  ctsv_result <- list()
  for (k in 1:K) {
    cell_result <- data.frame(matrix(NA, ncol = 6, nrow = nrow(ctsv_obj$qval)))
    rownames(cell_result) <- rownames(ctsv_obj$qval)
    colnames(cell_result) <- c("p_value", "q_value", "p_value_s1", "p_value_s2", "q_value_s1", "q_value_s2")
    
    cell_result$p_value <- apply(ctsv_obj$pval[, c(k, k + K)], 1, min)
    cell_result$q_value <- apply(ctsv_obj$qval[, c(k, k + K)], 1, min)
    cell_result$p_value_s1 <- ctsv_obj$pval[, k]
    cell_result$p_value_s2 <- ctsv_obj$pval[, k + K]
    cell_result$q_value_s1 <- ctsv_obj$qval[, k]
    cell_result$q_value_s2 <- ctsv_obj$qval[, k + K]
    
    ctsv_result[[colnames(sim_obj$spot_composition)[k]]] <- cell_result
  }
  
  # Identify significant genes per cell type (matches your original)
  ctsv_sig_genes <- CTSV::svGene(ctsv_obj$qval)
  ctsv_sig_genes <- ctsv_sig_genes$SVGene
  names(ctsv_sig_genes) <- colnames(sim_obj$spot_composition)
  
  cell_types <- unique(sim_obj$cell_metadata$type)
  results <- list()
  
  for (ct in cell_types) {
    # sometimes cell types are not tested
    if (ct %in% names(ctsv_result)) {
      result <- ctsv_result[[ct]]
      gene_tested <- rownames(result)
      
      # double check the sig result fits with the model (exactly as your original)
      sig_genes_user <- sort(rownames(result)[which(result$q_value < 0.05)])
      sig_genes_method <- sort(ctsv_sig_genes[[ct]])
      if (any(sig_genes_user != sig_genes_method)) {
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
    
    # allocate output (exactly 15 columns as your original)
    result_ct <- as.data.frame(matrix(NA, nrow = length(all_gene), ncol = 15))
    colnames(result_ct) <- c(
      "simulation_name", "test_method", "deconv_method", "seed", "cell_type", "cell_proportion",
      "gene_name", "dispersion", "gene_type", "p_value", "p_adj",
      "p_value_s1", "p_value_s2", "q_value_s1", "q_value_s2"
    )
    
    result_ct$simulation_name <- sim_name
    result_ct$test_method <- "ctsv"
    result_ct$deconv_method <- deconv_method
    result_ct$seed <- sim_obj$seed
    result_ct$cell_type <- ct
    result_ct$cell_proportion <- sim_obj$cell_type_proportion[ct]
    result_ct$gene_name <- all_gene
    
    if (length(sim_obj$dispersion) == 1) {
      result_ct$dispersion <- sim_obj$dispersion
    } else {
      result_ct$dispersion <- rowMeans(sim_obj$dispersion[all_gene, ])
    }
    
    # fill gene_type and stats exactly as your original logic
    for (gene in all_gene) {
      gene_labels <- names(check_list)[sapply(check_list, function(lst) gene %in% lst)]
      gene_label <- if (length(gene_labels) == 0) "Null" else paste(gene_labels, collapse = ",")
      
      idx <- which(result_ct$gene_name == gene)
      result_ct[idx, "gene_type"] <- gene_label
      
      if (gene %in% gene_tested) {
        result_ct[idx, "p_value"] <- result[gene, "p_value"]
        result_ct[idx, "p_adj"] <- result[gene, "q_value"]
        result_ct[idx, "p_value_s1"] <- result[gene, "p_value_s1"]
        result_ct[idx, "p_value_s2"] <- result[gene, "p_value_s2"]
        result_ct[idx, "q_value_s1"] <- result[gene, "q_value_s1"]
        result_ct[idx, "q_value_s2"] <- result[gene, "q_value_s2"]
      } else {
        result_ct[idx, "p_value"] <- -1
        result_ct[idx, "p_adj"] <- -1
        result_ct[idx, "p_value_s1"] <- -1
        result_ct[idx, "p_value_s2"] <- -1
        result_ct[idx, "q_value_s1"] <- -1
        result_ct[idx, "q_value_s2"] <- -1
      }
    }
    
    results[[ct]] <- result_ct
  }
  
  long_df <- do.call(rbind, results)
  long_df
}
