# src/methods/formatters/stance_format.R

stance_details_to_df <- function(
    sim_obj,
    stance_obj,
    sim_name,
    deconv_method,
    test_method = "stance"
) {
  all_gene <- sim_obj$gene_names
  
  # recreate what has been tested in stance test 2 (matches your original logic)
  gene.list <- rownames(stance_obj@gene_expression)
  utsvgs <- gene.list[stance_obj@Test_1$p_value_adj < 0.05]  # uses fixed 0.05 as in your code
  
  cell_types <- unique(sim_obj$cell_metadata$type)
  results <- list()
  
  for (ct in cell_types) {
    result <- stance_obj@Test_2[[ct]]
    
    # gene categories (same as your original)
    marker <- sim_obj$special_genes[[ct]]$marker
    ctsvg <- sim_obj$special_genes[[ct]]$ct_svg
    other_marker <- c()
    other_ctsvg <- c()
    for (oct in setdiff(cell_types, ct)) {
      other_marker <- c(other_marker, sim_obj$special_genes[[oct]]$marker)
      other_ctsvg <- c(other_ctsvg, sim_obj$special_genes[[oct]]$ct_svg)
    }
    check_list <- list(
      marker = marker,
      ctsvg = ctsvg,
      other_marker = other_marker,
      other_ctsvg = other_ctsvg
    )
    
    # pre-allocate (exactly 10 columns as your original)
    result_ct <- as.data.frame(matrix(NA, nrow = length(all_gene), ncol = 10))
    colnames(result_ct) <- c(
      "simulation_name", "test_method", "deconv_method", "seed", "cell_type", "cell_proportion",
      "gene_name", "dispersion", "gene_type", "p_value"
    )
    
    result_ct$simulation_name <- sim_name
    result_ct$test_method <- test_method
    result_ct$deconv_method <- deconv_method
    result_ct$seed <- sim_obj$seed
    result_ct$cell_type <- ct
    result_ct$cell_proportion <- sim_obj$cell_type_proportion[ct]
    result_ct$gene_name <- all_gene
    
    if (length(sim_obj$dispersion) == 1) {
      result_ct$dispersion <- sim_obj$dispersion
    } else {
      result_ct$dispersion <- rowMeans(sim_obj$dispersion[all_gene, ])  # same TODO comment not needed here
    }
    
    # fill gene_type and p_value exactly as your original logic
    for (gene in all_gene) {
      gene_labels <- names(check_list)[sapply(check_list, function(lst) gene %in% lst)]
      gene_label <- if (length(gene_labels) == 0) "Null" else paste(gene_labels, collapse = ",")
      
      idx <- which(result_ct$gene_name == gene)
      result_ct[idx, "gene_type"] <- gene_label
      
      if (gene %in% utsvgs) {
        if (!is.matrix(result) & !is.data.frame(result)) {
          result_ct[idx, "p_value"] <- result["p_value"]
        } else {
          result_ct[idx, "p_value"] <- result[gene, "p_value"]
        }
      } else {
        result_ct[idx, "p_value"] <- -1
      }
    }
    
    results[[ct]] <- result_ct
  }
  
  long_df <- do.call(rbind, results)
  
  # add adjusted p value column exactly as your original
  tested_idx <- which(long_df$p_value >= 0)
  p_adj <- p.adjust(long_df$p_value[tested_idx], method = "BH")
  long_df$p_adj <- -1
  long_df$p_adj[tested_idx] <- p_adj
  
  long_df
}
