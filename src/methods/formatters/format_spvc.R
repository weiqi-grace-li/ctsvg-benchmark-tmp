# src/methods/formatters/spvc_format.R

spvc_details_to_df <- function(
    sim_obj,
    spvc_obj,
    sim_name,
    deconv_method,
    test_method = "spvc",
    two_step = TRUE
) {
  all_gene <- sim_obj$gene_names
  cell_types <- unique(sim_obj$cell_metadata$type)
  
  if (two_step) {
    spvc_raw <- spvc_obj$results.varying
  } else {
    spvc_raw <- spvc_obj$results.full
  }
  
  # get all covariates that appear
  all_covariates <- unique(unlist(lapply(spvc_raw, function(genelist) {
    if (length(genelist) > 1) {
      names(genelist$p.value)
    }
  })))
  
  gamma_covariates <- setdiff(grep("^gamma", all_covariates, value = TRUE), "gamma_0")
  print(gamma_covariates)
  
  spvc_results <- list()
  spvc_gam_fail <- list() # kept (unused) to match your original structure
  
  for (i in seq_along(gamma_covariates)) {
    cov <- gamma_covariates[i]
    
    temp_result <- sapply(spvc_raw, function(genelist) {
      if (length(genelist) > 1) {
        genelist$p.value[cov]
      } else {
        -1
      }
    })
    
    temp_fail <- sapply(spvc_raw, function(genelist) {
      if (length(genelist) > 1) {
        genelist$p.value.fail[cov]
      } else {
        2
      }
    })
    
    temp_result <- temp_result[!is.na(temp_result)]
    temp_fail <- temp_fail[!is.na(temp_fail)]
    
    split_names <- do.call(rbind, strsplit(names(temp_result), "\\."))
    
    temp_result_df <- data.frame(
      p_value = as.numeric(temp_result),
      fail = as.numeric(temp_fail),
      row.names = split_names[, 1]
    )
    
    spvc_results[[colnames(sim_obj$spot_composition)[i]]] <- temp_result_df
  }
  
  results <- list()
  
  for (ct in cell_types) {
    # sometimes cell types are not tested
    if (ct %in% names(spvc_results)) {
      result <- spvc_results[[ct]]
      gene_tested <- rownames(result)
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
    
    # create result table and prefill items (exactly 13 columns as your original)
    result_ct <- as.data.frame(matrix(NA, nrow = length(all_gene), ncol = 13))
    colnames(result_ct) <- c(
      "simulation_name", "test_method", "deconv_method", "seed", "cell_type", "cell_proportion",
      "gene_name", "dispersion", "twostep", "gene_type", "p_value", "p_adj", "p_fail"
    )
    
    result_ct$simulation_name <- sim_name
    result_ct$test_method <- test_method
    result_ct$deconv_method <- deconv_method
    result_ct$seed <- sim_obj$seed
    result_ct$cell_type <- ct
    result_ct$cell_proportion <- sim_obj$cell_type_proportion[ct]
    result_ct$gene_name <- all_gene
    result_ct$twostep <- ifelse(two_step, "TRUE", "FALSE")
    
    if (length(sim_obj$dispersion) == 1) {
      result_ct$dispersion <- sim_obj$dispersion
    } else {
      result_ct$dispersion <- rowMeans(sim_obj$dispersion[all_gene, ])
    }
    
    for (gene in all_gene) {
      gene_labels <- names(check_list)[sapply(check_list, function(lst) gene %in% lst)]
      gene_label <- if (length(gene_labels) == 0) "Null" else paste(gene_labels, collapse = ",")
      
      idx <- which(result_ct$gene_name == gene)
      result_ct[idx, "gene_type"] <- gene_label
      
      if (gene %in% gene_tested) {
        result_ct[idx, "p_value"] <- result[gene, "p_value"]
        result_ct[idx, "p_fail"] <- result[gene, "fail"]
      } else {
        result_ct[idx, "p_value"] <- -1
        result_ct[idx, "p_fail"] <- -1
      }
    }
    
    results[[ct]] <- result_ct
  }
  
  long_df <- do.call(rbind, results)
  
  tested <- which(long_df$p_value >= 0)
  adjusted_p <- p.adjust(long_df$p_value[tested], method = "BH")
  long_df[tested, "p_adj"] <- adjusted_p
  long_df[-tested, "p_adj"] <- -1
  
  long_df
}
