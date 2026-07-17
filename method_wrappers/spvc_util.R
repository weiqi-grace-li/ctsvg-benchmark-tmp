library(spVC)
library(sp)
library(BPST)
library(Triangulation)
library(MGLM)

run_spvc_full_square = function(sp_counts, sp_coords, sp_comp, Tr.cell, ncores = 4){
  V = as.matrix(Tr.cell$V)
  Tr = as.matrix(Tr.cell$Tr)
  mySPVC = test.spVC(Y = sp_counts, X = sp_comp, S = sp_coords, V = V, Tr = Tr, 
                     para.cores = ncores, filter.min.nonzero = 100, twostep = FALSE)
  return(mySPVC)
}



spvc_details_save = function(sim_obj, spvc_obj, sim_name, deconv_method, test_method = "spvc", two_step = TRUE, file_path = NA){
  all_gene = sim_obj$gene_names
  cell_types = unique(sim_obj$cell_metadata$type)
  if (two_step){
    spvc_raw = spvc_obj$results.varying
  } else {
    spvc_raw = spvc_obj$results.full 
  }
  
  ## need to transform spvc results to common form 
  all_covariates = unique(unlist(lapply(spvc_raw, function(genelist) {
    if (length(genelist)>1) {
      names(genelist$p.value)
    }
  })))
  # gamma_covariates = sort(setdiff(grep("^gamma", all_covariates, value = TRUE), "gamma_0"))
  gamma_covariates = setdiff(grep("^gamma", all_covariates, value = TRUE), "gamma_0")
  print(gamma_covariates)
  spvc_results = list()
  for (i in seq_along(gamma_covariates)){
    cov = gamma_covariates[i]
    temp_result = sapply(spvc_raw, function(genelist) {
      if (length(genelist) > 1){
        genelist$p.value[cov] 
      } else {
        -1 
      }
    })
    
    temp_result = temp_result[!is.na(temp_result)]
    
    # Split names into parts
    split_names <- do.call(rbind, strsplit(names(temp_result), "\\.")) # names of those p values are genename.gamma_*
    
    # Create data frame with p-values and set row names to gene names
    temp_result_df <- data.frame(
      p_value = as.numeric(temp_result),
      row.names = split_names[, 1]
    )
    spvc_results[[colnames(sim_obj$spot_composition)[i]]] = temp_result_df
  }
  
  results = list()
  ## now for each cell type record the detailed p value results
  for (ct in cell_types){
    ## sometimes cell types are not tested 
    if (ct %in% names(spvc_results)){
      result = spvc_results[[ct]]  
      gene_tested = rownames(result)
    } else{
      result = NULL
      gene_tested = character(0)
    }
    
    gene_tested = row.names(result)

    # create result table and prefil items
    result_ct = as.data.frame(matrix(NA, nrow = length(all_gene), ncol = 12))
    colnames(result_ct) = c("simulation_name", "test_method", "deconv_method", "seed", "cell_type", "cell_proportion",
                            "gene_name", "dispersion","twostep","gene_type", "p_value", "p_adj")
    result_ct$simulation_name = sim_name
    result_ct$test_method = test_method
    result_ct$deconv_method = deconv_method
    result_ct$seed = sim_obj$seed
    result_ct$cell_type = ct 
    result_ct$cell_proportion = sim_obj$cell_type_proportion[ct]
    result_ct$gene_name = all_gene
    result_ct$twostep = ifelse(two_step, "TRUE", "FALSE")
    if (length(sim_obj$dispersion)==1){
      result_ct$dispersion = sim_obj$dispersion
    } else {
      result_ct$dispersion = rowMeans(sim_obj$dispersion[all_gene, ])
    }
    result_ct$gene_type = label_gene_types(all_gene, sim_obj, ct, cell_types)

    # generate the p_value for each gene
    for (gene in all_gene){
      if (gene %in% gene_tested) {
        result_ct[which(result_ct$gene_name == gene), "p_value"] = result[gene, "p_value"] 
      } else {
        result_ct[which(result_ct$gene_name == gene), "p_value"] =  -1
      }
    }
    results[[ct]] = result_ct
  }
  
  long_df = do.call(rbind, results)
  ## put adjusted p value 
  tested = which(long_df$p_value >= 0)
  adjusted_p = p.adjust(long_df$p_value[tested], method = "BH")
  long_df[tested, "p_adj"] = adjusted_p
  long_df[-tested, "p_adj"] = -1
  
  
  save_result_with_file_switching(long_df, file_path, prefix = "spvc", max_rows = 300000)

  return(long_df)
}
