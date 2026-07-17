library(CTSV)
library(SpatialExperiment)
library(BiocParallel)

run_ctsv = function(sp_counts, sp_coords, sp_comp, ncores = 4){
  spe = SpatialExperiment(
    assays = list(counts = sp_counts), 
    spatialCoords = sp_coords
  )
  BPPARAM = MulticoreParam(workers = ncores, progressbar = TRUE)
  myCTSV = CTSV(spe, sp_comp, num_core = ncores, BPPARAM)
  return(myCTSV)
}

ctsv_details_save = function(sim_obj, ctsv_obj, sim_name, deconv_method, file_path = NA){
  all_gene = sim_obj$gene_names
  
  ## need to assume the order of CTSV results match the column of spot_composition
  ## first let's fix the ctsv results, a list of cell types, with column p value, q value, pvalue s1, p value s2, q value s1, q value s2
  K = ncol(sim_obj$spot_composition)
  ctsv_result = list()
  for (k in 1:K){
    cell_result = data.frame(matrix(NA, ncol = 6, nrow = nrow(ctsv_obj$qval)))
    rownames(cell_result) = rownames(ctsv_obj$qval)
    colnames(cell_result) = c("p_value", "q_value", "p_value_s1", "p_value_s2", "q_value_s1", "q_value_s2")
    cell_result$p_value = apply(ctsv_obj$pval[,c(k, k+K)], 1, min)
    cell_result$q_value = apply(ctsv_obj$qval[,c(k, k+K)], 1, min)
    cell_result$p_value_s1 = ctsv_obj$pval[,k]
    cell_result$p_value_s2 = ctsv_obj$pval[,k+K]
    cell_result$q_value_s1 = ctsv_obj$qval[,k]
    cell_result$q_value_s2 = ctsv_obj$qval[,k+K]
    ctsv_result[[colnames(sim_obj$spot_composition)[k]]] = cell_result
  }
  
  ## now for each cell type record the detailed p value results 
  cell_types = unique(sim_obj$cell_metadata$type)
  results = list()
  ctsv_sig_genes = svGene(ctsv_obj$qval)
  ctsv_sig_genes = ctsv_sig_genes$SVGene
  names(ctsv_sig_genes) = colnames(sim_obj$spot_composition)
  for (ct in cell_types){
    ## sometimes cell types are not tested 
    if (ct %in% names(ctsv_result)){
      result = ctsv_result[[ct]]  
      gene_tested = rownames(result)
      
      # double check the sig result fits with the model 
      sig_genes_user = sort(rownames(result)[which(result$q_value<0.05)])
      sig_genes_method = sort(ctsv_sig_genes[[ct]])
      if (any(sig_genes_user != sig_genes_method)){
        stop("BH adjusted pvalue result doesn't match CSIDE")
      }
      
    } else{
      result = NULL
      gene_tested = character(0)
    }
    
    gene_tested = row.names(result)
    
    # create result table and prefil items 
    result_ct = as.data.frame(matrix(NA, nrow = length(all_gene), ncol = 15))
    colnames(result_ct) = c("simulation_name", "test_method", "deconv_method", "seed", "cell_type", "cell_proportion",
                            "gene_name", "dispersion","gene_type", "p_value", "p_adj", "p_value_s1", "p_value_s2", "q_value_s1", "q_value_s2")
    result_ct$simulation_name = sim_name
    result_ct$test_method = "ctsv"
    result_ct$deconv_method = deconv_method
    result_ct$seed = sim_obj$seed
    result_ct$cell_type = ct 
    result_ct$cell_proportion = sim_obj$cell_type_proportion[ct]
    result_ct$gene_name = all_gene
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
        result_ct[which(result_ct$gene_name == gene), "p_adj"] = result[gene, "q_value"]
        result_ct[which(result_ct$gene_name == gene), "p_value_s1"] = result[gene, "p_value_s1"]
        result_ct[which(result_ct$gene_name == gene), "p_value_s2"] = result[gene, "p_value_s2"]
        result_ct[which(result_ct$gene_name == gene), "q_value_s1"] = result[gene, "q_value_s1"]
        result_ct[which(result_ct$gene_name == gene), "q_value_s2"] = result[gene, "q_value_s2"]
      } else {
        result_ct[which(result_ct$gene_name == gene), "p_value"] =  -1
        result_ct[which(result_ct$gene_name == gene), "p_adj"] = -1
        result_ct[which(result_ct$gene_name == gene), "p_value_s1"] = -1
        result_ct[which(result_ct$gene_name == gene), "p_value_s2"] = -1
        result_ct[which(result_ct$gene_name == gene), "q_value_s1"] = -1
        result_ct[which(result_ct$gene_name == gene), "q_value_s2"] = -1
      }
    }
    results[[ct]] = result_ct
  }
  
  long_df = do.call(rbind, results)
  save_result_with_file_switching(long_df, file_path, prefix = "ctsv", max_rows = 300000)

  return(long_df)
}
