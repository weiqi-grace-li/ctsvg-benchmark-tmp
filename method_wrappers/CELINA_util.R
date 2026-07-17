run_celina_spot = function(count, pos, comp, sc_count, sc_type, ncores = 4){
  library(CELINA)

  if (is.null(names(sc_type))){
    names(sc_type) = colnames(sc_count)
  }
  
  myCELINA = Create_Celina_Object(celltype_mat = t(comp),
                                  gene_expression_mat = count,
                                  location = pos, 
                                  covariates = NULL,
                                  project = "Sim"
  )

  filtered_cell_types = unique(colnames(comp))
  
  myCELINA = preprocess_input(myCELINA,
                              cell_types_to_test = filtered_cell_types,
                              scRNA_count = sc_count,
                              sc_cell_type_labels = as.factor(sc_type),
                              threshold = 5e-5
  )

  if (dim(count)[2] <= 8000){
    myCELINA = Calculate_Kernel(myCELINA, approximation = FALSE)  
  } else {
    myCELINA = Calculate_Kernel(myCELINA, approximation = TRUE, sparseKernel = TRUE, sparseKernel_tol = 1e-4)  
  }
  
  myCELINA = Testing_interaction_all(
    myCELINA, celltype_to_test = filtered_cell_types,
    num_cores = ncores
  )
  return(myCELINA)
}

celina_details_save = function(sim_obj, celina_obj, sim_name, deconv_method, file_path = NA){
  all_gene = sim_obj$gene_names
  celina_results = celina_obj@result
  ## now for each cell type record the detailed p value results 
  cell_types = unique(sim_obj$cell_metadata$type)
  results = list()
  for (ct in cell_types){
    result = celina_results[[ct]]
    gene_tested = row.names(result)
    
    # create result table and prefil items 
    result_ct = as.data.frame(matrix(NA, nrow = length(all_gene), ncol = 10))
    colnames(result_ct) = c("simulation_name", "test_method", "deconv_method", "seed", "cell_type", "cell_proportion",
                            "gene_name", "dispersion","gene_type", 
                            "p_value")
    result_ct$simulation_name = sim_name
    result_ct$test_method = "celina"
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
        result_ct[which(result_ct$gene_name == gene), "p_value"] = result[gene, "CombinedPvals"]
      } else {
        result_ct[which(result_ct$gene_name == gene), "p_value"] = -1
      }
    }
    results[[ct]] = result_ct
  }
  
  long_df = do.call(rbind, results)
  
  ## --- add adjusted p value 
  tested_idx = which(long_df$p_value >= 0)
  p_adj = p.adjust(long_df$p_value[tested_idx], method = "BH") 
  long_df$p_adj = -1 
  long_df$p_adj[tested_idx] = p_adj
  
  ## ------------------------
  
  save_result_with_file_switching(long_df, file_path, prefix = "celina", max_rows = 1000000)
  return(long_df)
}