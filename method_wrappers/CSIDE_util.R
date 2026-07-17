suppressPackageStartupMessages(library(spacexr))

run_cside_import = function(
    sp_counts, sp_coords, sp_comp, sc_counts, sc_metadata, ncores = 4, # data
    CELL_MIN_INSTANCE = 25, counts_MIN = 10, fc_cutoff_reg = 0.75, fc_cutoff = 0.5, # for create.RCTD
    gene_cutoff_reg = 0.0002, gene_cutoff = 0.000125,                                # for create.RCTD
    df = 15, cell_type_threshold = 125, weight_threshold = 0.8                       # for run.CSIDE.nonparam
  ){
  # create spatial data
  puck = SpatialRNA(as.data.frame(sp_coords), sp_counts)
  sc_types = factor(sc_metadata$type)
  names(sc_types) = rownames(sc_metadata)

  # create reference data, it is NOT used
  reference = Reference(sc_counts, sc_types, colSums(sc_counts), min_UMI = 25) # this is different than default

  myRCTD = create.RCTD(
    spatialRNA = puck, reference, max_cores = ncores,
    CELL_MIN_INSTANCE = CELL_MIN_INSTANCE, counts_MIN = counts_MIN, fc_cutoff_reg = fc_cutoff_reg, fc_cutoff = fc_cutoff,
    gene_cutoff_reg = gene_cutoff_reg, gene_cutoff = gene_cutoff
  )
  
  # directly import user given weight
  myRCTD = import_weights(myRCTD, sp_comp)
  myRCTD@config$RCTDmode <- 'full' ## should be in the most updated package 

  myRCTD = run.CSIDE.nonparam(
    myRCTD, df = df, cell_types = NULL, gene_threshold = 5e-5, 
    cell_type_threshold = cell_type_threshold, fdr = 0.05, doublet_mode = FALSE, 
    weight_threshold = weight_threshold
  )
    return(myRCTD) 
}


cside_details_save = function(
    sim_obj, cside_obj, sim_name, deconv_method, 
    parameter_setting = "CELL_MIN_INSTANCE = 25 | counts_MIN = 10 | fc_cutoff_reg = 0.75 | fc_cutoff = 0.5 | df = 15 | cell_type_threshold = 125 | weight_threshold = 0.8", # default parameters 
    file_path = NA){
  
  all_gene = sim_obj$gene_names
  cside_results = cside_obj@de_results$all_gene_list
  cside_convmat = cside_obj@de_results$gene_fits$con_mat
  ## now for each cell type record the detailed p value results 
  cell_types = unique(sim_obj$cell_metadata$type)
  results = list()
  for (ct in cell_types){
    ## sometimes cell types are not tested 
    if (ct %in% names(cside_results)){
      result = cside_results[[ct]]  
      gene_tested = rownames(result)
      # create BH adjusted p value 
      p_adj = p.adjust(result$p_val, method = "BH")
      result = cbind(result, p_adj)
      
      # make sure we have found the same significant genes as cside 
      sig_genes_user = sort(rownames(result)[p_adj<0.05])
      sig_genes_cside = sort(rownames(cside_obj@de_results$sig_gene_list[[ct]]))
      if (any(sig_genes_user != sig_genes_cside)){
        stop("BH adjusted pvalue result doesn't match CSIDE")
      }
    } else{
      result = NULL
      gene_tested = character(0)
    }
    
    gene_tested = row.names(result)
    
    # create result table and prefil items 
    result_ct = as.data.frame(matrix(NA, nrow = length(all_gene), ncol = 13))
    colnames(result_ct) = c("simulation_name", "test_method", "deconv_method", "parameter_setting", "seed", "cell_type", "cell_proportion",
                            "gene_name", "dispersion","gene_type", "p_value", "p_adj", "convergence")
    result_ct$simulation_name = sim_name
    result_ct$test_method = "cside"
    result_ct$deconv_method = deconv_method
    result_ct$parameter_setting = parameter_setting
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
        result_ct[which(result_ct$gene_name == gene), "p_value"] = result[gene, "p_val"] 
        result_ct[which(result_ct$gene_name == gene), "p_adj"] = result[gene, "p_adj"]
      } else {
        result_ct[which(result_ct$gene_name == gene), "p_value"] = -1
        result_ct[which(result_ct$gene_name == gene), "p_adj"] = -1
      }
      if ((gene %in% rownames(cside_convmat)) & (ct %in% colnames(cside_convmat))){
        result_ct[which(result_ct$gene_name == gene), "convergence"] = cside_convmat[gene, ct]
      } else {
        result_ct[which(result_ct$gene_name == gene), "convergence"] = "unknown"
      }
      
    }
    results[[ct]] = result_ct
  }
  
  long_df = do.call(rbind, results)
  save_result_with_file_switching(long_df, file_path, prefix = "cside", max_rows = 300000)

  return(long_df)
}
