library(Matrix)
library(spacexr)
library(MCube)

# create RCTD object 

run_mmm = function(
    sp_counts, sp_coords, sp_comp, sc_counts, sc_metadata, ncores = 4,
    proportion_threshold = 0.1, celltype_threshold = 100, reference_threshold = 0.5,
    CELL_MIN_INSTANCE = 25, counts_MIN = 10, fc_cutoff_reg = 0.75, fc_cutoff = 0.5, # for create.RCTD
    gene_cutoff_reg = 0.0002, gene_cutoff = 0.000125,                                # for create.RCTD
    RCTD_obj = NULL
){
  
  message(sprintf("Proportion_threshold = %.2f, celltype_threshold = %.2f, reference_threshold = %.2f", proportion_threshold, celltype_threshold, reference_threshold))

  if (!is.null(RCTD_obj)){
    message("Using the given RCTD object...")
    myRCTD = RCTD_obj
    weights_RCTD <- as.matrix(myRCTD@results$weights)
    spot_effects_RCTD <- log(rowSums(weights_RCTD))
    names(spot_effects_RCTD) <- rownames(weights_RCTD)
  } else {
    message("Initializing RCTD object...")
    cell_types = factor(sc_metadata$type, levels = colnames(sp_comp))
    names(cell_types) = rownames(sc_metadata)
    
    reference = Reference(
      sc_counts, 
      cell_types, 
      colSums(sc_counts), 
      min_UMI = 25 # this is a different choice compared to CSIDE 
    )
    
    puck = SpatialRNA(as.data.frame(sp_coords), sp_counts, colSums(sp_counts))
    
    myRCTD = create.RCTD(
      puck, reference, max_cores = 8,
      CELL_MIN_INSTANCE = CELL_MIN_INSTANCE, counts_MIN = counts_MIN, fc_cutoff_reg = fc_cutoff_reg,
      fc_cutoff = fc_cutoff, gene_cutoff_reg = gene_cutoff_reg, gene_cutoff = gene_cutoff
    )
    
    spot_effects_RCTD <- log(rowSums(sp_comp))
  }
  
  mcube_object <- createMCube(
    counts = t(sp_counts),
    coordinates = sp_coords,
    proportions = sp_comp,
    library_sizes = colSums(sp_counts),
    reference = t(myRCTD@cell_type_info$info[[1]]),
    used_for_deconvolution = rownames(myRCTD@spatialRNA@counts),
    spot_effects = spot_effects_RCTD,
    platform_effect = NULL, 
    proportion_threshold = proportion_threshold, 
    celltype_threshold = celltype_threshold, 
    reference_threshold = reference_threshold
  )
  
  mcube_object <- mcubeFitNull(
    mcube_object,
    num_workers = ncores, num_threads = 1
  )
  
  mcube_object <- mcubeTest(
    mcube_object,
    num_workers = ncores, num_threads = 1, shared_memory = TRUE
  )
  
  return(mcube_object)
}  

mmm_details_save = function(
    sim_obj, mcube_obj, sim_name, deconv_method, file_path = NA,
    parameter_setting = "proportion_threshold = 0.1 | celltype_threshold = 100 | reference_threshold = 0.5 | CELL_MIN_INSTANCE = 25 | counts_MIN = 10 | fc_cutoff_reg = 0.75 | fc_cutoff = 0.5"
  ){
  all_gene = sim_obj$gene_names
  mmm_results = mcube_obj@pvalues
  ## now for each cell type record the detailed p value results 
  cell_types = unique(sim_obj$cell_metadata$type)
  results = list()
  for (ct in cell_types){
    result = mmm_results[[ct]]
    gene_tested = row.names(result)
    
    # create result table and prefil items 
    result_ct = as.data.frame(matrix(NA, nrow = length(all_gene), ncol = 11))
    rownames(result_ct) = all_gene
    colnames(result_ct) = c("simulation_name", "test_method", "deconv_method", "parameter_setting", "seed", "cell_type", "cell_proportion",
                            "gene_name", "dispersion","gene_type", "p_value")
    result_ct$simulation_name = sim_name
    result_ct$test_method = "mmm"
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
        result_ct[which(result_ct$gene_name == gene), "p_value"] = result[gene, "combined_pvalue"] 
      } else {
        result_ct[which(result_ct$gene_name == gene), "p_value"] = -1
      }
    }
    
    results[[ct]] = result_ct
  }

  long_df = do.call(rbind, results)

  tested_idx = which(long_df$p_value >= 0)
  long_df$p_adj = -1
  long_df$p_adj[tested_idx] = p.adjust(long_df$p_value[tested_idx], method = "BH")
  
  ## ------------------------
  
  save_result_with_file_switching(long_df, file_path, prefix = "mmm", max_rows = 1000000)
  return(long_df)
}





