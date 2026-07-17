## require 
# 1. count matrix, [gene x spots].  Have to have column names and row names 
# 2. position matrix, [spots x 2]. 
# 3. composition matrix, [spots x n_cell_types]
# 4. gene threshold = 0.05
# 5. spot.threshold = 10
# 6. normalized = FALSE 
# 7. correction = FALSE 
# 8. ncores = 4

run_stance = function(count, pos, comp, gene_thres = 0.05, spot_thres = 10, normalized = F, correction = F, ncores = 4, utsvg_thres = 0.05, topn_utsvg = NULL, dist_preserve = FALSE){
  library(STANCE)
  # create stance object
  stance_obj = creatSTANCEobject(counts = count, pos = pos, prop = comp)
  stance_obj = data_preprocess(object = stance_obj, gene.threshold = gene_thres, spot.threshold = spot_thres, normalized = normalized)

  if (dist_preserve) stance_obj@location = dist_perserve_trans(pos)

  # run STANCE individual test
  stance_obj = build_kernelMatrix(object = stance_obj)
  stance_obj = runTest1(object = stance_obj, correction = correction, pv.adjust = "BY")
  gene.list = rownames(stance_obj@gene_expression)
  if (is.null(topn_utsvg)){ # 20260120 change for real data 
    utSVG.list = gene.list[stance_obj@Test_1$p_value_adj < utsvg_thres] # use adjusted p value  
  } else { 
    topn_utsvg <- min(topn_utsvg, length(gene.list))
    ord = order(stance_obj@Test_1$p_value_adj)
    utSVG.list = gene.list[ord[1:topn_utsvg]]
    best_i <- ord[1]
    message(paste0("Selecting only top ", topn_utsvg,
                   " | best gene ", gene.list[best_i],
                   " | best p_adj ", stance_obj@Test_1$p_value_adj[best_i]))
  }
  
  message("##-----")
  message(paste0("utSVG has selected ", length(utSVG.list), " | ",  length(gene.list)," genes for stage 2 test."))
  message("##-----")
  stance_obj = runTest2(object = stance_obj, Genes_to_test = utSVG.list, Cell_types_to_test = NULL, # test for all cell types 
                        correction = correction, ncores = ncores)
  
  stance_obj
}

# "Modified" STANCE run used in the scdesign3_fitting realistic simulations
# (breast_small.R / ovarian_small.R / lymph_small.R): same pipeline as run_stance,
# but with spot_thres = 0 and distance-preserving location normalization
# (dist_perserve_trans) applied instead of the raw spot coordinates.
run_stance_mod = function(count, pos, comp, gene_thres = 0.05, spot_thres = 0, normalized = F, correction = F, ncores = 4, utsvg_thres = 0.05, topn_utsvg = NULL){
  run_stance(
    count = count, pos = pos, comp = comp, gene_thres = gene_thres, spot_thres = spot_thres,
    normalized = normalized, correction = correction, ncores = ncores,
    utsvg_thres = utsvg_thres, topn_utsvg = topn_utsvg, dist_preserve = TRUE
  )
}


stance_details_save = function(sim_obj, stance_obj, sim_name, deconv_method, file_path = NA, test_method = "stance"){
  all_gene = sim_obj$gene_names
  ## recreate what has been tested in stance test 2
  gene.list = rownames(stance_obj@gene_expression)
  utsvgs = gene.list[stance_obj@Test_1$p_value_adj < 0.05] # use adjusted p value
  
  ## now for each cell type record the detailed p value results 
  cell_types = unique(sim_obj$cell_metadata$type)
  results = list()
  for (ct in cell_types){
    result = stance_obj@Test_2[[ct]]

    # create result table and prefil items
    result_ct = as.data.frame(matrix(NA, nrow = length(all_gene), ncol = 10))
    colnames(result_ct) = c("simulation_name", "test_method", "deconv_method", "seed", "cell_type", "cell_proportion",
                            "gene_name", "dispersion","gene_type", 
                            "p_value")
    result_ct$simulation_name = sim_name
    result_ct$test_method = test_method
    result_ct$deconv_method = deconv_method
    result_ct$seed = sim_obj$seed
    result_ct$cell_type = ct 
    result_ct$cell_proportion = sim_obj$cell_type_proportion[ct]
    result_ct$gene_name = all_gene
    if (length(sim_obj$dispersion)==1){
      result_ct$dispersion = sim_obj$dispersion
    } else {
      result_ct$dispersion = rowMeans(sim_obj$dispersion[all_gene, ]) ### TODO: need to specify more once scDesign3 are formalized
    }
    result_ct$gene_type = label_gene_types(all_gene, sim_obj, ct, cell_types)

    # generate the p_value for each gene
    for (gene in all_gene){
      if (gene %in% utsvgs) {
        if (!is.matrix(result) & !is.data.frame(result)){
          result_ct[which(result_ct$gene_name == gene), "p_value"] = result["p_value"]
        } else{
           result_ct[which(result_ct$gene_name == gene), "p_value"] = result[gene, "p_value"]  
        }
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
  save_result_with_file_switching(long_df, file_path, prefix = "stance", max_rows = 1000000)

  return(long_df)
}