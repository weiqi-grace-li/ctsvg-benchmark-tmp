## --- collect scDesign3 copula and para 
scdesign3_fit_raw = function(data_ori, save_path_raw, save_path_dev, seed = 42, ncores = 4, ncores_small = NULL, cell_count_threshold = 100){
  library(scDesign3)
  library(SingleCellExperiment)
  
  sorted_cell_types = sort(table(data_ori$cell_metadata$type), decreasing = TRUE)
  cell_types = names(sorted_cell_types)[sorted_cell_types>cell_count_threshold]
  
  # -- initialize list for copula and para 
  copula_list = vector("list", length = length(cell_types))
  para_list = vector("list", length = length(cell_types))
  data_list = vector("list", length = length(cell_types))
  sp_list = vector("list", length = length(cell_types))
  names(copula_list) = names(para_list) = names(data_list) = names(sp_list) = cell_types
  
  dev_table = matrix(NA, ncol = length(cell_types), nrow = nrow(data_ori$cell_counts_ori))
  colnames(dev_table) = cell_types
  rownames(dev_table) = rownames(data_ori$cell_counts_ori)
  
  for (ct in cell_types){
    cell_names = rownames(data_ori$cell_metadata)[which(data_ori$cell_metadata$type==ct)]
    
    cell_sp = SingleCellExperiment(assays = list(counts = data_ori$cell_counts_ori[,cell_names]))
    sp_list[[ct]] = cell_sp
    
    
    colData(cell_sp) = DataFrame(cell_type = data_ori$cell_metadata[cell_names,"type"], spatial1 = data_ori$cell_metadata[cell_names,"x"],
                                 spatial2 = data_ori$cell_metadata[cell_names,"y"],
                                 row.names = cell_names)
    
    message(paste0("Fitting cell type number ", which(cell_types==ct)," out of ", length(cell_types), " | ", ct))
    
    set.seed(seed)
    cell_data = construct_data(
      sce = cell_sp,
      assay_use = "counts",
      celltype = "cell_type",
      pseudotime = NULL,
      spatial = c("spatial1", "spatial2"),
      other_covariates = NULL,
      corr_by = "1"
    )
    
    data_list[[ct]] = cell_data
    
    # here we allow a gene to have just one intercept, ignoring the cell type issue
    cell_marginal = fit_marginal(
      data = cell_data,
      predictor = "gene",
      mu_formula = "s(spatial1, spatial2, bs ='gp', k = 50)",
      sigma_formula = "1",
      family_use = "nb",
      n_cores = ncores,
      parallelization = "pbmcmapply",
      usebam = FALSE
    )
    
    for (gene in rownames(dev_table)){
      if (!(is.na(cell_marginal[[gene]]$fit)[1])){
        dev = summary(cell_marginal[[gene]]$fit)$dev.expl
      } else{
        dev = 0
      }
      dev_table[gene, ct] = dev
    }
    message("The Top 10 deviation explained by fitting is noted below: ")
    print(sort(dev_table[,ct], decreasing = TRUE)[1:10])
    
    set.seed(seed)
    cell_copula = fit_copula(
      sce = cell_sp,
      assay_use = "counts",
      marginal_list = cell_marginal,
      family_use = "nb",
      copula = "gaussian",
      n_cores = ifelse(is.null(ncores_small), ncores, ncores_small),
      parallelization = "pbmcmapply",
      input_data = cell_data$dat
    )
    
    copula_list[[ct]] = cell_copula
    
    cell_para = extract_para(
      sce = cell_sp,
      marginal_list = cell_marginal,
      family_use = "nb",
      new_covariate = cell_data$newCovariate,
      n_cores = ifelse(is.null(ncores_small), ncores, ncores_small),
      parallelization = "pbmcmapply",
      data = cell_data$dat
    )
    
    para_list[[ct]] = cell_para
  }
  
  scdesign_raw = list(copula = copula_list, para = para_list, data = data_list, sp = sp_list)
  
  save(dev_table, file = save_path_dev)
  save(scdesign_raw, file = save_path_raw)
  return(scdesign_raw)
}

## --- generate null data with special_genes and scdesign_raw 
generate_scdesign_newdata = function(gene_names, special_genes, scdesign_raw_path, save_path, seed = 42, ncores = 4){
  library(scDesign3)
  library(SingleCellExperiment)
  load(scdesign_raw_path)
  cell_types = names(scdesign_raw$sp)
  newcount_list = vector("list", length(cell_types))
  names(newcount_list) = cell_types
  
  for (ct in cell_types){
    message(paste0("Generating data for cell type number ", which(cell_types==ct)," out of ", length(cell_types), " | ", ct))
    se_genes = special_genes[[ct]]$ct_svg
    mu_matrix = scdesign_raw$para[[ct]]$mean_mat
    non_se_genes = setdiff(gene_names, se_genes)
    for (gene in non_se_genes){
      new_mu = mean(mu_matrix[, gene])
      mu_matrix[, gene] = new_mu
    }
    
    set.seed(seed)
    
    scdesign_raw$para[[ct]]$mean_mat = mu_matrix
    cell_newcount <- simu_new(
      sce = scdesign_raw$sp[[ct]],
      mean_mat = scdesign_raw$para[[ct]]$mean_mat,
      sigma_mat = scdesign_raw$para[[ct]]$sigma_mat,
      zero_mat = scdesign_raw$para[[ct]]$zero_mat,
      quantile_mat = NULL,
      copula_list = scdesign_raw$copula[[ct]]$copula_list,
      n_cores = ncores,
      parallelization = "pbmcmapply",
      family_use = "nb",
      input_data = scdesign_raw$data[[ct]]$dat,
      new_covariate = scdesign_raw$data[[ct]]$newCovariate,
      important_feature = rep(TRUE, dim(scdesign_raw$sp[[ct]])[1]),
      filtered_gene = scdesign_raw$data[[ct]]$filtered_gene
    )
    
    newcount_list[[ct]] = cell_newcount
  }
  
  scdesign_new = list(para = scdesign_raw$para, newcount = newcount_list)
  save(scdesign_new, file = save_path)
  return(scdesign_new)
}

generate_simulator = function(data_ori, scdesign_new_path, save_path, seed, special_genes, strat = "A"){
  source("./scDesign_simulators/pseudo_spot_simulator.R")
  load(scdesign_new_path)
  data_new = pseudo_simulator$new(
    data_ori, scdesign_new$newcount,scdesign_new$para, seed = seed, special_genes = special_genes, strat = strat
  )
  save(data_new, file = save_path)
  message("Simulator generated...")
  return(data_new)
}
