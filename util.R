## --- 1) a method that runs all methods (spvc let's use original and gam version, both remove 0)
run_all_tests = function(
    sp_count, sp_comp, sp_coords, sc_count, sc_metadata, spvc_tri, 
    sim_name, save_dir, 
    sim_obj, 
    time_save_name = "runtime_all", failure_save_name = "runtime_fail", 
    ncores = 4, ncores_cside = NULL,
    run_method = c("celina", "stance", "spvc", "ctsv", "cside", "mmm"), # 20260303 add mmm 
    timeout_sec = Inf, 
    
    # stance parameters 
    utsvg_thres = 0.05, # utsvg_thres added 20260119 (for STANCE)
    topn_utsvg = NULL, # utsvg_thres added 20260120 (for STANCE)
    stance.spot_thres = 10, 
    stance.correction = F, 
    
    # CSIDE parameters
    cside.CELL_MIN_INSTANCE = 25, cside.counts_MIN = 10, cside.fc_cutoff_reg = 0.75, cside.fc_cutoff = 0.5, # for create.RCTD
    cside.gene_cutoff_reg = 0.0002, cside.gene_cutoff = 0.000125,                                           # for create.RCTD
    cside.df = 15, cside.cell_type_threshold = c(125, 100, 75, 50, 25, 0), cside.weight_threshold = 0.8,   # for run.CSIDE.nonparam, CSIDE cell_type_threshold become a vector for auto-switch
    
    # MMM parameters
    mmm.proportion_threshold = 0.1, mmm.celltype_threshold = c(100, 75, 50, 25, 0), mmm.reference_threshold = 0.5,
    mmm.CELL_MIN_INSTANCE = 25, mmm.counts_MIN = 10, mmm.fc_cutoff_reg = 0.75, mmm.fc_cutoff = 0.5, # for create.RCTD
    mmm.gene_cutoff_reg = 0.0002, mmm.gene_cutoff = 0.000125,                                        # for create.RCTD
    RCTD_obj = NULL 
 ){
  source("./method_wrappers/CELINA_util.R")
  source("./method_wrappers/stance_util.R")
  source("./method_wrappers/CSIDE_util.R")
  source("./method_wrappers/CTSV_util.R")
  source("./method_wrappers/spvc_util.R")
  source("./method_wrappers/mmm_util.R")
  
  if ("celina" %in% run_method){
    message("Running celina ...")
    
    # - run celina
    myCELINA = record_run(
      paste0(sim_name, "-", sim_obj$seed), "celina", ncores, paste0(save_dir, time_save_name, "_celina.csv"), function(){
        run_celina_spot(
          count = sp_count,
          pos = sp_coords,
          comp = sp_comp,
          sc_count = sc_count,
          sc_type = sc_metadata$type,
          # approximation = FALSE,  # revision 20260421 respect CELINA's 8000 rules 
          ncores = ncores
        )
      }, 
      timeout_sec = timeout_sec, 
      failure_outfile = paste0(save_dir, failure_save_name, "_celina.csv")
    )
    
    if (!is.null(myCELINA)){
      save_file = retrieve_file_start_with(save_dir, "celina")
      save_file = ifelse(is.na(save_file), "celina.xlsx", save_file)
      celina_details_save(sim_obj = sim_obj, celina_obj = myCELINA, sim_name = sim_name, deconv_method = "oracle", file_path = paste0(save_dir, save_file))  
    }
  }
  
  if ("cside" %in% run_method){
    
    for (thres in cside.cell_type_threshold){
      message(paste0("Running cside with cell type threshold ", thres, "..."))
      cside_parameter_setting = sprintf(
        "CELL_MIN_INSTANCE = %d | counts_MIN = %d | fc_cutoff_reg = %g | fc_cutoff = %g | df = %d | cell_type_threshold = %d | weight_threshold = %g",
        cside.CELL_MIN_INSTANCE, cside.counts_MIN, cside.fc_cutoff_reg, cside.fc_cutoff, cside.df, thres, cside.weight_threshold
      )
      # - run cside 
      myCSIDE = record_run(
        paste0(sim_name, "-", sim_obj$seed), "cside", ifelse(is.null(ncores_cside), ncores, ncores_cside), paste0(save_dir, time_save_name, "_cside.csv"), function(){
          run_cside_import(
            sp_counts = round(sp_count, 0), sp_coords = sp_coords, sp_comp = sp_comp,
            sc_counts = round(sc_count, 0), sc_metadata = sc_metadata,
            ncores = ifelse(is.null(ncores_cside), ncores, ncores_cside),
            # parameters
            CELL_MIN_INSTANCE = cside.CELL_MIN_INSTANCE, counts_MIN = cside.counts_MIN, fc_cutoff_reg = cside.fc_cutoff_reg, fc_cutoff = cside.fc_cutoff,
            gene_cutoff_reg = cside.gene_cutoff_reg, gene_cutoff = cside.gene_cutoff,
            df = cside.df, cell_type_threshold = thres, weight_threshold = cside.weight_threshold
          )
        }, 
        timeout_sec = timeout_sec, 
        failure_outfile = paste0(save_dir, failure_save_name, "_cside.csv"), 
        parameter_setting = cside_parameter_setting
      )
      
      if (!is.null(myCSIDE)){
        save_file = retrieve_file_start_with(save_dir, "cside", format = "csv")
        save_file = ifelse(is.na(save_file), "cside.csv", save_file)
        cside_details_save(
          sim_obj, myCSIDE, sim_name = sim_name, deconv_method = "oracle", 
          parameter_setting = cside_parameter_setting, 
          file_path = paste0(save_dir, save_file)
        )
        break 
      } 
    }
  }
  
  if ("ctsv" %in% run_method){
    message("Running ctsv ...")
    
    # - ctsv 
    my_ctsv = record_run(
      paste0(sim_name, "-", sim_obj$seed), "ctsv", ncores, paste0(save_dir, time_save_name, "_ctsv.csv"), function(){
        run_ctsv(
          sp_counts = sp_count,
          sp_coords = sp_coords,
          sp_comp = sp_comp,
          ncores = ncores
        )
      }, 
      timeout_sec = timeout_sec, 
      failure_outfile = paste0(save_dir, failure_save_name, "_ctsv.csv") 
    )
    
    if (!is.null(my_ctsv)){
      save_file = retrieve_file_start_with(save_dir, "ctsv", format = "csv")
      save_file = ifelse(is.na(save_file), "ctsv.csv", save_file)
      ctsv_details_save(
        sim_obj, my_ctsv, sim_name = sim_name, deconv_method = "oracle", 
        file_path = paste0(save_dir, save_file)
      )  
    }
  }
  
  if ("stance" %in% run_method){
    message("Running stance ...")
    
    # - stance 
    mySTANCE = record_run(
      paste0(sim_name, "-", sim_obj$seed), "stance", ncores, paste0(save_dir, time_save_name, "_stance.csv"), function(){
        run_stance(
          count = sp_count, pos = sp_coords, comp = sp_comp, ncores = ncores, 
          spot_thres = stance.spot_thres, utsvg_thres = utsvg_thres, topn_utsvg = topn_utsvg, 
          correction = stance.correction
        )   
      }, 
      timeout_sec = timeout_sec, 
      failure_outfile = paste0(save_dir, failure_save_name, "_stance.csv")
    )
    
    if (!is.null(mySTANCE)){
      save_file = retrieve_file_start_with(save_dir, "stance")
      save_file = ifelse(is.na(save_file), "stance.xlsx", save_file)
      stance_details_save(
        sim_obj, mySTANCE, sim_name, "oracle", 
        file_path = paste0(save_dir, save_file)
        , test_method = "stance"
      )   
    }
  }
  
  
  if ("stance-mod" %in% run_method){
    message("Running stance-mod ...")

    mySTANCE = record_run(
      paste0(sim_name, "-", sim_obj$seed), "stance-mod", ncores, paste0(save_dir, time_save_name, "_stance.csv"), function(){
        run_stance_mod(
          count = sp_count, pos = sp_coords, comp = sp_comp, ncores = ncores,
          spot_thres = stance.spot_thres, utsvg_thres = utsvg_thres, topn_utsvg = topn_utsvg,
          correction = stance.correction
        )
      },
      timeout_sec = timeout_sec,
      failure_outfile = paste0(save_dir, failure_save_name, "_stance.csv")
    )

    if (!is.null(mySTANCE)){
      save_file = retrieve_file_start_with(save_dir, "stance")
      save_file = ifelse(is.na(save_file), "stance.xlsx", save_file)
      stance_details_save(
        sim_obj, mySTANCE, sim_name, "oracle",
        file_path = paste0(save_dir, save_file),
        test_method = "stance-mod"
      )
    }
  }

  if ("spvc" %in% run_method){
    message("Running spvc-noalt ...")
    
    # - spvc 
    myspVC = record_run(
      paste0(sim_name, "-", sim_obj$seed), "spvc-noalt", ncores, paste0(save_dir, time_save_name, "_spvc.csv"), function(){
        run_spvc_full_square(
          sp_counts = sp_count,
          sp_coords = sp_coords,
          sp_comp = sp_comp,
          Tr.cell = spvc_tri,
          ncores = ncores
        )
      },
      timeout_sec = timeout_sec,
      failure_outfile = paste0(save_dir, failure_save_name, "_spvc.csv")
    )

    if (!is.null(myspVC)){
      save_file = retrieve_file_start_with(save_dir, "spvc", format = "csv")
      save_file = ifelse(is.na(save_file), "spvc.csv", save_file)
      spvc_details_save(
        sim_obj, myspVC, sim_name, deconv_method = "oracle", test_method = "spvc-noalt", two_step = FALSE,
        file_path = paste0(save_dir, save_file)
      )
    }
  }
  
  if ("mmm" %in% run_method){

    for (thres in mmm.celltype_threshold){
      message(paste0("Running mmm with cell type threshold ", thres, "..."))
      mmm_parameter_setting = sprintf(
        "proportion_threshold = %g | celltype_threshold = %d | reference_threshold = %g | CELL_MIN_INSTANCE = %d | counts_MIN = %d | fc_cutoff_reg = %g | fc_cutoff = %g",
        mmm.proportion_threshold, thres, mmm.reference_threshold, mmm.CELL_MIN_INSTANCE, mmm.counts_MIN, mmm.fc_cutoff_reg, mmm.fc_cutoff
      )

      # - run mmm
      mcube_obj = record_run(
        paste0(sim_name, "-", sim_obj$seed), "mmm", ncores, paste0(save_dir, time_save_name, "_mmm.csv"), function(){
          run_mmm(
            sp_counts = round(sp_count, 0), sp_coords = sp_coords, sp_comp = sp_comp,
            sc_counts = round(sc_count, 0), sc_metadata = sc_metadata, ncores = ncores,

            # parameters
            proportion_threshold = mmm.proportion_threshold, celltype_threshold = thres, reference_threshold = mmm.reference_threshold,
            CELL_MIN_INSTANCE = mmm.CELL_MIN_INSTANCE, counts_MIN = mmm.counts_MIN, fc_cutoff_reg = mmm.fc_cutoff_reg, fc_cutoff = mmm.fc_cutoff,
            gene_cutoff_reg = mmm.gene_cutoff_reg, gene_cutoff = mmm.gene_cutoff,
            RCTD_obj = RCTD_obj
          )
        },
        timeout_sec = timeout_sec,
        failure_outfile = paste0(save_dir, failure_save_name, "_mmm.csv"),
        parameter_setting = mmm_parameter_setting
      )

      if (!is.null(mcube_obj)){
        save_file = retrieve_file_start_with(save_dir, "mmm")
        save_file = ifelse(is.na(save_file), "mmm.xlsx", save_file)
        mmm_details_save(
          sim_obj = sim_obj,
          mcube_obj = mcube_obj,
          sim_name = sim_name,
          deconv_method = "oracle",
          file_path = paste0(save_dir, save_file),
          parameter_setting = mmm_parameter_setting
        )
        break
      }
    }
  }
}



record_run <- function(
    data_name, method, ncores, outfile, fun, timeout_sec = Inf, 
    failure_outfile = NULL, parameter_setting = ""
) {
  # reset GC counters so "max used (Mb)" reflects this run
  gc(reset = TRUE)
  
  # run and time
  err <- NULL; err_type = NULL; result <- NULL
  
  tm <- system.time({
    result <- tryCatch({
      if (is.finite(timeout_sec)) {
        R.utils::withTimeout(
          fun(), 
          onTimeout = "error",
          timeout = timeout_sec,
          cpu = Inf, 
          elapsed = timeout_sec # limit on wall time not cpu time 
        )
      } else {
        fun()
      }
    }, error = function(e) {
      # Distinguish timeout vs other errors if possible
      if (inherits(e, "TimeoutException")) {
        err_type <<- "timeout"
        err <<- sprintf("Timed out after %s seconds", timeout_sec)
      } else {
        err_type <<- "error"
        err <<- conditionMessage(e)
      }
      NULL
    })
  })
  
  # peak memory (sum of Ncells/Vcells peaks, in MB)
  g <- gc()
  bpn <- if (.Machine$sizeof.pointer == 8) 56L else 28L  # bytes per Ncell
  
  used_mb <- (g["Vcells","used"]*8 + g["Ncells","used"]*bpn) / 1024^2
  peak_mb <- (g["Vcells","max used"]*8 + g["Ncells","max used"]*bpn) / 1024^2
  
  # write one row immediately (append if file exists)
  dir.create(dirname(outfile), showWarnings = FALSE, recursive = TRUE)
  row <- data.frame(
    dataset   = data_name,
    method    = method,
    elapsed   = unname(tm["elapsed"]),
    user      = unname(tm["user.self"]),
    sys       = unname(tm["sys.self"]),
    peak_mb   = peak_mb,
    used_mb   = used_mb,
    ncores    = ncores,
    status    = if (is.null(err)) "ok" else err_type,
    error     = if (is.null(err)) "" else err,
    timestamp = as.character(Sys.time()),
    parameter_setting = parameter_setting, 
    stringsAsFactors = FALSE
  )
  write.table(
    row, file = outfile, sep = ",",
    row.names = FALSE,
    col.names = !file.exists(outfile),
    append = file.exists(outfile),
    quote = TRUE
  )
  
  if (!is.null(err) && identical(err_type, "error") && !is.null(failure_outfile)) {
    dir.create(dirname(failure_outfile), showWarnings = FALSE, recursive = TRUE)
    write.table(
      row, file = failure_outfile, sep = ",",
      row.names = FALSE,
      col.names = !file.exists(failure_outfile),
      append = file.exists(failure_outfile),
      quote = TRUE
    )
  }
  
  result
}


# Reads a realistic-simulator preprocessed data directory (cell_meta.parquet,
# cell_counts_10x/, gene_dropout_frac.*, cell_spot_prop.*, pseudo_meta.parquet),
# aligns cells/genes/spots, quality-filters empty cells/genes, and builds the
# "_ori" object used to fit scDesign3 and generate simulated data.
build_realistic_ori = function(preprocessed_path) {
  cell_meta = as.data.frame(read_parquet(paste0(preprocessed_path, "cell_meta.parquet")))
  rownames(cell_meta) = cell_meta[, "__index_level_0__"]
  cell_meta = cell_meta[, c("x", "y", "type")]

  cell_count = Read10X(paste0(preprocessed_path, "cell_counts_10x/"))
  stopifnot(length(intersect(rownames(cell_meta), colnames(cell_count))) == length(rownames(cell_meta)))
  stopifnot(length(intersect(rownames(cell_meta), colnames(cell_count))) == length(colnames(cell_count)))

  cell_names = intersect(rownames(cell_meta), colnames(cell_count))
  cell_meta = cell_meta[cell_names, ]
  cell_count = cell_count[, cell_names]

  # quality filter: SELECT non-zero rows/cols (never `-which(...)`, which drops everything
  # when rowSums(...) == 0 has no hits, since `x[-integer(0)]` selects nothing, not everything)
  cell_count = cell_count[rowSums(cell_count) > 0, , drop = FALSE]
  cell_count = cell_count[, colSums(cell_count) > 0, drop = FALSE]
  cell_names = colnames(cell_count)
  cell_meta = cell_meta[cell_names, ]

  cell_gene_drop = readMM(paste0(preprocessed_path, "gene_dropout_frac.mtx.gz"))
  rownames(cell_gene_drop) = read_lines(paste0(preprocessed_path, "gene_dropout_frac.rows.tsv.gz"))
  colnames(cell_gene_drop) = read_lines(paste0(preprocessed_path, "gene_dropout_frac.cols.tsv.gz"))
  gene_cell_drop = t(as(cell_gene_drop, "CsparseMatrix"))

  cell_spot_prop = readMM(paste0(preprocessed_path, "cell_spot_prop.mtx.gz"))
  rownames(cell_spot_prop) = read_lines(paste0(preprocessed_path, "cell_spot_prop.rows.tsv.gz"))
  colnames(cell_spot_prop) = read_lines(paste0(preprocessed_path, "cell_spot_prop.cols.tsv.gz"))
  cell_spot_prop = as(cell_spot_prop, "CsparseMatrix")

  pseudo_meta = as.data.frame(read_parquet(paste0(preprocessed_path, "pseudo_meta.parquet")))
  rownames(pseudo_meta) = pseudo_meta[, "__index_level_0__"]
  pseudo_meta = pseudo_meta[, c("x_transformed", "y_transformed")]
  colnames(pseudo_meta) = c("x", "y")

  # align cell order
  cell_spot_prop = cell_spot_prop[cell_names, ]
  gene_cell_drop = gene_cell_drop[, cell_names]
  rownames(cell_meta) = colnames(cell_count) = rownames(cell_spot_prop) = paste0("Cell", cell_names)
  colnames(gene_cell_drop) = paste0("Cell", cell_names)

  # align gene order
  genes = intersect(rownames(cell_count), rownames(gene_cell_drop))
  cell_count = cell_count[genes, ]
  gene_cell_drop = gene_cell_drop[genes, ]

  # align spot order
  spot = intersect(colnames(cell_spot_prop), rownames(pseudo_meta))
  cell_spot_prop = cell_spot_prop[, spot]
  pseudo_meta = pseudo_meta[spot, ]

  pseudo_data_ori$new(cell_count, cell_meta, pseudo_meta, cell_spot_prop, gene_cell_drop)
}


# Selects cell-type-specific spatially variable genes (ctSVGs) as ground truth for the
# realistic simulators: filters genes by per-cell-type prevalence, filters by Moran's I
# spatial autocorrelation (via Seurat), then ranks the remaining (gene, cell type) pairs
# by scDesign3 deviance-explained and keeps the global top `ctsvg_top_n`.
select_ctsvg_by_deviance = function(
  data_ori, dev_table,
  ct_prop_threshold = 0.01, non_zero_prop_min = 0.01, non_zero_count_min = 100,
  moran_top_n = 200, nfeatures = 200, ctsvg_top_n = 50,
  save_dir = NULL, plot_dir = NULL, width = 4, height = 3.5, dpi = 300
) {
  library(Seurat)

  non_zero_prop = matrix(NA, ncol = ncol(dev_table), nrow = nrow(dev_table),
                          dimnames = dimnames(dev_table))
  non_zero_count = non_zero_prop
  moran = non_zero_prop

  for (ct in colnames(non_zero_prop)) {
    idx = which(data_ori$cell_metadata$type == ct)
    for (gene in rownames(non_zero_prop)) {
      non_zero_prop[gene, ct]  = length(which(data_ori$cell_counts_ori[gene, idx] > 0)) / length(idx)
      non_zero_count[gene, ct] = length(which(data_ori$cell_counts_ori[gene, idx] > 0))
    }
    features = FindSpatiallyVariableFeatures(
      data_ori$cell_counts_ori[, idx], spatial.location = data_ori$cell_metadata[idx, ],
      selection.method = "moransi", nfeatures = nfeatures
    )
    moran[, ct] = features[rownames(moran), "p.value"]
  }

  selected_cells = names(which(data_ori$cell_type_proportion > ct_prop_threshold))
  dev_table_1 = dev_table[, selected_cells]
  dev_table_1[which((non_zero_prop[, selected_cells] < non_zero_prop_min) |
                     (non_zero_count[, selected_cells] < non_zero_count_min))] = 0
  moran_thres = sort(moran)[moran_top_n]
  dev_table_1[which(moran[, selected_cells] > moran_thres)] = 0

  ord = order(dev_table_1, decreasing = TRUE)
  top_idx = ord[1:ctsvg_top_n]
  coords = arrayInd(top_idx, dim(dev_table_1)); rownames(coords) = NULL

  res = vector("list", length(data_ori$cell_type_proportion))
  names(res) = names(data_ori$cell_type_proportion)
  for (i in seq_len(ctsvg_top_n)) {
    col = colnames(dev_table_1)[coords[i, 2]]; row = rownames(dev_table_1)[coords[i, 1]]
    res[[col]] = c(res[[col]], row)
  }

  special_genes = res
  for (name in names(special_genes)) {
    special_genes[[name]] = list(ct_svg = special_genes[[name]], marker = NULL)
  }
  if (!is.null(save_dir)) save(special_genes, file = paste0(save_dir, "special_genes.RData"))

  if (!is.null(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
    col_map = viridis_pal(option = "magma")(10)
    for (ct in names(res)) {
      idx_ct = which(data_ori$cell_metadata$type == ct)
      for (gene in res[[ct]]) {
        base_df = data_ori$cell_metadata[idx_ct, c("x","y"), drop = FALSE]
        z = as.numeric(scale(as.numeric(data_ori$cell_counts_ori[gene, idx_ct])))
        df = cbind(base_df, z = z)
        p = ggplot(df, aes(x, y, colour = z)) + geom_point(size = 2) +
          stat_density_2d(aes(x = x, y = y), colour = "black", linewidth = 0.25, alpha = 0.4) +
          scale_colour_gradientn(colors = col_map, name = "z") + coord_equal() + theme_bw(base_size = 11) +
          theme(panel.grid = element_blank(), legend.position = "right") + labs(title = paste0(ct, " — ", gene), x = NULL, y = NULL)
        fn = file.path(plot_dir, sprintf("[%.2f_%.2f_%d]_%s_%s.png", dev_table[gene, ct], non_zero_prop[gene, ct], non_zero_count[gene, ct],
                                          gsub("[^[:alnum:]_]+", "-", ct), gsub("[^[:alnum:]_]+", "-", gene)))
        ggsave(fn, p, width = width, height = height, dpi = dpi, bg = "white")
      }
    }
  }
  special_genes
}


# Collects, across the given seeds' packaged data_new_v3_<seed>.RData, genes whose max
# simulated cell count exceeds `max_count_threshold` -- extreme outlier genes to exclude
# before running test methods on the simulated data. Ground-truth ctSVGs are protected:
# they're excluded from the returned list rather than the whole run failing on overlap.
compute_genes_over300 = function(save_base_dir, seeds, max_count_threshold = 300, save_dir = NULL) {
  genes_over300 = c()
  ct_svg_all = c()
  for (seed in seeds) {
    load(paste0(save_base_dir, "data_new_v3_", seed, ".RData"))  # -> data_new
    temp_genes = names(which(apply(data_new$cell_counts, MARGIN = 1, FUN = max) > max_count_threshold))
    ct_svg_all = union(ct_svg_all, unique(unlist(lapply(data_new$special_genes, FUN = function(x) x$ct_svg))))
    genes_new = setdiff(temp_genes, genes_over300)
    if (length(genes_new) > 0) {
      message(paste0("New genes detected in dataset, ", paste0(genes_new, collapse = ", ")))
    }
    genes_over300 = union(genes_over300, temp_genes)
  }
  genes_over300 = setdiff(genes_over300, ct_svg_all)
  if (!is.null(save_dir)) save(genes_over300, file = paste0(save_dir, "extreme_genes_over300.RData"))
  genes_over300
}


# Loads the four common inputs for a named dataset.
# Returns list(save_base_dir, data_new, data_ori, boundary, genes_over300).
load_sim_inputs = function(data_name, seed) {
  data_name = match.arg(data_name, c("breast", "ovarian", "lymph"))

  cfg = switch(data_name,
    breast  = list(save_base_dir = "./data/realistic_packaged/breast/",
                   boundary_file = "breast_boundary.RData",
                   ori_file      = "breast_small_ori.RData"),
    ovarian = list(save_base_dir = "./data/realistic_packaged/ovarian/",
                   boundary_file = "ovarian_boundary.RData",
                   ori_file      = "ovarian_small_ori.RData"),
    lymph   = list(save_base_dir = "./data/realistic_packaged/lymph/",
                   boundary_file = "lymph_small_boundary.RData",
                   ori_file      = "lymph_small_ori.RData")
  )
  load(paste0(cfg$save_base_dir, "data_new_v3_", seed, ".RData"))  # -> data_new
  load(paste0(cfg$save_base_dir, "extreme_genes_over300.RData"))                          # -> genes_over300

  tmp = new.env()
  load(paste0(cfg$save_base_dir, cfg$ori_file), envir = tmp)
  data_ori = tmp[[ls(tmp)[1]]]

  tmp2 = new.env()
  load(paste0(cfg$save_base_dir, cfg$boundary_file), envir = tmp2)
  boundary = tmp2[[ls(tmp2)[1]]]

  list(
    save_base_dir = cfg$save_base_dir,
    data_new      = data_new,
    data_ori      = data_ori,
    boundary      = boundary,
    genes_over300 = genes_over300
  )
}


# Shared filtering logic for any sim-data object (pseudo_simulator or decompose_simulator).
# obj must expose: cell_counts, cell_counts_predrop, spot_counts, cell_metadata,
#   cell_type_proportion, gene_names, special_genes, seed, dispersion, spot_composition.
preprocess_sim_data = function(obj, genes_over300, ct_prop_threshold = 0.01) {
  if (any(rownames(obj$cell_counts) != rownames(obj$spot_counts))) {
    stop("Order of genes are different between cell count matrix and spot count matrix!")
  }
  genes_to_save = names(which(
    (!rownames(obj$cell_counts) %in% genes_over300) & (rowSums(obj$spot_counts) > 0)
  ))
  genes_to_remove = setdiff(obj$gene_names, genes_to_save)
  ct_svgs = unique(unlist(lapply(obj$special_genes, FUN = function(x) x$ct_svg)))
  if (any(genes_to_remove %in% ct_svgs)) {
    stop("ct_svg gene has usually high new count but cannot be removed!")
  } else {
    message(paste0("Removing gene due to high usually high new cell count: ",
                   paste(genes_to_remove, collapse = ", ")))
  }

  spot_comp = obj$spot_composition
  save_cell_types = names(which(obj$cell_type_proportion > ct_prop_threshold))

  if (any(rownames(obj$cell_metadata) != colnames(obj$cell_counts))) {
    stop("Cell order doesn't match for cell metadata and cell count data")
  }
  save_cells = which(
    (obj$cell_metadata$type %in% save_cell_types) &
    colSums(obj$cell_counts_predrop[genes_to_save, ]) > 0
  )

  spot_comp = spot_comp[, save_cell_types]

  if (any(colnames(obj$spot_counts) != rownames(spot_comp))) {
    stop("Spot counts and spot comp spot order don't match!")
  }

  spot_to_save = which(
    colSums(obj$spot_counts[genes_to_save, ]) > 0 & rowSums(spot_comp) > 0
  )

  spot_comp = t(apply(spot_comp[spot_to_save, ], MARGIN = 1, FUN = function(x) x / sum(x)))
  if (any(round(rowSums(spot_comp), 2) != 1)) {
    stop("Spot composition doesn't add up to 1 for some rows!")
  }

  sim_obj = list()
  sim_obj[["gene_names"]]           = obj$gene_names[which(obj$gene_names %in% genes_to_save)]
  sim_obj[["cell_metadata"]]        = obj$cell_metadata[save_cells, ]
  sim_obj[["special_genes"]]        = obj$special_genes
  sim_obj[["seed"]]                 = obj$seed
  sim_obj[["cell_type_proportion"]] = obj$cell_type_proportion[save_cell_types]
  sim_obj[["dispersion"]]           = obj$dispersion[genes_to_save, save_cells]
  sim_obj[["spot_composition"]]     = spot_comp

  list(
    genes_to_save   = genes_to_save,
    save_cells      = save_cells,
    save_cell_types = save_cell_types,
    spot_to_save    = spot_to_save,
    spot_comp       = spot_comp,
    sim_obj         = sim_obj
  )
}


load_process_data_new = function(data_name, seed, rotate = 0, ct_prop_threshold  = 0.01) {
  library(Triangulation)
  inputs = load_sim_inputs(data_name, seed)

  pp = preprocess_sim_data(inputs$data_new, inputs$genes_over300, ct_prop_threshold)

  spot_coords_rotated = rotate_translate_points(
    inputs$data_new$spot_coords[pp$spot_to_save, ], angle_degrees = rotate, type = "C")
  boundary_rotated = rotate_translate_points(inputs$boundary, angle_degrees = rotate, type = "C")
  Tr_cell_rotated  = remove_colinear(TriMesh(boundary_rotated, n = 2))

  list(
    spot_count          = inputs$data_new$spot_counts[pp$genes_to_save, pp$spot_to_save],
    spot_composition    = pp$spot_comp,
    spot_coords         = inputs$data_new$spot_coords[pp$spot_to_save, ],
    spot_coords_rotated = spot_coords_rotated,
    cell_counts         = inputs$data_new$cell_counts[pp$genes_to_save, pp$save_cells],
    cell_counts_predrop = inputs$data_new$cell_counts_predrop[pp$genes_to_save, pp$save_cells],
    cell_metadata       = inputs$data_new$cell_metadata[pp$save_cells, ],
    boundary            = inputs$boundary,
    boundary_rotated    = boundary_rotated,
    Tr_cell_rotated     = Tr_cell_rotated,
    sim_obj             = pp$sim_obj
  )
}

# backward-compatible alias
data_new_fast_loader = load_process_data_new


load_process_decompose = function(
  data_name, seed, # data_new, data_ori parameter 
  sim_type,
  dispersion           = 0.7,
  # n_markers            = 25, # this follows previous runs (after update marker rule, let's not pass marker)
  cell_type_proportion = c("Type1" = 0.1, "Type2" = 0.3, "Type3" = 0.6),
  ct_prop_threshold    = 0.01
) {
  library(Triangulation)
  inputs = load_sim_inputs(data_name, seed)


  decompose = decompose_simulator$new(
    inputs$data_new, inputs$data_ori,
    dispersion           = dispersion, seed = seed,
    sim_type             = sim_type, 
    cell_type_proportion = cell_type_proportion
  )

  pp = preprocess_sim_data(decompose, inputs$genes_over300, ct_prop_threshold)

  # use real spatial triangulation when sim_type includes 2, else bounding-box
  if (2 %in% sim_type) {
    Tr_cell  = remove_colinear(TriMesh(inputs$boundary, n = 2))
  } else {
    coords   = decompose$cell_metadata
    boundary = matrix(
      c(min(coords$x), min(coords$y),
        max(coords$x), min(coords$y),
        max(coords$x), max(coords$y),
        min(coords$x), max(coords$y)),
      ncol = 2, byrow = TRUE
    )
    Tr_cell = TriMesh(boundary, n = 2)
  }

  list(
    spot_count          = decompose$spot_counts[pp$genes_to_save, pp$spot_to_save],
    spot_composition    = pp$spot_comp,
    spot_coords         = decompose$spot_coords[pp$spot_to_save, ],
    cell_counts         = decompose$cell_counts[pp$genes_to_save, pp$save_cells],
    cell_counts_predrop = decompose$cell_counts_predrop[pp$genes_to_save, pp$save_cells],
    cell_metadata       = decompose$cell_metadata[pp$save_cells, ],
    Tr_cell             = Tr_cell,
    sim_obj             = pp$sim_obj,
    decompose           = decompose
  )
}



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

generate_simulator = function(data_ori, scdesign_new_path, save_path, seed, special_genes){
  source("./simulators/realistic/pseudo_spot_simulator.R")
  load(scdesign_new_path)
  data_new = pseudo_simulator$new(
    data_ori, scdesign_new$newcount,scdesign_new$para, seed = seed, special_genes = special_genes
  )
  save(data_new, file = save_path)
  message("Simulator generated...")
  return(data_new)
}


collect_details = function(input_dir, analyze_method = c("celina", "stance", "spvc", "ctsv", "cside", "mmm")){
  detail_result = vector("list", length(analyze_method))
  names(detail_result) = analyze_method

  collect_files = function(pattern, full_names = TRUE) {
    unlist(lapply(input_dir, function(dir) {
      list.files(path = dir, pattern = pattern, full.names = full_names)
    }))
  }

  # - read results
  if ("celina" %in% analyze_method){
    celina_files <- collect_files("^celina_?([0-9]*)\\.xlsx$")

    message(paste0("CELINA | Reading ", length(celina_files), " result: ", paste0(basename(celina_files), collapse = ", ")))

    if (length(celina_files) > 0) {
      celina_result <- do.call(
        rbind,
        lapply(celina_files, function(file) {
          readxl::read_excel(file)
        })
      )

      # if (!("p_adj" %in% colnames(celina_result))){
      #   stop("Celina result lack p adjusted values, add them first")
      # }
      # celina_result$p_fail = 0

      celina_result$parameter_setting = ""
      detail_result[["celina"]] = celina_result
    }
  }

  if ("stance" %in% analyze_method){
    stance_files <- collect_files("^stance_?([0-9]*)\\.xlsx$")

    message(paste0("STANCE | Reading ", length(stance_files), " result: ", paste0(basename(stance_files), collapse = ", ")))

    if (length(stance_files) > 0) {
      stance_result <- do.call(
        rbind,
        lapply(stance_files, function(file) {
          readxl::read_excel(file)
        })
      )

      # if (!("p_adj" %in% colnames(stance_result))){
      #   stop("Stance result lack p adjusted values, add them first")
      # }
      # stance_result$p_fail = 0
      stance_result$parameter_setting = ""
      detail_result[["stance"]] = stance_result
    }
  }

  if ("ctsv" %in% analyze_method){
    ctsv_files <- collect_files("^ctsv_?([0-9]*)\\.csv$")

    message(paste0("CTSV | Reading ", length(ctsv_files), " result: ", paste0(basename(ctsv_files), collapse = ", ")))

    if (length(ctsv_files) > 0) {
      ctsv_result = do.call(
        rbind,
        lapply(ctsv_files, function(file) {
          read.csv(file)
        })
      )

      # ctsv_result$p_fail = 0
      ctsv_result$parameter_setting = ""
      detail_result[["ctsv"]] = ctsv_result
    }
  }

  if ("cside" %in% analyze_method){
    cside_files <- collect_files("^cside_?([0-9]*)\\.csv$")

    message(paste0("CSIDE | Reading ", length(cside_files), " result: ", paste0(basename(cside_files), collapse = ", ")))

    if (length(cside_files) > 0) {
      cside_result = do.call(
        rbind,
        lapply(cside_files, function(file) {
          read.csv(file)
        })
      )
      # cside_result$p_fail = 0
      detail_result[["cside"]] = cside_result
    }
  }

  if ("spvc" %in% analyze_method){
    spvc_files <- collect_files("^spvc_?([0-9]*)\\.csv$")

    message(paste0("SPVC | Reading ", length(spvc_files), " result: ", paste0(basename(spvc_files), collapse = ", ")))

    if (length(spvc_files) > 0) {
      spvc_result = do.call(
        rbind,
        lapply(spvc_files, function(file) {
          read.csv(file)
        })
      )
      spvc_result$parameter_setting = ""
      detail_result[["spvc"]] = spvc_result
    }
  }

  if ("mmm" %in% analyze_method){
    mmm_files <- collect_files("^mmm_?([0-9]*)\\.xlsx$")

    message(paste0("MMM | Reading ", length(mmm_files), " result: ", paste0(basename(mmm_files), collapse = ", ")))

    if (length(mmm_files) > 0) {
      mmm_result <- do.call(
        rbind,
        lapply(mmm_files, function(file) {
          readxl::read_excel(file)
        })
      )

      # mmm_result$p_fail = 0
      detail_result[["mmm"]] = mmm_result
    }
  }


  return(detail_result)
}

### --- analyze results
write_sheet = function(write_df, save_path, sheet_name){
  if (file.exists(save_path)){
    wb = openxlsx::loadWorkbook(save_path)
    if (sheet_name %in% names(wb)) {
      openxlsx::removeWorksheet(wb, sheet = sheet_name)
    }
    openxlsx::addWorksheet(wb, sheet_name)
    openxlsx::writeData(wb, sheet = sheet_name, x = write_df)
    openxlsx::saveWorkbook(wb, save_path, overwrite = TRUE)
  } else {
    writexl::write_xlsx(setNames(list(write_df), sheet_name), path = save_path)
  }
}

analyze_result = function(input_dir = NULL, save_path, save_sheet, detail_df = NULL, analyze_method = c("celina", "stance", "spvc", "ctsv", "cside", "mmm"), threshold = c(0.01, 0.03, 0.05, 0.07, 0.09, 0.1, 0.2)){
  if (!is.null(detail_df)) {
    keep = sapply(detail_df$test_method, function(m) any(startsWith(m, analyze_method)))
    detail_result = split(detail_df[keep, ], detail_df$test_method[keep])
  } else {
    detail_result = collect_details(input_dir, analyze_method)
  }
  summary_result = list()
  for (thres in threshold){
    for (method in names(Filter(Negate(is.null), detail_result))){
      summary_result[[paste0(method, "-",thres)]] = detail_result[[method]]%>%
        mutate(
          has_ctsvg = sapply(strsplit(gene_type, ","), function(x) "ctsvg" %in% x),
          has_marker = sapply(strsplit(gene_type, ","), function(x) "marker" %in% x), 
          tested = (p_value >= 0)
        ) %>%
        ## --- first group operation for ep calculation 
        group_by(simulation_name, test_method, deconv_method, seed, parameter_setting) %>%
        mutate(
          ep_denom = sum(has_ctsvg, na.rm = TRUE),
          padj_rank = NA_real_
        ) %>%
        mutate(
          padj_rank = replace(
            padj_rank,
            tested,
            rank(p_adj[tested], ties.method = "first")
          ),
          in_top_ep = tested & (padj_rank <= ep_denom)
        ) %>%
        ungroup()%>%
        ## --- second group operation 
        group_by(simulation_name, test_method, deconv_method, seed, cell_type, cell_proportion, gene_type, parameter_setting)%>%
        summarise(
          # --- fdp & power 
          rejected_adj = sum(p_value>=0 & p_adj < thres), 
          rejected_raw = sum(p_value>=0 & p_value < thres),
          tp = sum(p_value>=0 & p_adj < thres & has_ctsvg),
          total_true = sum(has_ctsvg),
          dispersion = mean(dispersion),
          
          # --- fdp without marker 
          rejected_adj_no_marker = sum(p_value>=0 & p_adj < thres & ifelse(has_ctsvg, TRUE, !has_marker)), 
          
          true_tested = sum(p_value >= 0 & has_ctsvg),
          tested = sum(p_value >= 0), 
          total = n(),
          threshold = thres,
          
          # --- type 1 & power 
          type1_num = sum((p_value >= 0) & (p_value < thres) & (!has_ctsvg)),
          type1_den = sum(!has_ctsvg),
          power_num = sum((p_value >= 0) & (p_value < thres) & has_ctsvg),
          power_den = sum(has_ctsvg),
          
          # convergence
          convergence_ctsvg = sum(p_value == 1 & has_ctsvg),
          convergence_null = sum(p_value == 1 & !has_ctsvg),
          
          # --- ep 
          ep_denom_check = mean(ep_denom), 
          in_top = sum(in_top_ep), 
          ep_num = sum(in_top_ep & has_ctsvg),
          
          .groups = "drop"
        )
    }
  }
  
  summary = do.call(rbind, summary_result)

  # --- simulation-level metadata fields
  sim_setup_path <- "./simulators/idealized/simulator_setup.csv"
  if (file.exists(sim_setup_path)) {
    sim_setup_df <- read.csv(sim_setup_path, stringsAsFactors = FALSE)
    sim_setup_df$pattern <- trimws(sim_setup_df$pattern)
    pure_names <- vapply(seq_len(nrow(sim_setup_df)), function(i) {
      row <- sim_setup_df[i, ]
      paste0(
        row$sim_name,
        ifelse(!is.na(row$phi)     & row$phi     != "", paste0("_", row$phi),     ""),
        ifelse(!is.na(row$scene)   & row$scene   != "", paste0("_", row$scene),   ""),
        ifelse(!is.na(row$pattern) & row$pattern != "", paste0("_", row$pattern), "")
      )
    }, character(1))
  } else {
    pure_names <- character(0)
  }

  summary <- summary %>%
    group_by(simulation_name, seed) %>%
    mutate(n_cell_type = n_distinct(cell_type)) %>%
    ungroup() %>%
    mutate(
      sim_type = case_when(
        simulation_name %in% pure_names                                  ~ "idealized",
        grepl("^(breast|ovarian|lymph)_decompose_exp_", simulation_name) ~ "decomposed",
        grepl("^(breast|ovarian|lymph)_rotate_[0-9]+$", simulation_name) ~ "realistic",
        TRUE                                                             ~ NA_character_
      ),
      dataset_name = case_when(
        sim_type %in% c("decomposed", "realistic") ~
          sub("^(breast|ovarian|lymph)_.*$", "\\1", simulation_name),
        TRUE ~ NA_character_
      ),
      additional_feature = case_when(
        sim_type == "decomposed" ~
          sub("^(breast|ovarian|lymph)_decompose_exp_", "", simulation_name),
        sim_type == "realistic" ~
          sub("^(breast|ovarian|lymph)_(rotate_[0-9]+)$", "\\2", simulation_name),
        TRUE ~ NA_character_
      )
    )

  if (is.na(save_path)){
    return(summary)
  } else {
    write_sheet(summary, save_path = save_path, sheet_name = save_sheet)
  }
}




library(ggplot2)
library(scales)   # for viridis_pal

# Append `new_df` to `file_path` (format inferred from its .csv/.xlsx extension),
# combining with any existing data already there. Once the combined row count
# reaches `max_rows`, instead switch to a new numbered file (e.g. mmm_1.xlsx)
# holding just `new_df`, rather than growing the existing file forever.
save_result_with_file_switching = function(new_df, file_path, prefix, max_rows = 300000){
  if (is.na(file_path)) return(invisible(NULL))

  ext = tools::file_ext(file_path)
  dir_name = dirname(file_path)
  read_fn = if (ext == "xlsx") readxl::read_xlsx else read.csv
  write_fn = if (ext == "xlsx") {
    function(df, path) writexl::write_xlsx(list(Results = df), path = path)
  } else {
    function(df, path) write.csv(df, file = path, row.names = FALSE)
  }

  if (file.exists(file_path)) {
    old_df = read_fn(file_path)
    combined_df = rbind(old_df, new_df)
  } else {
    combined_df = new_df
  }

  if (nrow(combined_df) >= max_rows){
    current_num = sub(paste0(".*", prefix, "_?([0-9]*)\\.", ext, "$"), "\\1", file_path)
    current_num = if (current_num == "") 0 else as.integer(current_num)

    if (is.na(current_num)){
      stop(sprintf("File observations exceed %d, automatic new file logic failed!", max_rows))
    }

    new_path = paste0(dir_name, "/", prefix, "_", current_num + 1, ".", ext)
    message(paste0("Switch ", prefix, " results to ", new_path))
    write_fn(new_df, new_path)
  } else{
    write_fn(combined_df, file_path)
  }

  invisible(NULL)
}

# Label each gene in `all_gene` as "marker"/"ctsvg"/"other_marker"/"other_ctsvg" (comma-joined if
# multiple apply, "Null" if none) based on sim_obj$special_genes for cell type `ct` vs. the rest.
label_gene_types = function(all_gene, sim_obj, ct, cell_types){
  marker = sim_obj$special_genes[[ct]]$marker
  ctsvg = sim_obj$special_genes[[ct]]$ct_svg
  other_marker = c()
  other_ctsvg = c()
  for (oct in setdiff(cell_types, ct)){
    other_marker = c(other_marker, sim_obj$special_genes[[oct]]$marker)
    other_ctsvg = c(other_ctsvg, sim_obj$special_genes[[oct]]$ct_svg)
  }
  check_list = list(marker = marker, ctsvg = ctsvg, other_marker = other_marker, other_ctsvg = other_ctsvg)

  vapply(all_gene, function(gene){
    gene_labels = names(check_list)[sapply(check_list, function(l) gene %in% l)]
    if (length(gene_labels) == 0) "Null" else paste(gene_labels, collapse = ",")
  }, character(1), USE.NAMES = FALSE)
}

retrieve_file_start_with = function(save_dir, test_method, format = "xlsx"){
  all_files <- list.files(
    path = save_dir,
    pattern = paste0("^", test_method, "_?([0-9]*)\\.", format, "$"),
    full.names = FALSE
  )

  message(paste0("All files for ", test_method, ": ", paste0(all_files, collapse = ", ")))
  
  if (length(all_files) == 0){
    return(NA)
  } else if (length(all_files) == 1){
    return(all_files[1])
  } else{
    numbered_files = setdiff(all_files, paste0(test_method, ".", format))
    nums = as.integer(sub(paste0(".*", test_method, "_([0-9]+)\\.", format, "$"), "\\1", numbered_files))
    last_file = numbered_files[which.max(nums)]
    message(paste0("Saving new result to ", last_file))
    return(last_file)
  }
}


# rotate_translate_points always ensures x, y
# coordinates are positive
# add type parameter 20260119:
#   A - adjust x, y separately and only adjust if min less than 0 
#   B - adjust x, y together, and only adjust the total min less than 0 
#   c - adjust unconditionally and bring x, y to 0, 0 
rotate_translate_points <- function(original_points, angle_degrees, type = "A") {

  if (!type %in% c("A", "B", "C")){
    stop("Choose translate type A, B, C")
  }

  if (angle_degrees == 0) {
    output = as.matrix(original_points)
    colnames(output) = c("x", "y")
    return(output)
  }
  
  angle_radians <- angle_degrees * (pi / 180) # Convert degrees to radians
  rotation_matrix <- matrix(c(cos(angle_radians), -sin(angle_radians),
                              sin(angle_radians), cos(angle_radians)),
                            nrow = 2, byrow = TRUE)
  rotated_points <- t(rotation_matrix %*% t(original_points))
  output <- as.matrix(rotated_points)
  # perform translation, the lowest x, y can be is (0, 0)
  min.x = min(output[,1]) 
  min.y = min(output[,2])
  
  # different translate strategies 
  if (type == "A"){
    if (min.x < 0) output[,1] = output[,1] + abs(min.x)
    if (min.y < 0) output[,2] = output[,2] + abs(min.y)
  } else if (type == "B"){
    adjust = abs(min(min(min.x, min.y),0))
    output[,1] = output[,1]+adjust
    output[,2] = output[,2]+adjust
  } else if (type == "C"){ # bring to (0, 0) regardless 
    output[,1] = output[,1] - min.x
    output[,2] = output[,2] - min.y 
  }
  

  colnames(output) <- c('x', 'y')
  return(output)
}

dist_perserve_trans = function(X){
  Xc = sweep(X, 2, colMeans(X), "-")
  
  # step 2: calculate pairwise distance 
  scaler = max(as.matrix(dist(Xc)))
  X_normal = Xc/scaler
  return(X_normal)
}

are_colinear = function(A, B, C, tol = 1e-6){
  abs(det(rbind(B-A, C-A))) < tol
}

remove_colinear = function(Tr.cell){
  save = c()
  for (i in 1:nrow(Tr.cell$Tr)){
    A = Tr.cell$V[Tr.cell$Tr[i, 1], ]
    B = Tr.cell$V[Tr.cell$Tr[i, 2], ]
    C = Tr.cell$V[Tr.cell$Tr[i, 3], ]
    if (!are_colinear(A, B, C)){
      save = c(save, i)
    }
  }
  Tr.cell$Tr = Tr.cell$Tr[save, ]
  return(Tr.cell)
}

## --- Real-data configuration lookup -------------------------------------------
## Returns a list with data_path, tr_cell_path, sim_name, save_dir for the given
## (data_name, filter_type) pair.  Run parameters are intentionally excluded.
get_real_data_config = function(data_name, filter_type) {
  configs = list(
    lung = list(
      filter1 = list(
        data_path    = "./data/real_packaged/lung/filter1.RData",
        tr_cell_path = "./data/real_packaged/lung/Tr.cell.RData",
        sim_name     = "lung",
        save_dir     = "./data/detection_results/real_data/lung/"
      ),
      filter2 = list(
        data_path    = "./data/real_packaged/lung/filter2.RData",
        tr_cell_path = "./data/real_packaged/lung/Tr.cell.RData",
        sim_name     = "lung_filtered_0.1",
        save_dir     = "./data/detection_results/real_data/lung/"
      )
    ),
    kidney_core = list(
      filter1 = list(
        data_path    = "./data/real_packaged/kidney_core/filter1.RData",
        tr_cell_path = "./data/real_packaged/kidney_core/Tr.cell.RData",
        sim_name     = "kidney_core",
        save_dir     = "./data/detection_results/real_data/kidney/"
      ),
      filter2 = list(
        data_path    = "./data/real_packaged/kidney_core/filter2.RData",
        tr_cell_path = "./data/real_packaged/kidney_core/Tr.cell.RData",
        sim_name     = "kidney_core_filtered_0.1",
        save_dir     = "./data/detection_results/real_data/kidney/"
      )
    ),
    kidney_inter = list(
      filter1 = list(
        data_path    = "./data/real_packaged/kidney_inter/filter1.RData",
        tr_cell_path = "./data/real_packaged/kidney_inter/Tr.cell.RData",
        sim_name     = "kidney_inter",
        save_dir     = "./data/detection_results/real_data/kidney/"
      ),
      filter2 = list(
        data_path    = "./data/real_packaged/kidney_inter/filter2.RData",
        tr_cell_path = "./data/real_packaged/kidney_inter/Tr.cell.RData",
        sim_name     = "kidney_inter_filtered_0.1",
        save_dir     = "./data/detection_results/real_data/kidney/"
      )
    ),
    breast = list(
      filter1 = list(
        data_path    = "./data/real_packaged/breast/filter1.RData",
        tr_cell_path = "./data/real_packaged/breast/Tr.cell.RData",
        sim_name     = "breast",
        save_dir     = "./data/detection_results/real_data/breast/"
      )
      # breast has no filter2
    ),
    brain = list(
      filter1 = list(
        data_path    = "./data/real_packaged/dlpfc/filter1.RData",
        tr_cell_path = "./data/real_packaged/dlpfc/Tr.cell.RData",
        sim_name     = "brain",
        save_dir     = "./data/detection_results/real_data/dlpfc/"
      ),
      filter2 = list(
        data_path    = "./data/real_packaged/dlpfc/filter2.RData",
        tr_cell_path = "./data/real_packaged/dlpfc/Tr.cell.RData",
        sim_name     = "brain_filtered_0.1",
        save_dir     = "./data/detection_results/real_data/dlpfc/"
      )
    )
  )

  valid_datasets = names(configs)
  if (!data_name %in% valid_datasets)
    stop(sprintf("Unknown data_name '%s'. Valid options: %s",
                 data_name, paste(valid_datasets, collapse = ", ")))

  valid_filters = names(configs[[data_name]])
  if (!filter_type %in% valid_filters)
    stop(sprintf("No config for %s / %s. Valid filter types: %s",
                 data_name, filter_type, paste(valid_filters, collapse = ", ")))

  configs[[data_name]][[filter_type]]
}


source("./simulators/celina_simulator_alt.R")
source("./simulators/celina_simulator_null.R")
source("./simulators/stance_simulator_1alt.R")
source("./simulators/stance_simulator_alt.R")

run_simulator = function(name, seed, phi = 0.95, scene = "I", pattern = "streak", control_UMI = FALSE){
  if (name == "stance_simulator_alt"){
    temp_sim = stance_simulator_alt$new(dispersion = phi, seed = seed, cell_type_proportion = c(0.1, 0.3, 0.6))
  } else if (name == "stance_simulator_1alt"){
    temp_sim = stance_simulator_1alt$new(dispersion = phi, seed = seed)
  } else if (name == "celina_simulator_alt"){
    if (!control_UMI){
      fold_change = c(ifelse(pattern == "hotspot", 2, 1.5), ifelse(pattern == "hotspot", 0.4, 0.5))
    } else{
      fold_change = c(ifelse(pattern == "hotspot", 47/36, 1.5), ifelse(pattern == "hotspot", 1/18, 0.5))
    }
    pi = ifelse(pattern == "hotspot", 0.4, ifelse(pattern == "streak", 0.25, 0.15))
    temp_sim = celina_simulator_alt$new(
      dispersion = 0.95, seed = seed, scenario = scene, resolution = "spot", pattern_type = pattern,
      pi = pi, fold_change_high = fold_change[1], fold_change_low = fold_change[2]
    )
  } else if (name == "celina_simulator_null"){
    temp_sim = celina_simulator_null$new(dispersion = 0.95, seed= seed, scenario = scene, resolution = "spot")
  }
  return(temp_sim)
}

