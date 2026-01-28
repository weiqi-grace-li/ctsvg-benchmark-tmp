# src/evaluation/run_all_tests.R
# Experiment runner: runs selected methods, formats outputs, and saves detail tables.
# IMPORTANT:
# - Assumes wrappers + formatters + utils are already loaded (via src/load_all.R)
# - Does NOT define wrappers/formatters itself.
# - Uses save_results() from src/utils/io_results.R for chunked saving.

run_all_tests <- function(
    sp_count, sp_comp, sp_coords,
    sc_count, sc_metadata,
    spvc_tri,
    sim_name,
    save_dir,
    time_save_path,
    sim_obj,
    ncores = 4,
    ncores_cside = NULL,
    run_method = c("celina", "stance", "spvc", "ctsv", "cside", "stance-alt", "spvc-gam"),
    timeout_sec = Inf,
    utsvg_thres = 0.05,
    topn_utsvg = NULL
) {
  # --- basic checks
  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
  
  # --- required helper (you already have this elsewhere)
  if (!exists("record_run", mode = "function", inherits = TRUE)) {
    stop("record_run() not found. Load it before calling run_all_tests().")
  }
  
  # normalize run_method
  run_method <- unique(run_method)
  
  # -----------------------
  # 1) CELINA
  # -----------------------
  if ("celina" %in% run_method) {
    message("Running celina ...")
    
    myCELINA <- record_run(
      sim_name, "celina", ncores, time_save_path,
      function() {
        run_celina_spot(
          count = sp_count,
          pos = sp_coords,
          comp = sp_comp,
          sc_count = sc_count,
          sc_type = sc_metadata$type,
          approximation = FALSE,
          ncores = ncores
        )
      },
      timeout_sec = timeout_sec
    )
    
    df <- celina_details_to_df(
      sim_obj = sim_obj,
      celina_obj = myCELINA,
      sim_name = sim_name,
      deconv_method = "oracle"
    )
    
    save_results(df, save_dir, method = "celina", format = "xlsx", max_rows = 1e6)
  }
  
  # -----------------------
  # 2) CSIDE
  # -----------------------
  if ("cside" %in% run_method) {
    message("Running cside ...")
    
    set.seed(42)
    if (nrow(sc_metadata) > 10000) {
      selected <- sample(rownames(sc_metadata), 10000, replace = FALSE)
    } else {
      selected <- 1:nrow(sc_metadata)
    }
    
    cs_cores <- ifelse(is.null(ncores_cside), ncores, ncores_cside)
    
    myCSIDE <- record_run(
      sim_name, "cside", cs_cores, time_save_path,
      function() {
        run_cside_import(
          sp_counts = round(sp_count, 0),
          sp_coords = sp_coords,
          sp_comp = sp_comp,
          sc_counts = round(sc_count[, selected], 0),
          sc_metadata = sc_metadata[selected, , drop = FALSE],
          df = 5,
          cell_type_threshold = 0,
          ncores = cs_cores,
          cell_min = 0,
          counts_min = 0
        )
      },
      timeout_sec = timeout_sec
    )
    
    df <- cside_details_to_df(
      sim_obj = sim_obj,
      cside_obj = myCSIDE,
      sim_name = sim_name,
      deconv_method = "oracle",
      cell_threshold = 0,
      df = 5,
      weight_threshold = 0.8
    )
    
    save_results(df, save_dir, method = "cside", format = "csv", max_rows = 5e5)
  }
  
  # -----------------------
  # 3) CTSV
  # -----------------------
  if ("ctsv" %in% run_method) {
    message("Running ctsv ...")
    
    my_ctsv <- record_run(
      sim_name, "ctsv", ncores, time_save_path,
      function() {
        run_ctsv(
          sp_counts = sp_count,
          sp_coords = sp_coords,
          sp_comp = sp_comp,
          ncores = ncores
        )
      },
      timeout_sec = timeout_sec
    )
    
    df <- ctsv_details_to_df(
      sim_obj = sim_obj,
      ctsv_obj = my_ctsv,
      sim_name = sim_name,
      deconv_method = "oracle"
    )
    
    save_results(df, save_dir, method = "ctsv", format = "csv", max_rows = 5e5)
  }
  
  # -----------------------
  # 4) STANCE
  # -----------------------
  if ("stance" %in% run_method) {
    message("Running stance ...")
    
    mySTANCE <- record_run(
      sim_name, "stance", ncores, time_save_path,
      function() {
        run_stance(
          count = sp_count,
          pos = sp_coords,
          comp = sp_comp,
          ncores = ncores,
          spot_thres = 0,
          utsvg_thres = utsvg_thres,
          topn_utsvg = topn_utsvg
        )
      },
      timeout_sec = timeout_sec
    )
    
    df <- stance_details_to_df(
      sim_obj = sim_obj,
      stance_obj = mySTANCE,
      sim_name = sim_name,
      deconv_method = "oracle",
      test_method = "stance"
    )
    
    save_results(df, save_dir, method = "stance", format = "xlsx", max_rows = 1e6)
  }
  
  # -----------------------
  # 5) STANCE-ALT
  # -----------------------
  if ("stance-alt" %in% run_method) {
    message("Running stance-alt ...")
    
    if (!exists("run_stance_alt", mode = "function", inherits = TRUE)) {
      stop("run_stance_alt() not found. Load it before calling stance-alt.")
    }
    
    mySTANCE_alt <- record_run(
      sim_name, "stance-alt", ncores, time_save_path,
      function() {
        run_stance_alt(
          count = sp_count,
          pos = sp_coords,
          comp = sp_comp,
          ncores = ncores,
          spot_thres = 0
        )
      },
      timeout_sec = timeout_sec
    )
    
    df <- stance_details_to_df(
      sim_obj = sim_obj,
      stance_obj = mySTANCE_alt,
      sim_name = sim_name,
      deconv_method = "oracle",
      test_method = "stance-alt"
    )
    
    # keep saving under "stance" prefix (matches your previous behavior)
    save_results(df, save_dir, method = "stance", format = "xlsx", max_rows = 1e6)
  }
  
  # -----------------------
  # 6) SPVC (original)
  # -----------------------
  if ("spvc" %in% run_method) {
    message("Running spvc-original ...")
    
    # prefer your existing name if present; otherwise fall back to run_spvc_alt
    spvc_runner <- NULL
    if (exists("run_spvc_full_square_alt", mode = "function", inherits = TRUE)) {
      spvc_runner <- get("run_spvc_full_square_alt", mode = "function", inherits = TRUE)
    } else if (exists("run_spvc_alt", mode = "function", inherits = TRUE)) {
      spvc_runner <- get("run_spvc_alt", mode = "function", inherits = TRUE)
    } else {
      stop("No spVC runner found (run_spvc_full_square_alt or run_spvc_alt).")
    }
    
    myspVC <- record_run(
      sim_name, "spvc-original", ncores, time_save_path,
      function() {
        spvc_runner(
          sp_counts = sp_count,
          sp_coords = sp_coords,
          sp_comp = sp_comp,
          Tr.cell = spvc_tri,
          ncores = ncores,
          ori = FALSE
        )
      },
      timeout_sec = timeout_sec
    )
    
    df <- spvc_details_to_df(
      sim_obj = sim_obj,
      spvc_obj = myspVC,
      sim_name = sim_name,
      deconv_method = "oracle",
      test_method = "spvc-original",
      two_step = FALSE
    )
    
    save_results(df, save_dir, method = "spvc", format = "csv", max_rows = 5e5)
  }
  
  # -----------------------
  # 7) SPVC-GAM
  # -----------------------
  if ("spvc-gam" %in% run_method) {
    message("Running spvc-gam ...")
    
    if (!exists("run_spvc_full_square_gam", mode = "function", inherits = TRUE)) {
      stop("run_spvc_full_square_gam() not found. Load it before calling spvc-gam.")
    }
    
    myspVC_gam <- record_run(
      sim_name, "spvc-gam", ncores, time_save_path,
      function() {
        run_spvc_full_square_gam(
          sp_counts = sp_count,
          sp_coords = sp_coords,
          sp_comp = sp_comp,
          Tr.cell = spvc_tri,
          ncores = ncores,
          ori = FALSE
        )
      },
      timeout_sec = timeout_sec
    )
    
    df <- spvc_details_to_df(
      sim_obj = sim_obj,
      spvc_obj = myspVC_gam,
      sim_name = sim_name,
      deconv_method = "oracle",
      test_method = "spvc-gam",
      two_step = FALSE
    )
    
    save_results(df, save_dir, method = "spvc", format = "csv", max_rows = 5e5)
  }
  
  invisible(TRUE)
}
