library(ggplot2)
library(arrow)
library(Seurat)
library(Matrix)
library(readr)
library(Triangulation)
source("./util.R")

## =============================================================
## User inputs — edit this section before each run
## =============================================================

# --- which (dataset, filter) pairs to run ---
# valid data_names : "lung", "kidney_core", "kidney_inter", "breast", "brain"
# valid filter_types: "filter1", "filter2"  (breast has filter1 only)

data_names   <- c("breast")
filter_types <- c("filter1")


# --- hardware ---
ncores       <- 40

# --- time limit per method call ---
timeout_sec  <- 48*3600  

# --- which methods to run ---
# options: "celina", "stance", "spvc", "ctsv", "cside", "mmm"
run_method   <- c("celina", "stance", "spvc", "ctsv", "cside", "mmm")

# --- STANCE-specific threshold ---
utsvg_thres  <- 0.01   # use 0.01 for filter2 runs

## =============================================================
## Run ctsvg detection 
## =============================================================

for (i in seq_along(data_names)) {
  cfg = get_real_data_config(data_names[i], filter_types[i])
  message(sprintf("=== %s / %s  [sim_name: %s] ===",
                  data_names[i], filter_types[i], cfg$sim_name))

  load(cfg$data_path)     # loads 'filter1' or 'filter2' into global env
  load(cfg$tr_cell_path)  # loads 'Tr.cell' into global env

  dat = get(filter_types[i])

  run_all_tests(
    sp_count     = as.matrix(dat$spot_counts),
    sp_comp      = as.matrix(dat$spot_composition),
    sp_coords    = as.matrix(dat$spot_coords),
    sc_count     = as.matrix(dat$cell_counts),
    sc_metadata  = dat$cell_metadata,
    spvc_tri     = Tr.cell,
    sim_name     = cfg$sim_name,
    save_dir     = cfg$save_dir,
    sim_obj      = dat,
    ncores       = ncores,
    run_method   = run_method,
    # utsvg_thres  = utsvg_thres,
    timeout_sec  = timeout_sec
  )
}
