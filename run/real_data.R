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

# processed breast data is included in the repo
# make you have either processed the other real data or 
# downloaded processed real data to run them 
data_names   <- c("breast")
filter_types <- c("filter1")


# --- hardware ---
ncores       <- 20

# --- time limit per method call ---
timeout_sec  <- 48*3600  

# --- which methods to run ---
# options: "celina", "stance", "spvc", "ctsv", "cside", "mmm"
# run_method   <- c("celina", "stance", "spvc", "ctsv", "cside", "mmm") # uncomment to run all methods 
run_method = c("celina")

# --- STANCE-specific threshold ---
utsvg_thres  <- 0.01   # use 0.01 for filter2 runs

## =============================================================
## Run ctsvg detection 
## =============================================================

for (i in seq_along(data_names)) {
  cfg = get_real_data_config(data_names[i], filter_types[i])
  message(sprintf("=== %s / %s  [sim_name: %s] ===",
                  data_names[i], filter_types[i], cfg$sim_name))

  data_env  = new.env()
  data_obj  = load(cfg$data_path, envir = data_env)     # actual object name varies by dataset (e.g. 'breast', not 'filter1')
  dat       = get(data_obj[1], envir = data_env)

  tr_env    = new.env()
  tr_obj    = load(cfg$tr_cell_path, envir = tr_env)    # detect actual object name rather than assuming 'Tr.cell'
  Tr.cell   = get(tr_obj[1], envir = tr_env)

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
    utsvg_thres  = utsvg_thres,
    timeout_sec  = timeout_sec
  )
}
