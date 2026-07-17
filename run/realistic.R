library(ggplot2)
library(arrow)
library(Seurat)
library(Matrix)
library(readr)
library(Triangulation)
source("./util.R")
source("./simulators/realistic/pseudo_spot_simulator.R")

ncores = 96
set.seed(253)
seeds = c(42, sample(1:1000, 5, replace = FALSE)) # 42, 734, 60, 67, 798, 419
seeds = seeds[1:5]
data_names = c("breast", "ovarian", "lymph")

for (seed in seeds){
  for (data_name in data_names){
    for (rotate in c(0, 30, 60, 90)) {
      
      message(paste0("Running ", data_name, " | replicate ", which(seeds == seed), " | rotate = ", rotate))
      processed_data_new = data_new_fast_loader(data_name = data_name, seed = seed, rotate = rotate)
      
      sim_name = paste0(data_name, "_rotate_", rotate)
      save_dir  = paste0("./data/detection_results/realistic/", data_name, "/")
      
      str(processed_data_new$sim_obj)
      
      run_all_tests(
        sp_count  = processed_data_new$spot_count,
        sp_comp   = processed_data_new$spot_composition,
        sp_coords = processed_data_new$spot_coords_rotated,
        sc_count  = processed_data_new$cell_counts_predrop, # use predrop as reference 
        sc_metadata = processed_data_new$cell_metadata,
        spvc_tri  = processed_data_new$Tr_cell_rotated,
        sim_name  = sim_name,
        save_dir  = save_dir,
        sim_obj   = processed_data_new$sim_obj,
        ncores    = ncores,
        run_method = c("stance")
      )
    }
  }
}