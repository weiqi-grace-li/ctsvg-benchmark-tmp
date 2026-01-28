load("./real_data/lung/processed/lung_filtered_v3.RData")
load("./real_data/lung/processed/Tr.cell.RData")
# run the methods 
source("./scDesign_simulators/simulation_util.R")
run_all_tests(
  sp_count = as.matrix(lung_filtered$spot_counts),
  sp_comp = as.matrix(lung_filtered$spot_composition),
  sp_coords = as.matrix(lung_filtered$spot_coords),
  
  # feed newly generated cell counts, not original cells 
  sc_count = as.matrix(lung_filtered$cell_counts),
  sc_metadata = lung_filtered$cell_metadata,
  spvc_tri = Tr.cell,
  sim_name = "lung_filtered_0.1",
  save_dir = "./real_data/lung/results/",
  time_save_path = "./real_data/lung/results/runtime.csv",
  sim_obj = lung_filtered,
  ncores = 32,
  ncores_cside = 30, 
  run_method = c("celina", "cside")
)



