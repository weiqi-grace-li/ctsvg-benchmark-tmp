library(ggplot2)
library(arrow)
library(Seurat)
library(Matrix)
library(readr)
library(Triangulation)
source("./util.R")
source("./simulators/realistic/pseudo_spot_simulator.R")
source("./simulators/decomposed/scDesign_decompose.R")

ncores = 1
seeds = c(42) # one seed for demonstration 

# breast datasets are included in this repo 
# to run decomposed simulation on the other datasets 
# make sure you either process and package them or 
# download packaged data from zenodo 
data_names = c("breast")

# run_method <- c("celina", "stance", "spvc", "ctsv", "cside", "mmm") #uncomment to run all methods 
run_method = c("celina")

sim_types = list(
  c(1, 2, 3, 4, 5),
  c(1), c(3),
  c(1, 2),
  c(1, 2, 4),
  c(1, 2, 5),
  c(1, 2, 3),
  c(0), c(2), c(2, 4),
  c()
)

for (data_name in data_names) {
  for (sim in seq_along(sim_types)) {
    for (seed in seeds) {

      message(paste0("Running ", data_name,
                     " | replicate ", which(seeds == seed),
                     " | sim_type: ", paste(sim_types[[sim]], collapse = ",")))

      processed = load_process_decompose(
        data_name = data_name,
        seed      = seed,
        sim_type  = sim_types[[sim]]
      )

      sim_name = paste0(data_name, "_decompose_exp_", paste(sim_types[[sim]], collapse = ","))
      save_dir = paste0("./data/detection_results/decomposed/", data_name, "/")

      str(processed$sim_obj)

      run_all_tests(
        sp_count    = processed$spot_count,
        sp_comp     = processed$spot_composition,
        sp_coords   = processed$spot_coords,
        sc_count    = processed$cell_counts_predrop,
        sc_metadata = processed$cell_metadata,
        spvc_tri    = processed$Tr_cell,
        sim_name    = sim_name,
        save_dir    = save_dir,
        sim_obj     = processed$sim_obj,
        ncores      = ncores,
        run_method  = run_method
      )
    }
  }
}
