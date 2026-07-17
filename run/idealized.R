library(ggplot2)
library(arrow)
library(Seurat)
library(Matrix)
library(readr)
library(Triangulation)
source("./util.R")
source("./simulators/idealized/celina_simulator_alt.R")
source("./simulators/idealized/celina_simulator_null.R")
source("./simulators/idealized/stance_simulator_alt.R")
source("./simulators/idealized/stance_simulator_1alt.R")

ncores = 1

seeds = c(42) # one seed for demonstration 

# run_method = c("celina", "stance", "spvc", "ctsv", "cside", "mmm") # uncomment to run all methods 
run_method = c("celina")

sim_setup = read.csv("./simulators/idealized/simulator_setup.csv")

boundary = matrix(c(0, 0, 1, 0, 1, 1, 0, 1), ncol = 2, byrow = TRUE)

for (seed in seeds){
  for (i in 1:nrow(sim_setup)){
    temp_sim <- run_simulator(sim_setup$sim_name[i], seed = seed, phi = sim_setup$phi[i], scene = sim_setup$scene[i], pattern = sim_setup$pattern[i], 
                              control_UMI = FALSE)
      sim_name = paste0(
        sim_setup$sim_name[i], 
        ifelse(((!is.na(sim_setup$phi[i]))& (sim_setup$phi[i]!="")), paste0("_",sim_setup$phi[i]), ""),
        ifelse(((!is.na(sim_setup$scene[i]))& (sim_setup$scene[i]!="")), paste0("_",sim_setup$scene[i]), ""),
        ifelse(((!is.na(sim_setup$pattern[i]))& (sim_setup$pattern[i]!="")), paste0("_",sim_setup$pattern[i]), "")
      )
      
      print(paste0("Now running replicate: ", which(seeds == seed), " | sim_type: ", sim_name))
      
      Tr.cell = TriMesh(boundary, n = 2)
      
      # --- run the tests 
      run_all_tests(
        sp_count = temp_sim$spot_counts,
        sp_comp = temp_sim$spot_composition,
        sp_coords = temp_sim$spot_coords,
        
        sc_count = temp_sim$cell_counts,
        sc_metadata = temp_sim$cell_metadata,
        
        spvc_tri = Tr.cell,
        sim_name = sim_name,
        save_dir = "./data/detection_results/idealized/",
        sim_obj = temp_sim,
        ncores = ncores, 
        run_method = run_method
      )    
  }
}

