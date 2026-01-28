library(ggplot2)
library(arrow)
library(Seurat)
library(Matrix)
library(readr)
library(Triangulation)
source("./scDesign_simulators/simulation_util.R")
source("./utility/simulator_util.R")

ncores = 48
ncores_small = 20
ncores_cside = 20

set.seed(253)
seeds = c(42, sample(1:1000, 5, replace = FALSE)) # 42, 734, 60, 67, 798
seeds = seeds[1:5]
sim_setup = read.csv("./utility/simulator_setup.csv")

boundary = matrix(c(0, 0, 1, 0, 1, 1, 0, 1), ncol = 2, byrow = TRUE)

for (seed in c(42, 60, 798)){
  
  if (seed == 42) {
    sim_number = c(20)
  } else if (seed == 60) {
    sim_number = c(16:20)
  } else if (seed == 798) {
    sim_number = c(8:20)
  }
  
  for (i in sim_number){
    temp_sim <- run_simulator(sim_setup$sim_name[i], seed = seed, phi = sim_setup$phi[i], scene = sim_setup$scene[i], pattern = sim_setup$pattern[i], 
                              control_UMI = FALSE)
    # for (rotate in rotations){
    sim_name = paste0(
      sim_setup$sim_name[i], 
      ifelse(((!is.na(sim_setup$phi[i]))& (sim_setup$phi[i]!="")), paste0("_",sim_setup$phi[i]), ""),
      ifelse(((!is.na(sim_setup$scene[i]))& (sim_setup$scene[i]!="")), paste0("_",sim_setup$scene[i]), ""),
      ifelse(((!is.na(sim_setup$pattern[i]))& (sim_setup$pattern[i]!="")), paste0("_",sim_setup$pattern[i]), "")
      # ,"_rotate", rotate
    )
    
    print(paste0("Now running replicate: ", which(seeds == seed), " | sim_type: ", sim_name))
    
    # pos.rotated <- rotate_points(temp_sim$spot_coords, angle_degrees = rotate)
    # boundary.rotated = rotate_points(boundary, angle_degrees = rotate)
    # Tr.cell.rotated = TriMesh(boundary.rotated, n = 2)
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
      save_dir = "./Results/pure_simulation/",
      time_save_path = "./Results/pure_simulation/runtime.csv",
      sim_obj = temp_sim,
      ncores = ncores, 
      ncores_cside =  ncores_cside,
      run_method = c("stance")
    )    
  }
}

### --------------------------------------------------------------------- ###
### analyze results  ###
### --------------------------------------------------------------------- ###
analyze_result(
  input_dir = "./Results/pure_simulation/",
  save_path = "./Results/pure_simulation/pure_sim_summary.xlsx",
  save_sheet = "aggregate",
  threshold = c(0.05)
)
