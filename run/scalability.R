source("./simulators/scalability/scalable_sim.R")
source("./util.R")
library(Triangulation)

seed         <- 42
ncores       <- 16
ncores_cside <- 16
run_method  <- c("ctsv", "stance", "celina", "spvc", "cside", "mmm")
save_dir     <- "./data/detection_results/scalability/"

cside_gene_cutoff_reg <- 5e-5
cside_gene_cutoff     <- 5e-5
cside_fc_cutoff_reg   <- 0.5  
mmm_gene_cutoff_reg   <- 5e-5
mmm_gene_cutoff       <- 5e-5
mmm_fc_cutoff_reg     <- 0.5

boundary <- matrix(c(0, 0, 1, 0, 1, 1, 0, 1), ncol = 2, byrow = TRUE)
Tr.cell  <- TriMesh(boundary, n = 2)

# ── Experiment 1: Variable genes (n_side fixed at 30) ────────────────────────
n_side_fixed <- 30
n_genes_grid <- c(100, 500, 1000, 2000, 5000, 10000, 15000, 20000)

for (n_gene in n_genes_grid) {
  print(paste0("VarGene | Now running n_gene: ", n_gene, " | n_side: ", n_side_fixed))
  temp_sim <- scalable_sim$new(dispersion = 0.7, seed = seed,
                               n_cell_types = 3, n_side = n_side_fixed, n_gene = n_gene)
  run_all_tests(
    sp_count     = temp_sim$spot_counts,
    sp_comp      = temp_sim$spot_composition,
    sp_coords    = temp_sim$spot_coords,
    sc_count     = temp_sim$cell_counts,
    sc_metadata  = temp_sim$cell_metadata,
    spvc_tri     = Tr.cell,
    sim_name             = paste0("Cell: 3 | n_side: ", n_side_fixed, " | n_gene: ", n_gene),
    save_dir             = save_dir,
    sim_obj              = temp_sim,
    ncores               = ncores,
    ncores_cside         = ncores_cside,
    run_method           = run_method,
    timeout_sec          = 48 * 3600,
    cside.fc_cutoff_reg   = cside_fc_cutoff_reg,
    cside.gene_cutoff_reg = cside_gene_cutoff_reg,
    cside.gene_cutoff     = cside_gene_cutoff,
    mmm.fc_cutoff_reg     = mmm_fc_cutoff_reg,
    mmm.gene_cutoff_reg   = mmm_gene_cutoff_reg,
    mmm.gene_cutoff       = mmm_gene_cutoff
  )
}

# ── Experiment 2: Variable spots (n_gene fixed at 100) ───────────────────────
n_gene_fixed <- 100
n_sides_grid <- seq(10, 200, by = 10)

for (n_side in n_sides_grid) {
  print(paste0("VarSpot | Now running n_side: ", n_side, " | n_gene: ", n_gene_fixed))
  temp_sim <- scalable_sim$new(dispersion = 0.7, seed = seed,
                               n_cell_types = 3, n_side = n_side, n_gene = n_gene_fixed)
  run_all_tests(
    sp_count     = temp_sim$spot_counts,
    sp_comp      = temp_sim$spot_composition,
    sp_coords    = temp_sim$spot_coords,
    sc_count     = temp_sim$cell_counts,
    sc_metadata  = temp_sim$cell_metadata,
    spvc_tri     = Tr.cell,
    sim_name             = paste0("Cell: 3 | n_side: ", n_side, " | n_gene: ", n_gene_fixed),
    save_dir             = save_dir,
    sim_obj              = temp_sim,
    ncores               = ncores,
    ncores_cside         = ncores_cside,
    run_method           = run_method,
    timeout_sec          = 48 * 3600,
    cside.fc_cutoff_reg   = cside_fc_cutoff_reg,
    cside.gene_cutoff_reg = cside_gene_cutoff_reg,
    cside.gene_cutoff     = cside_gene_cutoff,
    mmm.fc_cutoff_reg     = mmm_fc_cutoff_reg,
    mmm.gene_cutoff_reg   = mmm_gene_cutoff_reg,
    mmm.gene_cutoff       = mmm_gene_cutoff
  )
}
