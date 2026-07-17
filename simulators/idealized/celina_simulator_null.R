library(spatstat)
library(dplyr)
library(ggplot2)
library(ggforce)
library(patchwork)

celina_simulator_null <- setRefClass("celina_simulator_null",
                                fields = list(
                                  dispersion = "numeric",
                                  seed = "numeric",
                                  cell_type_proportion = "numeric",
                                  scenario = "character",
                                  cell_counts = "matrix",
                                  cell_metadata = "data.frame",
                                  spot_counts = "matrix",
                                  spot_coords = "matrix",
                                  spot_composition = "matrix",
                                  gene_names = "character",
                                  mu_matrix = "matrix",
                                  special_genes = "list", 
                                  resolution = "character"
                                ),
                                
                                methods = list(
                                  #' @description Initialization method. Sets simulation parameters and runs simulation.
                                  #' @param dispersion Dispersion parameter for NB sampling.
                                  #' @param seed Random seed.
                                  #' @param scenario Simulation scenario label ("I", "II", "III", or "IV").
                                  #' @param resolution If cell, 5000 cells, if spot, 20000 cells. 
                                  initialize = function(dispersion = 0.95, seed = 42, scenario = "I", resolution = "cell") {
                                    .self$dispersion <- dispersion
                                    .self$seed <- seed
                                    .self$scenario <- scenario
                                    .self$resolution <- resolution 
                                    .self$simulate()
                                  },
                                  
                                  #' @description Run full simulation.
                                  simulate = function() {
                                    set.seed(.self$seed)
                                    n_genes <- 1000
                                    n_marker_genes <- 200
                                    cell_types <- 1:4
                                    n_cells <- ifelse(.self$resolution == "cell", 5000, 20000)
                                    fold_changes <- list(
                                      Type1 = list(up = 4, down = 0.25),
                                      Type2 = list(up = 2, down = 0.5),
                                      Type3 = list(up = 2, down = 0.5),
                                      Type4 = list(up = 4, down = 0.25)
                                    )
                                    
                                    win <- owin(c(0, 1), c(0, 1))
                                    ppp_obj <- rpoispp(lambda = n_cells, win = win)
                                    x_coords <- ppp_obj$x
                                    y_coords <- ppp_obj$y
                                    n_cells <- ppp_obj$n
                                    
                                    assign_spatial_domain <- function(x) {
                                      if (x < 0.25) return(1)
                                      else if (x < 0.5) return(2)
                                      else if (x < 0.75) return(3)
                                      else return(4)
                                    }
                                    spatial_domain <- sapply(x_coords, assign_spatial_domain)
                                    
                                    scenario_compositions <- list(
                                      I = list(c("1"), c("2"), c("3"), c("4")),
                                      II = list(c("1", "2", "3"), c("2", "3", "4"), c("3", "4", "1"), c("4", "1", "2")),
                                      III = list(c("1", "2", "3"), c("2", "3", "4"), c("3", "4", "1"), c("4", "1", "2")),
                                      IV = list(c("1", "2", "3"), c("2", "3", "4"), c("3", "4", "1"), c("4", "1", "2"))
                                    )
                                    get_composition_probs <- function(scenario, domain_idx) {
                                      if (scenario == "I") return(c(1))
                                      if (scenario == "II") return(c(0.9, 0.05, 0.05))
                                      if (scenario == "III") return(c(0.5, 0.25, 0.25))
                                      if (scenario == "IV") return(rep(1/3, 3))
                                    }
                                    assign_cell_types <- function() {
                                      ct_vec <- character(n_cells)
                                      for (i in seq_len(n_cells)) {
                                        d <- spatial_domain[i]
                                        comp <- scenario_compositions[[.self$scenario]][[d]]
                                        probs <- get_composition_probs(.self$scenario, d)
                                        ct_vec[i] <- sample(comp, 1, prob = probs)
                                      }
                                      as.integer(ct_vec)
                                    }
                                    cell_type_vec <- assign_cell_types()
                                    
                                    mu <- matrix(1, nrow = n_cells, ncol = n_genes)
                                    markers <- sample(1:n_genes, n_marker_genes)
                                    special_gene_list <- list()
                                    for (ct in 1:4) {
                                      ct_idx <- which(cell_type_vec == ct)
                                      gene_ids <- markers[((ct - 1) * 50 + 1):(ct * 50)]
                                      fc_vals <- rep(c(fold_changes[[ct]]$up, fold_changes[[ct]]$down), each = 25)
                                      for (j in 1:50) {
                                        mu[ct_idx, gene_ids[j]] <- fc_vals[j]
                                      }
                                      special_gene_list[[paste0("Type",ct)]] = list(marker = paste0("Gene", gene_ids), ct_svg = "")
                                    }
                                    
                                    counts <- matrix(rnbinom(length(mu), mu = mu, size = .self$dispersion), nrow = nrow(mu))
                                    mu_matrix <- mu
                                    
                                    coords_df <- data.frame(x = x_coords, y = y_coords, type = paste0("Type", cell_type_vec), domain = spatial_domain)
                                    colnames(mu_matrix) <- colnames(counts) <- paste0("Gene", 1:n_genes)
                                    rownames(coords_df) <- rownames(mu_matrix) <- rownames(counts) <- paste0("Cell", 1:n_cells)
                                    
                                    bin_size_x <- 1 / 70
                                    bin_size_y <- 1 / 71
                                    coords_df$ix <- cut(coords_df$x, breaks = seq(0, 1, by = bin_size_x), labels = FALSE)
                                    coords_df$iy <- cut(coords_df$y, breaks = seq(0, 1, by = bin_size_y), labels = FALSE)
                                    coords_df$spot <- paste0("Spot", interaction(coords_df$ix, coords_df$iy))
                                    
                                    df <- cbind(coords_df, as.data.frame(counts))
                                    df_spot <- df %>%
                                      group_by(spot) %>%
                                      summarise(
                                        x = mean(x), y = mean(y), # CELINA didn't specify
                                        counts = list(colSums(across(starts_with("Gene")))),
                                        types = list(table(type)), .groups = "drop"
                                      )
                                    
                                    gene_names <- paste0("Gene", 1:n_genes)
                                    count_mat <- do.call(rbind, lapply(df_spot$counts, function(v) {
                                      v[setdiff(gene_names, names(v))] <- 0
                                      v[gene_names]
                                    }))
                                    rownames(count_mat) <- df_spot$spot
                                    count_mat <- t(count_mat)
                                    
                                    comp_mat <- matrix(0, nrow = nrow(df_spot), ncol = 4)
                                    colnames(comp_mat) <- paste0("Type", 1:4)
                                    rownames(comp_mat) <- df_spot$spot
                                    for (i in 1:nrow(df_spot)) {
                                      type_counts <- df_spot$types[[i]]
                                      comp_mat[i, names(type_counts)] <- type_counts / sum(type_counts)
                                    }
                                    .self$spot_composition <- comp_mat
                                    
                                    pos <- df_spot %>% dplyr::select(x, y)
                                    rownames(comp_mat) <- rownames(pos) <- df_spot$spot
                                    
                                    .self$cell_counts <- t(counts)
                                    .self$cell_metadata <- coords_df[, c("x", "y", "type", "spot", "domain")]
                                    .self$spot_counts <- count_mat
                                    .self$spot_coords <- as.matrix(pos)
                                    .self$spot_composition <- comp_mat
                                    .self$gene_names <- gene_names
                                    .self$mu_matrix <- t(mu_matrix)
                                    .self$special_genes <- special_gene_list
                                    .self$cell_type_proportion = as.vector(prop.table(table(.self$cell_metadata$type)[paste0("Type", 1:4)]))
                                    names(.self$cell_type_proportion) = paste0("Type", 1:4)
                                  },
                                  
                                  sim_plot = function(gene_name, desc = "") {
                                    by_cell_mu = cbind(.self$cell_metadata, expr = .self$mu_matrix[gene_name,])
                                    all_mu = cbind(.self$cell_metadata, expr = .self$mu_matrix[gene_name,])
                                    all_mu$type = "all"
                                    plot_mu = rbind(by_cell_mu, all_mu)
                                    
                                    by_cell_raw = cbind(.self$cell_metadata, expr = .self$cell_counts[gene_name,])
                                    all_raw = cbind(.self$cell_metadata, expr = .self$cell_counts[gene_name,])
                                    all_raw$type = "all"
                                    plot_raw = rbind(by_cell_raw, all_raw)
                                    
                                    types = unique(plot_mu$type)
                                    mu = ggplot(plot_mu, aes(x = x, y= y, color = expr)) + 
                                      geom_point(size = 0.6) + 
                                      geom_vline(xintercept = c(0.25, 0.5, 0.75), linetype = "dashed", color = "gray80") + 
                                      scale_color_gradient(low = "lightyellow", high = "darkgreen", na.value = "gray80") + 
                                      facet_wrap(~type, nrow = 1) + 
                                      coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) + 
                                      labs(title = "mu", color = "Expression")
                                    raw = ggplot(plot_raw, aes(x = x, y= y, color = expr)) + 
                                      geom_point(size = 0.6) + 
                                      geom_vline(xintercept = c(0.25, 0.5, 0.75), linetype = "dashed", color = "gray80") + 
                                      scale_color_gradient(low = "lightyellow", high = "darkgreen", na.value = "gray80") + 
                                      facet_wrap(~type, nrow = 1) + 
                                      coord_fixed(xlim = c(0, 1), ylim = c(0, 1))+ 
                                      labs(title = "raw expression", color = "Expression")
                                    mu / raw + plot_annotation(title = paste(gene_name, desc))
                                  }
                                )
)