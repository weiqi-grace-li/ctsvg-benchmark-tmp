# Load required libraries
library(spatstat)
library(STANCE)
library(dplyr)
library(purrr)
library(ggplot2)
library(ggforce)
library(patchwork)



stance_simulator_alt <- setRefClass("stance_simulator_alt",
                                    fields = list(
                                      dispersion = "numeric",
                                      seed = "numeric",
                                      cell_type_proportion = "numeric",
                                      cell_counts = "matrix",
                                      cell_metadata = "data.frame",
                                      spot_counts = "matrix",
                                      spot_coords = "matrix",
                                      spot_composition = "matrix",
                                      special_genes = "list",
                                      gene_names = "character",
                                      center = "data.frame",
                                      radius = "numeric",
                                      mu_matrix = "matrix"
                                    ),
                                    
                                    methods = list(
                                      initialize = function(dispersion = 0.7, seed = 42, cell_type_proportion = c(0.1, 0.3, 0.6)) {
                                        .self$dispersion <- dispersion
                                        .self$seed <- seed
                                        .self$cell_type_proportion <- cell_type_proportion
                                        names(.self$cell_type_proportion) = paste0("Type", 1:3)
                                        .self$simulate()
                                      },
                                      
                                      simulate = function() {
                                        set.seed(seed)
                                        
                                        # Step 1: Cell coordinates
                                        pp <- rpoispp(lambda = 4000, win = owin(c(0, 1), c(0, 1)))
                                        coords <- data.frame(x = pp$x, y = pp$y)
                                        
                                        # Step 2: Define domain center and radius
                                        center_idx <- sample(1:nrow(coords), 1)
                                        center <- coords[center_idx, , drop = FALSE]
                                        radius <- runif(1, 0.2, 0.4)
                                        dists <- sqrt((coords$x - center$x)^2 + (coords$y - center$y)^2)
                                        coords$domain <- ifelse(dists <= radius, "SD", "D")
                                        
                                        # Step 3: Cell types
                                        coords$type <- sample(1:3, size = nrow(coords), replace = TRUE, prob = cell_type_proportion)
                                        coords$type <- paste0("Type", coords$type)
                                        
                                        # Step 4: Gene expression
                                        G <- 1000
                                        mu <- matrix(1, nrow = nrow(coords), ncol = G)
                                        
                                        marker_list <- list()
                                        ct_svg_list <- list()
                                        
                                        for (k in 1:3) {
                                          marker_idx <- ((k - 1) * 100 + 1):(k * 100)
                                          mu[coords$type == paste0("Type", k), marker_idx] <- 2
                                          marker_list[[k]] <- marker_idx
                                        }
                                        
                                        for (k in 1:3) {
                                          base_idx <- 300 + (k - 1) * 200
                                          type_k <- coords$type == paste0("Type", k)
                                          sd_k <- coords$domain == "SD"
                                          svg_idx <- c()
                                          for (i in 0:3) {
                                            gene_idx <- base_idx + (1:50) + i * 50
                                            fc <- c(4, 2, 0.5, 0.25)[i + 1] # revised back to 4, 2, 0.5, 0.25 Dec 8, 2025.  Was 40, xxx for Figure 1
                                            mu[type_k & sd_k, gene_idx] <- fc
                                            svg_idx <- c(svg_idx, gene_idx)
                                          }
                                          ct_svg_list[[k]] <- svg_idx
                                        }
                                        
                                        counts <- matrix(rnbinom(length(mu), mu = mu, size = dispersion), nrow = nrow(mu))
                                        colnames(mu) = colnames(counts) <- paste0("Gene", 1:G)
                                        rownames(mu) = rownames(counts) <- paste0("Cell", 1:nrow(counts))
                                        
                                        
                                        # Step 5: Assign to spots
                                        bin_size <- 0.05
                                        coords$ix <- cut(coords$x, breaks = seq(0, 1, by = bin_size), labels = FALSE)
                                        coords$iy <- cut(coords$y, breaks = seq(0, 1, by = bin_size), labels = FALSE)
                                        coords$spot <- paste0("Spot", interaction(coords$ix, coords$iy))
                                        
                                        df <- cbind(coords, counts)
                                        df_spot <- df %>%
                                          group_by(spot) %>%
                                          summarise(
                                            x = mean(x), y = mean(y),
                                            counts = list(colSums(across(starts_with("Gene")))),
                                            types = list(table(type)), .groups = "drop"
                                          ) %>% ungroup()
                                        
                                        count_mat <- do.call(rbind, df_spot$counts)
                                        rownames(count_mat) <- df_spot$spot
                                        colnames(count_mat) <- paste0("Gene", 1:G)
                                        
                                        # Build cell and spot metadata
                                        .self$cell_counts <- t(counts)
                                        .self$cell_metadata <- coords[, c("x", "y", "type", "spot", "domain")]
                                        rownames(.self$cell_metadata) <- paste0("Cell", 1:nrow(.self$cell_metadata))
                                        
                                        .self$spot_counts <- t(count_mat)
                                        .self$spot_coords <- as.matrix(df_spot[, c("x", "y")])
                                        rownames(.self$spot_coords) <- df_spot$spot
                                        
                                        # Spot composition
                                        comp_mat <- matrix(0, nrow = nrow(df_spot), ncol = 3)
                                        colnames(comp_mat) <- paste0("Type", 1:3)
                                        rownames(comp_mat) <- df_spot$spot
                                        for (i in 1:nrow(df_spot)) {
                                          type_counts <- df_spot$types[[i]]
                                          comp_mat[i, names(type_counts)] <- type_counts / sum(type_counts)
                                        }
                                        .self$spot_composition <- comp_mat
                                        
                                        # Gene lists
                                        .self$special_genes <- list()
                                        for (k in 1:3) {
                                          .self$special_genes[[paste0("Type", k)]] <- list(marker = paste0("Gene", marker_list[[k]]), ct_svg = paste0("Gene", ct_svg_list[[k]]))
                                        }
                                        .self$gene_names <- paste0("Gene", 1:G)
                                        .self$center <- center
                                        .self$radius <- radius
                                        # .self$mu_matrix <- compute_mu_per_celltype(mu, coords$type)
                                        .self$mu_matrix <- t(mu)                                   
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
                                        circle_df = data.frame(
                                          type = types, 
                                          x0 = rep(.self$center[1, "x"], length(types)), 
                                          y0 = rep(.self$center[1, "y"], length(types)),
                                          radius = rep(.self$radius, length(types))
                                        )
                                        mu = ggplot(plot_mu, aes(x = x, y= y, color = expr)) + 
                                          geom_point(size = 0.6) + 
                                          scale_color_gradient(low = "lightyellow", high = "darkgreen", na.value = "gray80") + 
                                          geom_circle(data = circle_df, aes(x0=x0, y0=y0, r = radius), color = "red", size = 0.8, 
                                                      inherit.aes = FALSE)+
                                          facet_wrap(~type, nrow = 1) + 
                                          coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) + 
                                          labs(title = "mu", color = "Expression")
                                        raw = ggplot(plot_raw, aes(x = x, y= y, color = expr)) + 
                                          geom_point(size = 0.6) + 
                                          scale_color_gradient(low = "lightyellow", high = "darkgreen", na.value = "gray80") + 
                                          geom_circle(data = circle_df, aes(x0=x0, y0=y0, r = radius), color = "red", size = 0.8, 
                                                      inherit.aes = FALSE)+
                                          facet_wrap(~type, nrow = 1) + 
                                          coord_fixed(xlim = c(0, 1), ylim = c(0, 1))+ 
                                          labs(title = "raw expression", color = "Expression")
                                        mu / raw + plot_annotation(title = paste(gene_name, desc))
                                      }
                                    )
)