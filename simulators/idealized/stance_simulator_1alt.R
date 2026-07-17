library(spatstat)
library(STANCE)
library(dplyr)
library(purrr)
library(ggplot2)
library(ggforce)
library(patchwork)

# This simulation follows the simulation description in STANCE

stance_simulator_1alt <- setRefClass("stance_simulator_1alt",
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
                                   initialize = function(dispersion = 0.7, seed = 42) {
                                     .self$dispersion <- dispersion
                                     .self$seed <- seed
                                     .self$simulate()
                                   },
                                   
                                   simulate = function() {
                                     set.seed(seed)
                                     
                                     # Generate coordinates and domains
                                     pp <- rpoispp(lambda = 10000, win = owin(c(0, 1), c(0, 1)))
                                     coords <- data.frame(x = pp$x, y = pp$y, domain = "D1")
                                     center_ids <- sample(1:nrow(coords), 2)
                                     centers <- coords[center_ids, ]
                                     radii <- runif(2, 0.1, 0.5)
                                     
                                     dist1 <- sqrt((coords$x - centers$x[1])^2 + (coords$y - centers$y[1])^2)
                                     dist2 <- sqrt((coords$x - centers$x[2])^2 + (coords$y - centers$y[2])^2)
                                     
                                     # coords$domain <- ifelse(dist1 <= radii[1] & (dist1 <= dist2 | radii[1] < radii[2]), "D2", coords$domain)
                                     # coords$domain <- ifelse(dist2 <= radii[2] & (dist2 < dist1 | radii[2] < radii[1]), "D3", coords$domain)
                                     coords$domain <- ifelse(dist1 <= radii[1], "D2", coords$domain)
                                     coords$domain <- ifelse(((coords$domain == "D1" & dist2 <= radii[2]) 
                                                             | (coords$domain == "D2" & dist2 <= radii[2] & radii[2] < radii[1])), "D3", coords$domain)
                                     
                                     # Assign domain-specific cell types
                                     coords$type <- NA
                                     coords$type[coords$domain == "D1"] <- sample(c(1, 2), sum(coords$domain == "D1"), replace = TRUE, prob = c(0.5, 0.5))
                                     coords$type[coords$domain == "D2"] <- sample(c(2, 3), sum(coords$domain == "D2"), replace = TRUE, prob = c(0.5, 0.5))
                                     coords$type[coords$domain == "D3"] <- sample(c(1, 3), sum(coords$domain == "D3"), replace = TRUE, prob = c(0.5, 0.5))
                                     
                                     coords$type = paste0("Type", coords$type)
                                     
                                     # Simulate gene expression
                                     G <- 1000
                                     mu <- matrix(1, nrow = nrow(coords), ncol = G)
                                     special_genes = list()
                                     for (k in 1:3) {
                                       base <- (k - 1) * 200
                                       fcs <- c(4, 2, 0.5, 0.25)
                                       for (i in 0:3) {
                                         idx <- base + (1:50) + 50 * i
                                         mu[coords$type == paste0("Type", k), idx] <- fcs[i + 1]
                                       }
                                       special_genes[[paste0("Type", k)]] = list(marker = paste0("Gene", (base+1):(base+200)), ct_svg = "")
                                     }
                                     
                                     
                                     counts <- matrix(rnbinom(length(mu), mu = mu, size = dispersion), nrow = nrow(mu))
                                     rownames(mu) <- rownames(counts) <- paste0("Cell", 1:nrow(counts))
                                     colnames(mu) <- colnames(counts) <- paste0("Gene", 1:G)
                                     
                                     bin_size <- 1 / 32
                                     coords$ix <- cut(coords$x, breaks = seq(0, 1, by = bin_size), labels = FALSE)
                                     coords$iy <- cut(coords$y, breaks = seq(0, 1, by = bin_size), labels = FALSE)
                                     coords$spot <- paste0("Spot", interaction(coords$ix, coords$iy))
                                     
                                     gene_names <- paste0("Gene", 1:G)
                                     df <- cbind(coords, counts)
                                     
                                     df_spot <- df %>%
                                       group_by(spot) %>%
                                       summarise(
                                         x = mean(x), y = mean(y),
                                         counts = list(colSums(across(starts_with("Gene")))),
                                         types = list(table(type)),
                                         .groups = "drop"
                                       )
                                     
                                     counts_mat <- do.call(rbind, lapply(df_spot$counts, function(v) {
                                       v[setdiff(gene_names, names(v))] <- 0
                                       v[gene_names]
                                     }))
                                     rownames(counts_mat) <- df_spot$spot
                                     counts_mat <- t(counts_mat)
                                     
                                     .self$cell_counts <- t(counts)
                                     .self$cell_metadata <- coords[, c("x", "y", "type", "spot", "domain")]
                                     rownames(.self$cell_metadata) <- paste0("Cell", 1:nrow(coords))
                                     .self$spot_counts <- counts_mat
                                     .self$spot_coords <- as.matrix(df_spot[, c("x", "y")])
                                     rownames(.self$spot_coords) <- df_spot$spot
                                     
                                     comp_mat <- matrix(0, nrow = nrow(df_spot), ncol = 3)
                                     colnames(comp_mat) <- paste0("Type", 1:3)
                                     rownames(comp_mat) <- df_spot$spot
                                     for (i in 1:nrow(df_spot)) {
                                       type_counts <- df_spot$types[[i]]
                                       comp_mat[i, names(type_counts)] <- type_counts / sum(type_counts)
                                     }
                                     .self$spot_composition <- comp_mat
                                     
                                     .self$special_genes <- special_genes
                                     .self$gene_names <- gene_names
                                     .self$center <- centers
                                     .self$radius <- radii
                                     .self$mu_matrix <- t(mu)
                                     .self$cell_type_proportion = as.vector(prop.table(table(.self$cell_metadata$type)[paste0("Type", 1:3)]))
                                     names(.self$cell_type_proportion) = paste0("Type", 1:3)
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
                                       type = rep(types, 2), 
                                       x0 = c(rep(temp_sim$center[1, "x"], length(types)), rep(temp_sim$center[2, "x"], length(types))), 
                                       y0 = c(rep(temp_sim$center[1, "y"], length(types)), rep(temp_sim$center[2, "y"], length(types))),
                                       radius = c(rep(temp_sim$radius[1], length(types)), rep(temp_sim$radius[2], length(types)))
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
