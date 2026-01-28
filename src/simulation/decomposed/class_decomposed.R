# decompose_simulator.R
# Assumes SimDataset is already defined (e.g., source("SimDataset.R"))

decompose_simulator <- setRefClass(
  "decompose_simulator",
  contains = "SimDataset",
  
  fields = list(
    sim_type = "ANY",
    center = "data.frame",
    radius = "numeric",
    cell_spot_prop = "matrix"
  ),
  
  methods = list(
    
    initialize = function(
    pseudo_simulator,
    data_ori,
    dispersion = 0.7,
    seed = 42,
    sim_type = c(),
    n_markers = 50,
    cell_type_proportion = c("Type1" = 0.1, "Type2" = 0.3, "Type3" = 0.6)
    ) {
      .self$seed <- seed
      .self$sim_type <- sim_type
      
      # ---- checks (same logic you had)
      if (!all(sim_type %in% c(-1, 0, 1, 2, 3, 4, 5))) {
        stop("Wrong decomposition type, only combinations of -1, 0, 1, 2, 3, 4, 5 are accepted!")
      }
      if (5 %in% sim_type) stopifnot(all(c(1, 2) %in% sim_type))
      if (4 %in% sim_type) stopifnot(2 %in% sim_type)
      
      # ---- cell type proportions
      if (1 %in% sim_type) {
        .self$cell_type_proportion <- pseudo_simulator$cell_type_proportion
      } else {
        .self$cell_type_proportion <- cell_type_proportion
      }
      
      # ---- gene names
      if (any(c(-1, 0, 1, 2, 3, 4, 5) %in% sim_type)) {
        .self$gene_names <- pseudo_simulator$gene_names
      } else {
        .self$gene_names <- paste0("Gene", 1:1000)
      }
      
      # ---- assumption helper for mode 3 if no mode 1
      assume_3 <- NULL
      if ((3 %in% sim_type) && !(1 %in% sim_type)) {
        assume_3 <- names(sort(pseudo_simulator$cell_type_proportion, decreasing = TRUE))[1:length(.self$cell_type_proportion)]
      }
      
      .self$simulate(
        sim_type = sim_type,
        pseudo_simulator = pseudo_simulator,
        data_ori = data_ori,
        n_markers = n_markers,
        assume_3 = assume_3,
        dispersion = dispersion
      )
      invisible(.self)
    },
    
    simulate = function(sim_type, pseudo_simulator, data_ori, n_markers, assume_3, dispersion) {
      set.seed(.self$seed)
      
      # ---- 1) cell locations
      if (2 %in% sim_type) {
        coords <- pseudo_simulator$cell_metadata[, c("x", "y")]
      } else if (any(c(0, 1, 3, 4, 5) %in% sim_type)) {
        pp <- spatstat.random::runifpoint(
          n = nrow(pseudo_simulator$cell_metadata),
          win = spatstat.geom::owin(
            c(min(pseudo_simulator$cell_metadata$x), max(pseudo_simulator$cell_metadata$x)),
            c(min(pseudo_simulator$cell_metadata$y), max(pseudo_simulator$cell_metadata$y))
          )
        )
        coords <- data.frame(x = pp$x, y = pp$y)
      } else {
        pp <- spatstat.random::rpoispp(lambda = 4000, win = spatstat.geom::owin(c(0, 1), c(0, 1)))
        coords <- data.frame(x = pp$x, y = pp$y)
      }
      
      # ---- 2) dispersion matrix (cell x gene here; later transpose like your original)
      .self$dispersion <- matrix(
        dispersion,
        nrow = nrow(coords),
        ncol = length(.self$gene_names),
        dimnames = list(NULL, .self$gene_names)
      )
      
      # ---- 3) domain center & radius
      center_idx <- sample.int(nrow(coords), 1)
      center <- coords[center_idx, , drop = FALSE]
      radius <- runif(1, 0.2, 0.4)
      
      if (any(c(0, 1, 2, 3, 4, 5) %in% sim_type)) {
        scale_xy <- c(max(coords$x) - min(coords$x), max(coords$y) - min(coords$y))
        axis <- which.min(scale_xy)
        radius <- radius * scale_xy[axis]
      }
      
      dists <- sqrt((coords$x - center$x)^2 + (coords$y - center$y)^2)
      coords$domain <- ifelse(dists <= radius, "SD", "D")
      
      # ---- 4) cell types
      if (all(c(1, 2) %in% sim_type)) {
        coords$type <- pseudo_simulator$cell_metadata$type
      } else {
        coords$type <- sample(
          names(.self$cell_type_proportion),
          size = nrow(coords),
          replace = TRUE,
          prob = .self$cell_type_proportion
        )
      }
      
      # ---- 5) special genes + mu
      special_genes <- vector("list", length(.self$cell_type_proportion))
      names(special_genes) <- names(.self$cell_type_proportion)
      
      mu <- matrix(1, nrow = nrow(coords), ncol = length(.self$gene_names))
      colnames(mu) <- .self$gene_names
      
      for (k in seq_along(special_genes)) {
        ct_name <- names(special_genes)[k]
        idx <- coords$type == ct_name
        sd_k <- coords$domain == "SD"
        
        # --- ctSVG set
        if (1 %in% sim_type) {
          if (length(pseudo_simulator$special_genes[[k]]$ct_svg) > 0) {
            special_genes[[k]]$ct_svg <- pseudo_simulator$special_genes[[k]]$ct_svg
          } else {
            special_genes[[k]]$ct_svg <- character(0)
          }
        } else if (any(c(-1, 0, 1, 2, 3, 4, 5) %in% sim_type)) {
          ct_svgs_all <- unlist(lapply(pseudo_simulator$special_genes, function(x) x$ct_svg))
          n_ctsvg <- floor(length(ct_svgs_all) / length(special_genes))
          if (n_ctsvg > 0) {
            special_genes[[k]]$ct_svg <- unique(ct_svgs_all[((k - 1) * n_ctsvg + 1):(k * n_ctsvg)])
          } else {
            special_genes[[k]]$ct_svg <- character(0)
          }
        } else {
          special_genes[[k]]$ct_svg <- paste0("Gene", (300 + (k - 1) * 200 + 1):(300 + k * 200))
        }
        
        # --- mu for ctSVG
        if (length(special_genes[[k]]$ct_svg) > 0) {
          if (5 %in% sim_type) {
            mu[idx, special_genes[[k]]$ct_svg] <-
              as.matrix(t(pseudo_simulator$mu_matrix[special_genes[[k]]$ct_svg, idx]), drop = FALSE)
            .self$dispersion[idx, special_genes[[k]]$ct_svg] <-
              as.matrix(t(pseudo_simulator$dispersion[special_genes[[k]]$ct_svg, idx]), drop = FALSE)
          } else if (any(c(-1, 0, 1, 2, 3, 4) %in% sim_type)) {
            mu[idx & sd_k, special_genes[[k]]$ct_svg] <- 4
          } else {
            base_idx <- 300 + (k - 1) * 200
            for (i in 0:3) {
              gene_idx <- base_idx + (1:50) + i * 50
              fc <- c(4, 2, 0.5, 0.25)[i + 1]
              mu[idx & sd_k, paste0("Gene", gene_idx)] <- fc
            }
          }
        }
        
        # --- marker genes
        if (3 %in% sim_type) {
          special_genes[[k]]$marker <- c("")
          if (1 %in% sim_type) {
            non_ctsvg <- which(!(.self$gene_names %in% special_genes[[k]]$ct_svg))
            mu[idx, non_ctsvg] <- as.matrix(t(pseudo_simulator$mu_matrix[non_ctsvg, idx]), drop = FALSE)
            .self$dispersion[idx, non_ctsvg] <- as.matrix(t(pseudo_simulator$dispersion[non_ctsvg, idx]), drop = FALSE)
          } else {
            replace_idx <- which(pseudo_simulator$cell_metadata$type == assume_3[k])
            replace_nonctsvg <- which(
              !(.self$gene_names %in% pseudo_simulator$special_genes[[assume_3[k]]]$ct_svg) &
                !(.self$gene_names %in% special_genes[[k]]$ct_svg)
            )
            mu[idx, replace_nonctsvg] <- matrix(
              rep(pseudo_simulator$mu_matrix[replace_nonctsvg, replace_idx[1]], sum(idx)),
              nrow = sum(idx),
              ncol = length(replace_nonctsvg),
              byrow = TRUE
            )
            .self$dispersion[idx, replace_nonctsvg] <- matrix(
              rep(pseudo_simulator$dispersion[replace_nonctsvg, replace_idx[1]], sum(idx)),
              nrow = sum(idx),
              ncol = length(replace_nonctsvg),
              byrow = TRUE
            )
          }
        } else if (any(c(-1, 0, 1, 2, 4, 5) %in% sim_type)) {
          special_genes[[k]]$marker <- sample(.self$gene_names, size = n_markers, replace = FALSE)
          marker_ctsvg <- intersect(special_genes[[k]]$marker, special_genes[[k]]$ct_svg)
          pure_marker <- setdiff(special_genes[[k]]$marker, marker_ctsvg)
          mu[idx, pure_marker] <- 2
          mu[idx & !sd_k, marker_ctsvg] <- 2
        } else {
          special_genes[[k]]$marker <- paste0("Gene", ((k - 1) * 100 + 1):(k * 100))
          mu[idx, special_genes[[k]]$marker] <- 2
        }
      }
      
      # ---- 6) counts with guard for mu==0 / NA dispersion
      counts <- matrix(0, nrow = nrow(mu), ncol = ncol(mu))
      nonzero_idx <- which(mu > 0 & !is.na(mu))
      counts[nonzero_idx] <- stats::rnbinom(
        n = length(nonzero_idx),
        mu = mu[nonzero_idx],
        size = .self$dispersion[nonzero_idx]
      )
      colnames(counts) <- colnames(mu)
      
      if (any(c(0, 1, 2, 3, 4, 5) %in% sim_type)) {
        rownames(mu) <- rownames(counts) <- rownames(.self$dispersion) <- rownames(pseudo_simulator$cell_metadata)
      } else {
        rownames(mu) <- rownames(counts) <- rownames(.self$dispersion) <- paste0("Cell", 1:nrow(counts))
      }
      
      # ---- 7) capture ratio adjustment (mode 4)
      if (4 %in% sim_type) {
        gene_cell_drop <- data_ori$gene_cell_drop[colnames(counts), rownames(counts)]
        counts <- counts * t(1 - gene_cell_drop)
      }
      
      # ---- 8) spots
      if (2 %in% sim_type) {
        count_mat <- round(as.matrix(t(counts) %*% data_ori$cell_spot_prop), 0)
        count_mat <- t(count_mat)
        colnames(count_mat) <- .self$gene_names
        coords$spot <- apply(data_ori$cell_spot_prop, 1, function(x) which.max(x))
        coords$spot <- paste0("spot", coords$spot)
        .self$spot_coords <- as.matrix(data_ori$pseudo_meta[, c("x", "y")])
        rownames(count_mat) <- rownames(.self$spot_coords) <- rownames(data_ori$pseudo_meta)
      } else {
        if (any(c(0, 1, 2, 3, 4, 5) %in% sim_type)) {
          n_bins <- floor(sqrt(nrow(pseudo_simulator$spot_coords)))
          bin_size_x <- floor((max(coords$x) - min(coords$x)) / n_bins)
          bin_size_y <- floor((max(coords$y) - min(coords$y)) / n_bins)
          coords$ix <- cut(coords$x, breaks = c(seq(min(coords$x), max(coords$x), by = bin_size_x), max(coords$x)), labels = FALSE)
          coords$iy <- cut(coords$y, breaks = c(seq(min(coords$y), max(coords$y), by = bin_size_y), max(coords$y)), labels = FALSE)
        } else {
          bin_size <- 0.05
          coords$ix <- cut(coords$x, breaks = seq(0, 1, by = bin_size), labels = FALSE)
          coords$iy <- cut(coords$y, breaks = seq(0, 1, by = bin_size), labels = FALSE)
        }
        
        coords$spot <- paste0("spot", interaction(coords$ix, coords$iy))
        
        df <- cbind(counts, coords)
        df_spot <- df %>%
          dplyr::group_by(spot) %>%
          dplyr::summarise(
            x = mean(x), y = mean(y),
            counts = list(colSums(dplyr::across(1:length(.self$gene_names)))),
            types = list(table(type)),
            .groups = "drop"
          ) %>% dplyr::ungroup()
        
        count_mat <- do.call(rbind, df_spot$counts)
        rownames(count_mat) <- df_spot$spot
        colnames(count_mat) <- .self$gene_names
        
        .self$spot_coords <- as.matrix(df_spot[, c("x", "y")])
        rownames(.self$spot_coords) <- df_spot$spot
      }
      
      # ---- 9) attach standardized outputs (SimDataset contract)
      .self$cell_metadata <- coords[, c("x", "y", "type", "spot", "domain")]
      if (any(c(0, 1, 2, 3, 4, 5) %in% sim_type)) {
        rownames(.self$cell_metadata) <- rownames(pseudo_simulator$cell_metadata)
      } else {
        rownames(.self$cell_metadata) <- paste0("Cell", 1:nrow(counts))
      }
      
      .self$cell_counts <- t(as.matrix(counts))
      .self$spot_counts <- t(count_mat)
      .self$special_genes <- special_genes
      .self$center <- center
      .self$radius <- radius
      .self$mu_matrix <- t(mu)
      .self$dispersion <- t(.self$dispersion)
      .self$cell_spot_prop <- as.matrix(data_ori$cell_spot_prop)
      
      # ---- 10) spot composition
      comp_mat <- matrix(0, nrow = nrow(.self$spot_coords), ncol = length(.self$cell_type_proportion))
      colnames(comp_mat) <- names(.self$cell_type_proportion)
      rownames(comp_mat) <- rownames(.self$spot_coords)
      
      if (all(c(1, 2) %in% sim_type)) {
        comp_mat <- pseudo_simulator$spot_composition
      } else if (2 %in% sim_type) {
        for (i in seq_len(nrow(.self$spot_coords))) {
          spot_idx <- which(coords$spot == rownames(.self$spot_coords)[i])
          type_counts <- table(coords$type[spot_idx])
          comp_mat[i, names(type_counts)] <- type_counts / sum(type_counts)
        }
      } else {
        for (i in seq_len(nrow(df_spot))) {
          type_counts <- df_spot$types[[i]]
          comp_mat[i, names(type_counts)] <- type_counts / sum(type_counts)
        }
      }
      .self$spot_composition <- comp_mat
      
      # ---- 11) drop all-zero genes/cells/spots (your original behavior)
      .self$cell_counts <- .self$cell_counts[.self$gene_names, , drop = FALSE]
      .self$spot_counts <- .self$spot_counts[.self$gene_names, , drop = FALSE]
      
      keep_genes <- which((rowSums(.self$cell_counts) > 0) & (rowSums(.self$spot_counts) > 0))
      keep_gene_names <- .self$gene_names[keep_genes]
      
      keep_cells <- which(colSums(.self$cell_counts) > 0)
      keep_cell_names <- colnames(.self$cell_counts)[keep_cells]
      
      keep_spots <- which(colSums(.self$spot_counts) > 0)
      keep_spot_names <- colnames(.self$spot_counts)[keep_spots]
      
      .self$gene_names <- keep_gene_names
      
      .self$cell_counts <- .self$cell_counts[keep_gene_names, keep_cell_names, drop = FALSE]
      .self$cell_metadata <- .self$cell_metadata[keep_cell_names, , drop = FALSE]
      .self$mu_matrix <- .self$mu_matrix[keep_gene_names, keep_cell_names, drop = FALSE]
      if (length(.self$dispersion) > 1) {
        .self$dispersion <- .self$dispersion[keep_gene_names, keep_cell_names, drop = FALSE]
      }
      
      .self$spot_counts <- .self$spot_counts[keep_gene_names, keep_spot_names, drop = FALSE]
      .self$spot_coords <- .self$spot_coords[keep_spot_names, , drop = FALSE]
      .self$spot_composition <- .self$spot_composition[keep_spot_names, , drop = FALSE]
      
      invisible(.self)
    },
    
    sim_plot = function(gene_name, desc = "") {
      # This is your existing sim_plot, unchanged in behavior (just wrapped inside method)
      ct_svg_types <- names(.self$special_genes)[
        vapply(.self$special_genes, function(x) gene_name %in% x$ct_svg, logical(1))
      ]
      marker_types <- names(.self$special_genes)[
        vapply(.self$special_genes, function(x) gene_name %in% x$marker, logical(1))
      ]
      
      by_cell_mu <- cbind(.self$cell_metadata, expr = .self$mu_matrix[gene_name, , drop = TRUE])
      all_mu <- cbind(.self$cell_metadata, expr = .self$mu_matrix[gene_name, , drop = TRUE])
      all_mu$type <- "all"
      plot_mu <- rbind(by_cell_mu, all_mu)
      
      by_cell_raw <- cbind(.self$cell_metadata, expr = .self$cell_counts[gene_name, , drop = TRUE])
      all_raw <- cbind(.self$cell_metadata, expr = .self$cell_counts[gene_name, , drop = TRUE])
      all_raw$type <- "all"
      plot_raw <- rbind(by_cell_raw, all_raw)
      
      types <- unique(plot_mu$type)
      
      mu_summary <- stats::aggregate(expr ~ type, data = plot_mu, FUN = mean)
      label_data <- merge(data.frame(type = types), mu_summary, by = "type", all.x = TRUE)
      
      label_map_mu <- setNames(
        paste0(
          label_data$type,
          ifelse(label_data$type %in% ct_svg_types, " · ct_svg", ""),
          ifelse(label_data$type %in% marker_types, " · marker", ""),
          " · mu = ",
          round(label_data$expr, 2)
        ),
        label_data$type
      )
      
      label_map_raw <- setNames(
        paste0(
          types,
          ifelse(types %in% ct_svg_types, " · ct_svg", ""),
          ifelse(types %in% marker_types, " · marker", "")
        ),
        types
      )
      
      circle_df <- data.frame(
        type = types,
        x0 = rep(.self$center[1, "x"], length(types)),
        y0 = rep(.self$center[1, "y"], length(types)),
        radius = rep(.self$radius, length(types))
      )
      
      mu_p <- ggplot2::ggplot(plot_mu, ggplot2::aes(x = x, y = y, color = expr)) +
        ggplot2::geom_point(size = 0.6) +
        ggplot2::scale_color_gradient(low = "lightyellow", high = "darkgreen", na.value = "gray80") +
        ggforce::geom_circle(
          data = circle_df,
          ggplot2::aes(x0 = x0, y0 = y0, r = radius),
          color = "red", size = 0.8, inherit.aes = FALSE
        ) +
        ggplot2::facet_wrap(~type, labeller = ggplot2::as_labeller(label_map_mu)) +
        ggplot2::labs(title = paste0(gene_name, "- mu"), color = "Expression")
      
      raw_p <- ggplot2::ggplot(plot_raw, ggplot2::aes(x = x, y = y, color = expr)) +
        ggplot2::geom_point(size = 0.6) +
        ggplot2::scale_color_gradient(low = "lightyellow", high = "darkgreen", na.value = "gray80") +
        ggforce::geom_circle(
          data = circle_df,
          ggplot2::aes(x0 = x0, y0 = y0, r = radius),
          color = "red", size = 0.8, inherit.aes = FALSE
        ) +
        ggplot2::facet_wrap(~type, labeller = ggplot2::as_labeller(label_map_raw)) +
        ggplot2::labs(title = paste0(gene_name, "- raw expression"), color = "Expression")
      
      # spot plots
      by_spot_all <- cbind(.self$spot_coords, expr = colSums(.self$spot_counts), type = "all")
      by_spot_gene <- cbind(.self$spot_coords, expr = .self$spot_counts[gene_name, , drop = TRUE], type = "gene")
      plot_spot <- as.data.frame(rbind(by_spot_all, by_spot_gene))
      plot_spot$expr <- as.numeric(plot_spot$expr)
      plot_spot$x <- as.numeric(plot_spot$x)
      plot_spot$y <- as.numeric(plot_spot$y)
      
      circle_df_spot <- data.frame(
        type = c("all", "gene"),
        x0 = .self$center[1, "x"],
        y0 = .self$center[1, "y"],
        radius = .self$radius
      )
      
      spot_gene <- ggplot2::ggplot(plot_spot[plot_spot$type == "gene", ], ggplot2::aes(x = x, y = y, color = expr)) +
        ggplot2::geom_point(size = 2) +
        ggplot2::scale_color_gradient(low = "lightyellow", high = "darkgreen", na.value = "gray80") +
        ggforce::geom_circle(
          data = circle_df_spot,
          ggplot2::aes(x0 = x0, y0 = y0, r = radius),
          color = "red", size = 0.8, inherit.aes = FALSE
        ) +
        ggplot2::labs(title = paste0(gene_name, "- Spot Raw Expression"), color = "Expression")
      
      if (2 %in% .self$sim_type) {
        spot_cell_prop <- t(.self$cell_spot_prop)[rownames(.self$spot_coords), , drop = FALSE]
      } else {
        spot_cell_prop <- as.matrix(stats::xtabs(~.self$cell_metadata$spot + rownames(.self$cell_metadata)))
      }
      
      spot_density_df <- data.frame(
        x = .self$spot_coords[, "x"],
        y = .self$spot_coords[, "y"],
        density = rowSums(spot_cell_prop)
      )
      
      spot_density <- ggplot2::ggplot(spot_density_df, ggplot2::aes(x = x, y = y, color = density)) +
        ggplot2::geom_point(size = 2) +
        ggplot2::scale_color_gradient(low = "lightyellow", high = "darkgreen", na.value = "gray80") +
        ggforce::geom_circle(
          data = circle_df_spot,
          ggplot2::aes(x0 = x0, y0 = y0, r = radius),
          color = "red", size = 0.8, inherit.aes = FALSE
        ) +
        ggplot2::labs(title = "Spot Cell Density", color = "Cell Counts")
      
      (mu_p / raw_p) | (spot_density / spot_gene) +
        patchwork::plot_annotation(title = paste(gene_name, desc))
    }
  )
)
