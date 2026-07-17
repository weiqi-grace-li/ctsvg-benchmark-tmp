# Load required libraries
library(spatstat)
library(STANCE)
library(dplyr)
library(purrr)
library(ggplot2)
library(ggforce)
library(patchwork)


# ── helper functions ───────────────────────────────────────────────────────────

decomp_generate_coords <- function(sim_type, pseudo_simulator) {
  # 2, 4, 5, would always fall into first basket, so 4, 5 not separated triggered. 
  if (2 %in% sim_type) {
    pseudo_simulator$cell_metadata[, c("x", "y")]
  } else if (any(c(0, 1, 3) %in% sim_type)) {
    pp <- runifpoint(
      n   = nrow(pseudo_simulator$cell_metadata),
      win = owin(c(min(pseudo_simulator$cell_metadata$x), max(pseudo_simulator$cell_metadata$x)),
                 c(min(pseudo_simulator$cell_metadata$y), max(pseudo_simulator$cell_metadata$y)))
    )
    data.frame(x = pp$x, y = pp$y)
  } else {
    pp <- rpoispp(lambda = 4000, win = owin(c(0, 1), c(0, 1)))
    data.frame(x = pp$x, y = pp$y)
  }
}

decomp_define_domain <- function(coords, sim_type) {
  center_idx <- sample(1:nrow(coords), 1)
  center     <- coords[center_idx, , drop = FALSE]
  radius     <- runif(1, 0.2, 0.4)
  if (any(c(0, 1, 2, 3, 4, 5) %in% sim_type)) {
    scale  <- c(max(coords$x) - min(coords$x), max(coords$y) - min(coords$y))
    radius <- radius * scale[which.min(scale)]
  }
  dists         <- sqrt((coords$x - center$x)^2 + (coords$y - center$y)^2)
  coords$domain <- ifelse(dists <= radius, "SD", "D")
  list(coords = coords, center = center, radius = radius)
}

decomp_assign_cell_types <- function(coords, sim_type, cell_type_proportion, pseudo_simulator) {
  if (all(c(1, 2) %in% sim_type)) {
    coords$type <- pseudo_simulator$cell_metadata$type
  } else {
    coords$type <- sample(names(cell_type_proportion), size = nrow(coords),
                          replace = TRUE, prob = cell_type_proportion)
  }
  coords
}

decomp_build_mu_dispersion <- function(coords, sim_type, pseudo_simulator, dispersion,
                                       gene_names, cell_type_proportion, n_markers) {
  K             <- length(cell_type_proportion)
  special_genes <- vector("list", K)
  names(special_genes) <- names(cell_type_proportion)

  disp_mat <- matrix(dispersion, nrow = nrow(coords), ncol = length(gene_names))
  colnames(disp_mat) <- gene_names

  mu <- matrix(1, nrow = nrow(coords), ncol = length(gene_names))
  colnames(mu) <- gene_names

  assume_3 <- if ((3 %in% sim_type) & !(1 %in% sim_type))
    names(sort(pseudo_simulator$cell_type_proportion, decreasing = TRUE))[1:K]
  else NULL

  # per-cell-type marker counts (even split; remainder distributed to first cell types)
  base_ct       <- floor(n_markers / K)
  extra_ct      <- n_markers %% K
  per_ct_counts <- c(rep(base_ct + 1L, extra_ct), rep(base_ct, K - extra_ct))
  used_markers  <- character(0)

  for (k in seq_along(special_genes)) {
    idx  <- coords$type == names(special_genes)[k]
    sd_k <- coords$domain == "SD"

    # set ct_svgs
    # 1, 5, would always fall into first basket, so 5 not separated triggered. 
    if (1 %in% sim_type) {
      if (length(pseudo_simulator$special_genes[[k]]$ct_svg) > 0)
        special_genes[[k]]$ct_svg <- pseudo_simulator$special_genes[[k]]$ct_svg
    } else if (any(c(-1, 0, 2, 3, 4) %in% sim_type)) {
      ct_svgs <- unlist(lapply(pseudo_simulator$special_genes, function(x) x$ct_svg))
      n_ctsvg <- floor(length(ct_svgs) / length(special_genes))
      special_genes[[k]]$ct_svg <- unique(ct_svgs[((k - 1) * n_ctsvg + 1):(k * n_ctsvg)])
    } else {
      special_genes[[k]]$ct_svg <- paste0("Gene", (300 + (k - 1) * 200 + 1):(300 + k * 200))
    }

    # set mu and dispersion for ct_svg
    if (length(special_genes[[k]]$ct_svg) > 0) {
      if (5 %in% sim_type) {
        mu[idx, special_genes[[k]]$ct_svg] <- as.matrix(
          t(pseudo_simulator$mu_matrix[special_genes[[k]]$ct_svg, idx]), drop = FALSE)
        disp_mat[idx, special_genes[[k]]$ct_svg] <- as.matrix(
          t(pseudo_simulator$dispersion[special_genes[[k]]$ct_svg, idx]), drop = FALSE)
      } else if (any(c(-1, 0, 1, 2, 3, 4) %in% sim_type)) {
        mu[idx & sd_k, special_genes[[k]]$ct_svg] <- 4
      } else {
        base_idx <- 300 + (k - 1) * 200
        for (i in 0:3) {
          gene_idx <- base_idx + (1:50) + i * 50
          mu[idx & sd_k, gene_idx] <- c(4, 2, 0.5, 0.25)[i + 1]
        }
      }
    }

    # set marker genes
    if (3 %in% sim_type) {
      special_genes[[k]]$marker <- c("")
      if (1 %in% sim_type) {
        non_ctsvg <- which(!(gene_names %in% special_genes[[k]]$ct_svg))
        mu[idx, non_ctsvg] <- as.matrix(
          t(pseudo_simulator$mu_matrix[non_ctsvg, idx]), drop = FALSE)
        disp_mat[idx, non_ctsvg] <- as.matrix(
          t(pseudo_simulator$dispersion[non_ctsvg, idx]), drop = FALSE)
      } else {
        replace_idx      <- which(pseudo_simulator$cell_metadata$type == assume_3[k])
        replace_nonctsvg <- which(
          !(gene_names %in% pseudo_simulator$special_genes[[assume_3[k]]]$ct_svg) &
          !(gene_names %in% special_genes[[k]]$ct_svg)
        )
        mu[idx, replace_nonctsvg] <- matrix(
          rep(pseudo_simulator$mu_matrix[replace_nonctsvg, replace_idx[1]], sum(idx)),
          nrow = sum(idx), ncol = length(replace_nonctsvg), byrow = TRUE)
        disp_mat[idx, replace_nonctsvg] <- matrix(
          rep(pseudo_simulator$dispersion[replace_nonctsvg, replace_idx[1]], sum(idx)),
          nrow = sum(idx), ncol = length(replace_nonctsvg), byrow = TRUE)
      }
    } else if (5 %in% sim_type) {
      pool                      <- setdiff(gene_names, c(used_markers, special_genes[[k]]$ct_svg))
      special_genes[[k]]$marker <- sample(pool, size = per_ct_counts[k], replace = FALSE)
      used_markers              <- c(used_markers, special_genes[[k]]$marker)
      mu[idx, special_genes[[k]]$marker] <- 3
    } else if (any(c(-1, 0, 1, 2, 4) %in% sim_type)) {
      pool                      <- setdiff(gene_names, used_markers)
      special_genes[[k]]$marker <- sample(pool, size = per_ct_counts[k], replace = FALSE)
      used_markers              <- c(used_markers, special_genes[[k]]$marker)
      marker_ctsvg              <- intersect(special_genes[[k]]$marker, special_genes[[k]]$ct_svg)
      pure_marker               <- setdiff(special_genes[[k]]$marker, marker_ctsvg)
      mu[idx, pure_marker]           <- 3
      mu[idx & !sd_k, marker_ctsvg] <- 3
    } else {
      special_genes[[k]]$marker <- paste0("Gene", ((k - 1) * 100 + 1):(k * 100))
      mu[idx, special_genes[[k]]$marker] <- 2
    }
  }

  list(mu = mu, disp_mat = disp_mat, special_genes = special_genes)
}

decomp_sample_counts <- function(mu, disp_mat) {
  counts      <- matrix(0, nrow = nrow(mu), ncol = ncol(mu))
  nonzero_idx <- which(mu > 0 & !is.na(mu))
  counts[nonzero_idx] <- rnbinom(
    n    = length(nonzero_idx),
    mu   = mu[nonzero_idx],
    size = disp_mat[nonzero_idx]
  )
  colnames(counts) <- colnames(mu)
  counts
}

decomp_assign_spots <- function(counts, coords, sim_type, pseudo_simulator, data_ori, gene_names) {
  # 2, 4, 5, would always fall into first basket, so 4, 5 not separated triggered.
  if (2 %in% sim_type) {
    count_mat           <- round(as.matrix(t(counts) %*% data_ori$cell_spot_prop), 0)
    count_mat           <- t(count_mat)
    colnames(count_mat) <- gene_names
    coords$spot         <- paste0("spot", apply(data_ori$cell_spot_prop, MARGIN = 1,
                                                 FUN = function(x) which.max(x)))
    spot_coords         <- as.matrix(data_ori$pseudo_meta[, c("x", "y")])
    rownames(count_mat) <- rownames(spot_coords) <- rownames(data_ori$pseudo_meta)
    df_spot             <- NULL
  } else {
    if (any(c(0, 1, 3) %in% sim_type)) {
      n_bins    <- floor(sqrt(nrow(pseudo_simulator$spot_coords)))
      coords$ix <- cut(coords$x, breaks = n_bins, labels = FALSE)
      coords$iy <- cut(coords$y, breaks = n_bins, labels = FALSE)
    } else {
      bin_size  <- 0.05
      coords$ix <- cut(coords$x, breaks = seq(0, 1, by = bin_size), labels = FALSE)
      coords$iy <- cut(coords$y, breaks = seq(0, 1, by = bin_size), labels = FALSE)
    }
    coords$spot <- paste0("spot", interaction(coords$ix, coords$iy))
    df          <- cbind(counts, coords)
    df_spot     <- df %>%
      group_by(spot) %>%
      summarise(
        x      = mean(x), y = mean(y),
        counts = list(colSums(across(seq_len(length(gene_names))))),
          # code explanation: seq_len(length(gene_names)) picks out the first G (n_genes) column;
          # ---  each of those column are sumed across rows in the group, output a vector of gene counts for this spot across all genes;
          # ---  list packs the vector into list to store in one cell. 
        types  = list(table(type)), .groups = "drop"
      ) %>% ungroup()
    count_mat             <- do.call(rbind, df_spot$counts)
    rownames(count_mat)   <- df_spot$spot
    colnames(count_mat)   <- gene_names
    spot_coords           <- as.matrix(df_spot[, c("x", "y")])
    rownames(spot_coords) <- df_spot$spot
  }
  list(count_mat = count_mat, coords = coords, spot_coords = spot_coords, df_spot = df_spot)
}

# spot_composition: proportion of cells of each type in a spot (count-based)
decomp_compute_spot_composition <- function(sim_type, coords, spot_coords,
                                             pseudo_simulator, df_spot, cell_type_proportion, data_ori) {
  comp_mat <- matrix(0, nrow = nrow(spot_coords), ncol = length(cell_type_proportion),
                     dimnames = list(rownames(spot_coords), names(cell_type_proportion)))
  if (all(c(1, 2) %in% sim_type)) {
    comp_mat <- pseudo_simulator$spot_composition
  } else if (2 %in% sim_type) {
    for (ct in colnames(comp_mat)) {
      ct_mask <- coords$type == ct
      if (any(ct_mask))
        comp_mat[, ct] <- colSums(data_ori$cell_spot_prop[ct_mask, , drop = FALSE])
    }
    comp_mat <- sweep(comp_mat, 1, pmax(rowSums(comp_mat), .Machine$double.eps), "/")
  } else {
    for (i in seq_len(nrow(df_spot))) {
      type_counts <- df_spot$types[[i]]
      comp_mat[i, names(type_counts)] <- type_counts / sum(type_counts)
    }
  }
  comp_mat
}



# ── class definition ───────────────────────────────────────────────────────────

decompose_simulator <- setRefClass("decompose_simulator",
                                    fields = list(
                                      dispersion          = "matrix",
                                      sim_type            = "ANY",
                                      seed                = "numeric",
                                      cell_type_proportion = "numeric",
                                      cell_counts         = "matrix",
                                      cell_counts_predrop = "matrix",
                                      cell_metadata       = "data.frame",
                                      spot_counts         = "matrix",
                                      spot_coords         = "matrix",
                                      spot_composition    = "matrix",
                                      special_genes       = "list",
                                      gene_names          = "character",
                                      center              = "data.frame",
                                      radius              = "numeric",
                                      mu_matrix           = "matrix",
                                      cell_spot_prop      = "matrix"
                                    ),

                                   # if there is nothing in sim_type, then just STANCE simulation
                                   # if there is -1, then STANCE simulation with adjusted number of total genes, ctsvgs, markers
                                   # if there is any of 0, 1, 2, 3, 4, 5 in, the baseline adjustments is assume done.

                                    methods = list(
                                      initialize = function(pseudo_simulator, data_ori, dispersion = 0.7, seed = 42, sim_type = c(), n_markers = NULL, cell_type_proportion = c("Type1" = 0.1, "Type2" = 0.3, "Type3" = 0.6)) {
                                        .self$sim_type = sim_type

                                        # --- error checking
                                        if (!all(sim_type %in% c(-1, 0, 1, 2, 3, 4, 5))){
                                          stop("Wrong decomposition type, only combinations of -1, 0, 1, 2, 3, 4, 5 are accepted!")
                                        }

                                        # --- enforce solid dependencies
                                        if (5 %in% sim_type){
                                          stopifnot(all(c(1, 2) %in% sim_type))
                                        }
                                        if (4 %in% sim_type){
                                          stopifnot(2 %in% sim_type)
                                        }

                                        # --- set cell type proportions
                                        if (1 %in% sim_type){
                                          .self$cell_type_proportion = pseudo_simulator$cell_type_proportion
                                        } else{
                                          .self$cell_type_proportion = cell_type_proportion
                                        }

                                        # --- set all gene names
                                        if (any(c(-1, 0, 1, 2, 3, 4, 5) %in% sim_type)){
                                          .self$gene_names = pseudo_simulator$gene_names
                                        } else {
                                          # only when nothing is in sim_type would go back to STANCE simulation here
                                          .self$gene_names = paste0("Gene", 1:1000)
                                        }

                                        if (is.null(n_markers)) n_markers <- min(max(300L, floor(0.3 * length(.self$gene_names))), length(.self$gene_names))

                                        .self$seed <- seed

                                        .self$simulate(sim_type, pseudo_simulator, data_ori, n_markers, dispersion)
                                      },

                                      simulate = function(sim_type, pseudo_simulator, data_ori, n_markers, dispersion) {
                                        set.seed(seed)

                                        # --- generate cell locations
                                        coords <- decomp_generate_coords(sim_type, pseudo_simulator) # with cell x, y 

                                        # --- define domain center and radius
                                        domain_res <- decomp_define_domain(coords, sim_type)
                                        coords     <- domain_res$coords # with cell x, y + domain
                                        center     <- domain_res$center
                                        radius     <- domain_res$radius

                                        # --- assign cell types
                                        coords <- decomp_assign_cell_types(coords, sim_type, .self$cell_type_proportion, pseudo_simulator)
                                        # with cell x, y, domain + cell type

                                        # --- build mu, dispersion, and special genes per cell type
                                        mu_disp_res   <- decomp_build_mu_dispersion(
                                          coords, sim_type, pseudo_simulator, dispersion,
                                          .self$gene_names, .self$cell_type_proportion, n_markers)
                                        mu            <- mu_disp_res$mu
                                        disp_mat      <- mu_disp_res$disp_mat
                                        special_genes <- mu_disp_res$special_genes

                                        # --- sample counts from NB
                                        counts <- decomp_sample_counts(mu, disp_mat)

                                        # --- assign row names
                                        if (any(c(0, 1, 2, 3, 4, 5) %in% sim_type)) {
                                          cell_names <- rownames(pseudo_simulator$cell_metadata)
                                        } else {
                                          cell_names <- paste0("Cell", 1:nrow(counts))
                                        }
                                        rownames(mu) <- rownames(counts) <- rownames(disp_mat) <- cell_names

                                        # --- save pre-capture-ratio counts
                                        counts_predrop <- counts

                                        # --- capture ratio adjustment
                                        if (4 %in% sim_type) {
                                          gene_cell_drop <- data_ori$gene_cell_drop[colnames(counts), rownames(counts)]
                                          counts         <- counts * t(1 - gene_cell_drop)
                                        }

                                        # --- assign cells to spots
                                        spot_res          <- decomp_assign_spots(counts, coords, sim_type, pseudo_simulator, data_ori, .self$gene_names)
                                        count_mat         <- spot_res$count_mat
                                        coords            <- spot_res$coords
                                        .self$spot_coords <- spot_res$spot_coords
                                        df_spot           <- spot_res$df_spot

                                        # --- cell metadata
                                        .self$cell_metadata <- coords[, c("x", "y", "type", "spot", "domain")]
                                        rownames(.self$cell_metadata) <- cell_names

                                        # --- store counts
                                        colnames(counts)         <- .self$gene_names
                                        colnames(counts_predrop) <- .self$gene_names
                                        .self$cell_counts         <- t(as.matrix(counts))
                                        .self$cell_counts_predrop <- t(as.matrix(counts_predrop))
                                        .self$spot_counts         <- t(count_mat)

                                        # --- spot composition (proportion of cells per type)
                                        .self$spot_composition <- decomp_compute_spot_composition(
                                          sim_type, coords, .self$spot_coords,
                                          pseudo_simulator, df_spot, .self$cell_type_proportion, data_ori)

                                        # --- other fields
                                        .self$special_genes  <- special_genes
                                        .self$center         <- center
                                        .self$radius         <- radius
                                        .self$mu_matrix      <- t(mu)
                                        .self$dispersion     <- t(disp_mat)
                                        .self$cell_spot_prop <- as.matrix(data_ori$cell_spot_prop)

                                        # --- drop zero genes / cells / spots
                                        .self$cell_counts         <- .self$cell_counts[.self$gene_names, ]
                                        .self$cell_counts_predrop <- .self$cell_counts_predrop[.self$gene_names, ]
                                        .self$spot_counts         <- .self$spot_counts[.self$gene_names, ]

                                        save_gene_names <- .self$gene_names[
                                          (rowSums(.self$cell_counts_predrop) > 0) & (rowSums(.self$spot_counts) > 0)]
                                        save_cell_names <- colnames(.self$cell_counts_predrop)[colSums(.self$cell_counts_predrop) > 0]
                                        save_spot_names <- colnames(.self$spot_counts)[colSums(.self$spot_counts) > 0]

                                        .self$gene_names          <- save_gene_names
                                        
                                        # note it is possible and very likely cell_counts have UMI = 0 cells.  Use cell_counts_predrop as reference
                                        .self$cell_counts         <- .self$cell_counts[save_gene_names, save_cell_names]
                                        .self$cell_counts_predrop <- .self$cell_counts_predrop[save_gene_names, save_cell_names]
                                        .self$cell_metadata       <- .self$cell_metadata[save_cell_names, ]
                                        .self$mu_matrix           <- .self$mu_matrix[save_gene_names, save_cell_names]
                                        if (length(.self$dispersion) > 1)
                                          .self$dispersion <- .self$dispersion[save_gene_names, save_cell_names]
                                        .self$spot_counts         <- .self$spot_counts[save_gene_names, save_spot_names]
                                        .self$spot_coords         <- .self$spot_coords[save_spot_names, ]
                                        .self$spot_composition    <- .self$spot_composition[save_spot_names, ]
                                      },

                                      sim_plot = function(gene_name, desc = "") {
                                        ct_svg_types <- names(.self$special_genes)[
                                          vapply(.self$special_genes, function(x) gene_name %in% x$ct_svg, logical(1))
                                        ]
                                        marker_types <- names(.self$special_genes)[
                                          vapply(.self$special_genes, function(x) gene_name %in% x$marker, logical(1))
                                        ]

                                        by_cell_mu = cbind(.self$cell_metadata, expr = .self$mu_matrix[gene_name,])
                                        all_mu = cbind(.self$cell_metadata, expr = .self$mu_matrix[gene_name,])
                                        all_mu$type = "all"
                                        plot_mu = rbind(by_cell_mu, all_mu)

                                        by_cell_raw = cbind(.self$cell_metadata, expr = .self$cell_counts[gene_name,])
                                        all_raw = cbind(.self$cell_metadata, expr = .self$cell_counts[gene_name,])
                                        all_raw$type = "all"
                                        plot_raw = rbind(by_cell_raw, all_raw)

                                        types = unique(plot_mu$type)

                                        mu_summary <- aggregate(expr ~ type, data = plot_mu, FUN = mean)
                                        label_data <- merge(data.frame(type = types),
                                                            mu_summary,
                                                            by = "type",
                                                            all.x = TRUE)

                                        label_map_mu <- setNames(
                                          paste0(
                                            label_data$type,
                                            ifelse(label_data$type %in% ct_svg_types,  " · ct_svg",  ""),
                                            ifelse(label_data$type %in% marker_types,  " · marker",""),
                                            " · mu = ",
                                            round(label_data$expr, 2)
                                          ),
                                         label_data$type
                                        )

                                        label_map_raw <- setNames(
                                          paste0(
                                            types,
                                            ifelse(types %in% ct_svg_types,  " · ct_svg",  ""),
                                            ifelse(types %in% marker_types,  " · marker","")
                                          ),
                                          types
                                        )

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
                                          facet_wrap(~type,  labeller = as_labeller(label_map_mu)) +
                                          # coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
                                          labs(title = paste0(gene_name, "- mu"), color = "Expression")

                                        raw = ggplot(plot_raw, aes(x = x, y= y, color = expr)) +
                                          geom_point(size = 0.6) +
                                          scale_color_gradient(low = "lightyellow", high = "darkgreen", na.value = "gray80") +
                                          geom_circle(data = circle_df, aes(x0=x0, y0=y0, r = radius), color = "red", size = 0.8,
                                                      inherit.aes = FALSE)+
                                          facet_wrap(~type,  labeller = as_labeller(label_map_raw)) +
                                          # coord_fixed(xlim = c(0, 1), ylim = c(0, 1))+
                                          labs(title = paste0(gene_name, "- raw expression") , color = "Expression")

                                        # --- cell resolution plot
                                        mu / raw + plot_annotation(title = paste(gene_name, desc))

                                        # --- spot resolution plot
                                        by_spot_all = cbind(.self$spot_coords, expr = colSums(.self$spot_counts), type = "all")
                                        by_spot_gene = cbind(.self$spot_coords, expr = .self$spot_counts[gene_name,], type = "gene")
                                        plot_spot = as.data.frame(rbind(by_spot_all, by_spot_gene))
                                        plot_spot$expr = as.numeric(plot_spot$expr)
                                        plot_spot$x = as.numeric(plot_spot$x)
                                        plot_spot$y = as.numeric(plot_spot$y)

                                        circle_df_spot <- data.frame(
                                          type   = c("all","gene"),
                                          x0     = .self$center[1, "x"],
                                          y0     = .self$center[1, "y"],
                                          radius = .self$radius
                                        )

                                        spot_all = ggplot(plot_spot[which(plot_spot$type == "all"),], aes(x = x, y= y, color = expr)) +
                                          geom_point(size = 2) +
                                          scale_color_gradient(low = "lightyellow", high = "darkgreen", na.value = "gray80") +
                                          geom_circle(data = circle_df_spot, aes(x0=x0, y0=y0, r = radius), color = "red", size = 0.8,
                                                      inherit.aes = FALSE)+
                                          labs(title = "Spot All Raw Expression", color = "Expression")

                                        spot_gene = ggplot(plot_spot[which(plot_spot$type == "gene"),], aes(x = x, y= y, color = expr)) +
                                          geom_point(size = 2) +
                                          scale_color_gradient(low = "lightyellow", high = "darkgreen", na.value = "gray80") +
                                          geom_circle(data = circle_df_spot, aes(x0=x0, y0=y0, r = radius), color = "red", size = 0.8,
                                                      inherit.aes = FALSE)+
                                          labs(title = paste0(gene_name, "- Spot Raw Expression"), color = "Expression")

                                        if (2 %in% .self$sim_type){
                                          spot_cell_prop = t(.self$cell_spot_prop)[rownames(.self$spot_coords), ]
                                        } else{
                                          spot_cell_prop = as.matrix(xtabs(~.self$cell_metadata$spot+rownames(.self$cell_metadata)))
                                        }

                                        spot_density_df = data.frame(x = .self$spot_coords[, "x"], y = .self$spot_coords[, "y"], density = rowSums(spot_cell_prop))
                                        spot_density = ggplot(spot_density_df, aes(x = x, y= y, color = density)) +
                                          geom_point(size = 2) +
                                          scale_color_gradient(low = "lightyellow", high = "darkgreen", na.value = "gray80") +
                                          geom_circle(data = circle_df_spot, aes(x0=x0, y0=y0, r = radius), color = "red", size = 0.8,
                                                      inherit.aes = FALSE)+
                                          labs(title = "Spot Cell Density", color = "Cell Counts")


                                        (mu / raw) | (spot_density / spot_gene)  +
                                          plot_annotation(title = paste(gene_name, desc))

                                      }
                                    )
)
