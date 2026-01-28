library(dplyr)
library(ggplot2)
library(ggforce)
library(patchwork)
library(scales)
library(Matrix)
library(readr)
library(arrow)
source("./scDesign_simulators/simulation_util.R")

setClassUnion("dgCMatrixOrNULL", c("dgCMatrix", "NULL"))
setClassUnion("numericOrNULL",  c("numeric",  "NULL"))

pseudo_data_ori = setRefClass(
  "scDesign3_data_ori",
  fields = list(
    cell_counts_ori       = "dgCMatrix",
    cell_metadata         = "data.frame",
    cell_type_proportion  = "numeric",
    cell_names = "character", 
    gene_names = "character",
    spot_names = "character",
    
    # NEW: CSV-backed fields
    cell_spot_prop        = "dgCMatrix",      # from cell_spot_prop.csv     (cells x spots, as read)
    gene_cell_keep        = "dgCMatrixOrNULL",      # from gene_dropout_frac.csv (genes x cells, as read)
    gene_cell_drop        = "dgCMatrixOrNULL",      # from gene_dropout_frac.csv (genes x cells, as read)
    mask_keep             = "dgCMatrixOrNULL",
    mask_fallback         = "dgCMatrixOrNULL",
    cell_fallback         = "numericOrNULL",
    pseudo_meta           = "data.frame",    # from pseudo_meta.csv
    pseudo_count          = "matrix"
  ),
  methods = list(
    initialize = function(cell_count, cell_meta, pseudo_meta, pseudo_count, cell_spot_prop, gene_cell_keep = NULL, gene_cell_drop = NULL, 
                          mask_keep = NULL, mask_fallback= NULL, cell_fallback= NULL) {
      # original checks and fields
      if (any(rownames(cell_count) != rownames(gene_cell_keep))) {
        stop("Number or order of genes in cell_count and gene_cell_keep doesn't match!")
      }
      if (any(colnames(cell_count) != rownames(cell_meta)) || 
          any(colnames(cell_count) != rownames(cell_spot_prop)) ||
          any(colnames(cell_count) != colnames(gene_cell_keep)) ||
          any(colnames(cell_count) != colnames(mask_keep)) || 
          any(colnames(cell_count) != colnames(mask_fallback)) || 
          any(colnames(cell_count) != names(cell_fallback)) 
      ){
        stop("Number or order of cells in each input doesn't match!")
      }
      if (any(colnames(cell_spot_prop) != rownames(pseudo_meta))){
        stop("Number or order of spots in cell_spot_prop and pseudo_meta doesn't match!")
      }
      
      stopifnot(colnames(cell_meta) == c("x", "y", "type"))
      
      .self$gene_names = rownames(cell_count)
      .self$cell_names = colnames(cell_count)
      .self$spot_names = colnames(cell_spot_prop)
      stopifnot(length(.self$gene_names) > 0 & length(.self$cell_names) > 0 & length(.self$spot_names) > 0)
      message(paste0("Loaded ", length(.self$spot_names), " spots, ", length(.self$cell_names), " cells, ", length(.self$gene_names), " genes..."))
      
      .self$cell_counts_ori = cell_count
      .self$cell_metadata = cell_meta
      if (!is.null(gene_cell_keep)) .self$gene_cell_keep = gene_cell_keep
      if (!is.null( gene_cell_drop)) .self$gene_cell_drop = gene_cell_drop
      if (!is.null(mask_keep )) .self$mask_keep = mask_keep 
      if (!is.null( mask_fallback)) .self$mask_fallback = mask_fallback
      if (!is.null(cell_fallback )) .self$cell_fallback = cell_fallback 
      .self$cell_spot_prop = cell_spot_prop
      .self$pseudo_meta = pseudo_meta
      .self$pseudo_count = pseudo_count
      
      
      tp <- prop.table(table(.self$cell_metadata$type)[sort(unique(.self$cell_metadata$type))])
      .self$cell_type_proportion <- as.vector(tp)
      names(.self$cell_type_proportion) <- names(tp)
      
    },
    
    sim_plot = function(gene_name, desc = "") {
      # unchanged: original plot for raw expression
      expr <- log1p(.self$cell_counts_ori[gene_name, ])
      expr <- scales::rescale(expr)
      loc <- .self$cell_metadata[, c("x", "y")]
      cell_type <- .self$cell_metadata$type
      
      df_all <- tibble(Expression = expr, x = loc[, 1], y = loc[, 2], cell_type = "Overall")
      df_by_type <- tibble(Expression = expr, x = loc[, 1], y = loc[, 2], cell_type = cell_type)
      df_plot <- bind_rows(df_all, df_by_type)
      
      ordered_levels <- c("Overall", names(sort(.self$cell_type_proportion, decreasing = TRUE)))
      df_plot$cell_type <- factor(df_plot$cell_type, levels = ordered_levels)
      
      original <- ggplot(df_plot, aes(x = x, y = y, color = Expression)) +
        geom_point(size = 0.1) +
        facet_wrap(~cell_type) +
        scale_colour_gradientn(colors = viridis_pal(option = "magma")(10)) +
        coord_fixed(ratio = 1) +
        theme(axis.text.x = element_text(angle = 45)) +
        ggtitle(paste("Raw Expression of", gene_name, "(Overall + Cell Types)"))
      
      original
    }
  )
)
