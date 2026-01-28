# --- solicit all method choices 
ctsvg_by_var = function(scdesign_marginal_path, top_n = 50, plot = TRUE){
  # - first create the variance explained by spatial 1 and 2
  load(scdesign_marginal_path)
  dev_table = matrix(NA, ncol = length(scdesign_marginal), nrow = length(scdesign_marginal[[1]]))
  colnames(dev_table) = names(scdesign_marginal)
  rownames(dev_table) = names(scdesign_marginal[[1]])
  for (ct in names(scdesign_marginal)){
    for (gene in rownames(dev_table)){
      if (!(is.na(scdesign_marginal[[ct]][[gene]]$fit)[1])){
        dev = summary(scdesign_marginal[[ct]][[gene]]$fit)$dev.expl
      } else{
        dev = 0
      }
      dev_table[gene, ct] = dev
    }
  }
  
  return(dev_table)
}



# --- construct special_genes  
construct_special_genes = function(selected_list, cell_types, save_dir, data_ori, agreement = 2, cell_detect_threshold = 5, fallback = 0, plot_dir = NULL
                                   , width = 4, height = 3.5, dpi = 300){
  
  total_ctsvg = 0 
  special_genes = vector("list", length(cell_types))
  names(special_genes) = cell_types
  for (ct in names(special_genes)){
    cell_selected = c()
    ctsvg_fallback = c()
    for (method in names(selected_list)){
      cell_selected = c(cell_selected, selected_list[[method]]$gene_name[which(selected_list[[method]]$cell_type == ct)])
      ctsvg_fallback = c(ctsvg_fallback, selected_list[[method]]$gene_name[which((selected_list[[method]]$rank <= fallback) & (selected_list[[method]]$cell_type == ct))])
    }
    
    # - select by agreement
    ctsvg_agree = names(which(table(cell_selected)>=agreement))
    if (length(ctsvg_agree)>0){
      message(sprintf("selected %d genes for cell %s based on agreement...", length(ctsvg_agree), ct)) 
    }
    if (fallback > 0){
      ctsvg_agree = unique(c(ctsvg_agree, ctsvg_fallback))
      if (length(ctsvg_agree)>0){
        message(sprintf("Adding each method top %d genes, total %d genes for cell %s...", fallback, length(ctsvg_agree), ct)) 
      }
    }
    
    # - add fall back 
    
    # - double check whether the gene in the cell type has at least over cel_detect threshold non_zero count 
    idx_ct = which(data_ori$cell_metadata$type == ct)
    ctsvg = c()
    for (gene in ctsvg_agree){
      non_zero = length(which(data_ori$cell_counts_ori[gene, idx_ct]>0))
      if (non_zero > cell_detect_threshold){
        ctsvg = c(ctsvg, gene)
      }
    }
    
    total_ctsvg = total_ctsvg + length(ctsvg)
    # - plot 
    if (!is.null(plot_dir)){
      dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
      col_map <- viridis_pal(option = "magma")(10)
      
      for (gene in ctsvg){
        base_df <- data_ori$cell_metadata[idx_ct, c("x","y"), drop = FALSE]
        expr <- as.numeric(data_ori$cell_counts_ori[gene, idx_ct])
        z    <- as.numeric(scale(expr))
        df   <- cbind(base_df, z = z)
        
        p <- ggplot(df, aes(x, y, colour = z)) +
          geom_point(size = 2) +
          stat_density_2d(aes(x = x, y = y), colour = "black", linewidth = 0.25, alpha = 0.4) +
          scale_colour_gradientn(colors = col_map, name = "z") +
          coord_equal() + theme_bw(base_size = 11) +
          theme(panel.grid = element_blank(), legend.position = "right") +
          labs(title = paste0(ct, " — ", gene), x = NULL, y = NULL)
        
        safe_ct <- gsub("[^[:alnum:]_]+", "-", ct)
        safe_g  <- gsub("[^[:alnum:]_]+", "-", gene)
        fn <- file.path(plot_dir, sprintf("[%d]_%s_%s.png", sum(cell_selected==gene),safe_ct, safe_g))
        ggsave(fn, p, width = width, height = height, dpi = dpi, bg = "white")
      } 
    }
    
    
    if (length(ctsvg) == 0) {
      ctsvg = c("")
    }
    special_genes[[ct]] = list(
      ct_svg = ctsvg, 
      marker = c("")
    )
  }
  
  save(special_genes, file = paste0(save_dir, "special_genes.RData"))
  
  message(paste0("Generated ", total_ctsvg, " ct_svgs ..."))
  
  return(special_genes)
  
}

library(ggplot2)
library(scales)   # for viridis_pal

# ctsvg_by_type: named list like list("B_cell"=c("GENE1","GENE2"), "T_cell"=c(...))
# ovarian_ori$cell_counts_ori: genes x cells (matrix or dgCMatrix)
# ovarian_ori$cell_metadata: data.frame with columns type and coordinates (x_transformed/y_transformed or x/y)
save_ctsvg_plots <- function(data_ori, ctsvg_by_type, out_dir,
                             width = 4, height = 3.5, dpi = 300) {
  counts <- data_ori$cell_counts_ori
  meta   <- data_ori$cell_metadata
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  col_map <- viridis_pal(option = "magma")(10)
  
  for (ct in names(ctsvg_by_type)) {
    genes <- ctsvg_by_type[[ct]]$ctsvg
    if (length(genes) == 0) next
    idx_ct <- meta$type == ct
    if (!any(idx_ct)) next
    base_df <- meta[idx_ct, c("x","y"), drop = FALSE]
    
    for (g in genes) {
      expr <- as.numeric(counts[g, idx_ct])
      z    <- as.numeric(scale(expr))
      df   <- cbind(base_df, z = z)
      
      p <- ggplot(df, aes(x, y, colour = z)) +
        geom_point(size = 2) +
        stat_density_2d(aes(x = x, y = y), colour = "black", linewidth = 0.25, alpha = 0.4) +
        scale_colour_gradientn(colors = col_map, name = "z") +
        coord_equal() + theme_bw(base_size = 11) +
        theme(panel.grid = element_blank(), legend.position = "right") +
        labs(title = paste0(ct, " — ", g), x = NULL, y = NULL)
      
      safe_ct <- gsub("[^[:alnum:]_]+", "-", ct)
      safe_g  <- gsub("[^[:alnum:]_]+", "-", g)
      fn <- file.path(out_dir, sprintf("ctsvg_%s_%s.png", safe_ct, safe_g))
      ggsave(fn, p, width = width, height = height, dpi = dpi, bg = "white")
    }
  }
  invisible(TRUE)
}


