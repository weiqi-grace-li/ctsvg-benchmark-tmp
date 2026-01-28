# src/evaluation/evaluate_details
# Summary metrics from detail tables.

evaluate_details <- function(
    input_dir,
    save_path,
    save_sheet,
    analyze_method = c("celina", "stance", "spvc", "ctsv", "cside"),
    threshold = c(0.01, 0.03, 0.05, 0.07, 0.09, 0.1, 0.2)
) {
  if (!exists("collect_details", mode = "function", inherits = TRUE)) {
    stop("collect_details() not found. Load src/evaluation/collect_details.R first.")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package 'dplyr' is required.")
  if (!requireNamespace("magrittr", quietly = TRUE)) stop("Package 'magrittr' is required.")
  
  detail_result <- collect_details(input_dir, analyze_method)
  
  summary_result <- list()
  
  for (thres in threshold) {
    message(sprintf("ANALYZE | threshold = %s", thres))
    
    for (method in names(detail_result)) {
      df <- detail_result[[method]]
      if (is.null(df) || nrow(df) == 0) next
      
      summary_result[[paste0(method, "-", thres)]] <- df %>%
        dplyr::mutate(
          has_ctsvg = sapply(strsplit(gene_type, ","), function(x) "ctsvg" %in% x),
          has_marker = sapply(strsplit(gene_type, ","), function(x) "marker" %in% x),
          tested = (p_value >= 0)
        ) %>%
        # --- first group operation for ep calculation
        dplyr::group_by(simulation_name, test_method, deconv_method, seed) %>%
        dplyr::mutate(
          ep_denom = sum(has_ctsvg, na.rm = TRUE),
          padj_rank = NA_real_
        ) %>%
        dplyr::mutate(
          padj_rank = replace(
            padj_rank,
            tested,
            rank(p_adj[tested], ties.method = "first")
          ),
          in_top_ep = tested & (padj_rank <= ep_denom)
        ) %>%
        dplyr::ungroup() %>%
        # --- second group operation
        dplyr::group_by(simulation_name, test_method, deconv_method, seed, cell_type, cell_proportion, gene_type) %>%
        dplyr::summarise(
          # --- fdp & power
          rejected_adj = sum(p_value >= 0 & p_adj < thres),
          rejected_raw = sum(p_value >= 0 & p_value < thres),
          tp = sum(p_value >= 0 & p_adj < thres & has_ctsvg),
          total_true = sum(has_ctsvg),
          dispersion = mean(dispersion),
          
          # --- fdp without marker
          rejected_adj_no_marker = sum(p_value >= 0 & p_adj < thres & ifelse(has_ctsvg, TRUE, !has_marker)),
          
          true_tested = sum(p_value >= 0 & has_ctsvg),
          tested = sum(p_value >= 0),
          failed = sum(p_fail == 1),
          fail_rejected = sum(p_fail == 1 & p_value >= 0 & p_value < thres),
          not_fail_rejected = sum(p_fail == 0 & p_value >= 0 & p_value < thres),
          total = dplyr::n(),
          threshold = thres,
          
          # --- type 1 & power
          type1_num = sum((p_value >= 0) & (p_value < thres) & (!has_ctsvg)),
          type1_den = sum(!has_ctsvg),
          power_num = sum((p_value >= 0) & (p_value < thres) & has_ctsvg),
          power_den = sum(has_ctsvg),
          
          # convergence (kept exactly)
          convergence_ctsvg = sum(p_value == 1 & has_ctsvg),
          convergence_null = sum(p_value == 1 & !has_ctsvg),
          
          # --- ep
          ep_denom_check = mean(ep_denom),
          in_top = sum(in_top_ep),
          ep_num = sum(in_top_ep & has_ctsvg),
          
          .groups = "drop"
        )
    }
  }
  
  summary <- do.call(rbind, summary_result)
  
  if (is.na(save_path)) {
    return(summary)
  } else {
    if (!exists("write_sheet", mode = "function", inherits = TRUE)) {
      stop("write_sheet() not found. Load your sheet-writing utility before saving.")
    }
    write_sheet(summary, save_path = save_path, sheet_name = save_sheet)
    invisible(summary)
  }
}
