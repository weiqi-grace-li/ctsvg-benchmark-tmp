# src/evaluation/collect_details.R

collect_details <- function(
    input_dir,
    analyze_method = c("celina", "stance", "spvc", "ctsv", "cside")
) {
  if (!dir.exists(input_dir)) stop("input_dir does not exist: ", input_dir)
  
  detail_result <- vector("list", length(analyze_method))
  names(detail_result) <- analyze_method
  
  # --- CELINA (xlsx)
  if ("celina" %in% analyze_method) {
    files <- list.files(
      path = input_dir,
      pattern = "^[^$]*celina_?([0-9]*)\\.xlsx$",
      full.names = FALSE
    )
    
    message(sprintf(
      "CELINA | Reading %d file(s): %s",
      length(files), paste(files, collapse = ", ")
    ))
    
    res <- do.call(
      rbind,
      lapply(files, function(f) readxl::read_excel(file.path(input_dir, f)))
    )
    
    if (!("p_adj" %in% colnames(res))) {
      stop("CELINA result lacks p_adj column")
    }
    
    res$p_fail <- 0
    detail_result[["celina"]] <- res
  }
  
  # --- STANCE (xlsx)
  if ("stance" %in% analyze_method) {
    files <- list.files(
      path = input_dir,
      pattern = "^[^$]*stance_?([0-9]*)\\.xlsx$",
      full.names = FALSE
    )
    
    message(sprintf(
      "STANCE | Reading %d file(s): %s",
      length(files), paste(files, collapse = ", ")
    ))
    
    res <- do.call(
      rbind,
      lapply(files, function(f) readxl::read_excel(file.path(input_dir, f)))
    )
    
    if (!("p_adj" %in% colnames(res))) {
      stop("STANCE result lacks p_adj column")
    }
    
    res$p_fail <- 0
    detail_result[["stance"]] <- res
  }
  
  # --- CTSV (csv)
  if ("ctsv" %in% analyze_method) {
    files <- list.files(
      path = input_dir,
      pattern = ".*ctsv_?([0-9]*)\\.csv$",
      full.names = FALSE
    )
    files <- files[!grepl("runtime", files)]
    
    message(sprintf(
      "CTSV | Reading %d file(s): %s",
      length(files), paste(files, collapse = ", ")
    ))
    
    res <- do.call(
      rbind,
      lapply(files, function(f) read.csv(file.path(input_dir, f)))
    )
    
    res$p_fail <- 0
    detail_result[["ctsv"]] <- res
  }
  
  # --- CSIDE (csv)
  if ("cside" %in% analyze_method) {
    files <- list.files(
      path = input_dir,
      pattern = ".*cside_?([0-9]*)\\.csv$",
      full.names = FALSE
    )
    files <- files[!grepl("runtime", files)]
    
    message(sprintf(
      "CSIDE | Reading %d file(s): %s",
      length(files), paste(files, collapse = ", ")
    ))
    
    res <- do.call(
      rbind,
      lapply(files, function(f) read.csv(file.path(input_dir, f)))
    )
    
    res$p_fail <- 0
    detail_result[["cside"]] <- res
  }
  
  # --- SPVC (csv)
  if ("spvc" %in% analyze_method) {
    files <- list.files(
      path = input_dir,
      pattern = ".*spvc_?([0-9]*)\\.csv$",
      full.names = FALSE
    )
    files <- files[!grepl("runtime", files)]
    
    message(sprintf(
      "SPVC | Reading %d file(s): %s",
      length(files), paste(files, collapse = ", ")
    ))
    
    res <- do.call(
      rbind,
      lapply(files, function(f) read.csv(file.path(input_dir, f)))
    )
    
    detail_result[["spvc"]] <- res
  }
  
  detail_result
}
