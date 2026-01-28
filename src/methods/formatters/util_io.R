# src/utils/io_results.R

#' Find the most recent result file matching prefix
#'
#' @param dir character, directory to search
#' @param prefix character, e.g. "celina", "stance", "cside"
#' @param ext character, file extension without dot ("csv" or "xlsx")
#'
#' @return character path or NA if none found
find_latest_result_file <- function(dir, prefix, ext = "csv") {
  if (!dir.exists(dir)) return(NA)
  
  files <- list.files(
    dir,
    pattern = paste0("^", prefix, "(_[0-9]+)?\\.", ext, "$"),
    full.names = TRUE
  )
  
  if (length(files) == 0) return(NA)
  
  nums <- suppressWarnings(
    as.integer(sub(paste0("^.*_", prefix, "_?|\\.", ext, "$"), "", files))
  )
  nums[is.na(nums)] <- 0
  
  files[which.max(nums)]
}


#' Append a results table and split if too large
#'
#' @param df data.frame to append
#' @param dir output directory
#' @param prefix file prefix (e.g. "celina")
#' @param ext "csv" or "xlsx"
#' @param max_rows maximum rows per file
#'
#' @return invisible(path) of file written
append_results_chunked <- function(
    df,
    dir,
    prefix,
    ext = "csv",
    max_rows = 5e5
) {
  stopifnot(is.data.frame(df))
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  
  ext <- match.arg(ext, c("csv", "xlsx"))
  
  last_file <- find_latest_result_file(dir, prefix, ext)
  
  if (is.na(last_file)) {
    out_file <- file.path(dir, paste0(prefix, ".", ext))
    combined <- df
  } else {
    if (ext == "csv") {
      old <- utils::read.csv(last_file)
    } else {
      old <- readxl::read_xlsx(last_file)
    }
    combined <- rbind(old, df)
    out_file <- last_file
  }
  
  # split if too large
  if (nrow(combined) > max_rows) {
    idx <- sub(
      paste0("^.*", prefix, "_?|\\.", ext, "$"),
      "",
      out_file
    )
    idx <- ifelse(idx == "", 0, as.integer(idx))
    if (is.na(idx)) stop("Failed to infer chunk index for ", prefix)
    
    out_file <- file.path(dir, paste0(prefix, "_", idx + 1, ".", ext))
    combined <- df
  }
  
  if (ext == "csv") {
    utils::write.csv(combined, out_file, row.names = FALSE)
  } else {
    writexl::write_xlsx(list(Results = combined), out_file)
  }
  
  invisible(out_file)
}


#' Append a results table and split if too large
#'
#' @param df data.frame to append
#' @param dir output directory
#' @param prefix file prefix (e.g. "celina")
#' @param ext "csv" or "xlsx"
#' @param max_rows maximum rows per file
#'
#' @return invisible(path) of file written
append_results_chunked <- function(
    df,
    dir,
    prefix,
    ext = "csv",
    max_rows = 5e5
) {
  stopifnot(is.data.frame(df))
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  
  ext <- match.arg(ext, c("csv", "xlsx"))
  
  last_file <- find_latest_result_file(dir, prefix, ext)
  
  if (is.na(last_file)) {
    out_file <- file.path(dir, paste0(prefix, ".", ext))
    combined <- df
  } else {
    if (ext == "csv") {
      old <- utils::read.csv(last_file)
    } else {
      old <- readxl::read_xlsx(last_file)
    }
    combined <- rbind(old, df)
    out_file <- last_file
  }
  
  # split if too large
  if (nrow(combined) > max_rows) {
    idx <- sub(
      paste0("^.*", prefix, "_?|\\.", ext, "$"),
      "",
      out_file
    )
    idx <- ifelse(idx == "", 0, as.integer(idx))
    if (is.na(idx)) stop("Failed to infer chunk index for ", prefix)
    
    out_file <- file.path(dir, paste0(prefix, "_", idx + 1, ".", ext))
    combined <- df
  }
  
  if (ext == "csv") {
    utils::write.csv(combined, out_file, row.names = FALSE)
  } else {
    writexl::write_xlsx(list(Results = combined), out_file)
  }
  
  invisible(out_file)
}
