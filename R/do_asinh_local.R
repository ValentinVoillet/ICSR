#' Apply Inverse Hyperbolic Sine (Arcsinh) Transformation
#'
#' This function performs an arcsinh transformation on specific columns of a
#' \code{data.table}. This implementation is adapted from the \bold{Spectre} R package.
#' Arcsinh is the standard transformation for flow and mass cytometry data
#' to compress high-intensity signals and visualize data around zero.
#'
#' @param dat A \code{data.table} containing the raw events.
#' @param use.cols Character vector. The names of the columns to be transformed.
#' @param cofactor Numeric. The divisor used before applying \code{asinh}.
#'   Common values are 5 for CyTOF and 150-500 for Flow Cytometry. Default is 5.
#' @param append.cf Logical. If \code{TRUE}, appends the cofactor value to the
#'   new column names (e.g., "_asinh_cf5"). Default is \code{FALSE}.
#' @param reduce.noise Logical. An experimental step to reduce noise from negative
#'   values by shifting data and replacing negatives with a normal distribution
#'   around zero. Use with caution. Default is \code{FALSE}.
#'
#' @return The original \code{data.table} with additional columns containing the
#'   transformed values.
#'
#' @importFrom data.table as.data.table
#' @export
#'
#' @examples
#' # df <- do_asinh_local(my_data, use.cols = c("CD4", "CD8"), cofactor = 150)
do_asinh_local <- function(dat,
                           use.cols,
                           cofactor = 5,
                           append.cf = FALSE,
                           reduce.noise = FALSE) {

  # 1. Setup and Validation
  # ---------------------------------------------------------------------------
  # Extract only the columns of interest for transformation
  value <- dat[, use.cols, with = FALSE]
  # Ensure all targeted columns are numeric to prevent calculation errors
  is_numeric_col <- sapply(value, is.numeric)
  if (!all(is_numeric_col)) {
    message("Non-numeric columns detected in transformation list:")
    print(is_numeric_col)
    stop("Transformation aborted: All 'use.cols' must be numeric.")
  }

  # 2. Optional Experimental Noise Reduction
  # ---------------------------------------------------------------------------
  # This logic mimics cytofkit's approach to handling negative values
  # https://github.com/JinmiaoChenLab/cytofkit/issues/71
  if(reduce.noise == TRUE){
    message("This noise reduction function is experimental, and should be used with caution")
    value <- value-1
    loID <- which(value < 0)
    if(length(loID) > 0)
      value[loID] <- rnorm(length(loID), mean = 0, sd = 0.01)
  }

  # 3. Arcsinh Calculation
  # ---------------------------------------------------------------------------
  # Formula: f(x) = asinh(x / cofactor)
  value <- value / cofactor
  value <- asinh(value)

  # 4. Column Renaming Logic
  # ---------------------------------------------------------------------------
  # Determine the suffix based on whether the user wants the cofactor included
  if(append.cf == TRUE){
    if(length(use.cols) > 1){
      names(value) <- paste0(names(value), "_asinh_cf", cofactor)
    }
    if(length(use.cols) == 1){
      names(value) <- paste0(use.cols, "_asinh_cf", cofactor)
    }
  }
  if(append.cf == FALSE){
    if(length(use.cols) > 1){
      names(value) <- paste0(names(value), "_asinh")
    }
    if(length(use.cols) == 1){
      names(value) <- paste0(use.cols, "_asinh")
    }
  }

  # 5. Wrap up and Return
  # ---------------------------------------------------------------------------
  # Combine the original data with the new transformed columns
  dat <- cbind(dat, value)
  return(dat)
}

