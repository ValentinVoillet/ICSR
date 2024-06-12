#'
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Function
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' arcsinh transformation (from Spectre)
#'
#' @param dat data.table Input sample.
#' @param use.cols Vector of character column names – these columns will be transformed and added to the data.table as new columns.
#' @param cofactor Co-factor to use for arcsinh transformation.
#' @param append.cf Appends the co-factor used to the end of the name of the transformed columns.
#' @param reduce.noise This is an experimental calculation which should reduce noise from negative values. Use with caution.
#'
#' @return arcsinh transformed values
#' @export
#'
#' @examples
do_asinh_local <- function(dat,
                           use.cols,
                           cofactor = 5,
                           append.cf = FALSE,
                           reduce.noise = FALSE) {

  #- Setup data
  value <- dat[,use.cols,with = FALSE]

  #- Numeric checks
  if(isFALSE(all(sapply(value, is.numeric)))){
    message("It appears that one column in your dataset is non numeric")
    print(sapply(value, is.numeric))
    stop("do.asinh stopped")
  }

  #- Optional noise reduction
  # https://github.com/JinmiaoChenLab/cytofkit/issues/71
  if(reduce.noise == TRUE){
    message("This noise reduction function is experimental, and should be used with caution")
    value <- value-1
    loID <- which(value < 0)
    if(length(loID) > 0)
      value[loID] <- rnorm(length(loID), mean = 0, sd = 0.01)
  }

  #- Arcsinh calculation
  value <- value / cofactor
  value <- asinh(value)

  #- Options to append the CF used
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

  #- Wrap up
  dat <- cbind(dat, value)
  return(dat)
}

