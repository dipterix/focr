#' @docType package
NULL

#' @import stats
#' @import graphics
#' @importFrom utils read.table
#' @importFrom Rcpp sourceCpp
#' @useDynLib focr, .registration = TRUE
NULL

deparse_svec <- function(nums, connect = ':', concatenate = FALSE, collapse = ',', max_lag = 1){
  nums <- nums[is.finite(nums)]
  if(length(nums) == 0){
    return('')
  }
  alag <- seq_len(max(1, max_lag))
  nums <- sort(unique(nums))
  lg <- c(NA, nums)[seq_len(length(nums))]
  ind <- nums - lg
  ind[1] <- 0
  ind2 <- c(ind[-1], -1)
  apply(cbind(nums[!ind %in% alag], nums[!ind2 %in% alag]), 1,function(x){
    if(x[1] == x[2]){
      sprintf("%.0f", x[1])
    }else{
      paste(sprintf("%.0f", x[c(1,2)]), collapse = connect)
    }
  }) ->
    re
  if(concatenate){
    re <- paste(re, collapse = collapse)
  }
  re
}


debug_verbose <- function(..., appendLF = TRUE, verbose = TRUE){
  if(verbose){
    message(..., appendLF = appendLF)
  }
}


