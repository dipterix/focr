get_sapply <- function(){
  # if(system.file('', package = "future.apply") != ''){
  #   return(function(..., future.seed=TRUE){
  #     future.apply::future_sapply(..., future.seed=future.seed)
  #   })
  # } else {
    return(sapply)
  # }
}


fastcov2 <- function(data, col1, col2){

  if(missing(col1)){
    col1 <- seq_len(ncol(data))
  }
  if(missing(col2)){
    col2 <- col1
  }
  col1 <- col1[col1 >= 1 & col1 <= ncol(data) ]
  col2 <- col2[col2 >= 1 & col2 <= ncol(data) ]

  if(!length(col1) || !length(col2)){
    return(numeric(0))
  }

  # fastcov(data, nrow(data), ncol(data), col1, col2)
  setThreads(0, TRUE)
  fastcov(data, col1, col2)

}
