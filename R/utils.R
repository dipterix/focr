get_sapply <- function(){
  if(system.file('', package = "future.apply") != ''){
    return(function(..., future.seed=TRUE){
      future.apply::future_sapply(..., future.seed=future.seed)
    })
  } else {
    return(sapply)
  }
}
