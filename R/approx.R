gamma_approx <- function(cov){
  M <- nrow(cov)
  d <- sqrt(diag(cov))
  if(any(abs(d - 1) > 1e-6)){
    cor <- t(t(cov / d) / d)
  } else {
    cor <- cov
  }

  diag(cor) <- 0
  u <- (sum(cor^2) / M + 1) * 2
  list(
    M = M,
    u = u
  )
}



fdp <- function(rej, sig){
  if(!length(rej)){ return(0) }
  mean(!rej %in% sig)
}

fcr <- function(rej, sig){
  if(!length(rej)){ return(0) }
  x <- deparse_svec(rej, concatenate = FALSE, connect = ':')
  w <- sapply(x, function(s){
    s <- eval(parse(text = s))
    l <- length(s)
    if(any(s %in% sig)){
      rw <- 0L
    } else {
      rw <- l
    }
    c(rw, l)
  })
  sum(w[1,]) / sum(w[2,])
}

pwr <- function(rej, sig){
  if(!length(sig)){ return(1) }
  mean(sig %in% rej)
}

