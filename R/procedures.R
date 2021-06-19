BY <- function(pvals, alpha = 0.05) {
  BH(pvals, alpha = alpha / log(length(pvals)))
}

BH <- function (pvals, alpha = 0.05) {
  m <- length(pvals)
  pvals[is.na(pvals)] <- 1

  o <- order(pvals, decreasing = TRUE)
  ro <- order(o)
  i <- m:1L
  qvals <- pmin(1, cummin(m / i * pvals[o]))[ro]

  sel <- qvals <= alpha
  nrejs <- sum(sel)
  if(!length(nrejs) || !nrejs){
    rejs <- integer(0L)
  } else {
    rejs <- which(sel)
  }
  tau <- alpha * nrejs / m
  return(list(
    nrejs = nrejs, rejs = rejs, order = (m + 1L)-ro,
    qvals = qvals, tau = tau
  ))
}

BH2 <- function(pvals, alpha = 0.05){
  res <- BH(pvals, alpha)
  if(!res$nrejs){
    res$initial_rejs <- integer(0L)
    return(res)
  }
  m <- length(pvals)
  k <- res$order[res$rejs]
  o <- order(k, decreasing = TRUE, method = 'radix')
  kn <- k[o]
  mok <- (m + 1L - kn) / (1 - pvals[res$rejs][o])
  sel <- (mok[-res$nrejs] - mok[-1]) > 0
  if(any(sel)){
    k <- min(kn[sel])
  } else {
    k <- res$nrejs
  }
  m0 <- mok[ res$nrejs + 1L - k ]
  m0 <- ceiling(m0)

  res2 <- BH(pvals, alpha * m / m0)
  res2$initial_rejs <- res$rejs
  res2
}
