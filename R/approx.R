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

# simplified as data sd = 1
gamma_approx2 <- function(cov){
  M <- nrow(cov)
  # diag(cov) <- 0
  u <- (sumsquared(cov) / M) * 2
  list(
    M = M,
    u = u
  )
}

#' @name fdp-pwr
#' @title Calculates false-discovery proportions and statistical power
#' @param rej integer indices of rejected hypotheses
#' @param support integer indices of the underlying true hypotheses
#' (with true alternative hypotheses)
#' @return \code{fdp} returns the false-discovery proportions and \code{pwr}
#' returns statistical power
#' @examples
#'
#' # Underlying support is the first 100 hypotheses, but the rejection is
#' # 2 to 101
#' support <- 1:100
#' rejection <- c(2:101)
#'
#' # Total rejection (length(rejection)) is 100; false rejection is 1
#' # FDP = 1/100
#' fdp(rejection, support)
#'
#' # True positives = 99, total true hypotheses is 100
#' # power is 99/100
#' pwr(rejection, support)
#'
#' @export
fdp <- function(rej, support){
  if(!length(rej)){ return(0) }
  mean(!rej %in% support)
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

#' @rdname fdp-pwr
#' @export
pwr <- function(rej, support){
  if(!length(support)){ return(1) }
  mean(support %in% rej)
}

