
laws_sabha_gauss <- function(
  data, bandwidth, alpha = 0.05, method = c('laws', "sabha"),
  initial_filter = 0.9, side = c("two", "right", "left"),
  dimension = c("one", "two", "three")){

  method <- match.arg(method)
  side <- match.arg(side)
  dimension <- match.arg(dimension)
  x <- sqrt(nrow(data)) * colMeans(data) / apply(data, 2, sd)
  switch (side,
          'two' = {
            pv <- 2 * pnorm(-abs(x), 0, 1)
          },
          "left" = {
            pv <- pnorm(x, 0, 1)
          }, {
            pv <- pnorm(-x, 0, 1)
          }
  )

  if(method == 'laws'){
    res <- laws_pval(pv, bandwidth = bandwidth, dimension = dimension,
                     alpha = alpha, initial_filter = initial_filter)
  } else {
    res <- sabha_pval(pv, bandwidth = bandwidth, dimension = dimension,
                     alpha = alpha, initial_filter = initial_filter)
  }


  res$side <- side
  res
}

get_bandwidth <- function(l, d){
  if(system.file('', package = 'kedd') == ''){
    stop("Bandwidth is NA. Please install `kedd` package to automatically calculate the bandwidth for LAWS and SABHA methods")
  }
  if(length(d)){
    bw <- kedd::h.ccv(seq_len(max(d)))
  } else {
    bw <- kedd::h.ccv(seq_len(l))
  }
  bw$h
}

laws_pval <- function(pv, bandwidth, dimension = c("one", "two", "three"), alpha = 0.05, initial_filter = 0.9, verbose = FALSE){
  if(is.na(bandwidth)){
    bandwidth <- get_bandwidth(length(pv), dim(pv))
  }
  stopifnot(isTRUE(bandwidth > 0))
  stopifnot(isTRUE(initial_filter >= 0 && initial_filter <= 1))
  dimension <- match.arg(dimension)
  # Do not remove 1s
  bh.th<-BH(pv, initial_filter, filter = 2)$tau
  stopifnot(is.double(bh.th))
  switch (
    dimension,
    'one' = {
      pis.hat<-pis_1D(pv, tau=bh.th, h=bandwidth, verbose = verbose)
    },
    'two' = {
      dm <- dim(pv)
      pis.hat<-pis_2D(pv, dm[1], dm[2], tau=bh.th, h=bandwidth, verbose = verbose)
    },
    { pis.hat<-pis_3D.func(pv, tau=bh.th, h=bandwidth) }
  )

  res <- law.func(pvs=pv, pis.hat, alpha)
  pretty_list(
    nrejs = res$nr,
    rejs = which(as.logical(res$de)),
    pis_hat = pis.hat,
    alpha = alpha,
    bandwidth = bandwidth,
    dimension = dimension,
    initial_filter = initial_filter,
    method = "LAWS",
    details = res
  )
}


sabha_pval <- function(pv, bandwidth, dimension = c("one", "two", "three"), alpha = 0.05, initial_filter = 0.9, verbose = FALSE){
  if(is.na(bandwidth)){
    bandwidth <- get_bandwidth(length(pv), dim(pv))
  }
  stopifnot(isTRUE(bandwidth > 0))
  stopifnot(isTRUE(initial_filter >= 0 && initial_filter <= 1))
  dimension <- match.arg(dimension)
  bh.th<-BH(pv, initial_filter)$tau
  stopifnot(is.double(bh.th))
  switch (
    dimension,
    'one' = {
      pis.hat<-pis_1D(pv, tau=bh.th, h=bandwidth, verbose = verbose)
    },
    'two' = {
      dm <- dim(pv)
      pis.hat<-pis_2D(pv, dm[1], dm[2], tau=bh.th, h=bandwidth, verbose = verbose)
    },
    { pis.hat<-pis_3D.func(pv, tau=bh.th, h=bandwidth) }
  )

  res <- sab.func(pvs=pv, pis.hat, alpha)
  pretty_list(
    nrejs = res$nr,
    rejs = which(as.logical(res$de)),
    pis_hat = pis.hat,
    alpha = alpha,
    bandwidth = bandwidth,
    dimension = dimension,
    initial_filter = initial_filter,
    method = "SABHA",
    details = res
  )
}

