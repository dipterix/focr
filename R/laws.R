
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

laws_pval <- function(pv, bandwidth, dimension = c("one", "two", "three"), alpha = 0.05, initial_filter = 0.9){
  dimension <- match.arg(dimension)
  bh.th<-BH(pv, initial_filter)$tau
  if(bandwidth < 1){
    bandwidth <- bandwidth * length(pv)
  }
  switch (
    dimension,
    'one' = { pis.hat<-pis_1D.func(pv, tau=bh.th, h=bandwidth) },
    'two' = { pis.hat<-pis_2D.func(pv, tau=bh.th, h=bandwidth) },
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


sabha_pval <- function(pv, bandwidth, dimension = c("one", "two", "three"), alpha = 0.05, initial_filter = 0.9){
  dimension <- match.arg(dimension)
  bh.th<-BH(pv, initial_filter)$tau
  switch (
    dimension,
    'one' = { pis.hat<-pis_1D.func(pv, tau=bh.th, h=bandwidth) },
    'two' = { pis.hat<-pis_2D.func(pv, tau=bh.th, h=bandwidth) },
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

