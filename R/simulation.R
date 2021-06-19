matern <- function (u, phi, kappa) {
  if (is.vector(u)){
    names(u) <- NULL
  }
  if (is.matrix(u)) {
    dimnames(u) <- list(NULL, NULL)
  }
  uphi <- u/phi
  uphi <- ifelse(u > 0,
                 (((2 ^ (-(kappa - 1))) /
                     ifelse(0, Inf, gamma( kappa ))) *
                    (uphi ^ kappa) * besselK(x = uphi, nu = kappa)),
                 1)
  uphi[u > 600 * phi] <- 0
  return(uphi)
}

gen_error <- function(n_points, type = c("AR", "exponential", "matern", "any", "iid"), phi, kappa, ...){
  type <- match.arg(type)

  if(type == "AR"){
    message("AR(1) with rho=0.95")
    rho <- 0.95
    cor <- rho^(abs(outer(seq_len(n_points), seq_len(n_points), "-")))
    sd <- diag(cor)
  } else if (type == "iid") {
    message("iid error")
    sd <- rep(1, n_points)
    cor <- diag(sd)
  } else if(type %in% c("exponential", "matern")){
    s <- seq_len(n_points)
    d <- sapply(s, function(i){abs(s-i)}) / n_points
    diag(d) <- 0
    if(type == "exponential"){
      if(missing(phi)){ phi <- 0.01 }
      if(missing(kappa)){ kappa <- 0.5 }
    } else {
      if(missing(phi)){ phi <- 0.01 }
      if(missing(kappa)){ kappa <- 2.5 }
    }
    message(sprintf("%s with phi=%.4f & kappa=%.4f", type, phi, kappa))
    cor <- matern(d, phi, kappa)
    sd <- diag(cor)
  } else {
    message("Using an iEEG example")
    stopifnot(n_points <= 1000)
    tfile <- file.path(tempdir(), "efcr-power-decibel-simulation.rds")
    url <- "https://github.com/dipterix/ieeg-data-examples/blob/master/efcr-power-decibel/data.rds?raw=true"
    utils::download.file(url, destfile = tfile)
    power <- readRDS(tfile)

    mm <- ncol(power)
    idx <- seq(1, mm, by = floor(mm / n_points))
    idx <- idx[seq_len(n_points)]
    power <- power[,idx, drop = FALSE]
    cov <- stats::cov(power)
    # standardize
    sd <- sqrt(diag(cov))
    cor <- t(cov / sd) / sd
  }

  eigen <- eigen(cor)
  eigen$values = eigen$values - min(eigen$values, 0)
  A <- eigen$vectors %*% diag(sqrt(eigen$values))

  gen_data <- function(n_obs, mu){
    npt <- n_points * n_obs
    df <- 10
    # r <- rnorm(npt) / sqrt(rgamma(npt, shape = df/2, scale = 2) / df)
    r <- rnorm(npt)
    t((A %*% matrix(r, nrow = n_points) + mu) * sd)
  }

  pretty_list(
    n_points = n_points,
    type = type,
    sd = sd,
    cor = cor,
    gen_data = gen_data
  )

}

#' @export
simulation_data <- function(n_points, n_obs = 100, mu_type = c("sine", "step", "custom"),
                            cov_type = c("AR", "exponential", "matern", "any", "iid"),
                            custom = NULL, corrupt = FALSE, ...){
  mu_type <- match.arg(mu_type)
  cov_type <- match.arg(cov_type)
  dat <- gen_error(n_points, type = cov_type, ...)

  # generate data using cov of power

  if(mu_type == "custom"){
    mu <- custom(seq_len(n_points))
  } else {
    aa <- c(0.7, 0.8, 0.5, 0.8, 1, 0)
    # aa <- c(0.7, 0, 0, 0, 0, 0)
    if(mu_type == 'step'){
      sin <- function(x){rep(1, length(x))}
    }
    # aa <- c(3, 0.8, 0.5, 0.8, 1, 0)
    mu <- c(rep(0, 50), (1+sin(seq(-pi/2, pi*3/2, length.out = 30))) * aa[1],
            rep(0, 70), (1+sin(seq(-pi/2, pi*3/2, length.out = 46))) * aa[2],
            rep(0, 5), (1+sin(seq(-pi/2, pi*3/2, length.out = 60))) * aa[3],
            rep(0, 70 + 65), (1+sin(seq(-pi/2, pi*3/2, length.out = 100))) * aa[4],
            rep(0, 70), (1+sin(seq(-pi/2, pi*3/2, length.out = 333))) * aa[5],
            rep(0, 50), (1+sin(seq(-pi/2, pi*3/2, length.out = 3))) * aa[6],
            rep(0, 48)
    )
    if(mu_type == 'step'){
      mu <- mu * 0.5
    }
    idx <- seq(1, length(mu), by = floor(length(mu) / n_points))
    idx <- idx[seq_len(n_points)]
    mu <- mu[idx]
  }

  sig <- which(mu != 0)

  gen_data <- function(snr = 1){
    if(corrupt){
      if(snr <= 0 || snr > 1){
        stop("$gen_data(corrupt=TRUE) only accept snr within (0,1]")
      }
      d <- dat$gen_data(n_obs = n_obs, 0)
      mu <- matrix(mu, nrow = n_obs, ncol = n_points, byrow = TRUE)
      idx <- sample(n_obs * n_points, ceiling(n_obs * n_points * (1-snr)))
      mu[idx] <- 0
      # d[-idx] <- 0
      d <- mu + d
    } else {
      d <- dat$gen_data(n_obs = n_obs, 0) #* (1-snr)
      mu <- matrix(mu * snr, nrow = n_obs, ncol = n_points, byrow = TRUE)
      d <- mu + d
    }
    d
  }

  pretty_list(
    signal_type = mu_type,
    simulation_type = dat$type,
    mu = mu * dat$sd, support = sig,
    n_obs = n_obs, n_points = dat$n_points,
    sd = dat$sd, cor = dat$cor,
    gen_data = gen_data,
    .class = "focr_simulation_data"
  )
}

plot_clean <- function (xlim = c(0,1), ylim = c(0,1), x = 1, y = 1, type = "n",
                        xlab = "", ylab = "", ...) {
  plot(x, y, type = type, axes = F, ylab = ylab, xlab = xlab,
       xlim = range(xlim), ylim = range(ylim), ...)
}

#' @export
plot.focr_simulation_data <- function(x, snr = 1, which = 1,
                                      alpha = 0.05, ymin = -0.3, ymax = "auto",
                                      ylab = "", main = "", samples = TRUE,
                                      data, add = FALSE, ...){

  if(missing(data)){
    tmp <- x$gen_data(snr = snr)
  } else {
    tmp <- data
  }

  muhat <- colMeans(tmp)
  sds <- x$sd / sqrt(x$n_obs)

  if(which == 1){
    if(identical(ymax, "auto")){
      ymax <- max(muhat)
    }

    if(add){
      points(x$mu * snr, type = "l", col = 'red', ylab = ylab, lwd = 2, ...)
    } else {
      plot(x$mu * snr, type = "l", ylim = c(ymin, ymax), col = 'red',
           ylab = ylab, main = main, lwd = 2, ...)
    }


    if(samples){
      sds <- sds * 2
      polygon(x = c(seq_len(x$n_points), rev(seq_len(x$n_points))),
              y = c(muhat - sds, rev(muhat + sds)), border = NA, col = '#1874CD66')
      lines(muhat, type = 'l', col = "#1874CD")
    }
  } else {
    # 2 * (1 - pnorm(abs(sqrt(x$n_obs) * colMeans(y))))
    pvals <- 2 * (1 - pnorm(abs(muhat) / sds))
    if(identical(ymax, "auto")){
      ymax <- alpha * 3
    }
    if(!add){
      plot_clean(c(1, length(pvals)), ylim = c(ymin, ymax), ylab = 'P-value', main = main, ...)
    }

    lines(pvals, type = 'l', col = 'grey')
    axis(1, pretty(c(1, x$n_points)))
    axis(2, c(0, alpha, 1), labels = c('0', alpha, "1"), las = 1)
    abline(h = alpha, lty = 3, col = 'grey')
  }


  # support
  v <- sprintf("c(%s)", deparse_svec(x$support, connect = ",", concatenate = TRUE))
  v <- eval(parse(text = v))
  abline(v = v, lty = 2, col = 'grey40')
  # abseg(x$support, y = alpha, clear = FALSE, pch = ".", lwd = 0.5, col = "green", lty = 1)
}

