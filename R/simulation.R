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

gen_error_2D <- function (dim, sd, nobs, rho = 0.75, ...) {
  ans <- sapply(seq_len(nobs), function(z){
    start <- array(rnorm(prod(dim), 0, sd), dim = dim)
    noise <- array(0, dim = dim)
    noise[1, 1] <- start[1, 1]

    for (i in 2:dim[1]) {
      noise[i, 1] <- rho * noise[(i - 1), 1] + sqrt(1 - rho ^ 2) * start[i, 1]
    }
    for (j in 2:dim[2]) {
      noise[, j] <- rho * noise[, (j - 1)] + sqrt(1 - rho^2) * start[, j]
    }
    for (i in 2:dim[1]) {
      noise[i, 2:dim[2]] <- rho * noise[(i - 1), 2:dim[2]] + sqrt(1 - rho^2) * noise[i, 2:dim[2]]
    }
    noise
  })
  # obs x npoints
  return(t(ans))
}


#' @name simulation
#' @title Generate simulation data used in the paper and vignettes
#' @description \code{simulation_data_1D} generates one-dimensional data,
#' \code{simulation_data_2D} generates two-dimensional image data.
#' @param n_points number of hypotheses, less equal than 1000
#' @param n_obs number of observations to generate
#' @param mu_type underlying function type, choices are 'sine', 'step', and
#' 'custom'
#' @param cov_type covariance type, choices are \code{'AR'},
#' \code{'exponential'}, \code{'matern'},
#' \code{'any'} (arbitrary dependence), and \code{'iid'} (independent)
#' @param custom function to generate mean function if \code{mu_type="custom"}
#' @param corrupt whether to use corruption model to generate data instead of
#' addition model
#' @param generator list items returned by \code{simulation_data_1D}
#' @param snr positive number to control the signal-to-noise ratio
#' @param block_sizes integer vectors to control the block size;
#' see \code{\link{focr}}
#' @param initial_filter,bandwidth used by \code{LAWS} and \code{SABHA}; see
#' \code{\link{focr}}, \code{\link{fdr-controls}}
#' @param alpha 'FOCR' and 'FDR' level; default is \code{0.05}
#' @param ... passed to internal function; see 'Details'
#' @details
#' When \code{mu_type} is 'custom', parameter \code{custom} needs to be a
#' function that takes \code{1:n_points} as input and spit out the underlying
#' mean function.
#'
#' When \code{cov_type} is \code{'AR'}, the auto-correlation of adjacent column
#' will be 0.9; The \code{cov_type="exponential"} and \code{cov_type="matern"}
#' share the same \code{phi=0.01} (range parameter) but different \code{kappa}
#' (smoothness parameter). If you wish to change the range or smoothness
#' parameter, pass \code{phi} and \code{kappa} to \code{...} (see 'Examples').
#' For \code{'any'} \code{cov_type}, the underlying covariance will be
#' generated from a real data with arbitrary dependence. For \code{'iid'}
#' \code{cov_type}, the errors are independent standard normal distributed.
#'
#' By default, \code{corrupt} is false, then the generated data is an addition
#' of underlying signal plus random noises. When \code{corrupt} is true, the
#' underlying signal will be randomly corrupted for each observation. The amount
#' of corrupted points follows a binomial distribution.
#'
#' @examples
#'
#'
#' # -------------------- Basic usage ------------------------
#' generator <- simulation_data_1D(200, cov_type = 'matern')
#'
#' # generate date with signal-to-noise ratio = 0.4
#' data <- generator$gen_data(snr = 0.4)
#'
#' # Data is n_obs x n_points matrix
#' dim(data)
#'
#' image(cor(data), main = 'Matern correlation')
#'
#' # -------------------- Change Matern parameters ------------------------
#' # Control kappa/phi here
#' generator <- simulation_data_1D(200, cov_type = 'matern', kappa = 10)
#' data <- generator$gen_data(snr = 0.4)
#'
#' image(cor(data), main = 'Smooth Matern correlation with kappa=10')
#'
#' # -------------------- 2D data ------------------------
#' # generate a 2D triangle data
#' generator <- simulation_data_2D(cov_type = 'AR')
#' data <- generator$gen_data(snr = 0.6)
#'
#' par(mfrow = c(1,2))
#' image(matrix(colMeans(data), 32), main = '2D sample mean')
#' image(matrix(generator$mu, 32), main = '2D underlying mean')
#'
#' # -------------------- Simulation used by paper ----------
#' # might take a while to run
#' if(interactive()){
#'   generator <- simulation_data_1D(n_points = 1000, cov_type = 'AR')
#'   set.seed(1000)
#'   sim <- simulate_1D(generator, snr = 0.3)
#'   plot(sim)
#' }
#'
#'
#' @export
simulation_data_1D <- function(
  n_points, n_obs = 100, mu_type = c("step", "sine", "custom"),
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
    .class = "focr_simulation_data_1d"
  )
}

#' @rdname simulation
#' @export
simulation_data_2D <- function(
  n_obs = 100, corrupt = FALSE, rho = 0.3, ...
){
  # download.file(url, 'figure/tmp.png')
  # mu <- png::readPNG('inst/triangle.png')
  # write.table(1-mu[,,3], file = 'inst/triangle.txt', sep = '\t',
  #             row.names = FALSE, col.names = FALSE)
  # mu <- read.table(system.file('triangle.txt', package = 'focr'), sep = '\t', header = FALSE)
  mu <- 1 - t(png::readPNG(system.file('cards-sm.png', package = 'focr'))[,,1])
  mu <- as.matrix(mu)
  mu[mu < 0.1] <- 0
  mu[mu > 0.1 & mu < 0.6] <- 0.5
  dim <- dim(mu)
  n_points <- prod(dim)
  # dat <- gen_error(n_points, type = cov_type, ...)
  dat <- list(
    n_points = n_points,
    sd = 1,
    gen_data = function(n_obs, mu){
      gen_error_2D(dim, 1, n_obs, rho = rho)
    }
  )
  # generate data using cov of power

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
    simulation_type = "2D Spatial",
    mu = mu * dat$sd, support = sig,
    n_obs = n_obs, n_points = dat$n_points,
    sd = dat$sd, cor = dat$cor,
    gen_data = gen_data,
    .class = "focr_simulation_data_2d"
  )
}

#' @rdname simulation
#' @export
simulate_1D <- function(generator, snr = 0.2, block_sizes = 21,
                     initial_filter = 0.9, alpha = 0.05, bandwidth = 90){
  y <- generator$gen_data(snr = snr)
  FDPs <- list()
  PWRs <- list()
  FOCPs <- list()
  rejs <- list()
  purity_tau <- initial_filter
  support <- generator$support
  theta <- rep(0, generator$n_points)
  theta[generator$support] <- 1
  # Marginal p-values
  x <- sqrt(generator$n_obs) * colMeans(y) / apply(y, 2, sd)
  pv <- 2 * pnorm(-abs(x), 0, 1)
  ## BH
  bh.res<-BH(pv, alpha)
  FDPs[["BH"]] <- fdp(bh.res$rejs, support)
  PWRs[["BH"]] <- pwr(bh.res$rejs, support)
  rejs[["BH"]] <- bh.res$rejs
  laws_bw <- bandwidth
  ## laws
  law.dd.res <- LAWS(pv, laws_bw, dimension = "one", alpha = alpha,
                     initial_filter = purity_tau)
  FDPs[["LAWS"]] <- fdp(law.dd.res$rejs, support)
  PWRs[["LAWS"]] <- pwr(law.dd.res$rejs, support)
  rejs[["LAWS"]] <- law.dd.res$rejs
  ## SABHA
  sab.dd.res <- SABHA(pv, laws_bw, dimension = "one", alpha = alpha,
                      initial_filter = purity_tau)
  FDPs[["SABHA"]] <- fdp(sab.dd.res$rejs, support)
  PWRs[["SABHA"]] <- pwr(sab.dd.res$rejs, support)
  rejs[["SABHA"]] <- sab.dd.res$rejs
  for(blk_s in block_sizes){
    w <- blk_s
    ## focr-Overlap-XXXX
    focr.res <- focr(y, blk_s, alpha = alpha, fdr_method = "BH")
    # XXXX = none (no post-selection)
    tmp <- sapply(focr.res$rej_blocks, function(ctr){
      !any(focr.res$blocks(ctr) %in% generator$support)
    })
    focp <- ifelse(length(tmp), mean(tmp), 0)
    FOCPs[[sprintf("FOCR-O%d-RAW", w)]] <- focp
    FDPs[[sprintf("FOCR-O%d-RAW", w)]] <- fdp(focr.res$rej_hypotheses, support)
    PWRs[[sprintf("FOCR-O%d-RAW", w)]] <- pwr(focr.res$rej_hypotheses, support)
    rejs[[sprintf("FOCR-O%d-RAW", w)]] <- focr.res$rej_hypotheses
    # XXXX = BH
    FDPs[[sprintf("FOCR-O%d-BH", w)]] <- fdp(focr.res$post_selection$rejs, support)
    PWRs[[sprintf("FOCR-O%d-BH", w)]] <- pwr(focr.res$post_selection$rejs, support)
    rejs[[sprintf("FOCR-O%d-BH", w)]] <- focr.res$post_selection$rejs
    # XXXX = LAWS
    focr.laws.res <- focr(y, blk_s, alpha = alpha, fdr_method = "LAWS",
                          initial_filter = purity_tau)
    FDPs[[sprintf("FOCR-O%d-LAWS", w)]] <- fdp(focr.laws.res$post_selection$rejs, support)
    PWRs[[sprintf("FOCR-O%d-LAWS", w)]] <- pwr(focr.laws.res$post_selection$rejs, support)
    rejs[[sprintf("FOCR-O%d-LAWS", w)]] <- focr.laws.res$post_selection$rejs
    # XXXX = SABHA
    focr.sabha.res <- focr(y, blk_s, alpha = alpha, fdr_method = "SABHA",
                           initial_filter = purity_tau)
    FDPs[[sprintf("FOCR-O%d-SABHA", w)]] <- fdp(focr.sabha.res$post_selection$rejs, support)
    PWRs[[sprintf("FOCR-O%d-SABHA", w)]] <- pwr(focr.sabha.res$post_selection$rejs, support)
    rejs[[sprintf("FOCR-O%d-SABHA", w)]] <- focr.sabha.res$post_selection$rejs
    ## FOCR-Disjoint-XXXX
    blocks <- function(ii){
      floor((ii - 1) / w) * w + (1:w)
    }
    focr.res <- focr(y, alpha = alpha, fdr_method = 'BH', blocks = blocks)
    # XXXX = none (no post-selection)
    tmp <- sapply(focr.res$rej_blocks, function(ctr){
      !any(blocks(ctr) %in% generator$support)
    })
    focp <- ifelse(length(tmp), mean(tmp), 0)
    FOCPs[[sprintf("FOCR-D%d-RAW", w)]] <- focp
    FDPs[[sprintf("FOCR-D%d-RAW", w)]] <- fdp(focr.res$rej_hypotheses, support)
    PWRs[[sprintf("FOCR-D%d-RAW", w)]] <- pwr(focr.res$rej_hypotheses, support)
    rejs[[sprintf("FOCR-D%d-RAW", w)]] <- focr.res$rej_hypotheses
    # XXXX = BH
    FDPs[[sprintf("FOCR-D%d-BH", w)]] <- fdp(focr.res$post_selection$rejs, support)
    PWRs[[sprintf("FOCR-D%d-BH", w)]] <- pwr(focr.res$post_selection$rejs, support)
    rejs[[sprintf("FOCR-D%d-BH", w)]] <- focr.res$post_selection$rejs
    # XXXX = LAWS
    focr.laws.res <- focr(y, alpha = alpha, fdr_method = 'LAWS', blocks = blocks,
                          initial_filter = purity_tau)
    FDPs[[sprintf("FOCR-D%d-LAWS", w)]] <- fdp(focr.laws.res$post_selection$rejs, support)
    PWRs[[sprintf("FOCR-D%d-LAWS", w)]] <- pwr(focr.laws.res$post_selection$rejs, support)
    rejs[[sprintf("FOCR-D%d-LAWS", w)]] <- focr.laws.res$post_selection$rejs
    # XXXX = SABHA
    focr.sabha.res <- focr(y, alpha = alpha, fdr_method = 'SABHA', blocks = blocks,
                           initial_filter = purity_tau)
    FDPs[[sprintf("FOCR-D%d-SABHA", w)]] <- fdp(focr.sabha.res$post_selection$rejs, support)
    PWRs[[sprintf("FOCR-D%d-SABHA", w)]] <- pwr(focr.sabha.res$post_selection$rejs, support)
    rejs[[sprintf("FOCR-D%d-SABHA", w)]] <- focr.sabha.res$post_selection$rejs
  }
  pretty_list(
    FDPs = FDPs,
    PWRs = PWRs,
    FOCPs = FOCPs,
    rejs = rejs,
    generator = generator,
    snr = snr,
    block_sizes = block_sizes,
    initial_filter = initial_filter,
    .class = "focr_simulate_results_1d"
  )
}



plot_clean <- function (xlim = c(0,1), ylim = c(0,1), x = 1, y = 1, type = "n",
                        xlab = "", ylab = "", ...) {
  plot(x, y, type = type, axes = F, ylab = ylab, xlab = xlab,
       xlim = range(xlim), ylim = range(ylim), ...)
}

#' @export
plot.focr_simulate_results_1d <- function(x, ...){
  generator <- x$generator
  res <- x
  snr <- x$snr
  mu_type <- generator$signal_type
  cov_type <- generator$simulation_type
  save_file <- FALSE

  cex.text <- 0.5
  xtext <- -75;
  xfocr <- 1075; xfdr <- 1125; xpwr <- 1200
  yline <- -0.2
  sp_sm <- 0.02; sp_md <- 0.04; sp_lg <- 0.08
  ymin <- -1.06 + (sp_lg + sp_sm * 6 + sp_md) * (3 - length(x$block_sizes))

  title <- sprintf(
    "bquote('Rejections under '~r~'=%.2f, '~.('%s')~mu[s]~', and %s'~%s~' '~(alpha~'=0.05'))",
    snr,
    ifelse(mu_type=='step', "Step", "Smooth"),
    list("iid" = "iid ",
         "AR" = "AR(1) ",
         "matern" = "Matern correlation",
         "Any"="arbitrary dependence")[[cov_type]],
    ifelse(cov_type %in% c('iid', 'AR'), "epsilon[s]", "''")
  )

  plot.default(x = c(-100,1200), y = c(ymin,-0.2),#ylim = c(-1.06,-0.2),
       main = eval(parse(text=title)),
       xlab = "", ylab = "", type = 'n', axes = FALSE)
  mtext("Location", 1, line = 1.6)

  plot.focr_simulation_data_1d(generator, which = 2, snr = snr, alpha = NA, add = TRUE)

  abseg(generator$support, yline, clear = TRUE, col = 'gray')
  text(x = 65.5, y = yline, "1", cex = cex.text)
  text(x = 173.5, y = yline, "2", cex = cex.text)
  text(x = 231.5, y = yline, "3", cex = cex.text)
  text(x = 446.5, y = yline, "4", cex = cex.text)
  text(x = 733, y = yline, "5", cex = cex.text)
  text(x = xtext, y = yline, "Underlying", cex = cex.text)
  text(x = xfocr, y = yline, "FOCP", cex = cex.text)
  text(x = xfdr, y = yline, "FDP", cex = cex.text)
  text(x = xpwr, y = yline, "POWER", cex = cex.text)
  yline <- yline-sp_lg


  col <- 'dodgerblue3'
  abseg(res$rejs$BH, yline, clear = TRUE, col = col);
  text(x = xtext, y = yline, "BH", cex = cex.text, col = col)
  text(x = xfdr, y = yline, sprintf("%.1f%%", res$FDPs[["BH"]] * 100), cex = cex.text, col = col)
  text(x = xpwr, y = yline, sprintf("%.1f%%", res$PWRs[["BH"]] * 100), cex = cex.text, col = col)
  yline <- yline-sp_sm

  nm <- "SABHA"; col = 'brown'
  abseg(res$rejs[[nm]], yline, clear = TRUE, col = col)
  text(x = xtext, y = yline, nm, cex = cex.text, col = col)
  text(x = xfdr, y = yline, sprintf("%.1f%%", res$FDPs[[nm]] * 100), cex = cex.text, col = col)
  text(x = xpwr, y = yline, sprintf("%.1f%%", res$PWRs[[nm]] * 100), cex = cex.text, col = col)
  yline <- yline-sp_sm


  nm <- "LAWS"; col = 'purple3'
  abseg(res$rejs[[nm]], yline, clear = TRUE, col = col)
  text(x = xtext, y = yline, nm, cex = cex.text, col = col)
  text(x = xfdr, y = yline, sprintf("%.1f%%", res$FDPs[[nm]] * 100), cex = cex.text, col = col)
  text(x = xpwr, y = yline, sprintf("%.1f%%", res$PWRs[[nm]] * 100), cex = cex.text, col = col)
  yline <- yline-sp_lg

  for(w in x$block_sizes){
    radius <- floor(w/2)

    # --------------------------- FOCR overlap
    nm <- sprintf("FOCR-O%d-RAW", w); col = 'orange'
    abseg(res$rejs[[nm]], yline, clear = TRUE, col = col)
    text(x = xtext, y = yline, nm, cex = cex.text, col = col)
    text(x = xfocr, y = yline, sprintf("%.1f%%", res$FOCPs[[nm]] * 100), cex = cex.text, col = col)
    text(x = xfdr, y = yline, sprintf("%.1f%%", res$FDPs[[nm]] * 100), cex = cex.text, col = col)
    text(x = xpwr, y = yline, sprintf("%.1f%%", res$PWRs[[nm]] * 100), cex = cex.text, col = col)
    yline <- yline-sp_sm

    nm <- sprintf("FOCR-O%d-BH", w); col <- 'dodgerblue3'
    abseg(res$rejs[[nm]], yline, clear = TRUE, col = col)
    text(x = xtext, y = yline, nm, cex = cex.text, col = col)
    text(x = xfdr, y = yline, sprintf("%.1f%%", res$FDPs[[nm]] * 100), cex = cex.text, col = col)
    text(x = xpwr, y = yline, sprintf("%.1f%%", res$PWRs[[nm]] * 100), cex = cex.text, col = col)
    yline <- yline-sp_sm

    nm <- sprintf("FOCR-O%d-SABHA", w); col <- 'brown'
    abseg(res$rejs[[nm]], yline, clear = TRUE, col = col)
    text(x = xtext, y = yline, nm, cex = cex.text, col = col)
    text(x = xfdr, y = yline, sprintf("%.1f%%", res$FDPs[[nm]] * 100), cex = cex.text, col = col)
    text(x = xpwr, y = yline, sprintf("%.1f%%", res$PWRs[[nm]] * 100), cex = cex.text, col = col)
    yline <- yline-sp_sm

    nm <- sprintf("FOCR-O%d-LAWS", w); col <- 'purple3'
    abseg(res$rejs[[nm]], yline, clear = TRUE, col = col)
    text(x = xtext, y = yline, nm, cex = cex.text, col = col)
    text(x = xfdr, y = yline, sprintf("%.1f%%", res$FDPs[[nm]] * 100), cex = cex.text, col = col)
    text(x = xpwr, y = yline, sprintf("%.1f%%", res$PWRs[[nm]] * 100), cex = cex.text, col = col)
    yline <- yline-sp_md

    # --------------------------- FOCR disjoint
    nm <- sprintf("FOCR-D%d-RAW", w); col = 'orange'
    abseg(res$rejs[[nm]], yline, clear = TRUE, col = col)
    text(x = xtext, y = yline, nm, cex = cex.text, col = col)
    text(x = xfocr, y = yline, sprintf("%.1f%%", res$FOCPs[[nm]] * 100), cex = cex.text, col = col)
    text(x = xfdr, y = yline, sprintf("%.1f%%", res$FDPs[[nm]] * 100), cex = cex.text, col = col)
    text(x = xpwr, y = yline, sprintf("%.1f%%", res$PWRs[[nm]] * 100), cex = cex.text, col = col)
    yline <- yline-sp_sm

    nm <- sprintf("FOCR-D%d-BH", w); col <- 'dodgerblue3'
    abseg(res$rejs[[nm]], yline, clear = TRUE, col = col)
    text(x = xtext, y = yline, nm, cex = cex.text, col = col)
    text(x = xfdr, y = yline, sprintf("%.1f%%", res$FDPs[[nm]] * 100), cex = cex.text, col = col)
    text(x = xpwr, y = yline, sprintf("%.1f%%", res$PWRs[[nm]] * 100), cex = cex.text, col = col)
    yline <- yline-sp_sm

    nm <- sprintf("FOCR-D%d-SABHA", w); col <- 'brown'
    abseg(res$rejs[[nm]], yline, clear = TRUE, col = col)
    text(x = xtext, y = yline, nm, cex = cex.text, col = col)
    text(x = xfdr, y = yline, sprintf("%.1f%%", res$FDPs[[nm]] * 100), cex = cex.text, col = col)
    text(x = xpwr, y = yline, sprintf("%.1f%%", res$PWRs[[nm]] * 100), cex = cex.text, col = col)
    yline <- yline-sp_sm

    nm <- sprintf("FOCR-D%d-LAWS", w); col <- 'purple3'
    abseg(res$rejs[[nm]], yline, clear = TRUE, col = col)
    text(x = xtext, y = yline, nm, cex = cex.text, col = col)
    text(x = xfdr, y = yline, sprintf("%.1f%%", res$FDPs[[nm]] * 100), cex = cex.text, col = col)
    text(x = xpwr, y = yline, sprintf("%.1f%%", res$PWRs[[nm]] * 100), cex = cex.text, col = col)
    yline <- yline-sp_lg
  }



}

#' @export
plot.focr_simulation_data_1d <- function(x, snr = 1, which = 1,
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

