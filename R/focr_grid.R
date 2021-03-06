#' @name focr
#' @title False overlapped-cluster rate (FOCR) control procedures
#' @param data a n-by-p numerical matrix (no missing values) with \code{n} to
#' be the total number of observations and \code{p} is the total number of
#' hypotheses
#' @param data_corr the correlation matrix of \code{data}. If missing, then
#' the correlation will be calculated empirically
#' @param scale numerical vector of standard deviations by column; default is
#' missing (use empirical standard deviation)
#' @param blocks a list of indices or a function that returns indices
#' @param nblocks the total number of blocks, used when \code{blocks} is a
#' function
#' @param block_size block size of sliding window; used by \code{focr}.
#' @param mu the mean function value to compare with; see 'Details'
#' @param alpha FOCR level for stage-I, and FDR level for stage-II
#' @param verbose whether to print out information; default is false
#' @param side test type, \code{'two'} if alternative hypotheses are two-sided,
#' and \code{'left'} or \code{'right'} if one-sided.
#' @param fdr_method characters or function of post-selection FDR control
#' procedures. Built-in choices are \code{"BH"}, \code{"BY"}, \code{"SABHA"},
#' and \code{"LAWS"}. See vignette for details, see also
#' \code{\link{fdr-controls}}.
#' @param dimension the dimension information of input hypotheses. For
#' \code{\link{LAWS}} and \code{\link{LAWS}}, current implementation only
#' supports 1-3 dimensions.
#' @param bandwidth used by \code{\link{LAWS}} and \code{\link{LAWS}} as
#' smoothing parameters to estimate the underlying sparsity level. Default
#' is half of \code{block_size}. If \code{block_size} is missing,
#' \code{bandwidth} must be specified.
#' @param initial_filter used by \code{\link{LAWS}} and \code{\link{LAWS}} as
#' initial filters (purity) to remove large p-values
#' @param distance_measure distance measure used to form blocks; see 'Details'.
#' @param ... passed to \code{focr_initial} and \code{fdr_method}
#' @details
#' The function \code{focr} and \code{focr_initial} control the type-I error
#' for multiple testing problems with topological constraints:
#' \deqn{
#' H_{0}(s):f(s)=\mu(s), H_{1}(s):f(s)\neq \mu(s)
#' }{H_{0}(s):f(s)=\mu(s), H_{1}(s):f(s)\neq \mu(s)}
#'
#' The type-I error control procedure has two stages. In the first stage,
#' the FOCR is controlled at block (overlapped-cluster) level. This step is
#' to find regions of interests that respect the topological constraints. The
#' second stage further inspects the hypotheses rejected by the first stage.
#' During this stage, conditional p-values will be calculated in a
#' post-selection fashion. FDR control methods are further applied to these
#' conditional p-values to select significant hypotheses at individual level.
#'
#' Function \eqn{\mu(s)} is specified in \code{mu}. By default the alternative
#' hypothesis is two-sided. For one-sided tests, please change the parameter
#' \code{side} to either \code{"left"} or \code{"right"}.
#'
#' The function \code{focr_initial} controls the FOCR on the block level
#' (stage-I), and calculates the conditional p-values. The function \code{focr}
#' uses \code{focr_initial}, providing default block settings and built-in
#' post-selection inference on conditional p-values.
#'
#' By default, \code{focr} uses sliding window as blocks. Each block is a ball
#' with distance between the boundary and center point given
#' by \code{block_size/2}. The distance measure is specified by
#' \code{distance_measure}. The choices are \code{"euclidean"}, \code{"lmax"},
#' and \code{"manhattan"}. This default settings should work in many spatial
#' or temporal situations. However, in case the blocks are to be customized,
#' please specify \code{blocks} manually. The argument \code{blocks} can be
#' either a list of hypothesis indices, or a function that returns ones given
#' by locations of hypotheses. See 'vignette'
#' \href{../doc/false-overlapped-cluster-rate.html}{
#' \code{vignette('false-overlapped-cluster-rate', package='focr')}}.
#'
#' @return A list of results
#' \describe{
#' \item{\code{method}}{method name}
#' \item{\code{alpha}}{level of significance: FOCR in the stage-I and FDR in
#' the stage-II}
#' \item{\code{side}}{passed from input}
#' \item{\code{blocks}}{function that returns indices of blocks}
#' \item{\code{nblocks}}{number of total blocks}
#' \item{\code{rej_blocks}}{blocks being rejected}
#' \item{\code{rej_hypotheses}}{individual hypotheses rejected in the first
#' stage}
#' \item{\code{tau}}{p-value cutoff in the first stage}
#' \item{\code{cond_pvals}}{conditional p-values in the stage-II}
#' \item{\code{uncond_pvals}}{unconditional p-values}
#' \item{\code{details}}{details of initial rejections}
#' \item{\code{stats}}{block-level test statistics and p-values}
#' }
#' The following additional items are \code{focr} only.
#' \describe{
#' \item{\code{post_selection}}{a list returned by FDR controlling methods,
#' see also \code{\link{fdr-controls}}}
#' \item{\code{fdr_method}}{function used to control the FDR in stage-II}
#' \item{\code{block_size}}{block size if specified, passed from input}
#' }
#'
#' @examples
#'
#'
#' library(focr)
#' set.seed(100)
#' generator <- simulation_data_1D(n_points = 1000, mu_type = 'step',
#'                              cov_type = 'AR')
#' data <- generator$gen_data(snr = 0.34)
#' plot(generator, data = data, snr = 0.34)
#'
#' # -------------------- Basic usage -------------------------
#' # FOCR-BH procedure
#' res <- focr(data = data, block_size = 41,
#'             alpha = 0.05, fdr_method = 'BH')
#'
#' # False discovery proportion
#' fdp <- fdp(res$post_selection$rejs, generator$support)
#' fdp
#'
#' # Statistical power
#' power <- pwr(res$post_selection$rejs, generator$support)
#' power
#'
#' # Visualize
#' plot(generator$mu, type = 'l', col = 'red', ylim = c(-.5,1.5),
#'      main = sprintf('FOCR-BH, FDP=%.1f%%, Power=%.1f%%',
#'                     fdp*100, power * 100))
#' lines(res$cond_pvals, col = 'gray')
#' abseg(res$rej_hypotheses, y = -0.3, col = 'orange3', lwd = 2)
#' abseg(res$post_selection$rejs, y = -0.5, col = 'blue', lwd = 2)
#' legend('topleft', c("Underlying signal", "Conditional p-values",
#'                     "FOCR initial clusters", "FOCR-BH final rejections"),
#'        col = c('red', 'orange3', 'blue'), lty = 1, cex = 0.7)
#'
#' # ------------------------- Change FDR methods --------------------
#' # FOCR-LAWS
#' res <- focr(data = data, block_size = 41,
#'             alpha = 0.05, fdr_method = 'LAWS',
#'             initial_filter = 0.5)
#' fdp <- fdp(res$post_selection$rejs, generator$support)
#' fdp
#' power <- pwr(res$post_selection$rejs, generator$support)
#' power
#'
#' # Visualize
#' plot(generator$mu, type = 'l', col = 'red', ylim = c(-.5,1.5),
#'      main = sprintf('FOCR-LAWS, FDP=%.1f%%, Power=%.1f%%',
#'                     fdp*100, power * 100))
#' lines(res$cond_pvals, col = 'gray')
#' abseg(res$rej_hypotheses, y = -0.3, col = 'orange3', lwd = 2)
#' abseg(res$post_selection$rejs, y = -0.5, col = 'blue', lwd = 2)
#' legend('topleft', c("Underlying signal", "Conditional p-values",
#'                     "FOCR initial clusters", "FOCR-LAWS final rejections"),
#'        col = c('red', 'orange3', 'blue'), lty = 1, cex = 0.7)
#'
#' # ------------------------- Customized blocks --------------------
#'
#' # The following example uses disjoint blocks; each block has length of 40
#' res <- focr(data = data, alpha = 0.05, fdr_method = 'LAWS',
#'             initial_filter = 0.5, blocks = function(index){
#'               # Disjoint blocks with size 40
#'               floor((index -1)/40) * 40 + seq_len(40)
#'             }, bandwidth = 20)
#'
#'
#' # Compared to overlapped blocks, disjoint blocks are less powerful
#' # However, if this might be useful provided the underlying topological
#' # structure is disjoint
#' fdp <- fdp(res$post_selection$rejs, generator$support)
#' fdp
#' power <- pwr(res$post_selection$rejs, generator$support)
#' power
#'
#'
NULL

#' @rdname focr
#' @export
focr_initial <- function(data, data_corr, scale, blocks, nblocks = ncol(data),
                         mu = 0, alpha = 0.05, verbose = FALSE,
                         side = c('two', 'left', 'right'), ...){
  # debug
  debug <- function(...){ debug_verbose(..., verbose = verbose) }
  alpha_focr <- alpha
  side <- match.arg(side)
  if(side != 'two'){ alpha_focr <- alpha * 2 }

  # # prepare methods
  # if(is.character(fdr_method)){
  #   get0(fdr_method, mode = 'function')
  # }
  # if(!is.function(fdr_method)){ stop("FOCR: `fdr_method` must be function name such as `BH`, or a function itself.") }

  M <- ncol(data)
  n <- nrow(data)

  # blocks
  if(is.list(blocks)){
    nblocks <- length(blocks)
    block_idx <- function(i){
      if(i > nblocks){ stop("Block number overflow") }
      re <- blocks[[i]]
      re[re >= 1 & re <= M]
    }
  } else if(is.function(blocks)){
    force(nblocks)
    block_idx <- function(i){
      if(i > nblocks){ stop("Block number overflow") }
      re <- blocks(i)
      re[re >= 1 & re <= M]
    }
  } else {
    stop("`blocks` must be a list of column indices or a function that returns indices.")
  }
  debug("Total ", nblocks, " blocks")
  # if(nblocks != M){
  #   warning(sprintf("Current FOCR implementation requires the block count equaling to number of variates (duplicated blocks are allowed). However, the total number of blocks (%d) != # of variates (%d)", nblocks, M))
  # }

  # Normalize data
  debug("Re-scale data...")
  if(length(mu) > 1 || mu != 0){
    data <- sweep(data, 2L, mu, '-', check.margin = TRUE)
  }
  if(missing(scale)){
    scale <- apply(data, 2L, stats::sd, na.rm = TRUE)
  }
  if(length(scale) == 1){
    data <- data / scale
  } else {
    data <- sweep(data, 2L, scale, '/')
  }

  # data covariance/correlation (identical since data is rescaled)
  if(missing(data_corr)){
    # data_corr <- stats::cov(data)
    slice_cor <- function(rows, cols){
      if(getThreads() <= 1){
        cov(data[,rows, drop = FALSE], data[,cols, drop = FALSE])
        # stats::cor(data[,rows, drop = FALSE], data[,cols, drop = FALSE])
      } else {
        fastcov2(data, rows, cols)
      }
    }
  } else if(is.function(data_corr)){
    slice_cor <- data_corr
  } else {
    slice_cor <- function(rows, cols){
      data_corr[rows, cols, drop = FALSE]
    }
  }

  mean <- colMeans(data)

  slice_mean <- function(idx){
    mean[idx]
  }
  # slice_cor
  # sapply <- get_sapply()
  debug("ROCR initial rejection...")
  stats <- t(sapply(seq_len(nblocks), function(ii){
    if(ii %% 1000 == 0){
      debug(sprintf("Block %d \r", ii), appendLF = FALSE)
    }
    idx <- block_idx(ii)
    m <- slice_mean(idx)
    v <- slice_cor(idx, idx)
    tval <- n * sum(m^2)
    param <- gamma_approx2(v)
    shape <- param$M / param$u
    scale <- param$u
    # pval <- pgamma(tval, shape = shape, scale = scale, lower.tail = FALSE)
    df2 <- n - 1
    pval <- pf(tval / scale / shape, shape * 2, df2 = df2, lower.tail = FALSE)
    # for post-selection inference
    # mean_c <- slice_mean(ctr)
    c(tval, pval, shape, scale, df2)
  }))
  stats <- as.data.frame(stats)
  names(stats) <- c('tval', 'pval', 'shape', 'scale', 'df2')

  # Stage-I control the FOCR
  rej_stage_1 <- BH(stats$pval, alpha = alpha_focr)
  debug(sprintf("Initial FOCR rejection (level=%.3f): %d out of %d blocks (%.1f%%)",
                alpha_focr, rej_stage_1$nrejs, nblocks,
                rej_stage_1$nrejs / nblocks * 100))

  # get actual indices rejected
  rej_stage_1_locations <- sort(unique(unlist(lapply(rej_stage_1$rejs, block_idx))))
  debug(sprintf("  %d hypotheses need further inspections", length(rej_stage_1_locations)))

  # get p-value cutoffs
  tau <- min(rej_stage_1$tau, 1)

  # calculate p-value threshold
  p0 <- rep(NA, M)

  thresholds <- qf(tau, stats$shape * 2, df2 = stats$df2, lower.tail = FALSE) * stats$scale * stats$shape
  count <- 0
  for(ii in rej_stage_1$rejs){

    # get indices
    idx <- block_idx(ii)

    # get block test-stat
    tval <- stats$tval[[ii]]

    # get mean of idx
    mean_c <- slice_mean(idx)

    # threshold
    thred <- thresholds[[ii]]

    d <- thred - tval + n * mean_c^2
    d[d < 0] <- 0

    p <- pf(d, 1, df2 = n - 1, lower.tail = FALSE)
    p0[idx] <- pmax(p0[idx], p, na.rm = TRUE)

    count <- count + 1
    if(count %% 1000 == 0){
      debug(sprintf("Post-selection: processed %d\r", count), appendLF = FALSE)
    }
  }

  p1 <- pf(mean^2 * n, 1, df2 = n - 1, lower.tail = FALSE)
  post_pvals <- p1 / p0
  post_pvals[is.na(post_pvals)] <- 1

  if(side == 'left'){
    post_pvals[mean > 0] <- 1
  } else if(side == 'right'){
    post_pvals[mean < 0] <- 1
  }

  pretty_list(
    method = "FOCR Initial Rejection (Stage-I)",
    alpha = alpha,
    side = side,
    blocks = blocks,
    nblocks = nblocks,
    rej_blocks = rej_stage_1$rejs,
    rej_hypotheses = rej_stage_1_locations,
    tau = rej_stage_1$tau,
    cond_pvals = post_pvals,
    uncond_pvals = p1,
    details = rej_stage_1,
    stats = stats
  )
}



#' @rdname focr
#' @export
focr <- function(data, block_size, alpha = 0.05, fdr_method = c('BH', 'LAWS', 'SABHA', 'BY'),
                 bandwidth = if(missing(block_size)){NA}else{block_size/2},
                 initial_filter = 0.9, dimension = NULL,
                 distance_measure = c('euclidean', 'lmax', 'manhattan'),
                 side = c('two', 'left', 'right'), verbose = FALSE,
                 blocks, ...){
  debug <- function(...){ debug_verbose(..., verbose = verbose) }
  M <- ncol(data)

  if(is.na(bandwidth) && system.file('', package = 'kedd') == ''){
    warnings('bandwidth is NA. Please install `kedd` package if you are using "LAWS" or "SABHA" method')
  }

  if(length(dimension) <= 1){
    dimension <- NULL
    ndims <- 1
  } else {
    if(prod(dimension) != M){
      stop("Argument `dimension` does not match with total number of variates.")
    }
    ndims <- length(dimension)
  }
  side <- match.arg(side)
  distance_measure <- match.arg(distance_measure)
  debug('Distance measure: ', distance_measure)
  # prepare methods
  method_name <- 'customized'
  if(is.character(fdr_method)){
    if(length(fdr_method) > 1){
      fdr_method <- match.arg(fdr_method)
    }
    method_name <- fdr_method
    debug('FDR method: ', fdr_method)
    if(fdr_method %in% c('LAWS', 'SABHA') && length(dimension) > 3){
      stop("`LAWS` and `SABHA` implementations only support 1D, 2D, 3D spatial data.")
    }
    fdr_method <- get0(fdr_method, mode = 'function')
  } else {
    debug('FDR method: <customized>')
  }
  if(!is.function(fdr_method)){ stop("FOCR: `fdr_method` must be function name such as `BH`, or a function itself.") }

  debug('Variable dimension: ', ndims)
  if(!missing(blocks)){
    window_idx <- blocks
  } else if(ndims == 1){
    hw <- floor(block_size/2)
    window_idx <- function(ii){
      seq(max(ii - hw, 1), min(ii + hw, M))
    }
  } else {

    origin <- dimension * 0
    # Make sure block_size is odd number
    hw <- floor(block_size/2)

    block_size <- hw * 2 + 1
    idx_lookup <- array(seq_len(block_size^ndims), (origin + 1) * block_size)
    # center of cube
    origin_idx <- as.vector(arrayInd(mean(idx_lookup), dim(idx_lookup)))
    idx_lookup <- t(sweep(arrayInd(idx_lookup, dim(idx_lookup)), 2L,
                          origin_idx, '-', check.margin = FALSE))
    dimension_shift <- c(1, dimension)[seq_len(ndims)]

    if(distance_measure == 'euclidean'){
      # eucl-distance
      dist <- sqrt(colSums(idx_lookup^2))
      idx_lookup <- idx_lookup[,dist <= hw, drop = FALSE]
    } else if(distance_measure == 'manhattan'){
      # eucl-distance
      dist <- colSums(abs(idx_lookup))
      idx_lookup <- idx_lookup[,dist <= hw, drop = FALSE]
    }
    debug('# elements in blocks: ', ncol(idx_lookup))
    window_idx <- function(ii){
      ijk <- as.vector(arrayInd(ii, dimension))
      ijk <- idx_lookup + ijk
      ijk <- ijk[,colSums(ijk < 1 | ijk > dimension) == 0]
      as.vector(crossprod(dimension_shift, ijk - 1) + 1)
    }

  }
  args <- list(..., bandwidth = bandwidth,
               initial_filter = initial_filter,
               verbose = verbose)

  debug('Calculating FOCR initial rejections')
  if(length(args[['nblocks']]) == 1){
    res1 <- focr_initial(data, blocks = window_idx, side = side,
                         alpha = alpha, verbose = verbose, ...)
  } else {
    res1 <- focr_initial(data, blocks = window_idx, side = side,
                         nblocks = M, alpha = alpha, verbose = verbose, ...)
  }


  # Post-selection
  cond_pvals <- res1$cond_pvals

  # re-shape pvals
  debug('Re-shape conditional p-values\t\t\t\t')
  if(ndims > 1){
    dim(cond_pvals) <- dimension
    res1$dimension <- dimension
    dim_str <- ifelse(ndims >= 3, 'three', 'two')
  } else {
    dim_str <- 'one'
    res1$dimension <- NULL
  }
  args$dimension <- dim_str


  # get formals
  debug('Prepare for post-selection FDR control')
  fml <- names(formals(fdr_method))[-1]
  if(!'...' %in% fml){
    fml <- fml[fml %in% names(args)]
    args <- c(list(cond_pvals), args[fml])
  }

  debug('Call post-selection FDR control')
  res2 <- do.call(fdr_method, args)

  debug('Finalizing, post-selection results is included as "ret$post_selection".')
  res1$post_selection <- res2
  res1$fdr_method <- fdr_method
  if(!missing(block_size)){
    res1$block_size <- block_size
  }

  res1$method <- sprintf('FOCR-%s', method_name)
  res1

}

focr_sliding_window <- function(data, radius, alpha = 0.05,
                      cov = stats::cov(data), mu = 0,
                      method = c("BY", "BH", "BH2", "custom"),
                      method2 = c("BY", "BH", "BH2"),
                      custom = NULL, debug = FALSE,
                      alpha_focr = alpha, ...){

  # debug
  if(debug){
    debug <- function(..., appendLF = TRUE){
      message(..., appendLF = appendLF)
    }
  } else {
    debug <- function(..., appendLF = TRUE){}
  }

  method <- match.arg(method)
  if(method == "custom"){
    stopifnot(is.function(custom))
    proc <- custom
  } else {
    proc <- get(method, mode = 'function')
  }
  method2 <- match.arg(method2)
  proc2 <- get(method2, mode = 'function')
  debug("Using FDR control procedure: ", method)

  # generate grid
  M <- ncol(data)
  n <- nrow(data)
  mean <- colMeans(data) - mu
  var <- diag(cov)
  sd_inv <- 1 / sqrt(var)
  cor <- t(cov * sd_inv) * sd_inv

  res <- grid_pvals(mean, var, cor, n, radius = radius, window_type = 'slide')

  post_fdr <- function(pvals){
    idx <- which(!is.na(pvals))
    pvals <- pvals[idx]
    # idx[p.adjust(pvals, method = "holm") <= alpha]
    rej4 <- proc2(pvals, alpha = alpha)
    rej4 <- idx[rej4$rejs]
    rej4
  }

  rej <- proc(res$stats$pval, alpha = alpha_focr)
  debug(sprintf("Initial rejection: %d out of %d (%.1f%%)",
                rej$nrejs, res$total_points,
                rej$nrejs / res$total_points * 100))
  # estimate eta2
  rej3 <- unique(unlist(lapply(rej$rejs, res$window_idx)))
  n_rej <- rej$nrejs

  # estimate eta2
  tau <- rej$tau
  tau <- min(tau, 1)
  # calculate p-value threshold
  pvals <- grid_post_pvals(res, rej, tau)
  # post-selection
  rej4 <- post_fdr(pvals)

  # pvals <- NULL
  # post_pvals_blocked <- grid_post_pvals_blocked(res, rej, tau)

  pretty_list(
    method = method,
    alpha = alpha,
    post_pvals = pvals,
    # post_pvals_blocked = post_pvals_blocked,
    post_selection = rej4,
    rejection = rej3,
    rejection_center = rej$rejs,
    details = res,
    .class = c("focr_sliding_window", "focr_result")
  )
}


focr_partition <- function(data, radius, alpha = 0.05,
                           cov = stats::cov(data), mu = 0,
                           method = c("BY", "BH", "BH2", "custom"),
                           method2 = c("BY", "BH", "BH2"),
                           custom = NULL, debug = FALSE,
                           alpha_focr = alpha, ...){

  # debug
  if(debug){
    debug <- function(..., appendLF = TRUE){
      message(..., appendLF = appendLF)
    }
  } else {
    debug <- function(..., appendLF = TRUE){}
  }



  method <- match.arg(method)
  if(method == "custom"){
    stopifnot(is.function(custom))
    proc <- custom
  } else {
    proc <- get(method, mode = 'function')
  }
  method2 <- match.arg(method2)
  proc2 <- get(method2, mode = 'function')
  debug("Using FDR control procedure: ", method)

  # generate grid
  M <- ncol(data)
  n <- nrow(data)
  mean <- colMeans(data) - mu
  var <- diag(cov)
  sd_inv <- 1 / sqrt(var)
  cor <- t(cov * sd_inv) * sd_inv

  res <- grid_pvals(mean, var, cor, n, radius = radius, window_type = 'partition')

  post_fdr <- function(pvals){
    idx <- which(!is.na(pvals))
    pvals <- pvals[idx]
    # idx[p.adjust(pvals, method = "holm") <= alpha]
    rej4 <- proc2(pvals, alpha = alpha)
    rej4 <- idx[rej4$rejs]
    rej4
  }

  rej <- proc(res$stats$pval, alpha = alpha_focr)

  debug(sprintf("Initial rejection: %d out of %d (%.1f%%)",
                rej$nrejs, res$total_points,
                rej$nrejs / res$total_points * 100))
  # estimate eta2
  rej3 <- unique(unlist(lapply(rej$rejs, res$window_idx)))

  # estimate eta2
  tau <- rej$tau
  tau <- min(tau, 1)
  # calculate p-value threshold
  pvals <- grid_post_pvals(res, rej, tau)
  # post-selection
  rej4 <- post_fdr(pvals)

  # pvals <- NULL
  # post_pvals_blocked <- grid_post_pvals_blocked(res, rej, tau)

  pretty_list(
    method = method,
    alpha = alpha,
    post_pvals = pvals,
    # post_pvals_blocked = post_pvals_blocked,
    post_selection = rej4,
    rejection = rej3,
    rejection_center = rej$rejs,
    details = res,
    .class = c("focr_grid", "focr_result")
  )
}


# tmp <- slide_window(ncol(data), 10)
# res <- focr(data, blocks = tmp$window_idx, nblocks = length(tmp$center), verbose = TRUE)
# rej <- BH(res$cond_pvals)
# plot(res$cond_pvals, type='l')
# lines(generator$mu, col = 'green')
# lines(res$uncond_pvals, col = 'grey')
# abseg(res$rej_hypotheses, 0, col = 'blue')
# abseg(rej$rejs, 0, col = 'red')
#
# a <- focr_sliding_window(data, 10, debug = TRUE, method = 'BH', method2 = 'BH', alpha = 0.05)
# plot(res$cond_pvals, type='l')
# lines(a$post_pvals, col = 'purple')
# abseg(a$post_selection, 0.5, col = 'blue')
# abseg(rej$rejs, 0, col = 'red')
# lines(generator$mu, col = 'green')

# generator <- simulation_data_1D(n_points = 900, n_obs = 100, mu_type = 'step')
# data <- generator$gen_data(0.3)
# res <- focr_sliding(data, 21, fdr_method = "BH", alpha = 0.05, verbose = TRUE, dimension = c(30,30))
# plot(res$cond_pvals, type='l')
# lines(generator$mu, col = 'green')
# abseg(res$post_selection$rejs, 0.48, col = 'red')
# fdp(res$post_selection$rejs, generator$support)
#
# head(res$stats); head(a$details$stats)
# res$cond_pvals - a$post_pvals
# res <- focr_sliding_1D(data, 21, fdr_method = "laws_pval", alpha = 0.05, verbose = TRUE)
# res
# res$post_selection$rejs
