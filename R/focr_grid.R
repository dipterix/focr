#' @export
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


#' @export
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

