#' @title Add indices as ling segments or points to the figures
#' @description Draws isolated points as dots, and consecutive points as line
#' segments.
#' @param x integer indices, will be used to calculate locations in
#' \code{\link[graphics]{points}} or \code{\link[graphics]{arrows}}
#' @param y y-axis locations for points or line segments
#' @param code,length,... passed to \code{\link[graphics]{arrows}}
#' @param pch point styles, passed to \code{\link[graphics]{points}}
#' @param clear whether to clear the line before drawing
#' @param lwd used to calculate the line weight and point size; see \code{'lwd'}
#' and \code{'cex'} in \code{\link[graphics]{plot}}, Section 'Details'.
#' @return Nothing.
#' @examples
#' plot(1:10)
#' abseg(c(1,2,3,5,7,8,10), y = 5, lwd = 4, col = 'red')
#' @export
abseg <- function(x, y, code = 0, pch = '.', length = 0.05, clear = FALSE, lwd = 1, ...){
  if(!length(x)){ return(invisible()) }
  x <- deparse_svec(unlist(x), concatenate = FALSE, connect = ",")
  if(clear){
    abline(h = y, col = par('bg'), lwd = lwd + 1)
  }
  for(s in x){
    if(s == ""){ break }
    s <- eval(parse(text = sprintf("c(%s)", s)))
    if(length(s) == 1){
      points(s, y, pch = pch, ..., cex = lwd)
    } else {
      arrows(x0 = s[1], y0 = y, x1 = s[2], y1 = y, code = code, length = length, lwd = lwd, ...)
    }
  }
  invisible()
}

pretty_list <- function(..., .list = NULL, .class = NULL){
  structure(c(list(...), .list), class = c(.class, "focr.pretty_list"))
}

#' @export
print.focr.pretty_list <- function(x, ...){
  cat("<list of ", length(x), " items: ",
      paste(
        sapply(names(x), function(nm){
          sprintf("\033[34m%s\033[39m [%s]", nm, mode(x[[nm]]))
        }),
        collapse = ", "
      ), ">\n"
      , sep = "")
  invisible(x)
}

slide_window <- function(n, radius){
  stopifnot(radius > 0)
  stopifnot(abs(n - round(n)) < 1e-6)
  n <- round(n)
  width <- floor(radius * 2)
  if(width %% 2 == 0){
    width <- width + 1
  }
  radius <- width / 2

  hw <- (width - 1) / 2
  idx <- function(center, b_k = TRUE){
    stopifnot(center >= 1 && center <= n)
    center <- round(center)
    if(b_k){
      return(seq(max(center - hw, 1L), min(center + hw, n)))
    } else {
      # C_k
      center
    }
  }
  # calculate eta ratio
  eta_max <- width + 1

  pretty_list(
    center = seq_len(n),
    window_idx = idx,
    eta_max = eta_max,
    radius = radius,
    width = width,
    total_points = n
  )

}


partition <- function(n, radius) {
  stopifnot(radius > 0)
  stopifnot(abs(n - round(n)) < 1e-6)
  n <- round(n)
  width <- floor(radius * 2)
  radius <- width / 2

  hw <- floor(radius)
  width <- hw * 2 + 1
  radius <- width / 2
  idx <- function(center, b_k = TRUE){
    stopifnot(center >= 1 && center <= n)
    center <- round(center)
    if(b_k){
      center <- center - ((center - 1) %% width) + hw
      return(seq(max(center - hw, 1L), min(center + hw, n)))
    } else {
      # C_k
      center
    }
  }
  # calculate eta ratio
  eta_max <- width + 1

  pretty_list(
    center = seq_len(n),
    window_idx = idx,
    eta_max = eta_max,
    radius = radius,
    width = width,
    total_points = n
  )
}

grid_pvals <- function(mean, var, cor, radius, n, sigma = FALSE, window_type = c("slide", "partition")){
  stopifnot(is.matrix(cor) || is.function(cor))
  if(is.matrix(cor) && !is.function(mean)){
    stopifnot(nrow(cor) == length(mean) && ncol(cor) == length(mean))
  }
  M <- length(mean)

  if(!is.function(window_type)){
    window_type <- match.arg(window_type)
    if(window_type == 'slide'){
      gen_block <- slide_window
    } else {
      gen_block <- partition
    }
  } else {
    gen_block <- window_type
  }

  conf <- gen_block(M, radius)

  slice_mean <- function(idx){
    mean[idx] / sqrt(var[idx])
  }

  if(is.function(cor)){
    slice_cor <- function(idx, idx2 = idx){
      re <- cor(idx, idx2)
      if(!is.matrix(re)){
        re <- matrix(re, nrow = length(idx))
      }
      re
    }
  } else {
    slice_cor <- function(idx, idx2 = idx){
      cor[idx, idx2, drop = FALSE]
    }
  }


  stats <- sapply(conf$center, function(ctr){
    idx <- conf$window_idx(ctr)
    m <- slice_mean(idx)
    v <- slice_cor(idx)

    tval <- n * sum(m^2)
    param <- gamma_approx(v)
    shape <- param$M / param$u
    scale <- param$u
    # pval <- pgamma(tval, shape = shape, scale = scale, lower.tail = FALSE)
    df2 <- n - 1
    pval <- pf(tval / scale / shape, shape * 2, df2 = df2, lower.tail = FALSE)
    # print(c(ctr, pval))

    # shape = 10; scale = 4
    # pf(2 / scale / shape, shape * 2, df2 = 1)
    # pgamma(2, shape, scale = scale)

    # for post-selection inference
    mean_c <- slice_mean(ctr)
    c(tval, pval, shape, scale, mean_c, df2)
  })
  stats <- as.data.frame(t(stats))
  names(stats) <- c('tval', 'pval', 'shape', 'scale', 'mean_c', 'df2')

  # calculate covariance of the statistics
  if(isTRUE(sigma)){
    stop("sigma = TRUE not implemented")
    sigma <- sapply(conf$center, function(i1){
      if(i1 %% 10 ==0){print(i1)}
      idx1 <- conf$window_idx(i1)
      tmp <- slice_cor(idx1, seq_len(M))
      sapply(conf$center, function(i2){
        idx2 <- conf$window_idx(i2)
        v <- tmp[, idx2]
        2 * sum(v^2)
      })
    })
  }

  pretty_list(stats = stats, sigma = sigma, slice_mean = slice_mean, slice_cor = slice_cor, n = n, .list = conf)

}

grid_post_pvals_blocked <- function(res, rej, pval_cutoff){
  M <- res$total_points
  p0 <- rep(NA, M)

  rejpts <- unique(unlist(lapply(rej$rejs, function(idx){res$window_idx(idx)})))
  # eta2 <- length(rejpts) / rej$nrejs

  # calculate p-value threshold
  # pval_cutoff <- alpha2 * rej$nrejs / M

  thresholds <- qf(pval_cutoff, res$stats$shape * 2, df2 = res$stats$df2, lower.tail = FALSE) * res$stats$scale * res$stats$shape
  # thresholds <- qgamma(pval_cutoff, shape = res$stats$shape, scale = res$stats$scale, lower.tail = FALSE)
  # delta <- thresholds + res$n * (res$stats$mean_c)^2 - res$stats$tval
  # delta <= 0 means without C_k, B_k is significant

  rej2 <- rej$rejs

  m <- res$slice_mean(seq_len(M))
  p1 <- pf(m^2 * res$n, 1, df2 = res$stats$df2, lower.tail = FALSE)

  res <- lapply(rej$rejs, function(ctr){
    tval <- res$stats$tval[[ctr]]
    idx <- res$window_idx(ctr)
    m <- res$slice_mean(idx)
    thred <- thresholds[[ctr]]

    # print(c(thred, tval))
    d <- thred - (tval - res$n * m^2)
    # print(d)
    d[d < 0] <- 0
    # p <- pgamma(d, shape = 0.5, scale = 2, lower.tail = FALSE)
    p <- pf(d, 1, df2 = res$stats$df2[[ctr]], lower.tail = FALSE)
    cond_pvals = p1[idx] / p
    rej <- BH(cond_pvals)
    idx[rej$rejs]
  })
  unique(unlist(res))
}

grid_post_pvals <- function(res, rej, pval_cutoff){
  M <- res$total_points
  p0 <- rep(NA, M)

  rejpts <- unique(unlist(lapply(rej$rejs, function(idx){res$window_idx(idx)})))
  # eta2 <- length(rejpts) / rej$nrejs

  # calculate p-value threshold
  # pval_cutoff <- alpha2 * rej$nrejs / M

  thresholds <- qf(pval_cutoff, res$stats$shape * 2, df2 = res$stats$df2, lower.tail = FALSE) * res$stats$scale * res$stats$shape
  # thresholds <- qgamma(pval_cutoff, shape = res$stats$shape, scale = res$stats$scale, lower.tail = FALSE)
  # delta <- thresholds + res$n * (res$stats$mean_c)^2 - res$stats$tval
  # delta <= 0 means without C_k, B_k is significant

  rej2 <- rej$rejs

  for(ctr in rej$rejs){
    tval <- res$stats$tval[[ctr]]
    idx <- res$window_idx(ctr)
    m <- res$slice_mean(idx)
    thred <- thresholds[[ctr]]

    d <- thred - tval + res$n * m^2
    d[d < 0] <- 0
    # p <- pgamma(d, shape = 0.5, scale = 2, lower.tail = FALSE)
    p <- pf(d, 1, df2 = res$stats$df2[[ctr]], lower.tail = FALSE)

    # sub <- p0[idx]
    # sub[is.na(sub)] <- 1
    p0[idx] <- pmax(p0[idx], p, na.rm = TRUE)
    # print(p0[idx])
    # rej_sub <- step_up(pvals, alpha = alpha1)
    # s1_idx <- which(idx <= ctr)
    # s2_idx <- which(idx >= ctr)
    # c(rej_sub$nrejs / length(pvals),
    #   mean(s1_idx %in% rej_sub$rejs),
    #   mean(s2_idx %in% rej_sub$rejs)
    # )
  }
  # p0[is.na(p0)] <- 1
  m <- res$slice_mean(seq_len(M))
  # pgamma(m^2 * res$n, shape = 0.5, scale = 2, lower.tail = FALSE) / p0
  p1 <- pf(m^2 * res$n, 1, df2 = res$stats$df2, lower.tail = FALSE)

  # p1 / p0 should be stochastically above uniform. However, we want to keep the order
  p1 / p0
  # not_na <- !is.na(p0)
  # o <- order(p1[not_na])
  # p0v <- p0[not_na][o]
  # p1v <- p1[not_na][o]
  # p2 <- p1v / p0v
  # p2 <- cummax(p2)
  # p2 <- p2[order(o)]
  # p0[not_na] <- p2
  # p0
}
