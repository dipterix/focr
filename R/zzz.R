#' @name fdr-controls
#' @title External false discovery rate ('FDR') control methods
#' @description Procedures used by \code{\link{focr}} to control the
#' FDR in the second stage.
#' \describe{
#' \item{Benjamini-Hochberg procedure ('BH')}{
#' \doi{10.1111/j.2517-6161.1995.tb02031.x}}
#' \item{Benjamini-Yekutieli procedure ('BY')}{
#' \url{https://www.jstor.org/stable/2674075}}
#' \item{Structure-Adaptive Benjamini–Hochberg Algorithm ('SABHA')}{
#' \doi{10.1111/rssb.12298}}
#' \item{Locally Adaptive Weighting and Screening ('LAWS')}{
#' \doi{10.1080/01621459.2020.1859379}}
#' }
#' @param pv,pvals p-values
#' @param alpha 'FDR' level
#' @param filter,initial_filter initial p-value filters, the goal is to remove
#' large p-values when true signals are sparse. For \code{LAWS} and
#' \code{SABHA}, \code{initial_filter} helps determine sparsity levels.
#' @param bandwidth kernel smoothing bandwidth
#' @param dimension the spatial dimension of underlying data. Current
#' implementation only supports 1-3 dimensions
#' @return A list of rejection results:
#' \describe{
#' \item{\code{rejs}}{integer indices of rejected hypotheses;}
#' \item{\code{nrejs}}{total number of rejections;}
#' \item{\code{method}}{characters of method name;}
#' \item{\code{filter,initial_filter}}{passed from arguments, used to calculate "purity" values;}
#' \item{\code{pis_hat}}{estimated sparsity level (\code{LAWS} and
#' \code{SABHA} only);}
#' \item{\code{tau}}{p-value cutoff value (\code{BH} and
#' \code{BY} only);}
#' \item{\code{order,qvals}}{other variable from \code{BH} and
#' \code{BY};}
#' \item{\code{bandwidth,dimension,details}}{other variable from \code{LAWS} and
#' \code{SABHA}.}
#' }
#' @references
#' \cite{[1] Benjamini, Y. and Hochberg, Y., 1995. Controlling the false discovery
#' rate: a practical and powerful approach to multiple testing. Journal of the
#' Royal statistical society: series B (Methodological), 57(1), pp.289-300.}
#'
#' \cite{[2] Benjamini, Y. and Yekutieli, D., 2001. The control of the false
#' discovery rate in multiple testing under dependency. Annals of statistics,
#' pp.1165-1188.}
#'
#' \cite{[3] Li, A. and Barber, R.F., 2019. Multiple testing with the
#' structure‐adaptive Benjamini–Hochberg algorithm. Journal of the Royal
#' Statistical Society: Series B (Statistical Methodology), 81(1), pp.45-74.}
#'
#' \cite{[4] Cai, T.T., Sun, W. and Xia, Y., 2020. LAWS: A Locally Adaptive
#' Weighting and Screening Approach To Spatial Multiple Testing. Journal of the
#' American Statistical Association, pp.1-30.}
NULL

#' @rdname fdr-controls
#' @export
LAWS <- laws_pval
# function(pv, bandwidth, dimension = c("one", "two", "three"),
#                  alpha = 0.05, initial_filter = 0.9){
#   dimension <- match.arg(dimension)
#   tau <- BH(pv, initial_filter)$tau
#   laws_pval(pv = pv, bandwidth = bandwidth, dimension = dimension,
#             alpha = alpha, initial_filter = tau)
# }

#' @rdname fdr-controls
#' @export
SABHA <- sabha_pval

# function(pv, bandwidth, dimension = c("one", "two", "three"),
#                   alpha = 0.05, initial_filter = 0.9){
#   dimension <- match.arg(dimension)
#   tau <- BH(pv, initial_filter)$tau
#   sabha_pval(pv = pv, bandwidth = bandwidth, dimension = dimension,
#              alpha = alpha, initial_filter = tau)
# }

#' @rdname fdr-controls
#' @export
BH <- bh_pvals

#' @rdname fdr-controls
#' @export
BY <- by_pvals
