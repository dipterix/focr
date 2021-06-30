// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// getThreads
int getThreads();
RcppExport SEXP _focr_getThreads() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(getThreads());
    return rcpp_result_gen;
END_RCPP
}
// setThreads
int setThreads(int n, SEXP reset_after_fork);
RcppExport SEXP _focr_setThreads(SEXP nSEXP, SEXP reset_after_forkSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< SEXP >::type reset_after_fork(reset_after_forkSEXP);
    rcpp_result_gen = Rcpp::wrap(setThreads(n, reset_after_fork));
    return rcpp_result_gen;
END_RCPP
}
// pis_2D
SEXP pis_2D(NumericVector& pval, const R_xlen_t& dim1, const R_xlen_t& dim2, double tau, double h, int verbose);
RcppExport SEXP _focr_pis_2D(SEXP pvalSEXP, SEXP dim1SEXP, SEXP dim2SEXP, SEXP tauSEXP, SEXP hSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type pval(pvalSEXP);
    Rcpp::traits::input_parameter< const R_xlen_t& >::type dim1(dim1SEXP);
    Rcpp::traits::input_parameter< const R_xlen_t& >::type dim2(dim2SEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< int >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(pis_2D(pval, dim1, dim2, tau, h, verbose));
    return rcpp_result_gen;
END_RCPP
}
// sumsquared
double sumsquared(NumericVector& x);
RcppExport SEXP _focr_sumsquared(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(sumsquared(x));
    return rcpp_result_gen;
END_RCPP
}
// fastcov
arma::mat fastcov(arma::mat& x, const arma::uvec& col1, const arma::uvec& col2);
RcppExport SEXP _focr_fastcov(SEXP xSEXP, SEXP col1SEXP, SEXP col2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type col1(col1SEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type col2(col2SEXP);
    rcpp_result_gen = Rcpp::wrap(fastcov(x, col1, col2));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_focr_getThreads", (DL_FUNC) &_focr_getThreads, 0},
    {"_focr_setThreads", (DL_FUNC) &_focr_setThreads, 2},
    {"_focr_pis_2D", (DL_FUNC) &_focr_pis_2D, 6},
    {"_focr_sumsquared", (DL_FUNC) &_focr_sumsquared, 1},
    {"_focr_fastcov", (DL_FUNC) &_focr_fastcov, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_focr(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
