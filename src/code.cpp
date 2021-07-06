// #include <Rcpp.h>
#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

// -------- Set up OPENMP
#ifdef _OPENMP
#include <omp.h>
#include <pthread.h>
#else
#define omp_get_thread_num() 0
#define omp_get_max_threads() 1
#endif


static int openMPThreads = 0;

// stores n threads when fork occurs
static bool detectFork = false;
static int reset_forked = true;

void onForked(){
  detectFork = true;
}
void onLeaveFork(){
  if(!reset_forked){
    openMPThreads = 1;
  }
  detectFork = false;
}
int detectForked(DllInfo *dll){
  // To disable openmp if fork is detected
#ifdef _OPENMP
  return pthread_atfork(&onForked, &onLeaveFork, NULL);
#endif

  return 0;
}

// [[Rcpp::export]]
int getThreads(){
#ifdef _OPENMP
  if(detectFork){
    return 1;
  }
  int t = openMPThreads <= 0 ? omp_get_max_threads() : std::min(openMPThreads, omp_get_max_threads());
  t = std::max(t, 1);
  // no need to go beyond 8 threads
  t = std::min(t, 8);
  return t;
#else
  return 1;
#endif
}

// [[Rcpp::export]]
int setThreads(int n, SEXP reset_after_fork){
#ifdef _OPENMP
  if(!detectFork){
    openMPThreads = n;
  }
  if( reset_after_fork != R_NilValue ){
    if( Rf_asLogical(reset_after_fork) ){
      reset_forked = true;
    } else {
      reset_forked = false;
    }
  }

  return n;
#else
  return 1;
#endif
}

// ----- Util functions
// void checkInterruptFn(void* /*dummy*/) {
//   R_CheckUserInterrupt();
// }


// ------ Main functions ------
// [[Rcpp::export]]
SEXP pis_1D(NumericVector &pval, double tau = 0.1,
            double h = 50, int verbose = 0) {

  R_xlen_t m = pval.length();
  SEXP re = PROTECT(Rf_allocVector(REALSXP, m));
  double* re_ptr = REAL(re);

  double tmp, tmp1, tmp2;
  R_xlen_t count = 0;

#if defined(_OPENMP)

  R_xlen_t ncores = getThreads();
  if( ncores > m ){ ncores = m; }

#pragma omp parallel num_threads(ncores) private(tmp, tmp1, tmp2)
{
#pragma omp for
  for( R_xlen_t i = 0; i < m; i++ ){

    tmp1 = tmp2 = 0;

    for( R_xlen_t s = 0; s < m; s++ ){
      tmp = abs(s - i) / h;
      tmp = exp(-0.5 * tmp * tmp) / h * 0.398942280401432677939946059934;

      if(pval[s] >= tau){
        tmp1 += tmp;
      }
      tmp2 += tmp;

    }

    tmp = tmp1 / tmp2 / (1 - tau);
    if( tmp > 1 ){ tmp = 1; }
    *(re_ptr+i) = 1.0 - tmp;

    if(verbose){
#pragma omp critical
{
  // Rcpp::checkUserInterrupt();
  count += 1;
  // verbose
  if(count % 1000 == 0){
    Rprintf("Calculating sparsity-level (%d of %d)\r", count, m);
  }
}
    }

  }
}
#else
  for( R_xlen_t i = 0; i < m; i++ ){

    tmp1 = tmp2 = 0;

    for( R_xlen_t s = 0; s < m; s++ ){
      tmp = abs(s - i) / h;
      tmp = exp(-0.5 * tmp * tmp) / h * 0.398942280401432677939946059934;

      if(pval[s] >= tau){
        tmp1 += tmp;
      }
      tmp2 += tmp;

    }

    tmp = tmp1 / tmp2 / (1 - tau);
    if( tmp > 1 ){ tmp = 1; }
    *(re_ptr+i) = 1.0 - tmp;

    if(verbose){
      // No openmp, so we can do this
      Rcpp::checkUserInterrupt();
      count += 1;
      // verbose
      if(count % 1000 == 0){
        Rprintf("Calculating sparsity-level (%d of %d)\r", count, m);
      }
    }

  }
#endif

  UNPROTECT(1);
  return re;
}


void disvec(
    double* re_ptr,
    const R_xlen_t &dim1, const R_xlen_t &dim2,
    const R_xlen_t &i, const R_xlen_t &j){

  // double* re_ptr = REAL(re);
  for(R_xlen_t ii = 0; ii < dim1; ii++ ){

    for(R_xlen_t jj = 0; jj < dim2; jj++ ){
      *(re_ptr + ii + jj * dim1) = sqrt(pow(ii - i, 2.0) + pow( jj - j , 2.0));
    }

  }
}

// [[Rcpp::export]]
SEXP pis_2D(NumericVector &pval,
            const R_xlen_t &dim1, const R_xlen_t &dim2,
            double tau = 0.1, double h = 10.0, int verbose = 0) {
  /**
  * pis_2D.func calculates the conditional proportions pis
  * Arguments
  * x: a matrix of z-values
  * tau: the screening threshold, which can be prespecified or chosen adaptively
  * bdw: bandwidth
  *
  * Values
  * pis: conditional proportions
  **/

  R_xlen_t m = dim1 * dim2;
  NumericVector::iterator pval_ptr = pval.begin();

  SEXP re = PROTECT(Rf_allocVector(REALSXP, m));
  // SEXP dim = PROTECT(Rf_allocVector(INTSXP, 2));
  // INTEGER(dim)[0] = dim1;
  // INTEGER(dim)[1] = dim2;
  // Rf_setAttrib(re, wrap("dim"), dim);

  double* re_ptr = REAL(re);

  // Use OpenMP
  R_xlen_t ncores = getThreads();
  if(ncores > dim1){
    ncores = dim1;
  }

#if defined(_OPENMP)

  std::vector<SEXP> buffers = std::vector<SEXP>(ncores);
  for(R_xlen_t core = 0; core < ncores; core++){
    buffers[core] = PROTECT(Rf_allocVector(REALSXP, m));
  }
  R_xlen_t thread_blocks = ceil(((float)dim1) / ncores);
  R_xlen_t count = 0;

#pragma omp parallel num_threads(ncores)
{
#pragma omp for schedule(static, 1) nowait
  for(R_xlen_t core = 0; core < ncores; core++){

    double tmp1, tmp2, tmp;
    R_xlen_t start = core * thread_blocks;
    R_xlen_t end = (core + 1) * thread_blocks;
    if( end > dim1){ end = dim1; }

    double* dis_ptr = REAL(buffers[core]);

    for(R_xlen_t i = start; i < end; i++){
      for(R_xlen_t j = 0; j < dim2; j++){
        disvec(dis_ptr, dim1, dim2, i, j);
        tmp1 = 0;
        tmp2 = 0;
        int count = 0;

        for(R_xlen_t l = 0; l < m; l++){

          tmp = *( dis_ptr + l ) / h;
          tmp = exp(-0.5 * tmp * tmp) / h * 0.398942280401432677939946059934;

          if( *(pval_ptr + l) >= tau ){
            tmp1 += tmp;
            count++;
          }
          tmp2 += tmp;
        }

        tmp = 1.0 - tmp1 / tmp2 / (1.0 - tau);
        if( tmp < 1e-5 ){
          tmp = 1e-5;
        }

        *(re_ptr + j * dim1 + i) = tmp;
      }
      if(verbose){
#pragma omp critical
{
      // Rcpp::checkUserInterrupt();
        count += 1;
        // verbose
        if(count % 10 == 0){
          Rprintf("Calculating sparsity-level (%d of %d)\r", count, dim2);
        }
}
      }
    }
  }
}
  UNPROTECT(ncores);
  Rprintf("                                               \r");

#else
  SEXP dis = PROTECT(Rf_allocVector(REALSXP, m));
  double* dis_ptr = REAL(dis);
  double tmp1, tmp2, tmp;
  for(R_xlen_t i = 0; i < dim1; i++){
    for(R_xlen_t j = 0; j < dim2; j++){

      disvec(dis_ptr, dim1, dim2, i, j);
      tmp1 = 0;
      tmp2 = 0;
      int count = 0;

      for(R_xlen_t l = 0; l < m; l++){

        tmp = *( dis_ptr + l ) / h;
        tmp = exp(-0.5 * tmp * tmp) / h * 0.398942280401432677939946059934;

        if( *(pval_ptr + l) >= tau ){
          tmp1 += tmp;
          count++;
        }
        tmp2 += tmp;
      }

      tmp = 1.0 - tmp1 / tmp2 / (1.0 - tau);
      if( tmp < 1e-5 ){
        tmp = 1e-5;
      }

      *(re_ptr + j * dim1 + i) = tmp;

    }
    Rcpp::checkUserInterrupt();
  }
  UNPROTECT(1);
#endif


  UNPROTECT(1);
  return(re);
}


// [[Rcpp::export]]
double sumsquared(NumericVector &x){
  double re = 0;
  for(NumericVector::iterator ptr = x.begin(); ptr!= x.end(); ptr++ ){
    re += std::pow(*ptr, 2.0);
  }
  return re;
}


// [[Rcpp::export]]
arma::mat fastcov(
    arma::mat& x,
    const arma::uvec& col1,
    const arma::uvec& col2) {

  R_xlen_t ncores = getThreads();
  if( ncores > col2.n_elem ){ ncores = col2.n_elem; }
  R_xlen_t nrows = x.n_rows;

  // arma::mat cm1 = arma::mean(x.cols(col1 - 1), 0);
  // arma::mat cm2 = arma::mean(x.cols(col2 - 1), 0);
  // x.cols(col1 - 1).t() * x.col(0) - cm1.t() * cm2.at(1);

  if( ncores == 1){
    return arma::cov( x.cols(col1 - 1) , x.cols(col2 - 1) );
  } else {
#if defined(_OPENMP)
  arma::mat re = arma::mat(col1.size(), col2.size());
  arma::mat cm1 = arma::mean(x.cols(col1 - 1), 0);
  arma::mat cm2 = arma::mean(x.cols(col2 - 1), 0);

#pragma omp parallel num_threads(ncores)
{
#pragma omp for collapse(2)
  for( R_xlen_t c2 = 0 ; c2 < col2.size(); c2++ ){
    for( R_xlen_t c1 = 0 ; c1 < col1.size(); c1++ ){
      R_xlen_t
      i2 = col2.at(c2) - 1,
        i1 = col1.at(c1) - 1;

      re(c1, c2) = (arma::as_scalar(
        x.col(i1).t() * x.col(i2)
      ) - (cm1.at(c1) * cm2.at(c2) * nrows)) / (nrows - 1);

    }
  }
}
return re;
#else
// This code will never execute, but just in case...
return arma::cov( x.cols(col1 - 1) , x.cols(col2 - 1) );
#endif
  }



}


/***R
devtools::load_all()
d <- c(100,10)
pval = runif(1000); dim(pval) = d
tau = 0.5; h = 10

re1 <- focr:::pis_2D.func(pval, tau = tau, h = h)
re2 <- focr:::pis_2D(pval, d[1], d[2], tau = tau, h = h)
# range(re2)
range(re2-re1)


disvec.func(c(2,5), c(2,4))

# devtools::load_all()
#
# microbenchmark::microbenchmark(
#   {
#     disvec.func(c(100,10), c(2,4))
#   },{
#     re <- double(1000)
#     disvec(re, 100,10,2,4)
#   }
# )
# re <- double(1000)
# disvec(re, 100,10,2,4)
# re - disvec.func(c(100,10), c(2,4))

# x <- rnorm(10000)
# microbenchmark::microbenchmark(
#   native = {
#     sum(x^2)
#   },
#   rcpp = {
#     sumsquared(x)
#   },
#   check = function(v){
#     abs(v[[1]]-v[[2]]) < 1e-6
#   }
# )
#
# devtools::load_all()
#
# x <- matrix(1:16,4)
# cov(x[,1:3], x[,1:3])
# fastcov2(x, 1:3, 1:3)
#
# data <- matrix(rnorm(100* 23700),nrow = 100)
# # data <- matrix(1:(10000),100)
#
# # fastcov(data, col1, col2) - data[,col1]
#
# col1 <- sample(23700, 3)
# col2 <- sample(23700, 3)
#
# slice_cor <- function(rows, cols){
#   # cov(data[,rows, drop = FALSE], data[,cols, drop = FALSE])
#   # stats::cor(data[,rows, drop = FALSE], data[,cols, drop = FALSE])
#   fastcov2(data, rows, cols)
# }
#
# system.time({
#   replicate(100, {
#     # setThreads(4, TRUE)
#     slice_cor(col1, col2)
#   })
# })
#
devtools::load_all()
data <- matrix(rnorm(100000* 300),ncol = 30000); pryr::object_size(data)
col1 <- sample(ncol(data), 300)
col2 <- col1
microbenchmark::microbenchmark(
  native = {
    cov(data[,col1], data[,col2])
  },
  rcpp1 = {
    focr:::setThreads(1, TRUE)
    focr:::fastcov2(data, col1, col2)
  },
  rcpp2 = {
    focr:::setThreads(8, TRUE)
    focr:::fastcov2(data, col1, col2)
  },
  rcpp3 = {
    focr:::setThreads(5, TRUE)
    focr:::fastcov2(data, col1, col2)
  },
  check = function(v){
    max(abs(v[[1]]-v[[3]])) < 1e-6
  }, unit = 'ms', times = 10
)

*/
