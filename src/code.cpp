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
  return std::max(t, 1);
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
  // if( ncores > col2.n_elem ){
  //   ncores = col2.n_elem;
  // }
  arma::rowvec rowmean = arma::mean(x, 0);
  R_xlen_t nrows = x.n_rows;

#if defined(_OPENMP)
  if(ncores <= 1){
    return arma::cov( x.cols(col1 - 1) , x.cols(col2 - 1) );
  }
  arma::mat re = arma::mat(col1.size(), col2.size());
  R_xlen_t i1,i2;

#pragma omp parallel num_threads(ncores)
{
#pragma omp for private(i1, i2) collapse(2)
{
  for( R_xlen_t c2 = 0 ; c2 < col2.size(); c2++ ){
    for( R_xlen_t c1 = 0 ; c1 < col1.size(); c1++ ){
      i2 = col2.at(c2) - 1;
      i1 = col1.at(c1) - 1;
      re(c1, c2) = arma::as_scalar((x.col(i1) - rowmean[i1]).t() *
        (x.col(i2) - rowmean[i2])) / (nrows - 1);

    }
  }
}
}
  return re;
#else
  return arma::cov( x.cols(col1 - 1) , x.cols(col2 - 1) );
#endif



}


/***R
x <- rnorm(10000)
microbenchmark::microbenchmark(
  native = {
    sum(x^2)
  },
  rcpp = {
    sumsquared(x)
  },
  check = function(v){
    abs(v[[1]]-v[[2]]) < 1e-6
  }
)

devtools::load_all()
data <- matrix(rnorm(10000),ncol = 100)
# data <- matrix(1:(10000),100)

# fastcov(data, col1, col2) - data[,col1]

col1 <- 1:8
col2 <- col1

microbenchmark::microbenchmark(
  native = {
    cov(data[,col1], data[,col2])
  },
  rcpp1 = {
    setThreads(1, TRUE)
    fastcov2(data, col1, col2)
  },
  rcpp2 = {
    setThreads(8, TRUE)
    fastcov2(data, col1, col2)
  },
  check = function(v){
    max(abs(v[[1]]-v[[2]])) < 1e-6
  }
)

*/
