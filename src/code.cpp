#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double sumsquared(NumericVector x){
  double re = 0;
  for(NumericVector::iterator ptr = x.begin(); ptr!= x.end(); ptr++ ){
    re += std::pow(*ptr, 2.0);
  }
  return re;
}
