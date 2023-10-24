#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat matrix_mult(const arma::mat& X,
                      const arma::mat& Y){
  int m = X.n_rows;
  int n = Y.n_cols;
  arma::mat Z(m, n);
  Z = X * Y;
  return Z;
}


