#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix timesTwo(NumericVector ds, int B = 1000) {
  NumericMatrix boot_stat(B, 2);
  int n = ds.size();
  
  for(int i=0; i<B; i++) {
    
    NumericVector gen_data = ds[floor(runif(n, 0, n))];
    boot_stat(i, 0) = mean(gen_data);
    boot_stat(i, 1) = sd(gen_data);
  }
  
  return(boot_stat);
}
