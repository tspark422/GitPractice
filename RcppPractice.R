library(Rcpp)
library(microbenchmark)
sourceCpp("bootstat.cpp")
sourceCpp("arma_mm.cpp")
X = matrix(rnorm(30000), 300, 100)
Y = matrix(rnorm(20000), 100, 200)
prodCpp = matrix_mult(X, Y)
prodR = X%*%Y
all.equal(prodCpp, prodR)

microbenchmark(
  matrix_mult(X, Y),
  X%*%Y,
  times=50
)
