sim2011 <- function(n, beta, sigma, phi, M=100, max.iter=50, eps=1e-09, B=2000, method=c("FI", "CFI")){
  if (method=="FI") {
    source("myPFI.R")
    res <- function(n=n, beta=beta, sigma=sigma, phi=phi, M=M, max.iter=max.iter, eps=eps, B=B)
  }
}