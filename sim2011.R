sim2011 <- function(n, beta, sigma, phi, M=100, max.iter=50, eps=1e-09, B=2000, method=c("FI", "CFI")){
  # Check n
  if (n<=0 | n%%1!=0) stop("The number of observations, n, should be positive integer")
  
  # Check beta
  if (length(beta)!=2) stop("beta should be 2-dimensional vector")
  
  # Check sigma
  if (sigma<=0) stop("sigma should be positive")
  
  # Check phi
  if (length(phi)!=3) stop("phi should be 3-dimensional vector")
  
  # Check M
  if (M<=0 | M%%1!=0) stop("The number of imputed value, M, should be positive integer")
  
  # Check max.iter
  if (max.iter<=0 | max.iter%%1!=0) stop("The maximum number of iteration for EM algorithm, max.iter, should be positive integer")
  
  # Check eps
  if (eps <= 0) stop("The convergence criteria for EM algorithm, eps, should be positive")
  
  # Check B
  if (B<=0 | B%%1!=0) stop("The number of bootstrap process, B, should be positive integer")
  
  # Check method
  if (method=="FI"){
    source("myPFI.R")
    res <- myPFI(n=n, beta=beta, sigma=sigma, phi=phi, M=M, max.iter=max.iter, eps=eps, B=B)
  } else if (method=="CFI"){
    source("myPCFI.R")
    res <- myPCFI(n=n, beta=beta, sigma=sigma, phi=phi, M=M, max.iter=max.iter, eps=eps, B=B)
  } else {
    stop("method should be 'FI' or 'CFI'")
  }
}