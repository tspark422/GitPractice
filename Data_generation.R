# Data simulation
gen_data <- function(n, beta, sigma, phi){
  if (sigma <= 0) stop("sigma must be positive")
  # Generate data
  x <- rnorm(n, 2, 1)
  p <- exp(x/2)/(1+exp(x/2))
  delta1 <- rbinom(n, 1, p)
  delta2 <- rbinom(n, 1, 0.7)
  y1 <- rnorm(n, beta[1]+beta[2]*x, sd=sigma)
  p2 <- exp(phi[1]+phi[2]*x+phi[3]*y1)/(1+exp(phi[1]+phi[2]*x+phi[3]*y1))
  y2 <- rbinom(n, 1, p2)
  data <- data.frame(id=1:n, x=x, y1=y1, y1_mis=y1, y2=y2, y2_mis=y2, delta1=delta1, delta2=delta2)
  
  # Missing value processing
  idx1 <- which(delta1==0)
  data$y1_mis[idx1] <- NA
  idx2 <- which(delta2==0)
  data$y2_mis[idx2] <- NA
  return(data)
}