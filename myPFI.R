myPFI <- function(dat, M=100, max.iter=50, eps=1e-09){
  n <- nrow(dat)
  idx1 <- which(dat$delta1==0)
  idx2 <- which(dat$delta2==0)
  
  # Assign index set correspond to delta
  idx_A10 <- which(dat$delta2==0 & dat$delta1!=0)
  idx_A01 <- which(dat$delta1==0 & dat$delta2!=0)
  idx_A00 <- which(dat$delta1==0 & dat$delta2==0)
  
  # Obtain initial coefficients of beta and sigma using observed y1 only
  y1_obs <- dat$y1[!is.na(dat$y1_mis)]
  x1_obs <- dat$x[!is.na(dat$y1_mis)]
  res_lm <- lm(y1_obs ~ x1_obs)
  init_beta <- res_lm$coefficients
  init_sigma2 <- var(res_lm$residuals)*((n-1)/n)
  
  # generate y1_mis
  idat <- dat[rep(1:n, each=M),]           # idat = M times duplicated data
  idat$weight <- rep(1/M, dim(idat)[1])
  idat$y1_mis[is.na(idat$y1_mis)] <- rnorm(length(idx1)*M, init_beta[1] + init_beta[2]*idat$x[is.na(idat$y1_mis)], sd=sqrt(init_sigma2))
  
  # Obtain initial coefficients of phi
  res_logis <- glm(y2_mis ~ 1+x+y1_mis, family=binomial(link='logit'), data=idat, weights=weight)
  init_phi <- res_logis$coef
  
  # generate y2_mis
  idat$y2_mis[is.na(idat$y2_mis)] <- rbinom(length(idx2)*M,
                                            1,
                                            exp(init_phi[1]+init_phi[2]*idat$x[is.na(idat$y2_mis)]+init_phi[3]*idat$y1_mis[is.na(idat$y2_mis)])/(1+exp(init_phi[1]+init_phi[2]*idat$x[is.na(idat$y2_mis)]+init_phi[3]*idat$y1_mis[is.na(idat$y2_mis)])))
  
  # EM algorithm
  phi_old <- init_phi; beta_old <- init_beta; sigma2_old <- init_sigma2
  idx_A10_weight <- idat$id %in% idx_A10; idx_A01_weight <- idat$id %in% idx_A01; idx_A00_weight <- idat$id %in% idx_A00
  pi_denom <- exp(init_phi[1]+init_phi[2]*idat$x+init_phi[3]*idat$y1_mis)/(1+exp(init_phi[1]+init_phi[2]*idat$x+init_phi[3]*idat$y1_mis))
  for(iter in 1:max.iter){
    if(eps < 1e-09) break()
    # E-step(Weighting)
    pi_num <- exp(phi_old[1]+phi_old[2]*idat$x+phi_old[3]*idat$y1_mis)/(1+exp(phi_old[1]+phi_old[2]*idat$x+phi_old[3]*idat$y1_mis))
    idat$weight[idx_A10_weight] <- pi_num[idx_A10_weight]^(idat$y2_mis[idx_A10_weight])*(1-pi_num[idx_A10_weight])^(1-idat$y2_mis[idx_A10_weight])/
      (pi_denom[idx_A10_weight]^(idat$y2_mis[idx_A10_weight])*(1-pi_denom[idx_A10_weight])^(1-idat$y2_mis[idx_A10_weight]))
    idat$weight[idx_A01_weight] <- dnorm(idat$y1_mis[idx_A01_weight], mean=beta_old[1]+beta_old[2]*idat$x[idx_A01_weight], sd=sqrt(sigma2_old))*pi_num[idx_A01_weight]^(idat$y2_mis[idx_A01_weight])*(1-pi_num[idx_A01_weight])^(1-idat$y2_mis[idx_A01_weight])/
      (dnorm(idat$y1_mis[idx_A01_weight], mean=init_beta[1]+init_beta[2]*idat$x[idx_A01_weight], sd=sqrt(init_sigma2)))
    idat$weight[idx_A00_weight] <- dnorm(idat$y1_mis[idx_A00_weight], mean=beta_old[1]+beta_old[2]*idat$x[idx_A00_weight], sd=sqrt(sigma2_old))*pi_num[idx_A00_weight]^(idat$y2_mis[idx_A00_weight])*(1-pi_num[idx_A00_weight])^(1-idat$y2_mis[idx_A00_weight])/
      (dnorm(idat$y1_mis[idx_A00_weight], mean=init_beta[1]+init_beta[2]*idat$x[idx_A00_weight], sd=sqrt(init_sigma2))*pi_denom[idx_A00_weight]^(idat$y2_mis[idx_A00_weight])*(1-pi_denom[idx_A00_weight])^(1-idat$y2_mis[idx_A00_weight]))
    
    for (k in 1:n){
      idat$weight[idat$id==k] <- idat$weight[idat$id==k]/sum(idat$weight[idat$id==k])
    }
    
    # M-step(Maximizing)
    # beta
    score1 <- sum(idat$weight*(idat$y1_mis-beta_old[1]-beta_old[2]*idat$x))
    score2 <- sum(idat$weight*(idat$y1_mis-beta_old[1]-beta_old[2]*idat$x)*idat$x)
    U11 <- sum(idat$weight*(-1))
    U12 <- sum(idat$weight*(-idat$x))
    U22 <- sum(idat$weight*(-idat$x^2))
    information <- matrix(c(U11, U12, U12, U22), nrow=2)         # Information matrix for beta
    beta_new = beta_old - solve(information)%*%c(score1, score2) # Newton-Raphson method
    
    # Sigma^2
    sigma2_new <- sum(idat$weight*(idat$y1_mis-beta_new[1]-beta_new[2]*idat$x)^2)/n
    
    # phi
    eps_logis = 1
    phi_old_NR=phi_old
    for(l in 1:max.iter){
      if(eps_logis < 1e-06) break()
      pi_logis_old <- exp(phi_old_NR[1]+phi_old_NR[2]*idat$x+phi_old_NR[3]*idat$y1_mis)/(1+exp(phi_old_NR[1]+phi_old_NR[2]*idat$x+phi_old_NR[3]*idat$y1_mis))
      score_logis_1 <- sum(idat$weight*(idat$y2_mis-pi_logis_old))
      score_logis_2 <- sum(idat$weight*(idat$y2_mis-pi_logis_old)*idat$x)
      score_logis_3 <- sum(idat$weight*(idat$y2_mis-pi_logis_old)*idat$y1_mis)
      U_logis_11 <- sum(idat$weight*pi_logis_old*(1-pi_logis_old))
      U_logis_12 <- sum(idat$weight*pi_logis_old*(1-pi_logis_old)*idat$x)
      U_logis_13 <- sum(idat$weight*pi_logis_old*(1-pi_logis_old)*idat$y1_mis)
      U_logis_22 <- sum(idat$weight*pi_logis_old*(1-pi_logis_old)*idat$x^2)
      U_logis_23 <- sum(idat$weight*pi_logis_old*(1-pi_logis_old)*idat$x*idat$y1_mis)
      U_logis_33 <- sum(idat$weight*pi_logis_old*(1-pi_logis_old)*idat$y1_mis^2)
      information = matrix(c(U_logis_11, U_logis_12, U_logis_13, U_logis_12, U_logis_22, U_logis_23, U_logis_13, U_logis_23, U_logis_33), nrow=3) # Information matrix
      pi_new = phi_old_NR + solve(information)%*%c(score_logis_1, score_logis_2, score_logis_3)   # Newton-Raphson method
      eps_logis = norm(pi_new-phi_old_NR, type='2')
      phi_old_NR=phi_new
    }
    eps=norm(c(beta_new-beta_old, sigma2_new-sigma2_old, phi_new-phi_old), type='2')
    beta_old=beta_new; sigma2_old=sigma2_new; phi_old=phi_new
  }
  return(list(beta=beta_new, sigma=sqrt(sigma2_new), phi=phi_new))
}
