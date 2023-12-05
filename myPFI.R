myPFI <- function(data, M=100, max.iter=50, eps=1e-09, B=2000){
  
  # Data generation
  source('Data_generation.R')
  rst.beta <- matrix(0, nrow=B, ncol=length(beta))
  rst.sigma <- rep(0, B)
  rst.phi <- matrix(0, nrow=B, ncol=length(phi))
  rst.eta1 <- matrix(0, nrow=B, ncol=3)
  rst.eta2 <- matrix(0, nrow=B, ncol=3)
  rst.eta4 <- matrix(0, nrow=B, ncol=3)
  
  for (a in 1:B){
    # Data generation
    data <- gen_data(n=n, beta=beta, sigma=sigma, phi=phi)
    
    # Define arguments in the function
    n <- nrow(data)
    idx1 <- which(data$delta1==0)
    idx2 <- which(data$delta2==0)
    
    # Assign index set correspond to delta
    idx_A11 <- which(data$delta1==1 & data$delta2==1)
    idx_A10 <- which(data$delta1==1 & data$delta2==0)
    idx_A01 <- which(data$delta1==0 & data$delta2==1)
    idx_A00 <- which(data$delta1==0 & data$delta2==0)
    
    # Obtain initial coefficients of beta and sigma using observed y1 only
    y1_obs <- data$y1[!is.na(data$y1_mis)]
    x1_obs <- data$x[!is.na(data$y1_mis)]
    res_lm <- lm(y1_obs ~ x1_obs)
    init_beta <- res_lm$coefficients
    init_sigma2 <- var(res_lm$residuals)*((n-1)/n)
    
    # generate y1_mis
    idat <- data[rep(1:n, each=M),]           # idat = M times duplicated data
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
    fully_obs_idx <- idat$id %in% idx_A11
    pi_denom <- exp(init_phi[1]+init_phi[2]*idat$x+init_phi[3]*idat$y1_mis)/(1+exp(init_phi[1]+init_phi[2]*idat$x+init_phi[3]*idat$y1_mis))
    eps_value = 1
    for(iter in 1:max.iter){
      if(eps_value < eps) break()
      # E-step(Weighting)
      pi_num <- exp(phi_old[1]+phi_old[2]*idat$x+phi_old[3]*idat$y1_mis)/(1+exp(phi_old[1]+phi_old[2]*idat$x+phi_old[3]*idat$y1_mis))
      mis_weight <- dnorm(idat$y1_mis[!fully_obs_idx], mean=beta_old[1]+beta_old[2]*idat$x[!fully_obs_idx], sd=sqrt(sigma2_old))*pi_num[!fully_obs_idx]^(idat$y2_mis[!fully_obs_idx])*(1-pi_num[!fully_obs_idx])^(1-idat$y2_mis[!fully_obs_idx])/
        (dnorm(idat$y1_mis[!fully_obs_idx], mean=init_beta[1]+init_beta[2]*idat$x[!fully_obs_idx], sd=sqrt(init_sigma2))*pi_denom[!fully_obs_idx]^(idat$y2_mis[!fully_obs_idx])*(1-pi_denom[!fully_obs_idx])^(1-idat$y2_mis[!fully_obs_idx]))
      
      # normalize weight
      mis_weight_matrix <- matrix(mis_weight, nrow=M)
      mis_weight_matrix_normalized <- mis_weight_matrix/rep(colSums(mis_weight_matrix), each=M)
      idat$weight[!fully_obs_idx] <- c(mis_weight_matrix_normalized)
      
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
        phi_new = phi_old_NR + solve(information)%*%c(score_logis_1, score_logis_2, score_logis_3)   # Newton-Raphson method
        eps_logis = norm(phi_new-phi_old_NR, type='2')
        phi_old_NR=phi_new
      }
      eps_value=norm(c(beta_new-beta_old, sigma2_new-sigma2_old, phi_new-phi_old), type='2')
      beta_old=beta_new; sigma2_old=sigma2_new; phi_old=phi_new
    }
    
    ### Variance estimate ###
    pi_cali_old <- exp(phi_new[1]+phi_new[2]*idat$x+phi_new[3]*idat$y1_mis)/(1+exp(phi_new[1]+phi_new[2]*idat$x+phi_new[3]*idat$y1_mis))
    score_1_ij <- idat$y1_mis-beta_new[1]-beta_new[2]*idat$x
    score_2_ij <- (idat$y1_mis-beta_new[1]-beta_new[2]*idat$x)*idat$x
    score_3_ij <- (idat$y1_mis-beta_new[1]-beta_new[2]*idat$x)^2-sigma2_new
    score_4_ij <- idat$y2_mis-pi_cali_old
    score_5_ij <- (idat$y2_mis-pi_cali_old)*idat$x
    score_6_ij <- (idat$y2_mis-pi_cali_old)*idat$y1_mis
    idat$score_1_ij <- score_1_ij; idat$score_2_ij <- score_2_ij; idat$score_3_ij <- score_3_ij; idat$score_4_ij <- score_4_ij; idat$score_5_ij <- score_5_ij; idat$score_6_ij <- score_6_ij
    
    score_mean_i <- matrix(NA, nrow=6, ncol=n)
    for (m in 1:n){
      score_mean_i[, m] <- matrix(c(sum(idat$weight[idat$id==m]*idat$score_1_ij[idat$id==m]),
                                    sum(idat$weight[idat$id==m]*idat$score_2_ij[idat$id==m]),
                                    sum(idat$weight[idat$id==m]*idat$score_3_ij[idat$id==m]),
                                    sum(idat$weight[idat$id==m]*idat$score_4_ij[idat$id==m]),
                                    sum(idat$weight[idat$id==m]*idat$score_5_ij[idat$id==m]),
                                    sum(idat$weight[idat$id==m]*idat$score_6_ij[idat$id==m])), nrow=6)
    }
    first_term <- solve(score_mean_i%*%t(score_mean_i)/n)
    
    double_sig_eta1 <- 0
    for (i in 1:n){
      score_ij <- matrix(c(idat$score_1_ij[idat$id==i], 
                           idat$score_2_ij[idat$id==i], 
                           idat$score_3_ij[idat$id==i], 
                           idat$score_4_ij[idat$id==i], 
                           idat$score_5_ij[idat$id==i], 
                           idat$score_6_ij[idat$id==i]), nrow=6, byrow=T)
      double_sig_eta1 = double_sig_eta1 + apply(rep(c(idat$weight[idat$id==i]*idat$y1_mis[idat$id==i]), each=6)*(score_ij-score_mean_i[, i]), 1, sum) ## g(y) ????
    }
    double_sig_eta1 <- double_sig_eta1/n
    
    K_1_eta1 <- first_term%*%double_sig_eta1
    K_1_eta1
    e_i_eta1 <- matrix(NA, nrow=n, ncol=1)
    for(l in 1:n){
      score_ij <- matrix(c(idat$score_1_ij[idat$id==l], 
                           idat$score_2_ij[idat$id==l], 
                           idat$score_3_ij[idat$id==l], 
                           idat$score_4_ij[idat$id==l], 
                           idat$score_5_ij[idat$id==l], 
                           idat$score_6_ij[idat$id==l]), nrow=6, byrow=T)
      e_i_eta1[l, ] <- sum(idat$weight[idat$id==l]*(idat$y1_mis[idat$id==l]+t(K_1_eta1)%*%score_ij)) # ó��?? g(y_ij)?ڸ?
    }
    
    ### estimate for eta1 ###
    eta1 <- sum(idat$weight*idat$y1_mis)/n  # point estimator of eta1
    est_impvar1 <- 0                        # variance estimator of eta1 using linearization
    for(i in 1:n){
      for(j in 1:n){
        if (i==j){
          est_impvar1 <- est_impvar1 + (e_i_eta1[i, ])^2/(n^2)
        } else {
          est_impvar1 <- est_impvar1 - e_i_eta1[i, ]*e_i_eta1[j, ]/(n^2*(n-1))
        }
      }
    }
    
    est_compvar1 <- 0                       # Variance estimator of eta1 using complete sample
    for(i in 1:n){
      for(j in 1:n){
        if (i==j){
          est_compvar1 <- est_compvar1 + (data$y1[data$id==i])^2/(n^2)
        } else {
          est_compvar1 <- est_compvar1 - data$y1[data$id==i]*data$y1[data$id==j]/(n^2*(n-1))
        }
      }
    }
    
    ### estimate eta2 ###
    double_sig_eta2 <- 0              
    for (i in 1:n){
      score_ij <- matrix(c(idat$score_1_ij[idat$id==i], 
                           idat$score_2_ij[idat$id==i], 
                           idat$score_3_ij[idat$id==i], 
                           idat$score_4_ij[idat$id==i], 
                           idat$score_5_ij[idat$id==i], 
                           idat$score_6_ij[idat$id==i]), nrow=6, byrow=T)
      double_sig_eta2 = double_sig_eta2 + apply(rep(c(idat$weight[idat$id==i]*idat$y2_mis[idat$id==i]), each=6)*(score_ij-score_mean_i[, i]), 1, sum) 
    }
    double_sig_eta2 <- double_sig_eta2/n
    
    K_1_eta2 <- first_term%*%double_sig_eta2
    K_1_eta2
    e_i_eta2 <- matrix(NA, nrow=n, ncol=1)
    for(l in 1:n){
      score_ij <- matrix(c(idat$score_1_ij[idat$id==l], 
                           idat$score_2_ij[idat$id==l], 
                           idat$score_3_ij[idat$id==l], 
                           idat$score_4_ij[idat$id==l], 
                           idat$score_5_ij[idat$id==l], 
                           idat$score_6_ij[idat$id==l]), nrow=6, byrow=T)
      e_i_eta2[l, ] <- sum(idat$weight[idat$id==l]*(idat$y2_mis[idat$id==l]+t(K_1_eta2)%*%score_ij)) 
    }
    
    eta2 <- sum(idat$weight*idat$y2_mis)/n     # point estimator of eta2
    est_impvar2 <- 0                           # variance estimator of eta2 using linearization
    for(i in 1:n){
      for(j in 1:n){
        if (i==j){
          est_impvar2 <- est_impvar2 + (e_i_eta2[i, ])^2/(n^2)
        } else {
          est_impvar2 <- est_impvar2 - e_i_eta2[i, ]*e_i_eta2[j, ]/(n^2*(n-1))
        }
      }
    }
    
    est_compvar2 <- 0                           # variance estimator of eta2 using complete sample
    for(i in 1:n){
      for(j in 1:n){
        if (i==j){
          est_compvar2 <- est_compvar2 + (data$y2[data$id==i])^2/(n^2)
        } else {
          est_compvar2 <- est_compvar2 - data$y2[data$id==i]*data$y2[data$id==j]/(n^2*(n-1))
        }
      }
    }
    
    ### estimate for eta4 ###
    double_sig_eta4 <- 0
    for (i in 1:n){
      score_ij <- matrix(c(idat$score_1_ij[idat$id==i], 
                           idat$score_2_ij[idat$id==i], 
                           idat$score_3_ij[idat$id==i], 
                           idat$score_4_ij[idat$id==i], 
                           idat$score_5_ij[idat$id==i], 
                           idat$score_6_ij[idat$id==i]), nrow=6, byrow=T)
      double_sig_eta4 = double_sig_eta4 + apply(rep(c(idat$weight[idat$id==i]*(idat$y1_mis[idat$id==i]<3)), each=6)*(score_ij-score_mean_i[, i]), 1, sum) 
    }
    double_sig_eta4 <- double_sig_eta4/n
    
    K_1_eta4 <- first_term%*%double_sig_eta4
    K_1_eta4
    e_i_eta4 <- matrix(NA, nrow=n, ncol=1)
    for(l in 1:n){
      score_ij <- matrix(c(idat$score_1_ij[idat$id==l], 
                           idat$score_2_ij[idat$id==l], 
                           idat$score_3_ij[idat$id==l], 
                           idat$score_4_ij[idat$id==l], 
                           idat$score_5_ij[idat$id==l], 
                           idat$score_6_ij[idat$id==l]), nrow=6, byrow=T)
      e_i_eta4[l, ] <- sum(idat$weight[idat$id==l]*((idat$y1_mis[idat$id==l]<3)+t(K_1_eta4)%*%score_ij))
    }
    
    eta4 <- sum(idat$weight*(idat$y1_mis<3))/n # point estimator of eta4
    est_impvar4 <- 0                           # variance estimator of eta4 using linearization
    for(i in 1:n){
      for(j in 1:n){
        if (i==j){
          est_impvar4 <- est_impvar4 + (e_i_eta4[i, ])^2/(n^2)
        } else {
          est_impvar4 <- est_impvar4 - e_i_eta4[i, ]*e_i_eta4[j, ]/(n^2*(n-1))
        }
      }
    }
    
    est_compvar4 <- 0                           # variance estimator of eta2 using complete sample
    for(i in 1:n){
      for(j in 1:n){
        if (i==j){
          est_compvar4 <- est_compvar4 + (data$y1[data$id==i]<3)^2/(n^2)
        } else {
          est_compvar4 <- est_compvar4 - (data$y1[data$id==i]<3)*(data$y1[data$id==j]<3)/(n^2*(n-1))
        }
      }
    }
    
    rst.beta[a, ] = beta_new
    rst.sigma[a] = sigma2_new
    rst.phi[a, ] = phi_new
    rst.eta1[a, ] = c(eta1, est_impvar1, est_compvar1)
    rst.eta2[a, ] = c(eta2, est_impvar2, est_compvar2)
    rst.eta4[a, ] = round(c(eta4, est_impvar4, est_compvar4), 5)
    
    if(a%%10==0){
      print(a)
      cat("______beta__________",'\n')
      print(apply(rst.beta[1:a,],2,mean))
      print(apply(rst.beta[1:a,],2,sd))
      cat("______sigma__________",'\n')
      print(mean(rst.sigma[1:a]))
      print(sd(rst.sigma[1:a]))
      cat("______pi__________",'\n')
      print(apply(rst.phi[1:a,],2,mean))
      print(apply(rst.phi[1:a,],2,sd))
      cat("______eta1__________",'\n')
      print(apply(rst.eta1[1:a,],2,mean))
      print(apply(rst.eta1[1:a,],2,sd))
      cat("______eta2__________",'\n')
      print(apply(rst.eta2[1:a,],2,mean))
      print(apply(rst.eta2[1:a,],2,sd))
      cat("______eta4__________",'\n')
      print(apply(rst.eta4[1:a,],2,mean))
      print(apply(rst.eta4[1:a,],2,sd))
    }
  }
}
