myPCFI_prac <- function(n, beta, sigma, phi, M=100, max.iter=50, eps=1e-09, B=2000){
  
  # Data generation
  source('Data_generation.R')
  rst.eta1 <- matrix(0, nrow=B, ncol=3)
  rst.eta2 <- matrix(0, nrow=B, ncol=3)
  rst.eta3 <- matrix(0, nrow=B, ncol=3)
  rst.eta4 <- matrix(0, nrow=B, ncol=3)
  
  for (a in 1:B){
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
      
      
      # Calibration
      exp_prob <- exp(phi_old[1]+phi_old[2]*idat$x+phi_old[3]*idat$y1_mis)
      pi_cali_old <- exp_prob/(1+exp_prob)
      score_1_ij <- idat$y1_mis-beta_old[1]-beta_old[2]*idat$x
      score_2_ij <- (idat$y1_mis-beta_old[1]-beta_old[2]*idat$x)*idat$x
      score_3_ij <- (idat$y1_mis-beta_old[1]-beta_old[2]*idat$x)^2-sigma2_old
      score_4_ij <- idat$y2_mis-pi_cali_old
      score_5_ij <- (idat$y2_mis-pi_cali_old)*idat$x
      score_6_ij <- (idat$y2_mis-pi_cali_old)*idat$y1_mis
      idat$score_1_ij <- score_1_ij; idat$score_2_ij <- score_2_ij; idat$score_3_ij <- score_3_ij; idat$score_4_ij <- score_4_ij; idat$score_5_ij <- score_5_ij; idat$score_6_ij <- score_6_ij
      

      # Each row represents score function
      score_i <- matrix(c(colSums(matrix(idat$weight*idat$score_1_ij, nrow=M)),
                          colSums(matrix(idat$weight*idat$score_2_ij, nrow=M)),
                          colSums(matrix(idat$weight*idat$score_3_ij, nrow=M)),
                          colSums(matrix(idat$weight*idat$score_4_ij, nrow=M)),
                          colSums(matrix(idat$weight*idat$score_5_ij, nrow=M)),
                          colSums(matrix(idat$weight*idat$score_6_ij, nrow=M))), nrow=6, byrow=T)
      
      score_i_rep <- matrix(rep(t(score_i), each=M), nrow = 6, byrow=T)
      idat_cali <- cbind(idat, t(score_i_rep))
      
      first_term <- rowMeans(score_i) #first term in eq (10).
      
      middle_mat <- matrix(0, nrow=6, ncol=6)
      for (o in 1:n*M){
        score_mat <- matrix(c(idat_cali$score_1_ij[o], 
                              idat_cali$score_2_ij[o], 
                              idat_cali$score_3_ij[o], 
                              idat_cali$score_4_ij[o], 
                              idat_cali$score_5_ij[o], 
                              idat_cali$score_6_ij[o]), nrow=6)
        middle_mat <- middle_mat + (idat_cali$weight[o]*(score_mat-score_i_rep[, o]))%*%t(score_mat-score_i_rep[, o])
      }
      
      full_score_mat <- matrix(c(idat_cali$score_1_ij, 
                                 idat_cali$score_2_ij, 
                                 idat_cali$score_3_ij, 
                                 idat_cali$score_4_ij, 
                                 idat_cali$score_5_ij, 
                                 idat_cali$score_6_ij), nrow=6, byrow=T)
      final_score_mat <- full_score_mat-score_i_rep
      
      res <- first_term%*%crossprod(t(solve(middle_mat)), final_score_mat)
      idat_cali$weight2 <- idat_cali$weight - c(idat_cali$weight*res)         # final weight: weight2
      
      
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
    
    # Define score function (S(theta))
    pi_cali_old <- exp(phi_new[1]+phi_new[2]*idat$x+phi_new[3]*idat$y1_mis)/(1+exp(phi_new[1]+phi_new[2]*idat$x+phi_new[3]*idat$y1_mis))
    score_1_ij <- idat$y1_mis-beta_new[1]-beta_new[2]*idat$x
    score_2_ij <- (idat$y1_mis-beta_new[1]-beta_new[2]*idat$x)*idat$x
    score_3_ij <- (idat$y1_mis-beta_new[1]-beta_new[2]*idat$x)^2-sigma2_new
    score_4_ij <- idat$y2_mis-pi_cali_old
    score_5_ij <- (idat$y2_mis-pi_cali_old)*idat$x
    score_6_ij <- (idat$y2_mis-pi_cali_old)*idat$y1_mis
    idat$score_1_ij <- score_1_ij; idat$score_2_ij <- score_2_ij; idat$score_3_ij <- score_3_ij; idat$score_4_ij <- score_4_ij; idat$score_5_ij <- score_5_ij; idat$score_6_ij <- score_6_ij
    
    # Define \bar{s_{i}^{*}} by 6 by n matrix
    # Each column represents each sample
    # Each row represents score function
    score_mean_i <- matrix(c(colSums(matrix(idat$weight*idat$score_1_ij, nrow=M)),
                             colSums(matrix(idat$weight*idat$score_2_ij, nrow=M)),
                             colSums(matrix(idat$weight*idat$score_3_ij, nrow=M)),
                             colSums(matrix(idat$weight*idat$score_4_ij, nrow=M)),
                             colSums(matrix(idat$weight*idat$score_5_ij, nrow=M)),
                             colSums(matrix(idat$weight*idat$score_6_ij, nrow=M))), nrow=6, byrow=T)
    first_term <- solve(tcrossprod(score_mean_i))*n
    
    
    ### estimate for eta1 ###
    # Calculate double sigma w_iw_{ij}^{*}(s-s_{i})*g(y_{ij}^{*}) for eta1
    g1 <- idat$y1_mis # estimating function for eta1
    common_term <- idat$weight*g1
    c1 <- sum(common_term*(idat$score_1_ij-rep(score_mean_i[1, ], each=M)))
    c2 <- sum(common_term*(idat$score_2_ij-rep(score_mean_i[2, ], each=M)))
    c3 <- sum(common_term*(idat$score_3_ij-rep(score_mean_i[3, ], each=M)))
    c4 <- sum(common_term*(idat$score_4_ij-rep(score_mean_i[4, ], each=M)))
    c5 <- sum(common_term*(idat$score_5_ij-rep(score_mean_i[5, ], each=M)))
    c6 <- sum(common_term*(idat$score_6_ij-rep(score_mean_i[6, ], each=M)))
    double_sig_eta1 <- c(c1, c2, c3, c4, c5, c6)/n
    
    K_1_eta1 <- first_term%*%double_sig_eta1
    
    # Calculate e_i term for eta1
    score_mat <- matrix(c(idat$score_1_ij,
                          idat$score_2_ij,
                          idat$score_3_ij,
                          idat$score_4_ij,
                          idat$score_5_ij,
                          idat$score_6_ij), byrow=T, nrow=6)
    e_i_eta1 <- matrix(colSums(matrix(idat$weight*(idat$y1_mis+crossprod(K_1_eta1, score_mat)), nrow=M)), ncol=1)
    
    ### results for eta1 ###
    eta1 <- sum(1/n*idat$weight*idat$y1_mis)  # point estimator of eta1 
    est_impvar1 <- as.numeric(crossprod(e_i_eta1)/(n*(n-1)) - (sum(e_i_eta1))^2/(n^2*(n-1))) # variance estimator of eta1 using linearization (eq. (13))
    est_compvar1 <- as.numeric(crossprod(data$y1)/(n*(n-1)) - (sum(data$y1)^2/(n^2*(n-1))))  # Variance estimator of eta1 using complete sample
    
    
    ### estimate eta2 ###
    # Calculate double sigma in K_1
    g2 <- idat$y2_mis
    common_term <- idat$weight*g2
    c1 <- sum(common_term*(idat$score_1_ij-rep(score_mean_i[1, ], each=M)))
    c2 <- sum(common_term*(idat$score_2_ij-rep(score_mean_i[2, ], each=M)))
    c3 <- sum(common_term*(idat$score_3_ij-rep(score_mean_i[3, ], each=M)))
    c4 <- sum(common_term*(idat$score_4_ij-rep(score_mean_i[4, ], each=M)))
    c5 <- sum(common_term*(idat$score_5_ij-rep(score_mean_i[5, ], each=M)))
    c6 <- sum(common_term*(idat$score_6_ij-rep(score_mean_i[6, ], each=M)))
    double_sig_eta2 <- c(c1, c2, c3, c4, c5, c6)/n
    
    K_1_eta2 <- first_term%*%double_sig_eta2
    e_i_eta2 <- matrix(colSums(matrix(idat$weight*(idat$y2_mis+crossprod(K_1_eta2, score_mat)), nrow=M)), ncol=1)
    
    ### results for eta2 ###
    eta2 <- sum(1/n*idat$weight*idat$y2_mis)     # point estimator of eta2
    est_impvar2 <- as.numeric(crossprod(e_i_eta2)/(n*(n-1)) - (sum(e_i_eta2))^2/(n^2*(n-1))) # variance estimator of eta2 using linearization (eq. (13))
    est_compvar2 <- as.numeric(crossprod(data$y2)/(n*(n-1)) - (sum(data$y2)^2/(n^2*(n-1))))  # Variance estimator of eta2 using complete sample
    
    ### estimate eta3 ###
    # Calculate double sigma in K_1
    g3 <- (idat$x-mean(idat$x))*(idat$y1_mis-mean(idat$y1_mis))/sum((idat$x-mean(idat$x))^2)*n*M
    g_comp <- (data$x-mean(data$x))*(data$y1-mean(data$y1))/sum((data$x-mean(data$x))^2)*n
    common_term <- idat$weight*g3
    # common_term <- idat$weight*as.numeric(Sxy/Sxx)
    c1 <- sum(common_term*(idat$score_1_ij-rep(score_mean_i[1, ], each=M)))
    c2 <- sum(common_term*(idat$score_2_ij-rep(score_mean_i[2, ], each=M)))
    c3 <- sum(common_term*(idat$score_3_ij-rep(score_mean_i[3, ], each=M)))
    c4 <- sum(common_term*(idat$score_4_ij-rep(score_mean_i[4, ], each=M)))
    c5 <- sum(common_term*(idat$score_5_ij-rep(score_mean_i[5, ], each=M)))
    c6 <- sum(common_term*(idat$score_6_ij-rep(score_mean_i[6, ], each=M)))
    double_sig_eta3 <- c(c1, c2, c3, c4, c5, c6)/n
    
    K_1_eta3 <- first_term%*%double_sig_eta3
    e_i_eta3 <- matrix(colSums(matrix(idat$weight*(g3+crossprod(K_1_eta3, score_mat)), nrow=M)), ncol=1)
    
    ### results for eta3 ###
    eta3 <- sum(1/n*idat$weight*g3)     # point estimator of eta2
    est_impvar3 <- as.numeric(crossprod(e_i_eta3)/(n*(n-1)) - (sum(e_i_eta3))^2/(n^2*(n-1))) # variance estimator of eta3 using linearization (eq. (13))
    est_compvar3 <- as.numeric(crossprod(g_comp)/(n*(n-1)) - ((sum(g_comp))^2/(n^2*(n-1))))  # Variance estimator of eta3 using complete sample
    
    ### estimate for eta4 ###
    # Calculate double sigma in K_1
    g4 <- as.numeric(idat$y1_mis<3)
    common_term <- idat$weight*g4
    c1 <- sum(common_term*(idat$score_1_ij-rep(score_mean_i[1, ], each=M)))
    c2 <- sum(common_term*(idat$score_2_ij-rep(score_mean_i[2, ], each=M)))
    c3 <- sum(common_term*(idat$score_3_ij-rep(score_mean_i[3, ], each=M)))
    c4 <- sum(common_term*(idat$score_4_ij-rep(score_mean_i[4, ], each=M)))
    c5 <- sum(common_term*(idat$score_5_ij-rep(score_mean_i[5, ], each=M)))
    c6 <- sum(common_term*(idat$score_6_ij-rep(score_mean_i[6, ], each=M)))
    double_sig_eta4 <- c(c1, c2, c3, c4, c5, c6)/n
    
    K_1_eta4 <- first_term%*%double_sig_eta4
    e_i_eta4 <- matrix(colSums(matrix(idat$weight*(as.numeric(idat$y1_mis<3)+crossprod(K_1_eta4, score_mat)), nrow=M)), ncol=1)
    
    ### results for eta4 ###
    eta4 <- sum(1/n*idat$weight*(idat$y1_mis<3)) # point estimator of eta4
    est_impvar4 <- as.numeric(crossprod(e_i_eta2)/(n*(n-1)) - (sum(e_i_eta2))^2/(n^2*(n-1))) # variance estimator of eta1 using linearization (eq. (13))
    est_compvar4 <- as.numeric(crossprod(data$y1<3)/(n*(n-1)) - (sum(data$y1<3)^2/(n^2*(n-1))))  # Variance estimator of eta1 using complete sample
    
    rst.eta1[a, ] = c(eta1, est_impvar1, est_compvar1)
    rst.eta2[a, ] = c(eta2, est_impvar2, est_compvar2)
    rst.eta3[a, ] = c(eta3, est_impvar3, est_compvar3)
    rst.eta4[a, ] = c(eta4, est_impvar4, est_compvar4)
    
    if(a%%100==0){
      print(a)
      cat("______eta1__________",'\n')
      print(round(apply(rst.eta1[1:a,],2,mean), 5))
      print(round(apply(rst.eta1[1:a,],2,sd), 5))
      cat("______eta2__________",'\n')
      print(round(apply(rst.eta2[1:a,],2,mean), 5))
      print(round(apply(rst.eta2[1:a,],2,sd), 5))
      cat("______eta3__________",'\n')
      print(round(apply(rst.eta3[1:a,],2,mean), 5))
      print(round(apply(rst.eta3[1:a,],2,sd), 5))
      cat("______eta4__________",'\n')
      print(round(apply(rst.eta4[1:a,],2,mean), 5))
      print(round(apply(rst.eta4[1:a,],2,sd), 5))
    }
  }
}

