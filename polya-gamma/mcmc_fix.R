# REGRESSION MCMC - SAMPLER FUNCTION
mcmc_fix <- function (y, locs, K, message = 100,
                      adapt = K / 2, ## when to stop tuning the MCMC
                      ## default priors
                      alpha_tau = 1/2,
                      beta_tau = 10,
                      mean_nu = -1,
                      sd_nu = 1,
                      mean_range = 0,
                      sd_range = 10,
                      mu_sigma_prop=0.1) {
  require(pgdraw)
  require(mvnfast)
  ## define priors
  theta_mean <- c(mean_range, mean_nu)
  theta_sd <- c(sd_range, sd_nu)
  N <- nrow(y)
  J <- ncol(y)
  D <- fields::rdist(locs)
  ########### Calculate Mi ###################
  Mi <- data.frame(matrix(0, N, J-1))
  sumY <- rep(0, times = N)
  cumsumY <- data.frame(matrix(0, N, J-1))
  
  for(i in 1: N){
    sumY[i] <- sum(y[i, ])
  }
  
  for(i in 1: N){
    cumsumY[i,] <- c(0, cumsum(y[i,1:(J-2)]))
  }
  
  for(i in 1: N){
    Mi[i,] <- sumY[i] - cumsumY[i,]
  }
  
  Mi <- round(Mi)
  Mi <- apply(Mi, 2, function(x) as.integer(x))
  Mi <- data.frame(Mi)
  
  ####################initialize kappa###################
  kappa <- data.frame(matrix(0, N, J-1))
  for (i in 1: N) {
    kappa[i,] <- y[i, 1:(J-1)] - Mi[i, ] / 2 
  }
  
  ## initialize tau2 using an inverse-gamma(1, 1) prior
  tau2 <- rep(0, K)
  tau <- rep(0, K)
  tau2[1] <- rinvgamma(1, shape=alpha_tau, scale=beta_tau) #min(1 / rgamma(1, alpha_tau, beta_tau), 10)
  tau[1] <- sqrt(tau2[1])
  
  tau2_star_sd   = 0.6
  
  
  # ## initialize theta
  # theta <- matrix(0, K, 2)
  # theta[1, ] <- c(rmvn(1, theta_mean, diag(theta_sd)))
  
  ## initialize theta
  # theta <- matrix(0, K, 2)
  theta <- array(0, c(J-1, K, 2))
  for (j in 1:(J-1)){
    theta[j,1,] <- c(rmvn(1, theta_mean, diag(theta_sd)))
  }
  
  
  ## initialize mu
  mu <- matrix(0, K, J-1)
  mu0 = 0#rep(0, N)
  mu0_Sigma = 10
  
  mu[1,] <- rnorm(J-1, mu0, mu0_Sigma)
  
  ## construct Sigma covariance matrix
  # MAKE SURE IT'S INVERTIBLE
  diag(D) <- NA
  any(D == 0, na.rm = TRUE)   # check if there are off-diagonal zeros
  D <- ifelse(D == 0, 0.007, D)  # remove them
  diag(D) <- 0
  
  Sigma = array(0, c(J-1,N,N))
  Sigma_chol = array(0, c(J-1,N,N))
  Sigma_inv = array(0, c(J-1,N,N))
  for (j in 1:(J-1)){
    Sigma[j,,] <- tau2[1] * correlation_function(D, theta[j,1, ])
    Sigma_chol[j,,] <- chol(Sigma[j,,])
    Sigma_inv[j,,] <- chol2inv(Sigma_chol[j,,])
  }
  
  
  ## don't save eta due to memory constraints for larger datasets
  # eta <- matrix(0, N, J-1)
  eta <- array(0, dim = c(K, N, J-1))
  for (j in 1:(J-1)){
    eta[1, ,j] <- t(rmvn(1, rep(mu[1,j], N), Sigma_chol[j,,], isChol = TRUE))
  }
  
  ##
  ## initialize omega
  ##
  
  # omega <- matrix(0, N, J-1)
  omega <- array(0, dim=c(K, N, J-1))
  for (i in 1:N) {
    for (j in 1:(J-1)) {
      if (Mi[i, j] != 0) {
        omega[1, i, j] <- pgdraw(Mi[i, j], eta[1, i, j])
      }
    }
  }
  
  Omega <- array(0, dim=c(K, N, J-1))#vector(mode = "list", length = J-1)
  for (j in 1:(J - 1)) {
    Omega[1,,j] <- omega[1, , j] #diag(omega[1, , j])
  }
  
  ##
  ## tuning variables for adaptive MCMC
  ##
  theta_batch <- matrix(0, 50, 2)
  theta_accept <- 0
  theta_accept_batch <- 0
  lambda_theta <- 0.05
  Sigma_theta_tune <- 1.8 * diag(2) - .8
  Sigma_theta_tune_chol <- chol(Sigma_theta_tune)
  
  
  accept_mu  <- 0
  accept_tau <- 0
  
  ## function to update the parameter
  k = K
  accept = theta_accept_batch
  lambda = lambda_theta
  batch_samples = theta_batch
  Sigma_tune = Sigma_theta_tune
  Sigma_tune_chol = Sigma_theta_tune_chol
  
  message("Starting MCMC, fitting for ", K, " iterations")
  
  for (k in 2:K) {
    print(paste("Iteration", k, sep=" "))
    # if (k %% message == 0) {
    #   message("Fitting iteration ", k, " out of ", K)
    # }
    
    ##
    ## sample Omega 
    ## GIBBS STEP
    ##
    for (i in 1:N) {
      for (j in 1:(J-1)) {
        # if(!is.na(Mi[i, j])) {
        if(Mi[i, j] != 0){
          omega[k, i, j] <- pgdraw(Mi[i, j], eta[k-1, i, j])
        }
        else {
          omega[k, i, j] <- 0
        }
      }
    }
    
    for (j in 1:(J-1)) {
      Omega[k,,j] <- omega[k, , j]#diag(omega[k, , j])
    }
    
    ##
    ## sample eta
    ## GIBBS STEP
    ##
    theta_star = array(NA, c(J-1, 2))
    Sigma_star = array(0, c(J-1,N,N))
    Sigma_chol_star = array(0, c(J-1,N,N))
    Sigma_inv_star  = array(0, c(J-1,N,N))
    mu_tilde <- array(0, dim = c(K, N, J-1))
    Sigma_tilde = array(0, c(J-1,N,N))
    Sigma_tilde_star = array(0, c(J-1,N,N))
    mu_tilde_star = array(0, c(N, J-1))
    for (j in 1:(J-1)) {
      
      # not sure if we need these
      Sigma[j,,] <- tau[k-1]^2 * correlation_function(D, theta[j, k-1, ])
      Sigma_inv[j,,] <- solve(Sigma[j,,])
      
      ## can make this much more efficient
      Sigma_tilde[j,,] <- chol2inv(chol(Sigma_inv[j,,] + diag(Omega[k,,j])))
      mu_tilde[k,,j] <- c(Sigma_tilde[j,,] %*% (Sigma_inv[j,,] %*% rep(mu[k-1,j], N) + kappa[, j]))
      
      eta[k, , j] <- rmvn(1, mu_tilde[k,,j], Sigma_tilde[j,,])  # k
    }
    
    
    ##
    ## sample theta (spatial correlation parameters)
    ## MCMC STEP
    ##
    
    for (j in 1:(J-1)) {
      theta_star[j,] <- rmvn(
        n = 1,
        mu = theta[j, k-1, ],
        sigma = lambda_theta * Sigma_theta_tune_chol,
        isChol = TRUE
      )
      Sigma_star[j,,] <- tau[k-1]^2 * correlation_function(D, theta_star[j,])
      
      Sigma_chol_star[j,,] <- chol(Sigma_star[j,,])
      Sigma_inv_star[j,,]  <- chol2inv(Sigma_chol_star[j,,])
    }
    
    mh1 <- sum(
      sapply(
        1:(J-1),
        function(j) {
          dmvn(eta[k, , j],  rep(mu[k-1,j], N), Sigma_star[j,,], isChol = FALSE, log = TRUE)
        }
      )
    ) +
      ## prior
      sum(
        sapply(
          1:(J-1),
          function(j){ 
            dmvn(theta_star[j,], theta_mean, diag(theta_sd), log = TRUE)
          }
        )
      )
    
    mh2 <- sum(
      sapply(
        1:(J-1),
        function(j) {
          dmvn(eta[k, , j], rep(mu[k-1,j], N), Sigma[j,,], isChol = FALSE, log = TRUE)
        }
      )
    ) +
      ## prior
      sum(
        sapply(
          1:(J-1),
          function(j){ 
            dmvn(theta[j,k-1,], theta_mean, diag(theta_sd), log = TRUE) 
          }
        )
      )
    
    mh <- exp(mh1 - mh2)
    
    if(mh > runif(1, 0, 1)) {
      for (j in 1:(J-1)){
        theta[j,k, ] <- theta_star[j,]
        Sigma[j,,] <- Sigma_star[j,,]
        Sigma_chol[j,,] <- Sigma_chol_star[j,,]
        Sigma_inv[j,,] <- Sigma_inv_star[j,,]
        # if(k <= adapt) {
        #   theta_accept_batch <- theta_accept_batch + 1 / 50
        # } else {
        #   theta_accept <- theta_accept + 1 / (K - adapt)
        #   theta[j,k, ] <- theta[j,k-1, ]
        # }
      }
    } else {
      for (j in 1:(J-1)){
        theta[j,k, ] <- theta[j,k-1, ]
      }
    }
    
    # COMMENTED FOR TESTING
    # ## adapt the tuning
    # if (k <= adapt) {
    #   theta_batch[k %% 50, ] <- theta[k, ]#}
    #   if (k %% 50 == 0) {
    #     out_tuning <- updateTuningMV(
    #       k = k,
    #       accept = theta_accept_batch,
    #       lambda = lambda_theta,
    #       batch_samples = theta_batch,
    #       Sigma_tune = Sigma_theta_tune,
    #       Sigma_tune_chol = Sigma_theta_tune_chol
    #     )
    #     theta_accept_batch <- out_tuning$accept
    #     lambda_theta_tune <- out_tuning$lambda
    #     theta_batch <- out_tuning$batch_samples
    #     Sigma_theta_tune <- out_tuning$Sigma_tune
    #     Sigma_theta_tune_chol <- out_tuning$Sigma_tune_chol
    #   }
    # }
    
    ##
    ## sample spatial process variance tau^2
    ## MCMC STEP
    ##
    
    SS <- 0
    for (j in 1:(J-1)){
      devs <- eta[k, ,j] - mu[k-1,j]
      devs <- as.matrix(devs)
      # devs <- kappa - eta[k, , ]
      # doesn't actually depend on tau[k-1] because Sigma_Inv has a 1/tau[k-1] factor
      # SS <- sum(devs * (tau[k-1]^2 * Sigma_inv[j,,] %*% as.matrix(devs)))
      SS <- SS + t(devs) %*% (tau[k-1]^2 * Sigma_inv[j,,] %*% devs)
    }
    # tau2[k] <- 1 / rgamma(1, N * (J - 1) / 2 + alpha_tau,
    #                       SS / 2 + beta_tau)
    tau2[k] <- rinvgamma(1, N * (J - 1) / 2 + alpha_tau,
                         SS / 2 + beta_tau)
    tau[k] <- sqrt(tau2[k])
    for (j in 1:(J-1)){
      Sigma[j,,] <- tau[k]^2 * correlation_function(D, theta[j,k, ])
      Sigma_chol[j,,] <- chol(Sigma[j,,])
      Sigma_inv[j,,] <- chol2inv(Sigma_chol[j,,])
    }
    
    # tau2_star_mean = log(tau2[k-1])
    # tau2_star = rlnorm(1, meanlog=tau2_star_mean, sdlog=tau2_star_sd)
    # 
    # for (j in 1:(J-1)){
    #   Sigma_star[j,,] <- tau2_star * correlation_function(D, theta[j,k,])
    #   Sigma_chol_star[j,,] <- chol(Sigma_star[j,,])
    #   Sigma_inv_star[j,,]  <- chol2inv(Sigma_chol_star[j,,])
    # }
    # 
    # tau_mh1 <- sum(
    #   sapply(
    #     1:(J-1),
    #     function(j) {
    #       dmvn(eta[k, , j], rep(mu[k-1,j], N), Sigma_star[j,,], isChol = FALSE, log = TRUE)
    #     }
    #   )
    # ) +
    #   ## prior
    #   dinvgamma(tau2_star, shape=alpha_tau, scale=beta_tau, log = TRUE) +
    #   dlnorm(tau2[k-1], meanlog=log(tau2_star), sdlog=tau2_star_sd, log=TRUE)
    # 
    # tau_mh2 <- sum(
    #   sapply(
    #     1:(J-1),
    #     function(j) {
    #       dmvn(eta[k, , j], rep(mu[k-1,j], N), Sigma[j,,], isChol = FALSE, log = TRUE)
    #     }
    #   )
    # ) +
    #   ## prior
    #   dinvgamma(tau2[k-1], shape=alpha_tau, scale=beta_tau, log = TRUE) +
    #   dlnorm(tau2_star, meanlog=tau2_star_mean, sdlog=tau2_star_sd, log=TRUE)
    # 
    # tau_mh <- exp(tau_mh1 - tau_mh2)
    # 
    # if(tau_mh > runif(1, 0, 1)) {
    #   tau2[k]   <- tau2_star
    #   accept_tau <- accept_tau + 1
    #   print("Accept tau")
    # } else {
    #   tau2[k] <- tau2[k-1]
    # }
    # 
    # tau[k] <- sqrt(tau2[k])
    
    ##
    ## sample spatial process variance tau^2
    ## MCMC STEP (but can write down full conditional I think)
    ##
    
    mu_prop = rep(0,J-1)
    
    for (j in 1:(J-1)){
      mu_prop[j] = rnorm(1, mu[k-1,j], mu_sigma_prop)
      
      Sigma[j,,] <- tau2[k] * correlation_function(D, theta[j,k,])
      Sigma_chol[j,,] <- chol(Sigma[j,,])
      Sigma_inv[j,,] <- chol2inv(Sigma_chol[j,,])
   }  
    mu_mh1 = sum(sapply(1:(J-1), function(j){
      dmvn(eta[k,,j],  rep(mu_prop[j], N), Sigma[j,,], log=TRUE)})) + 
      sum(dnorm(mu_prop, mu0, mu0_Sigma, log=TRUE))
    mu_mh2 = sum(sapply(1:(J-1), function(j){
      dmvn(eta[k,,j], rep(mu[k-1,j], N), Sigma[j,,], log=TRUE)})) + 
      sum(dnorm(mu[k-1,], mu0, mu0_Sigma, log=TRUE))
    
    mu_mh <- exp(mu_mh1 - mu_mh2)
    
    if(mu_mh > runif(1, 0, 1)) {
      mu[k,]   <- mu_prop
      accept_mu <- accept_mu + 1
    } else {
      mu[k,] <- mu[k-1,]
    }
    
    ##
    ## update proposal variance for mu and tau
    ##
    
    if (k %% 50 == 0) {
      print(paste0("Acceptance rate for mu is: ", accept_mu/k))
      print(paste0("Acceptance rate for tau is: ", accept_tau/k))
      if ((accept_mu/k)<0.4){
        mu_sigma_prop = mu_sigma_prop*1.10
        # } else if ((accept_mu/k)>0.4) {
        #   
      }
      if ((accept_tau/k)<0.4){
        # tau2_star_sd = tau2_star_sd*1.10
        tau2_star_sd = tau2_star_sd*0.90
        # } else if ((accept_mu/k)>0.4) {
        #   
      }
    }
  }
  # }
  
  message("Acceptance rate for theta is ", theta_accept)
  
  return(
    list(
      tau = tau,
      theta = theta,
      eta = eta,
      omega = omega,
      mu = mu,
      mu_tilde=mu_tilde,
      Sigma_tilde=Sigma_tilde
    )
  )
}


updateTuningMV <- function(k, accept, lambda, batch_samples,
                           Sigma_tune, Sigma_tune_chol) {
  arr <- c(0.44, 0.35, 0.32, 0.25, 0.234)
  # std::vector<double> acceptance_rates (arr, arr + sizeof(arr) / sizeof(arr[0]))
  dimension <- nrow(batch_samples)
  if (dimension >= 5) {
    dimension <- 5
  }
  d <- ncol(batch_samples)
  batch_size <- nrow(batch_samples)
  optimal_accept <- arr[dimension]
  times_adapted <- floor(k / 50)
  gamma1 <- 1.0 / ((times_adapted + 3.0)^0.8)
  gamma2 <- 10.0 * gamma1
  adapt_factor <- exp(gamma2 * (accept - optimal_accept))
  ## update the MV scaling parameter
  lambda_out <- lambda * adapt_factor
  ## center the batch of MCMC samples
  batch_samples_tmp <- batch_samples
  for (j in 1:d) {
    mean_batch = mean(batch_samples[, j])
    for (i in 1:batch_size) {
      batch_samples_tmp[i, j] <- batch_samples[i, j] - mean_batch
    }
  }
  
  ## 50 is an MCMC batch size, can make this function more general later...
  Sigma_tune_out <- Sigma_tune + gamma1 *
    (t(batch_samples) %*% batch_samples / (50.0-1.0) - Sigma_tune)
  Sigma_tune_chol_out <- chol(Sigma_tune)
  accept_out <- 0.0
  batch_samples_out <- matrix(0, batch_size, d)
  return(list(
    batch_samples = batch_samples_out,
    Sigma_tune = Sigma_tune_out,
    Sigma_tune_chol = Sigma_tune_chol_out,
    lambda = lambda_out,
    accept = accept_out
  ))
}