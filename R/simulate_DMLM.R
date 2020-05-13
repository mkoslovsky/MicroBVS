# Code to simulate data for DMLMbvs model

simulate_DMLM <- function( subject_sim = 50, B_sim = 50, covariates_sim = 50, active_cov = 10, covar = 0.4, seed = 111 ){
  
  # Call libraries  
  library(mvtnorm) 
  library(MCMCpack)
  
  # Set seed 
  set.seed( seed )
  
  # Set covariance strucuture for covariates 
  sig <- diag(covariates_sim)
  for( i in 1:covariates_sim ){
    for( j in 1:covariates_sim ){
      if( i != j){
        sig[ i , j] = covar^abs(i - j)
      }
    }
  }
  
  X <- rmvnorm( subject_sim, rep( 0, covariates_sim ), sig )
  
  zeta_sim <- matrix( 0, B_sim, covariates_sim)
  true_cov <- cbind( sample( seq( 1, active_cov ), active_cov , replace = T ),sample( seq( 1, covariates_sim ), active_cov) )
  zeta_sim[  true_cov ] <- 1
  alpha_sim <- matrix( 1, nrow = subject_sim , ncol = 1)%*%( runif( n = B_sim, -2.3,2.3 ) )
  true_coeff <- runif( 10, 0.75, 1.25)
  phi_sim <- matrix( 0, B_sim, covariates_sim)
  phi_sim[ true_cov ] <- true_coeff*sample(c(1,-1), 10, replace = TRUE)
  
  inside_sim <- exp( alpha_sim + X%*%t(phi_sim) )
  psi_sim <- inside_sim/rowSums(inside_sim)
  psi_sim_overdispersed <- psi_sim*( 1 - 0.01 )/0.01
  
  # Simulate probabilities for each branch
  prob_sim <- apply( psi_sim_overdispersed, 1,  function(x){ rdirichlet(1,x) } )
  
  # Simulate count data for each subtree branch (simplified for bifurcating tree in this example)
  Z <- t( apply( t( prob_sim ), 1, function(x){ rmultinom(1, sample( seq(2500, 7500), 1 ), x) } ) )
  
  # adjust for zero counts
  
  Z_norm = Z
  Z_norm[Z_norm == 0] <- 0.5
  Z_norm <- Z_norm/rowSums(Z_norm)
  
  prob_sim [ prob_sim  < 0.5/7500 ] <- 0.5/7500
  prob_sim2 <-  t(prob_sim)/colSums(prob_sim)
  
  # Calculate Balances 
  Balances <- ilr_fun_cpp( prob_sim2  )
  
  # Select Balances
  xi_sim <- rep( 0, B_sim - 1 )
  true_cov_beta <-  seq(1:10)  
  xi_sim[  true_cov_beta ] <- 1
  true_coeff_beta <- runif(length( true_cov_beta ), 1.25, 1.75)
  beta_sim <- rep( 0, B_sim - 1)
  beta_sim[ true_cov_beta ] <- true_coeff_beta*sample(c(1 ,-1 ), length( true_cov_beta ) , replace = TRUE)
  
  Y <-  Balances%*%beta_sim + rnorm( subject_sim )
  
  return( list( Y = Y, Z = Z, X = X, true_cov = true_cov, true_coeff = true_coeff, true_coeff_beta = true_coeff_beta))
}