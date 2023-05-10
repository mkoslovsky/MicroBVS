Simulate_MedDM <- function( n_obs = 200, n_taxa = 50, bin_cov = 0, con_cov = 0, seed = 111, 
                           alpha_intercept_min = -2, alpha_intercept_max = 0.5, phi_min = 1, phi_max = 3, c0 = 0, trt_coeff = 1, 
                           covar_min = 1, covar_max = 3,  beta_coeff = NULL, active = NULL  ){
  library(mvtnorm) 
  library(MCMCpack)
  
  set.seed(seed)
  
  if( alpha_intercept_min > alpha_intercept_max ){
    stop("Bad input: Intercept minimumn should be less than maximum!")
  }
  
  if( bin_cov < 0 | con_cov < 0 | bin_cov %%1 !=0 | con_cov %%1 != 0 ){
    stop("Bad input: Number of covariate must be positive interger!")
  }
  
  if( length( beta_coeff ) != length( active ) ){
    stop("Bad input: Number of active mediators must match number of regression coefficients!")
  }
  
  if( sum( beta_coeff ) != 0 ){
    stop("Bad input: Regression coefficients must sum to one!")
  }
  
  I <- n_obs
  J <- n_taxa
  
  # Simulate intercept for DM model 
  alpha_intercept <- runif( J, alpha_intercept_min, alpha_intercept_max )
  
  # Simulate binary treatment
  sample_trt <- rbinom( I, 1, 0.5 ) 
  
  # If no active mediators are given, then set to 3rd, 4th, and 5th.
  if( is.null( active ) ){
    active <- c( 3, 4, 5 )
    beta_coeff <- c( 1.3, -0.7, -0.6 )
  }
  
  # Simulate treatment effect for each of the taxa in DM model 
  # If active, then specify non-zero effects in phi_min and phi_max range
  phi <- rep( 0 , J )
  phi[ active ] <- runif( length( active ), phi_min, phi_max )
  
  # Generate covariates for the DM and LM if number of terms > 0 
  # Specify the effects of the covariates (Note: same for each taxa j in DM)
  c2_bin  <- NULL
  phi_2_bin <- NULL
  if( bin_cov > 0 ){
    sample_bin_covar <- matrix( rbinom( bin_cov*I, 1, 0.5 ), nrow = I, ncol = bin_cov )
    phi_2_bin <- runif( bin_cov, covar_min, covar_max )
    c2_bin <- runif( bin_cov, covar_min, covar_max )
  }
  c2_con  <- NULL
  phi_2_con <- NULL
  if( con_cov > 0 ){
    sample_con_covar <-  matrix( rnorm( con_cov*I, 0, 1 ), nrow = I, ncol = con_cov )
    phi_2_con <- runif( con_cov, covar_min, covar_max )
    c2_con <- runif( con_cov, covar_min, covar_max )
  } 
  
  # Generate log concentration parameters for Dirichlet distribution 
  gamma_sim <- matrix( data = NA, nrow = I , ncol = J )
  for( i in 1:I ){
    for( j in 1:J ){
      # Incorporate intercept and treatment effect
      gamma_sim[ i , j ] <- alpha_intercept[ j ] + phi[ j ]*sample_trt[ i ] 
      
      # Incorporate covariates if necessary
      if( bin_cov > 0 ){
        gamma_sim[ i , j ] <- gamma_sim[ i , j ] + phi_2_bin%*%sample_bin_covar[ i, ] 
      }  
      if( con_cov > 0 ){
        gamma_sim[ i , j ] <- gamma_sim[ i , j ] + phi_2_con%*%sample_con_covar[ i ,  ]
      }
    }
  }
  
  # Calculate concentration parameter 
  gamma_sim <- exp( gamma_sim )
  
  # Based on concentration parameter , compute taxa proportions
  prob_sim <- apply( gamma_sim, 1,  function(x){ rdirichlet( 1, x ) } )
  prob_sim2 <-  t( prob_sim )/colSums( prob_sim )
  
  # Generate taxa counts
  Z <- t( apply( t( prob_sim ), 1, function(x){ rmultinom(1, sample( seq(2500, 7500), 1 ), x) } )) 
  
  # Set treatment effect in LM 
  c1 <- trt_coeff
  
  # Generate outcome
  y_1 <- c0 + c1*sample_trt + rnorm( I )
  
  # Add covariate effects
  if( bin_cov > 0 ){
    y_1 <- y_1 + sample_bin_covar%*%c2_bin
  }  
  if( con_cov > 0 ){
    y_1 <- y_1 + sample_con_covar%*%c2_con
  } 
  
  # Add the mediation effects to LM  
     y_1 <- y_1 + log( prob_sim2[ , active ] )%*%beta_coeff

  # Combine covariates for output
  sample_covar <- NULL
  if( bin_cov > 0 & con_cov == 0 ){
    sample_covar <- cbind( sample_bin_covar)
  }
  if( bin_cov == 0 & con_cov > 0 ){
    sample_covar <- cbind( sample_con_covar )
  }
  if( bin_cov > 0 & con_cov > 0 ){
    sample_covar <- cbind( sample_bin_covar, sample_con_covar )
  } 
    
  return( list( Y = y_1, Z = Z, proportion = prob_sim2, active = active, trt = sample_trt, covariate = sample_covar, beta_coeff = beta_coeff, c0 = c0, c1 = c1, c2_bin = c2_bin, c2_con = c2_con,  alpha_intercept = alpha_intercept, phi = phi, phi_2_bin = phi_2_bin, phi_2_con = phi_2_con, seed = seed ) )
}

