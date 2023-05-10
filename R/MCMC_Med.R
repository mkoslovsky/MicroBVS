MCMC_Med <- function( iterations = 5000, thin = 10, trt = NULL, Y = NULL, Z = NULL, covariate = NULL, sigma_alpha = sqrt( 10 ), 
                      sigma_beta = sqrt( 10 ), sigma_phi = sqrt( 10 ), h_alpha = 1, h_beta = 1, a = 1, b = 1,
                         a_0 = 1, b_0 = 1,  seed = 1, taxa = NULL ){
  library(mvtnorm)
  library(MCMCpack)
  library(Rcpp)
  library(ggplot2)
  library(RcppArmadillo)
  
  set.seed(seed)
  
  # Defense 
  if( is.null( trt ) | is.null( Y ) | is.null( Z ) | is.null( taxa ) ){
    stop("Bad input: Data must be supplied!")
  }
  
  if( taxa %%1 != 0){
    stop("Bad input: Taxon index should be an integer!")
  }
  
  B_sim <-  ncol( Z ) 
  
  if( taxa < 0 | taxa > B_sim ){
    stop("Bad input: Taxon index does not exist in the data provided!")
  }
  
  # Get number of samples
  samples <- floor( iterations/thin )
  
  if( samples < 100 ){
    stop("Bad input: Reduce the amount of thinning or increase number of iterations!")
  }
  
  # Prep the data for inference 
  # Set the covariate space (includes treatment and covariates)
  x <- as.matrix( trt, nrow = length( trt ) )
  
  if( !is.null( covariate ) ){
    covariate <- as.matrix( covariate, nrow = length(trt) )
    x <- cbind( x, covariate )
  }
 
  # Make y a matrix 
  Y <- as.matrix( Y, nrow = length( trt ) )
  
  # If zero count in Z, add 0.5. 
  Z[ Z == 0 ] <- 0.5
  
  # Rearrange so that the mediation term to investigate is in the first column of Z
  Z_insert <- Z[ , taxa ]
  Z_rest <- Z[ , -taxa ]
  Z_temp <- cbind( Z_insert, Z_rest)
  
  # Get dimensions of data to prep output memory allocation 
  covariates_sim <- ncol( x )
  subject_sim <- length( Y )
 
  # Initialize MCMC 
  alpha <- matrix( 0, nrow = B_sim, ncol = samples )
  phi <- array( 0, dim = c( B_sim, covariates_sim, samples ) )
  zeta <- array( 0, dim = c( B_sim, covariates_sim, samples ) )
  psi <- array( 0, dim = c( subject_sim, B_sim, samples ) )
  temp_cc <- Z_temp
  temp_uu <- rep( 0 , nrow( Z_temp ) )
  for( n in 1:nrow( x ) ){
    sum_Z <- sum( Z_temp[ n, ] )
    temp_uu[ n ] <- rgamma( 1, sum_Z, sum_Z );
  } 
  sigma2 <- rep( 0, samples )
  sigma2[ 1 ] <- 0.1
  xi <- matrix( 0, nrow = c( B_sim + covariates_sim ), ncol = samples )  
  xi[1,] <- 1 # Include intercept always
  beta2 <-  matrix( 0, nrow = c( B_sim + covariates_sim ), ncol = samples )
  a_m <- a
  b_m <- b
  
  output <- dm_lm_mediation( iterations, thin, alpha, Y, Z_temp, x, phi, psi, temp_cc, temp_uu, sigma2, sigma_alpha,
                             sigma_beta, zeta, xi, beta2, sigma_phi, a, b, a_0, b_0, h_alpha, h_beta, a_m, b_m)
  
  # Return the MCMC output and other inputs for reference
  
  return(list( alpha = output[[1]], zeta = output[[2]], phi = output[[3]], 
               psi = output[[4]], xi = output[[5]], beta = output[[6]], sigma2 = output[[7]], 
               sigma_alpha = sigma_alpha, sigma_beta = sigma_beta, a = a, b = b, a_m = a_m,
               b_m = b_m, h_alpha = h_alpha, h_beta = h_beta, a_0 = a_0, b_0 = b_0, Y = Y, Z = Z_temp, covariates = x,iterations = iterations, thin = thin ) )
  
}
