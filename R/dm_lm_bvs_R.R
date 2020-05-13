# Wrapper function for the Rcpp code to initiate with defaults and simulate data if requested
dm_lm_bvs_R <- function( iterations = 10000, thin = 10, y = NULL, z = NULL, x = NULL, alpha = NULL, phi = NULL,
                         psi = NULL, zeta = NULL, xi = NULL, sigma2_alpha = sqrt( 10 ), sigma2_phi = sqrt( 10 ),
                         h_alpha = 1, h_beta = 1, a_m = 1, b_m = 9, a = 1, b = 9,
                         a_0 = 2, b_0 = 2, seed = 1 ){
  library(mvtnorm)
  library(DMLMbvs)
  library(MCMCpack)
  library(Rcpp)
  library(ggplot2)
  library(RcppArmadillo)
  
  # iterations - Number of MCMC samples, Default = 20000
  # thin - Then MCMC by # thin, Default = 10
  # y - subject x 1 vector of continuous response
  # z - subject x part matrix of multivariate count data
  # x - subject x covariate matrix of measures
  # alpha - part X 1 vector of initial intercept values
  # phi - part X covariate matrix of initial regression coefficients
  # psi - subject X part matrix of initial compositional probabilities
  # zeta - part X covariate matrix of initial inclusion indicators for covariates ( 0 or 1 )
  # xi - part - 1 X 1 vector of initial inclusion indicators for balances
  # sigma2_alpha - prior value for alpha variance, Default = sqrt(10)
  # sigma2_phi - prior value for phi variance, Default = sqrt(10)
  # h_alpha - prior variance for intercept term a_0, Default = 1
  # h_beta - prior variance for beta terms for balances, Default = 1
  # a_m -  parameter for beta prior for balance inclusion probability, Default = 1
  # b_m -  parameter for beta prior for balance inclusion probability, Default = 1
  # a -  parameter for beta prior for covariate inclusion probability, Default = 1
  # b -  parameter for beta prior for covariate inclusion probability, Default = 1
  # a_0 - shape parameter for inverse-gamma prior for sigma, Default = 2
  # b_0 - scale parameter for inverse-gamme prior for sigma, Default = 2
  # seed - set random seed for simulated data, Default = 1
  
  
  # Set seed for replication
  set.seed( seed )
  
  # Simulate data 
  if( is.null( x ) | is.null( y ) | is.null( z ) ){
    stop("Bad input: Data must be supplied")
  }
  
  # Adjust inputs if x,y,z are provided
  Z_norm = z
  Z_norm[ Z_norm == 0 ] <- 0.5
  
  Z_norm <- Z_norm/rowSums( Z_norm )
  
  B_sim <-  ncol( z ) 
  covariates_sim <- ncol( x )
  subject_sim <- length( y )
  
  # Initiate starting values and allocate memory 
  samples <- floor( iterations/thin )
  
  # Intercept term alpha_j
  alpha. <- matrix( 0, nrow = B_sim, ncol = samples )
  
  # Inclusion indicators zeta_jp
  zeta. <- array( 0, dim = c( B_sim, covariates_sim, samples ) )
  
  # Regression Coefficients phi_jp
  phi. <- array( 0, dim = c( B_sim, covariates_sim, samples ) )
  
  # Psi for each i
  psi. <- array( 0, dim = c( subject_sim, B_sim, samples ) )
  
  # Inclusion indicators xi_m
  xi. <- matrix( 0, nrow = B_sim - 1, ncol = samples )
  
  # Adjust inital values for alpha, zeta, phi, psi, and xi if they are still NULL
  alpha.[ , 1] <- if( is.null( alpha ) ){ rnorm( n = 10 )}else{ alpha }  
  zeta.[ , , 1] <- ifelse( is.null( zeta ), 0, zeta )
  phi.[ , , 1]  <-  if( is.null( phi )){ rnorm( n = B_sim*covariates_sim )*zeta.[,,1]}else{ phi }
  psi.[ , , 1] <-  if( is.null( psi ) ){ as.matrix( Z_norm, nrow = subject_sim, ncol = B_sim) }else{ psi }
  xi.[ , 1] <- ifelse( is.null( xi ), 0, xi ) 
  
  
  # initialize latent variables
  cc <- z 
  uu <- rep( 0 , nrow( z ) ) 
  for( n in 1:nrow( x ) ){
    sum_Z <- sum( z[ n, ] )
    uu[ n ] <- rgamma( 1, sum_Z, sum_Z );
  }
  
  
  # Run model
  output <- dm_lm_bvs( iterations, thin, alpha., y, z, x, phi., psi., cc, uu, sigma2_alpha, zeta., xi., sigma2_phi, a, b, a_0, b_0, h_alpha, h_beta, a_m, b_m)
  
  
  return( output )
}