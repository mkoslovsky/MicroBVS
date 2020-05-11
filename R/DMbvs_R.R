# Wrapper function for the Rcpp code to initiate with defaults and simulate data if requested for DTMbvs package
# DTMbvs: Dirichlet Multinomial Regression Models with Bayesian Variable Selection for Microbiome Data - an R Package
# This runs a version of Wadsworth (2017) An integrative Bayesian Dirichlet-multinomial regression model for the analysis of taxonomic abundances in microbiome data
# This is much faster than DTMbvs for special case of #taxa = #branches in the tree

# Wrapper function for the Rcpp code to initiate with defaults and simulate data if requested
DMbvs_R <- function( iterations = 20000, thin = 10, z = NULL, x = NULL, alpha = NULL, phi = NULL, zeta = NULL,
                         sigma2_alpha = sqrt( 10 ), sigma2_phi = sqrt( 10 ), a = 1, b = 9, warmstart = T, seed = 1 ){
  library(mvtnorm)
  library(DMLMbvs)
  library(MCMCpack)
  library(Rcpp)
  library(ggplot2)
  library(RcppArmadillo)
  
  # iterations - Number of MCMC samples, Default = 20000
  # thin - Then MCMC by # thin, Default = 10
  # z - subject x part matrix of multivariate count data
  # x - subject x covariate matrix of measures
  # alpha - part X 1 vector of initial intercept values
  # phi - part X covariate matrix of initial regression coefficients
  # sigma2_alpha - prior value for alpha variance, Default = sqrt(10)
  # sigma2_phi - prior value for phi variance, Default = sqrt(10)
  # a -  parameter for beta prior for covariate inclusion probability, Default = 1
  # b -  parameter for beta prior for covariate inclusion probability, Default = 9
  # seed - set random seed for simulated data, Default = 1
  
  
  # Set seed for replication
  set.seed( seed )
  
  # Simulate data 
  if( is.null( x ) | is.null( z ) ){
    stop("Bad input: Data must be supplied")
  }
  
  # Adjust inputs if x,z are provided
  
  B_sim <- ncol( z ) 
  covariates_sim <- ncol( x ) 
  
  # Initiate starting values and allocate memory 
  samples <- floor( iterations/thin )
  
  # Intercept term alpha_j
  alpha. <- matrix( 0, nrow = B_sim, ncol = samples )
  
  # Inclusion indicators zeta_jp
  zeta. <- array( 0, dim = c( B_sim, covariates_sim, samples ) )
  
  # Regression Coefficients phi_jp
  phi. <- array( 0, dim = c( B_sim, covariates_sim, samples ) )
  
  # Adjust inital values for alpha, zeta, phi, psi, and xi if they are still NULL
  alpha.[ , 1] <- if( is.null( alpha ) ){ rnorm( n = 10 )}else{ alpha }  
  zeta.[ , , 1] <- ifelse( is.null( zeta ), 0, zeta )
  phi.[ , , 1]  <-  if( is.null( phi )){ rnorm( n = B_sim*covariates_sim )*zeta.[,,1]}else{ phi }  
  
  
  # initialize latent variables
  cc. <- z 
  uu <- rep( 0 , nrow( z ) ) 
  for( n in 1:nrow( x ) ){
    sum_Z <- sum( z[ n, ] )
    uu[ n ] <- rgamma( 1, sum_Z, sum_Z );
  }
  
  if( warmstart == TRUE ){
    alpha.[ , 1] <-   matrix( scale( log( colSums( z ) ) )  )
    
    cormat = matrix(0, B_sim, covariates_sim)
    pmat = matrix(0, B_sim, covariates_sim)
     zz = z/rowSums(z) # compositionalize

    for(rr in 1:B_sim){
      for(cc in 1:covariates_sim){
        pmat[rr, cc] = stats::cor.test( x[, cc], zz[, rr], method = "spearman", exact = F)$p.value
        cormat[rr, cc] = stats::cor( x[, cc], zz[, rr], method = "spearman")
      }
    }

    # defaults to 0.2 false discovery rate
    pm <- matrix((stats::p.adjust(c(pmat), method = "fdr") <= 0.2) + 0, B_sim, covariates_sim)
    betmat <- cormat * pm

    phi.[ , , 1] <- matrix( betmat , B_sim, covariates_sim)
    zeta.[ , , 1] <- matrix( (betmat != 0)*1, B_sim, covariates_sim)
    
  }
  
  # Run model
  output <- dm_bvs( iterations, thin, alpha., z, x, phi., cc., uu, sigma2_alpha, zeta., sigma2_phi, a, b )
  names( output ) <- c( "alpha", "zeta", "phi" )
  return( output )
}
