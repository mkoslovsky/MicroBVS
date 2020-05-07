# Wrapper function for the Rcpp code to initiate with defaults and simulate data if requested for DTMbvs package
# DTMbvs: Dirichlet Multinomial Regression Models with Bayesian Variable Selection for Microbiome Data - an R Package
# This runs a version of Wadsworth (2017) An integrative Bayesian Dirichlet-multinomial regression model for the analysis of taxonomic abundances in microbiome data
# This is much faster than DTMbvs for special case of #taxa = #branches in the tree

DMbvs_R <- function( iterations = 5000, thin = 10,  Y = NULL, X = NULL, MCMC = "Gibbs", alpha = NULL, beta = NULL, sigma2_alpha = 10, sigma2_beta = 10,
                      aa = 0.1, bb = 1.9, seed = 1212, warmstart = FALSE ){
  
  library(mvtnorm)
  library(MCMCpack) 
  
  # iterations - integer number of MCMC samples, default = 50000
  # thin - integer MCMC by # thin, default = 10
  # Y - N x branch matrix of count data.
  # X - N x P matrix of covariates
  # MCMC - name of MCMC type to use. Takes values "Gibbs" (default) and "SSVS" ( faster for each iteration, may need more iterations for convergence )
  # alpha branch x 1 vector of branch intercepts
  # beta branch x P vector of regression coefficients
  # sigma2_alpha - prior variance for alpha, default = 10
  # sigma2_beta - prior variance for phi, default = 10
  # aa - double beta-binomial hyperparameter
  # bb - double beta-binomial hyperparameter 
  # seed - set the seed, if you want 
  # warmstart - boolean if true, start model using informed initial values. 
  
  # Defense
  if( iterations%%1 != 0 | iterations <= 0){
    stop("Bad input: iterations should be a positive integer")
  }
  
  if( thin%%1 != 0 | thin < 0 ){
    stop("Bad input: thin should be a positive integer")
  }
  
  if( (MCMC != "Gibbs") + (MCMC != "SSVS")  == 2 ){
    stop("Bad input: prior is not in correct format")
  }
  
  # Set seed for replication
  set.seed( seed )

  n_cats = ncol( Y )
  n_vars = ncol( X )
  n_obs = nrow( Y )

  # Adjust inital values for alpha and beta if they are NULL
  alpha. <- if( is.null( alpha ) ){  rep( 0, n_cats )  }else{ alpha }
  beta. <- if( is.null( beta ) ){ rep(0, n_cats*n_vars ) }else{ beta }
  
  if( warmstart == TRUE ){
    alpha. <-   scale( log( colSums( Y ) ) )  
     
      cormat = matrix(0, n_cats, n_vars)
      pmat = matrix(0, n_cats, n_vars)
      yy = Y/rowSums(Y) # compositionalize
      for(rr in 1:n_cats){
        for(cc in 1:n_vars){
          pmat[rr, cc] = stats::cor.test(X[, cc], yy[, rr], method = "spearman",
                                         exact = F)$p.value
          cormat[rr, cc] = stats::cor(X[, cc], yy[, rr], method = "spearman")
        }
      }
      
      # defaults to 0.2 false discovery rate
      pm = matrix((stats::p.adjust(c(pmat), method = "fdr") <= 0.2) + 0, n_cats,
                  n_vars)
      betmat = cormat * pm
      betmat = c(t(betmat))
      beta. = betmat 
  }
  
  # Set everything else 
  mu_al = rep( 0, n_cats)
  sig_al = rep( sigma2_alpha, n_cats)
  mu_be = matrix( 0,n_cats, n_vars)
  sig_be = matrix( sigma2_beta, n_cats,n_vars)
  prop_per_alpha = rep( 0.5, n_cats)
  prop_per_beta = rep(0.5, (n_cats*n_vars))
  
  # Run DMbvs
  if( MCMC == "Gibbs" ){
    # Run model
    output <- dmbvs_gibbs( XX = X, YY = Y, alpha = alpha., beta = beta., mu_al = mu_al, sig_al = sig_al, mu_be = mu_be, sig_be = sig_be, aa_hp = aa, bb_hp = bb, prop_per_alpha = prop_per_alpha, prop_per_beta = prop_per_beta, iterations, thin = thin, 1 )
    names( output ) <- c("alpha", "beta")
    
    return( output )
  }else{
    output <- dmbvs_ss( XX = X, YY = Y, alpha = alpha., beta = beta., mu_al = mu_al, sig_al = sig_al, mu_be = mu_be, sig_be = sig_be, aa_hp = aa, bb_hp = bb, prop_per_alpha = prop_per_alpha, prop_per_beta = prop_per_beta, iterations, thin = thin, 1 )
    names( output ) <- c("alpha", "beta")
    
    return( output )
  }
}

