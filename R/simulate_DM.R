# This function can be used to simulate DM data.
# Taken from Wadsworth (2017) An integrative Bayesian Dirichlet-multinomial regression model for the analysis of taxonomic abundances in microbiome data
simulate_DM = function(                              n_obs = 100,
                                                     n_vars = 30,
                                                     n_taxa = 75,
                                                     n_relevant_vars = 4,
                                                     n_relevant_taxa = 4,
                                                     beta_min = 1,
                                                     beta_max = 1.25,
                                                     signoise = 1.0,
                                                     n_reads_min = 5000,
                                                     n_reads_max = 10000,
                                                     theta0 = 0.01,
                                                     rho = NULL,
                                                     Sigma = NULL ){
  
  # check for required packages
  if(!require(dirmult)){
    stop("dirmult package required")
  }
  if(!require(MASS)){
    stop("MASS package required")
  }
  if(!require(matrixcalc)){
    stop("matrixcalc package required")
  }
  
  # Defense
  if( !is.null( rho ) & !is.null( Sigma ) ){
    stop("Bad input: Please provide either rho or Sigma for covariate correlation structure.")
  }
  
  if( is.null( rho ) & is.null( Sigma ) ){
    stop("Bad input: Please provide either rho or Sigma for covariate correlation structure.")
  }
  
  if( !is.null( rho ) ){
    if( rho > 1 | rho < 0 ){
       stop("Bad input: Please provide rho between 0 and 1.")
    }
  }
  
  if( !is.null( Sigma ) ){
    if( !is.positive.definite( Sigma ) ){
      stop("Bad input: Please provide positive definite covariance matrix.")
    }
  }
  
  if( !is.null( Sigma ) ){
    if( ncol(Sigma) != n_vars ){
      stop("Bad input: Please provide covariance matrix to match the number of covariates")
    }
  }  
  
  # covariance matrix for predictors
  if( !is.null( rho ) ){
    Sigma <- matrix( 0, n_vars, n_vars )
    Sigma = rho^abs(row(Sigma) - col(Sigma))
  }
  
 
  # include the intercept
  XX = cbind(rep(1, n_obs),
             scale(MASS::mvrnorm(n = n_obs, mu = rep(0, n_vars), Sigma = Sigma)))
  # empties
  YY = matrix(0, n_obs, n_taxa)
  betas = matrix(0, n_taxa, n_vars)
  phi = matrix(0, n_obs, n_taxa)
  # parameters with signs alternating
  st = 0
  low_side = beta_min
  high_side = beta_max
  if(n_relevant_taxa != 1){
    # warning if the lengths don't match
    coef = suppressWarnings(seq(low_side, high_side, len = n_relevant_taxa) * c(1, -1))
  }else{
    coef = (low_side + high_side) / 2
  }
  coef_g = rep(1.0, len = n_relevant_vars)
  for(ii in 1:n_relevant_vars){
    # overlap species
    betas[(st:(st + n_relevant_taxa - 1)) %% n_taxa + 1, 3 * ii - 2] = coef_g[ii] * sample(coef)[((ii - 1):(ii + n_relevant_taxa - 2)) %% n_relevant_taxa + 1]
    st = st + 1
  }
  # -2.3 and 2.3 so that the intercept varies over three orders of magnitude
  intercept = runif(n_taxa, -2.3, 2.3)
  Beta = cbind(intercept, signoise * betas)
  # row totals
  ct0 = sample(n_reads_min:n_reads_max, n_obs, rep = T)
  for(ii in 1:n_obs){
    thisrow = as.vector(exp(Beta %*% XX[ii, ]))
    phi[ii, ] = thisrow/sum(thisrow)
    YY[ii, ] = dirmult::simPop(J = 1, n = ct0[ii], pi = phi[ii, ], theta = theta0)$data[1, ]
  }
  
  return(list(X = XX, Y = YY, alphas = intercept, betas = Beta,
              n_reads_min = n_reads_min, n_reads_max = n_reads_max,
              theta0 = theta0, phi = phi, rho = rho, signoise = signoise, Sigma = Sigma))
  
}
