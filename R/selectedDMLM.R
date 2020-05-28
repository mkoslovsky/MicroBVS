# Processes output from dm_lm_bvs_R for inference
# Returns selected xi and zeta at given threshold of inclusion
# Option to produce marginal PPI plots and numbers of selected zeta and xi versus MCMC iteration
selected_DMLM <- function( dmlm_obj = NULL, threshold = c(0.5, 0.5), plotting = FALSE, burnin = 0, G = FALSE ){
  # dmlm_obj - output from dm_lm_bvs or dm_lm_bvs_R
  # threshold - 2 X 1 vector for posterior probability of inclusion threshold for zeta and xi (in order)
  # plotting - logical (True or False) indicator to plot # selected xi and zeta as well as PPI for zeta and xi 
  # burnin - number of MCMC samples to drop before inference, Default = 0
  # G - boolean if true, provide posterior edge inclusion 
  
  library(ggplot2)
  
  if( burnin%%1 != 0){
    stop("Bad input: burn-in should be an integer")
  }
  
  if( burnin < 0 ){
    stop("Bad input: burn-in should be positive")
  }
  
  if( length( threshold ) != 2 ){
    stop("Bad input: threshold must be a two-dimensional vector")
  }
  
  if( threshold[ 1 ] > 1 | threshold[ 1 ] < 0 | threshold[ 2 ] > 1 | threshold[ 2 ] < 0){
    stop("Bad input: threshold should be between 0 and 1")
  }
  
  if( !is.logical(plotting) ){
    stop("Bad input: plotting should be a boolean")
  }
  
  if( !is.logical(G) ){
    stop("Bad input: G should be a boolean")
  }
  
  len <- dim(dmlm_obj[[1]])[2]
  zeta_means <- apply( dmlm_obj[[ 2 ]][ , , ( burnin + 1 ):len ],c( 1, 2), mean )
  xi_means <- apply( dmlm_obj[[ 5 ]][ , ( burnin + 1 ):len ], 1,mean )
  
  selected_zeta <-  which( zeta_means >= threshold[ 1 ], arr.ind = T )
  selected_xi <-  which( xi_means >= threshold[ 2 ] )
  
  if( G ){
    estimated_G = apply( dmlm_obj[[ 7 ]][ , , ( burnin + 1 ):len ],c( 1, 2), mean ) 
  }
  
  # Plots of number of selected indices and PPI
  
  if( plotting == TRUE ){
    plot( 1:len, apply( dmlm_obj[[ 2 ]], 3, sum ), xlab = "Samples", ylab = "Number of selected covariates", lty = 1, type = "l") 
    plot( 1:len, apply( dmlm_obj[[ 5 ]], 2, sum ), xlab = "Samples", ylab = "Number of selected parts", lty = 1, type = "l") 
    
    # Plot zeta PPI
    y <- c( zeta_means )
    x <- seq(1,length(y))
    data <- data.frame(cbind(y,x))
    print(
      ggplot(data, aes(x, y) ) +
        geom_segment(aes(xend = x, yend = 0), size = 1 , lineend = "butt") + 
        labs(x="Covariate Index", 
             y="PPI") + geom_abline(slope = 0, intercept = threshold[ 1 ], linetype = "dashed"))
    
    # Plot xi PPI
    y <- c( xi_means )
    x <- seq(1,length(y))
    data <- data.frame(cbind(y,x))
    
    print(
      ggplot(data, aes(x, y) ) +
        geom_segment(aes(xend = x, yend = 0), size = 1 , lineend = "butt") + 
        labs(x="Balance Index", 
             y="PPI") + geom_abline(slope = 0, intercept = threshold[ 2 ], linetype = "dashed"))
  }
  if( G ){
    output <- list( selected_zeta = selected_zeta, selected_xi = selected_xi, mppi_zeta = zeta_means, mppi_xi = xi_means, estimated_G = estimated_G )
  }else{
    output <- list( selected_zeta = selected_zeta, selected_xi = selected_xi, mppi_zeta = zeta_means, mppi_xi = xi_means )
  }
  return( output )
}

