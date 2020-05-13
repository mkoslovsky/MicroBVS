# Processes output from DMbvs_R for inference
# Returns selected zeta at given threshold of inclusion
# Option to produce marginal PPI plots and numbers of selected zeta versus MCMC iteration
selected_DM <- function( dm_obj = NULL, threshold = c(0.5), plotting = FALSE, burnin = 0, G = FALSE ){
  # dm_obj - output from DMbvs or DMbvs_R
  # threshold - double posterior probability of inclusion threshold for zeta 
  # plotting - logical (True or False) indicator to plot # selected zeta as well as PPI for zeta 
  # burnin - number of MCMC samples to drop before inference, Default = 0
  # G - return learned graph structure at threshold
  
  library(ggplot2)
  
  if( burnin%%1 != 0){
    stop("Bad input: burn-in should be an integer")
  }
  
  if( burnin < 0 ){
    stop("Bad input: burn-in should be positive")
  }
  
  
  if( threshold > 1 | threshold < 0 ){
    stop("Bad input: threshold should be between 0 and 1")
  }
  
  if( !is.logical(plotting) ){
    stop("Bad input: plotting should be a boolean")
  }
  
  if( !is.logical(G) ){
    stop("Bad input: G should be a boolean")
  }
  
  len <- dim(dm_obj[[1]])[2]
  zeta_means <- apply( dm_obj[[ 2 ]][ , , ( burnin + 1 ):len ],c( 1, 2), mean ) 
  selected_zeta <-  which( zeta_means >= threshold, arr.ind = T ) 
  
  if( G ){
    estimated_G = apply( dm_obj[[ 5 ]][ , , ( burnin + 1 ):len ],c( 1, 2), mean ) 
  }
  
  # Plots of number of selected indices and PPI
  
  if( plotting == TRUE ){
    plot( 1:len, apply( dmlm_obj[[ 2 ]], 3, sum ), xlab = "Samples", ylab = "Number of selected covariates", lty = 1, type = "l") 
 
    # Plot zeta PPI
    y <- c( zeta_means )
    x <- seq(1,length(y))
    data <- data.frame(cbind(y,x))
    print(
      ggplot(data, aes(x, y) ) +
        geom_segment(aes(xend = x, yend = 0), size = 1 , lineend = "butt") + 
        labs(x="Covariate Index", 
             y="PPI") + geom_abline(slope = 0, intercept = threshold, linetype = "dashed"))

  }
  if( G ){
    output <- list( selected_zeta = selected_zeta,mppi_zeta = zeta_means, estimated_G = estimated_G )
  }else{
    output <- list( selected_zeta = selected_zeta,mppi_zeta = zeta_means )
  }
  return( output )
}

