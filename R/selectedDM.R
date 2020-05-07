
selected_DM <- function( dm_obj = NULL, threshold = c( 0.5 ), plotting = FALSE, burnin = 0, cov_lab = NULL, edge_lab = NULL ){
  # dtm - output from DTMbvs_R
  # threshold - double for posterior probability of inclusion threshold
  # plotting - logical (True or False) indicator to plot # selected covarariates as well as MPPIs
  # burnin - number of MCMC samples to drop before inference, default = 0
  # cov_lab - P-dimensional vector of covariate names
  # edge_lab - Branch-dimensional vector of branch names
  
  library(ggplot2)
  
  if( burnin%%1 != 0){
    stop("Bad input: burn-in should be an integer")
  }
  
  if( burnin < 0 ){
    stop("Bad input: burn-in should be positive")
  }
  
  if( length( threshold ) > 1 ){
    stop("Bad input: threshold must be a single number")
  }
  
  if( threshold > 1 | threshold < 0 ){
    stop("Bad input: threshold should be between 0 and 1")
  }
  
  if( !is.logical( plotting ) ){
    stop("Bad input: plotting should be a boolean")
  }
  
  if( !is.null(edge_lab) & (length( edge_lab ) != dim( dm_obj[[ 1 ]] )[ 1 ]) ){
    stop("Bad input: length of edge labels does not match number of branches")
  }
  
  if( !is.null(cov_lab) & (length( cov_lab ) != dim( dm_obj[[ 2 ]] )[ 1 ]/dim( dm_obj[[ 1 ]] )[ 1 ]) ){
    stop("Bad input: length of covariate labels does not match number of covariates")
  }
  
  len <- dim( dm_obj[[ 1 ]] )[ 2 ]
  B <- dim( dm_obj[[ 1 ]] )[ 1 ]
  P <- dim( dm_obj[[ 2 ]] )[ 1 ]/dim( dm_obj[[ 1 ]] )[ 1 ]
  
  mppi <- rowMeans((dm_obj[[2]][, burnin:len] != 0) + 0)
  selected_zeta <- which(t(matrix(mppi, nrow = P, ncol = B)) > threshold)
  
  # Plots of number of selected indices and PPI
  
  if( plotting == TRUE ){
    plot( 1:len, apply( dm_obj[[ 2 ]] != 0 , 2, sum ), xlab = "Samples", ylab = "Number of selected covariates", lty = 1, type = "l")
    
    # Plot zeta PPI
    y <- c( mppi )
    x <- seq(1,length(y))
    data <- data.frame(cbind(y,x))
    print(
      ggplot(data, aes(x, y) ) +
        geom_segment(aes(xend = x, yend = 0), size = 1 , lineend = "butt") +
        labs(x="Branch/Covariate Index",
             y="PPI") + geom_abline(slope = 0, intercept = threshold[ 1 ], linetype = "dashed"))
    
  }
  
  covariate_list <-  ( selected_zeta - 1 ) %% B + 1
  branch_list <- floor( ( selected_zeta - 1 )/B ) + 1
  
  edge_lab. <- if( is.null( edge_lab ) ){ seq( 1:B ) }else{ edge_lab }
  cov_lab. <- if( is.null( cov_lab ) ){ seq( 1:P ) }else{ cov_lab }
  
  selected_name <- cbind( edge_lab.[ branch_list ], cov_lab.[ covariate_list ] )
  colnames(selected_name) <- c( "Branch", "Covariate" )
  
  return( list( selected = selected_zeta, mppi = mppi, selected_name = selected_name) )
}

