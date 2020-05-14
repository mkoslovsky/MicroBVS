# Wrapper function for the Rcpp code to simulate DTM data  
# DTMbvs: Dirichlet-tree Multinomial Regression Models with Bayesian Variable Selection for Microbiome Data - an R Package
simulate_DTM <- function( subject_sim = 100, tree = NULL, num_leaf = 5,
                      covariates_sim = 50, rho = NULL, Sigma = NULL, num_branch = 3, num_cov = 5, phi_min = 0.9, phi_max = 1.2, seed = 1212 ){
  
  library(mvtnorm)
  library(MCMCpack)
  library(ape)
  library(GGMselect)
  # iterations - integer number of MCMC samples, default = 50000
  # thin - integer MCMC by # thin, default = 10
  # tree - tree object .tre file to be sourced
  # Y - N x branch matrix of count data. Make sure branch order corresponds to edge list in tree
  # X - N x P matrix of covariates
  # prior - name of selection prior to use. Takes values "BB" beta-binomial (default), "MRF_fixed" known Markov random field, "MRF_unknown" unknown Markov random field
  # alpha branch x 1 vector of branch intercepts
  # phi branch x P matrix of regression coefficients
  # zeta - branch x P matrix of inclusion indicators
  # sigma2_alpha - prior variance for alpha, default = 10
  # sigma2_phi - prior variance for phi, default = 10
  # aa - double beta-binomial hyperparameter
  # bb - double beta-binomial hyperparameter
  # a_G - double MRF hyperparameter
  # b_G - double MRF hyperparameter
  # Omega - P x P concentration matrix
  # G - P x P adjacency matrix for relational graph
  # v0 - double variance of exclusion for the relational graph
  # v1- double variance of inclusion for the relational graph
  # pie - double prior probability of inclusion for the relational graph
  # lambda - double exponential prior hyperparameter
  # subject_sim - integer number of subjects to simulate, default = 100
  # num_leaf - integer number of leaves in tree to simulate, default = 5
  # covariates_sim - iteger number of covariates to simulate, default = 50
  # corr - double correlation structure, default = 0.20
  # num_branch = interger number of branches the associated covariates are in, default = 3
  # num_cov - integer number of associated covariates to simulate, default = 5
  # seed - integer seed #, default = 1234
  # eval - boolean if true, run model on simulated data. Otherwise, simply simulate
  # warmstart - boolean if true, start model using informed initial values. 
  
  # Defense
  
  # Set seed for replication
  set.seed( seed )
  
  # Simulate data if there is no tree, X, or Y
    
    if( subject_sim%%1 != 0 | subject_sim < 0 ){
      stop("Bad input: number of subjects should be a positive integer")
    }
    
    if( num_leaf%%1 != 0 | num_leaf < 0 ){
      stop("Bad input: number of leaves should be a positive integer")
    }
    
    if( covariates_sim%%1 != 0 | covariates_sim< 0 ){
      stop("Bad input: number of covariates should be a positive integer")
    }
    
    if( num_branch%%1 != 0 | num_branch < 0 ){
      stop("Bad input: number of branches for associated covariates should be a positive integer")
    }
    
    if( num_cov%%1 != 0 | num_cov < 0 ){
      stop("Bad input: number of associated covariates should be a positive integer")
    }
    
    if( num_cov >= covariates_sim ){
      stop("Bad input: number of associated covariates should be less than total number of covariates")
    }
  
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
      if( ncol(Sigma) != covariates_sim ){
        stop("Bad input: Please provide covariance matrix to match the number of covariates")
      }
    }  
  
    # covariance matrix for predictors
    if( !is.null( rho ) ){
      Sigma <- matrix( 0, covariates_sim, covariates_sim )
      Sigma = rho^abs(row(Sigma) - col(Sigma))
    }
    
    # Simulate random DTM with phylogenetic tree
    # Relies on 'ape' package
    tree.ex <- if( is.null( tree ) ){ rtree( n = num_leaf ) }else{ tree }
    
    # Get dimensions from tree
    # Set number of parent nodes = #subtrees = #parentheses sets
    V <- tree.ex$Nnode
    
    # Set number of child nodes for each parent node
    Cv <- table( tree.ex$edge[,1] )
    
    # Set number of leaves (tips) in the tree
    K <- length( tree.ex$tip.label )
    
    # Set parameters
    B_sim <- sum( Cv )
    
    # Simulate covariates for selection
    X <- scale( rmvnorm( subject_sim, rep( 0, covariates_sim ), Sigma ) ) 
    
    # Set true inclusion indicators
    zeta_sim <- matrix( 0 , B_sim, covariates_sim )
    for( i in sample( seq( 1, B_sim ), num_branch ) ){
      select <- sample( 1:covariates_sim, num_cov )
      zeta_sim[ i, select] <- 1
    }
    truth <<- which(zeta_sim == 1)
    
    # Simulate true alpha parameters and regression coefficients phi
    alpha_sim <- matrix( 1, nrow = subject_sim , ncol = 1)%*%( runif( n = B_sim, -1.3,1.3 ) )
    true_cov <- which( zeta_sim == 1)
    phi_sim <- matrix( 0 , B_sim , covariates_sim )
    phi_sim[ true_cov ] <- runif(sum(zeta_sim), phi_min, phi_max)*sample(c(-1,1), sum(zeta_sim), replace = TRUE)
    
    # Used for count probabilities
    inside_sim <- exp( alpha_sim + X%*%t(phi_sim) )
    
    # Look through the tree and seperate to get dirichlet parameters and then simulated
    
    node_counts <- matrix(0, nrow= subject_sim, ncol = ( V + K ) )
    node_counts[ , ( K + 1 ) ] <- sample( seq( 7500,10000 ) , subject_sim )
    
    for( b in ( K + 1 ):( sum( Cv ) + 1 )  ){
      node <- which( tree.ex$edge[ , 1] == b )
      
      # Split inside by each subtree
      inside_branches <- inside_sim[, node]
      
      # Simulate probabilities for each subtree
      prob_sim <- apply( inside_branches, 1,  function(x){ rdirichlet(1,x) } )
      
      # Simulate count data
      for( i in 1:subject_sim ){
        y <- t(rmultinom(1, node_counts[ i, b ], t( prob_sim )[ i, ] ))
        node_counts[ i , ( tree.ex$edge[ node, 2]  ) ] <- y
      }
    }
    
    Y <- node_counts[ , 1:K]
    X <- scale(X)
    
    return( list( Y = Y, X = X, phi_sim = phi_sim, zeta_sim = zeta_sim, tree = tree.ex, Sigma = Sigma ) )
  }


