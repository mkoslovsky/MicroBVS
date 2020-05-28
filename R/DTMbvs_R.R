# Wrapper function for the Rcpp code to initiate with defaults and simulate data if requested for DTM model 
DTMbvs_R <- function( iterations = 50000, thin = 10, tree = NULL, Y = NULL, X = NULL, prior = "BB", alpha = NULL, phi = NULL, zeta = NULL, sigma2_alpha = 10, sigma2_phi = 10,
                      aa = 1, bb = 1, a_G = log(0.1/0.9), b_G = 0.5, Omega = NULL, G = NULL, v0 = 0.1, v1 = 10, pie = NULL, lambda = 1, subject_sim = 100, num_leaf = 5,
                      covariates_sim = 50, corr = 0.2, num_branch = 3, num_cov = 5, seed = 1212, eval = TRUE, warmstart = FALSE ){
  
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
  if( iterations%%1 != 0 | iterations <= 0){
    stop("Bad input: iterations should be a positive integer")
  }
  
  if( thin%%1 != 0 | thin < 0 ){
    stop("Bad input: thin should be a positive integer")
  }
  
  if( prior == "MRF_fixed" & is.null( G ) ){
    stop("Bad input: MRF prior requires a known graphical structure")
  }
  
  if( (prior != "BB") + (prior != "MRF_fixed") + (prior != "MRF_unknown") == 3 ){
    stop("Bad input: prior is not in correct format")
  }
  
  
  # Set seed for replication
  set.seed( seed )
  
  # Simulate data if there is no tree, X, or Y
  if( is.null( tree ) | is.null( X ) | is.null( Y ) ){
    
    if( eval == TRUE ){
      input <- readline(prompt = "Tree, count, or covariate matrix is missing. Results will be based on simulated data. Proceed (Y/N)?")
      if( (input != "Y") + (input != "N") == 2 ){
        stop( "Please refer to vignette for help with this package." )
      }
      if( input == "N" ){
        stop( "Please refer to vignette for help with this package." )
      }
    }
    
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
    
    # Simulate random DTM with phylogenetic tree
    # Relies on 'ape' package
    tree.ex <- rtree( n = num_leaf )
    
    # Get dimensions from tree
    # Set number of parent nodes = #subtrees = #parentheses sets
    V <- tree.ex$Nnode
    
    # Set number of child nodes for each parent node
    Cv <- table( tree.ex$edge[,1] )
    
    # Set number of leaves (tips) in the tree
    K <- length( tree.ex$tip.label )
    
    # Set parameters
    B_sim <- sum( Cv )
    
    # Set correlation structure between covariates
    sig <- diag( covariates_sim )
    for( i in 1:covariates_sim ){
      for( j in 1:covariates_sim ){
        if( i != j){
          sig[ i , j] = corr^abs(i - j)
        }
      }
    }
    
    # Simulate covariates for selection
    X <- rmvnorm( subject_sim, rep( 0, nrow(sig) ), sig )
    
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
    phi_sim[ true_cov ] <- runif(sum(zeta_sim), 0.9, 1.2)*sample(c(-1,1), sum(zeta_sim), replace = TRUE)
    
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
  }
  
  # Adjust inputs if tree and X are provided
  if( is.null(tree) ){
    tree <- tree.ex
  }
  
  # Get dimensions from tree
  # Set number of parent nodes = #subtrees = #parentheses sets
  V <- tree$Nnode
  branches <-  sum( table( tree$edge[,1] ) )
  covariates <- ifelse( is.null( X ), covariates_sim, ncol( X ) )
  subject <- ifelse( is.null( Y ), subject_sim, nrow( Y ) )
  
  # Initiate starting values and allocate memory
  samples <- floor( iterations/thin )
  
  # Intercept term alpha_j
  alpha. <- matrix( 0, nrow = branches, ncol = samples )
  
  # Inclusion indicators zeta_jp
  zeta. <- array( 0, dim = c( branches, covariates, samples ) )
  
  # Regression Coefficients phi_jp
  phi. <- array( 0, dim = c( branches, covariates, samples ) )
  
  # Graph parameters
  Omega. <- array( 0, dim = c( covariates, covariates, samples ) )
  Var. <- array( 0, dim = c( covariates, covariates, samples ) )
  G. <- array( 0, dim = c( covariates, covariates, samples ) )
  
  # Adjust inital values for alpha, zeta, and phi if they are still NULL
  alpha.[ , 1] <- if( is.null( alpha ) ){ rnorm( n = branches ) }else{ alpha }
  zeta.[ , , 1] <- ifelse( is.null( zeta ), 0, zeta )
  phi.[ , , 1]  <-  ifelse( is.null( phi ), 0, phi )
  
  # Adjust inital values for graphical parameters if they are still NULL
  pie <- ifelse( is.null( pie ), 2/( covariates - 1), pie )
  
  set_V <- function( G. = G, v0. = v0, v1. = v1 ){
    G_in <- G.
    G_out <- 1 - G.
    
    V <- v1.*G_in + v0.*G_out
    diag( V ) <- 0
    return(V)
  }
  
  Omega.[ , , 1 ] <- if( is.null( Omega ) ){ diag( covariates ) }else{ Omega }
  G.[ , , 1 ] <- if( is.null( G ) ){ 0 }else{ G }
  
  Var.[ , , 1] <- set_V( G.[ , , 1 ], v0, v1 )
  S. <- t(X)%*%X
  
  # Prep tree information arguments
  # Function :: Get list of lineages for each leaf node and parent/children relationships
  get_lineage <- function( K. = K, tree. = tree ){
    
    # Note that it is standard for 'ape'package to label leaf nodes 1 to K
    # K is number of leave nodes in the tree
    # tree is a treeObject from 'ape' package
    
    # Initiate a list 'lineage' to collect all branches along the chain from
    # leaf i to the root node.
    lineage <- list()
    
    # Parent/child Matrix
    parent_children <- matrix( 0, nrow = max( tree.$edge ), ncol =  max( tree.$edge ) )
    
    for( i in 1:K. ){
      
      # Initate a stopping rule
      stopping <- 1
      
      # Initiate search term as the node we are looking for
      term <- i
      
      # Initiate place holder for current list i and set iteration
      line <- matrix(NA, nrow = 20, ncol = 2)
      l <- 1
      
      # Stop seaching each chain when we get to the root node
      # (i.e., there is no more parent nodes)
      while( stopping != 0 ){
        
        # Identify parent node and add branch to matrix, update parent/child matrix,
        # and increment index
        
        line[ l, ] <- tree.$edge[ which( tree.$edge[,2] == term ) , ]
        parent_children[ line[l,1], line[l,2] ] <- 1
        l = l + 1
        
        # Set up next term
        term <- tree.$edge[ which( tree.$edge[,2] == term ) , 1 ]
        
        # Update stopping rule
        stopping <- length( which( tree.$edge[,2] == term ) )
      }
      
      # Remove NA from 'line' matrix
      line <- line[ complete.cases( line ), ]
      
      # Add line to lineage
      lineage <- append( lineage, list( line ) )
    }
    # Return list of leaf lineage and parent matrix
    return( list( lineage = lineage, parent_children = parent_children ) )
  }
  
  # Note that for 'lineage' the first column of the kth item in the list are the parent nodes of leaf k
  # Each row in the kth item in the list represents a branch
  # For 'parent_children' row is parent and column is children
  # Siblings are derived by reading across each row
  leaf_lineage <- get_lineage( K. = length( tree$tip.label ), tree. = tree )
  
  # Function :: Get leaves that are in each branch, sum of counts for each branch (by subject),
  #             and sum of counts for each subtree (by subject)
  counter <- function( tree. = tree, subjects. = subjects, branches. = branches, V. = V , leaf_lineage. = leaf_lineage, Y = Y  ){
    # Initiate subtree and node children
    subtree <- list()
    node_children_pointer <- list()
    
    # Interate over each subtree (Note these are the unique items in 1st column of edges)
    for( v in unique( tree.$edge[ ,1 ] ) ){
      
      # Get pointer to child branches and the children of v
      children_pointer <- which( tree.$edge[ , 1 ] == v )
      node_children_pointer <- append( node_children_pointer, list( children_pointer ) )
      children <- tree.$edge[ children_pointer , 2 ]
      
      # Find leaves that connect to each of its branches
      for(c in 1:length( children ) ){
        leaves <- which( sapply( leaf_lineage.$lineage, function(x){ children[c] %in% x } ) )
        # Append to list of subtree leaves
        subtree[[ children_pointer[ c ] ]] <- leaves
      }
    }
    
    # Get number of counts for each subject in each branch (n x B)
    # Initiate matrix to hold summed counts
    branch_counts <- matrix( 0, nrow = subjects. , ncol = branches. )
    
    # Sum each row of 'Y' given a subset of columns
    for(b in 1:branches. ){
      # Don't use rowSums unless there is more than one column
      if( length( subtree[[b]] ) > 1  ){
        branch_counts[ , b ] <- rowSums( Y[, subtree[[b]] ] )  #
      }else{
        branch_counts[ , b ] <- Y[, subtree[[b]] ]
      }
    }
    
    # Sum the branches that are in the same subtree 'v'
    # Initiate matrix to hold summed counts for the subtree 'v'
    subtree_counts <- matrix( 0, nrow = subjects. , ncol = V. )
    
    # Sum the rows of the columns that are in the same subtree
    # Note these are in order of unique(tree.ex$edge[,1])
    for(v in 1:V.){
      
      if( length( node_children_pointer[[ v ]] ) > 1  ){
        subtree_counts[ , v ] <- rowSums( branch_counts[, node_children_pointer[[ v ]] ] )
      }else{
        subtree_counts[ , v ] <- branch_counts[ , node_children_pointer[[ v ]] ]
      }
    }
    # Return output
    return( list( node_children_pointer = node_children_pointer, branch_counts = branch_counts, subtree_counts = subtree_counts) )
    
  } # End of Counter Function
  
  Y_summary <- counter( subjects. = nrow(Y), tree. = tree, branches. = branches, V. = V, leaf_lineage. = leaf_lineage, Y = Y )
  
  # Get list of each branch's node
  branch_location = numeric()
  for(b in 1:branches){
    branch_location[ b ] <- which(  sapply(Y_summary$node_children_pointer, function(x){ b %in% x }) )
  }
  
  # Save data to environment
  X.model <<- X
  Y.model <<- Y
  tree.model <<- tree
  
  if( warmstart == TRUE ){
    # Warmstart
    cormat = matrix(0, branches, covariates)
    pmat = matrix(0, branches, covariates)
    
    ## Compositionalize by subtree 
    
    # Need subtree location for each branch 
    yy <- numeric()
    for( i in 1:length( branch_location ) ){
      yy <- cbind( yy, 0 + Y_summary$branch_counts[ , i]/Y_summary$subtree_counts[ , branch_location[ i ] ] )
    }
    yy[ is.nan( yy )] <- 0 
    
    for(rr in 1:branches){
      for(cc in 1:covariates){
        pmat[rr, cc] = stats::cor.test(X.model[, cc], yy[, rr], method = "spearman", exact = F)$p.value
        cormat[rr, cc] = stats::cor(X.model[, cc], yy[, rr], method = "spearman")
      }
    }
    
    # defaults to 0.2 false discovery rate
    pm = matrix((stats::p.adjust(c(pmat), method = "fdr") <= 0.2) + 0, branches, covariates )
    betmat = cormat * pm
    
    # alpha warmstart
    alpha.[ , 1] <- scale(log(colSums(Y_summary$branch_counts)))
    zeta.[ , , 1] <- pm
    phi.[ , , 1]  <- betmat
  }
  
  # Run DTMbvs
  if( eval == TRUE ){
    # Run model
    output <- DTMbvs( iterations, thin, prior, X, branch_location, Y_summary$node_children_pointer, Y_summary$branch_counts, Y_summary$subtree_counts, alpha., phi., zeta., sigma2_alpha, sigma2_phi, aa, bb, Omega., G., Var., S., v0, v1, a_G, b_G, pie, lambda )
    names( output ) <- c("alpha", "zeta", "psi", "Omega", "G")
    
    # Return based on prior
    if( prior == "BB"){
      output <- output[1:3]
    }
    return( output )
  }
}

