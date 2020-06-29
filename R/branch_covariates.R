# Generates a list of the covariates associated with a particular leaf node
branch_covariates <- function( tree = tree, dtm_obj = dtm_obj, covariate_name = NULL, branch_name = NULL, threshold = 0.50 ){
  # arg checking
  if( is.null( branch_name ) ){
    stop("Please provide a leaf name that matches a tip label of the tree.")
  }
  
  if( is.null( covariate_name ) ){
    stop("Please provide a vector of covariate names.")
  }
  
  if( length( threshold ) > 1 ){
    stop("Bad input: threshold must be a single number")
  }
  
  if( threshold > 1 | threshold < 0 ){
    stop("Bad input: threshold should be between 0 and 1")
  }
  
  # Note that it is standard for 'ape' package to label leaf nodes 1 to K
  # tree is a treeObject from 'ape' package
  # K is number of leaf nodes in the tree 
  
  get_lineage <- function( K. = K, tree. = tree ){
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
  
  leaf_lineage <- get_lineage( K. = length( tree$tip.label ), tree. = tree )
  
  hold <- numeric()
  
  for( item in 1:length( branch_name ) ){
    len <- leaf_lineage$lineage[  which( tree$tip.label == branch_name[ item ] ) ]  
    for( i in 1:dim(len[[1]])[1] ){
      for(j in 1:dim(tree$edge)[1]){
        if( sum( leaf_lineage$lineage[ which( tree$tip.label == branch_name[ item ] ) ][[1]][ i, ] == tree$edge[ j, ]) == 2 ){
          hold <- c( hold , j )
        }
      }
    }
  }
  
  selected <- selected_DTM( dtm_obj = dtm_obj, threshold = threshold )
  
  covariates <- covariate_name[ colSums( selected$mppi_zeta[ hold, ] > threshold ) > 0 ]
  
  return( covariates )
}
 
