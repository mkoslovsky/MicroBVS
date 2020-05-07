bfdr <- function( mppi_vector, threshold = 0.1 ){
  # Courtesy of D.Wadsworth
  # arg checking
  if( any( mppi_vector > 1 | mppi_vector < 0 ) ){
    stop("Bad input: mppi_vector should contain probabilities")
  }
  if( threshold > 1 | threshold < 0 ){
    stop( "Bad input: threshold should be a probability" )
  }
  # sorting the ppi's in decresing order
  sorted_mppi_vector <- sort( mppi_vector, decreasing = TRUE )
  
  # computing the fdr
  fdr <- cumsum( ( 1 - sorted_mppi_vector ) )/seq( 1:length( mppi_vector ) )
  
  # determine index of the largest fdr less than threshold
  thecut.index <- max( which( fdr < threshold ) )
  ppi_threshold <- sorted_mppi_vector[ thecut.index ]
  selected <- mppi_vector > ppi_threshold
  
  return( list( selected = selected, threshold = ppi_threshold ) )
}