// Performs Bayesian variable selection for multivariate count data that follow a tree stucture.
// Accommodates various inclusion priors
// M. Koslovsky 
// mkoslovsky12@gmail.com
// 2019

#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
#include <math.h>

// [[Rcpp::plugins("cpp11")]]

// Define helper functions
namespace help{

// Function :: Call sample function  from R
double sample_cpp( IntegerVector x ){
  // Calling sample()
  Function f( "sample" );
  IntegerVector sampled = f( x, Named( "size" ) = 1 );
  return sampled[0];
}

// Make a sequence
IntegerVector myseq( int first, int last) {
  IntegerVector y(abs(last - first) + 1);
  if (first < last) 
    std::iota(y.begin(), y.end(), first);
  else {
    std::iota(y.begin(), y.end(), last);
    std::reverse(y.begin(), y.end());
  }
  return y;
}

// Function :: Calculates log likelihood contribution of multivariate count data assuming each sub-tree
//             follows a Dirichlet-Multinomial distribution
double log_like_cpp( arma::mat x, int branch_loc, List node_children_pointer, arma::mat branch_counts, arma::mat subtree_counts, arma::mat loggamma ){
  
  // Initiate memory
  int B = loggamma.n_rows;
  int N = x.n_rows;
  arma::mat gamma( N, B );
  arma::vec expand( N );
  expand.ones();
  arma::vec gamma_subtree( N );
  gamma_subtree.zeros();
  double log_like = 0;
  
  // Calculate gamma
  gamma =  exp( loggamma );
  
  // Calculate sum of gammas within each subtree
  IntegerVector item = node_children_pointer[ branch_loc ];
  for( int n = 0; n < N; ++n ){
    for( int c = 0; c < item.size(); ++c ){
      gamma_subtree[ n ] += gamma( n, item( c ) - 1 );
    }
  }
  
  log_like = log_like + sum( lgamma( gamma_subtree ) ) ;
  log_like = log_like - sum( lgamma( subtree_counts.col( branch_loc ) + gamma_subtree ) ) ;
  
  for( int c = 0; c < item.size(); ++c ){
    log_like = log_like +  sum( lgamma( branch_counts.col( item( c ) - 1 ) + gamma.col( item( c ) - 1 ) ) ) ;
    log_like = log_like - sum( lgamma( gamma.col( item( c ) - 1 ) ) );
  }
  
  // Return output
  return log_like;
}

// Function :: Calculate log alpha contribution
double log_alpha_cpp( arma::vec alpha, double sigma2_alpha ){
  // sum the log p(alpha) across all branches
  // initiate log_alpha and set size of vector
  double log_alpha = 0;
  int B = alpha.size();
  
  // Sum over all alpha
  for( int b = 0; b < B; ++b ) {
    log_alpha += -0.50*log( 2*atan(1)*4*sigma2_alpha ) - 1/( 2*sigma2_alpha )*pow( alpha[ b ],2 );
  }
  
  // Return output
  return log_alpha;
}

// Function :: Propose a new alpha
double alpha_prop_cpp( double alpha ){
  // Simulate proposal with jitter around current alpha
  double alpha_prop = alpha + rnorm( 1 )[ 0 ];
  
  // Return output
  return alpha_prop ;
}

// Function :: Propose a new phi
double phi_prop_cpp( double phi ){
  // Simulate proposal with jitter around current phi
  double phi_prop = phi + rnorm( 1 )[ 0 ];
  
  // Return output
  return phi_prop ;
}

// Function :: Calculate log phi contribution
double log_phi_cpp( NumericMatrix phi, double sigma2_phi, NumericMatrix zeta ){
  // sum the log p(phi) across all branches and p
  
  // initiate log_phi and set size of vector
  double log_phi = 0;
  int B = phi.rows();
  int P = phi.cols();
  
  // Sum over all phi in B and P
  for( int b = 0; b < B; ++b ) {
    for( int p = 0; p < P; ++p ) {
      if( zeta( b, p ) == 1){
        log_phi += -0.50*log( 2*atan(1)*4*sigma2_phi ) - 1/( 2*sigma2_phi )*pow( phi( b, p ), 2 );
      }
    }
  }
  
  // Return output
  return log_phi;
}

// Function :: Calculate individual zeta
double log_zeta_pj_cpp( double t_pj, double aa, double bb ){
  
  double post_a = t_pj + aa;
  double post_b = 1 - t_pj + bb;
  double log_zeta_pj = lgamma( post_a ) + lgamma( post_b ) - lgamma( post_a + post_b ) - ( lgamma( aa ) + lgamma( bb ) - lgamma( aa + bb ) );
  
  // Return output
  return log_zeta_pj ;
}

// Function :: Calculate zeta probability given graph G
double log_zeta_G_cpp( double a_G, arma::mat zeta, double b_G, arma::mat G, int branch ){
  
  // Initiate memory
  int P = zeta.n_cols;
  arma::vec one( P );
  one.ones();
  
  arma::mat zeta_row = zeta.row( branch );
  
  arma::vec log_zeta_G = a_G*one.t()*zeta_row.t() + b_G*zeta_row*G*zeta_row.t();
  
  double log_zeta = log_zeta_G[ 0 ];
  
  // Return output
  return log_zeta ;
}

// Function :: Update alpha
List alpha_update_cpp(
    arma::mat x,
    arma::vec alpha,
    arma::mat phi,
    arma::vec branch_location,
    List node_children_pointer,
    arma::mat branch_counts,
    arma::mat subtree_counts,
    double sigma2_alpha,
    arma::mat loggamma
){
  // Initiate memory space
  int B = phi.n_rows;
  int N = x.n_rows;
  
  // Update each alpha
  for( int b = 0; b < B; ++b ){
    
    // Proposal alpha
    arma::vec alpha_proposal = alpha;
    double current_alpha = alpha[ b ];
    alpha_proposal[ b ] = help::alpha_prop_cpp( current_alpha );
    
    arma::mat loggamma_proposal = loggamma;
    for( int n = 0; n < N; ++n ){
      loggamma_proposal( n , b ) = loggamma_proposal( n , b ) + alpha_proposal[ b ] - current_alpha;
    }
    
    int branch_loc = branch_location[ b ] - 1 ;
    
    // Calculate ratio
    double r = help::log_like_cpp( x, branch_loc, node_children_pointer,  branch_counts, subtree_counts, loggamma_proposal ) + help::log_alpha_cpp( alpha_proposal, sigma2_alpha ) - ( help::log_like_cpp( x, branch_loc, node_children_pointer, branch_counts, subtree_counts, loggamma ) + help::log_alpha_cpp( alpha, sigma2_alpha ) );
    
    // Calculate acceptance probability
    double a  = log( runif( 1 )[ 0 ] );
    
    // Determine acceptance
    if(a < r){
      alpha[ b ] = alpha_proposal[ b ];
      loggamma = loggamma_proposal;
    }
  }
  // Return output
  List alpha_return( 2 );
  alpha_return[ 0 ] = alpha;
  alpha_return[ 1 ] = loggamma;
  return alpha_return;
}

// Function :: Between Step (jointly update phi and zeta) with add/delete or swap with MRF
List between_phi_zeta_swap_update_MRF_cpp(
    arma::mat x,
    arma::mat zeta,
    arma::vec alpha,
    arma::mat phi,
    arma::vec branch_location,
    List node_children_pointer,
    arma::mat branch_counts,
    arma::mat subtree_counts,
    double sigma2_phi,
    double a_G,
    double b_G,
    arma::mat G,
    arma::mat loggamma
){
  
  // Between Model Step
  // Initiate memory space
  int B = zeta.n_rows;
  int P = zeta.n_cols;
  int dim = B*P;
  int N = x.n_rows;
  
  // Set proposal to current
  arma::mat zeta_proposal( B, P );
  arma::mat phi_proposal( B, P );
  zeta_proposal = zeta;
  phi_proposal = phi;
  
  // Add/Delete or swap
  int swap = rbinom( 1, 1, 0.5 )[ 0 ];
  IntegerVector add_vcp;
  IntegerVector delete_vcp;
  
  // Swap
  if( swap == 1 ){
    
    // Get which included and excluded
    IntegerVector included( 0 );
    IntegerVector excluded( 0 );
    for( int bp = 0; bp < dim; ++bp ) {
      if( zeta[ bp ] == 1 ){
        included.push_back( bp );
      }else{
        excluded.push_back( bp );
      }
    }
    
    // Adjust missing indicator
    if( ( included.size() == 0 ) | ( excluded.size() == 0 ) ){
      swap = 0;
    }
    
    // If there is at least one included and excluded - perform a swap
    if( swap == 1 ){
      
      // Select covariate to include
      if( excluded.size() > 1 ){
        add_vcp =  sample_cpp( excluded );
      }else{
        add_vcp = excluded;
      }
      
      // Get a covariate to delete
      if( included.size() > 1 ){
        delete_vcp =  sample_cpp( included );
      }else{
        delete_vcp = included;
      }
      
      // Update proposal
      // For added
      zeta_proposal[ add_vcp[ 0 ] ] = 1;
      double current_phi = phi[ add_vcp[ 0 ] ];
      phi_proposal[ add_vcp[ 0 ] ] = help::phi_prop_cpp( current_phi );
      
      // For deleted
      zeta_proposal[ delete_vcp[ 0 ] ] = 0;
      phi_proposal[ delete_vcp[ 0 ] ] = 0;
      
      // Get branch location of added and deleted to simplify likelihood calculation.
      int add_mod = add_vcp[ 0 ] % B;
      int del_mod = delete_vcp[ 0 ] % B;
      int add_mod_P = floor( add_vcp[ 0 ]/B );
      int del_mod_P = floor( delete_vcp[ 0 ]/B );
      int branch_loc_add = branch_location[ add_mod ] - 1 ;
      int branch_loc_del = branch_location[ del_mod ] - 1 ;
      
      // Adjust loggamma for proposal
      arma::mat loggamma_proposal = loggamma;
      
      for( int n = 0; n < N; ++n ){
        loggamma_proposal( n , add_mod ) = loggamma( n , add_mod ) +  phi_proposal( add_mod, add_mod_P )*x( n, add_mod_P ) - phi( add_mod, add_mod_P )*x( n, add_mod_P ) ;
        loggamma_proposal( n , del_mod ) = loggamma_proposal( n , del_mod ) +  phi_proposal( del_mod, del_mod_P )*x( n, del_mod_P ) - phi( del_mod, del_mod_P )*x( n, del_mod_P ) ;
      }
      
      double r = 0;
      
      // Calculate ratio
      if( add_mod != del_mod){
        r = help::log_like_cpp( x, branch_loc_add, node_children_pointer, branch_counts, subtree_counts, loggamma_proposal ) + help::log_like_cpp( x, branch_loc_del, node_children_pointer, branch_counts, subtree_counts, loggamma_proposal ) + help::log_phi_cpp( wrap( phi_proposal ), sigma2_phi, wrap( zeta_proposal ) ) + help::log_zeta_G_cpp( a_G, zeta_proposal, b_G, G, add_mod ) + help::log_zeta_G_cpp( a_G, zeta_proposal, b_G, G, del_mod ) - ( help::log_like_cpp( x, branch_loc_add, node_children_pointer, branch_counts, subtree_counts, loggamma ) + help::log_like_cpp( x, branch_loc_del, node_children_pointer, branch_counts, subtree_counts, loggamma ) + help::log_phi_cpp( wrap( phi ), sigma2_phi, wrap( zeta ) ) + help::log_zeta_G_cpp( a_G, zeta, b_G, G, add_mod ) + help::log_zeta_G_cpp( a_G, zeta, b_G, G, del_mod ) );
      }else{
        r = help::log_like_cpp( x, branch_loc_add, node_children_pointer, branch_counts, subtree_counts, loggamma_proposal ) + help::log_phi_cpp( wrap( phi_proposal ), sigma2_phi, wrap( zeta_proposal ) ) + help::log_zeta_G_cpp( a_G, zeta_proposal, b_G, G, add_mod ) - ( help::log_like_cpp( x, branch_loc_add, node_children_pointer, branch_counts, subtree_counts, loggamma ) + help::log_phi_cpp( wrap( phi ), sigma2_phi, wrap( zeta ) ) + help::log_zeta_G_cpp( a_G, zeta, b_G, G, add_mod ) );
      }
      
      
      // Calculate acceptance probability
      double a  = log( runif( 1 )[ 0 ] );
      
      // Determine acceptance
      if( a < r ){
        zeta[ add_vcp[ 0 ] ] = zeta_proposal[ add_vcp[ 0 ] ];
        phi[ add_vcp[ 0 ] ] = phi_proposal[ add_vcp[ 0 ] ];
        zeta[ delete_vcp[ 0 ] ] = zeta_proposal[ delete_vcp[ 0 ] ];
        phi[ delete_vcp[ 0 ] ] = phi_proposal[ delete_vcp[ 0 ] ];
        
        loggamma = loggamma_proposal;
      }
      
    } // If not missing (still swap)
    
  } // If swap
  
  if( swap == 0 ){
    
    // Choose a branch 
    IntegerVector branches = seq( 0, B-1 );
    int branch_j = sample_cpp( branches );
    
    
    // Choose a covariate in that branch
    IntegerVector cov = seq( 0, P-1 );
    int p = sample_cpp( cov );
    
    // Set proposal to current
    arma::mat zeta_proposal( B, P );
    arma::mat phi_proposal( B, P );
    zeta_proposal = zeta;
    phi_proposal = phi;
    
    // If included propose delete
    if( zeta( branch_j, p ) == 1 ){
      
      // Update proposal
      zeta_proposal( branch_j, p ) = 0;
      phi_proposal( branch_j, p ) = 0;
      
      // Get branch location to simplify likelihood calculation.
      int branch_loc = branch_location[ branch_j ] - 1 ;
      
      // Adjust loggamma for proposal
      arma::mat loggamma_proposal = loggamma;
      
      for( int n = 0; n < N; ++n ){
        loggamma_proposal( n , branch_j ) = loggamma( n , branch_j ) +   phi_proposal( branch_j, p )*x( n, p )  -   phi( branch_j, p )*x( n, p ) ;
      }
      
      // Calculate ratio
      double r = help::log_like_cpp( x, branch_loc, node_children_pointer, branch_counts, subtree_counts, loggamma_proposal ) + help::log_phi_cpp( wrap( phi_proposal ), sigma2_phi, wrap( zeta_proposal ) ) + help::log_zeta_G_cpp( a_G, zeta_proposal, b_G, G , branch_j ) - ( help::log_like_cpp( x, branch_loc, node_children_pointer, branch_counts, subtree_counts, loggamma ) + help::log_phi_cpp( wrap( phi ), sigma2_phi, wrap( zeta ) ) + help::log_zeta_G_cpp( a_G, zeta, b_G, G , branch_j ) );
      
      // Calculate acceptance probability
      double a  = log( runif( 1 )[ 0 ] );
      
      // Determine acceptance
      if( a < r ){
        zeta( branch_j, p ) = zeta_proposal( branch_j, p );
        phi( branch_j, p ) = phi_proposal( branch_j, p );
        
        loggamma = loggamma_proposal;
      }
      
    }else{
      // If excluded propose inclusion
      
      // Update proposal
      zeta_proposal( branch_j, p ) = 1;
      double current_phi = phi( branch_j, p );
      phi_proposal( branch_j, p ) = help::phi_prop_cpp( current_phi );
      
      // Get branch location of added to simplify likelihood calculation.
      int branch_loc = branch_location[ branch_j ] - 1 ;   
      
      // Adjust loggamma for proposal
      arma::mat loggamma_proposal = loggamma;
      
      for( int n = 0; n < N; ++n ){
        loggamma_proposal( n , branch_j ) = loggamma( n , branch_j ) +  phi_proposal( branch_j, p )*x( n, p )  - phi( branch_j, p )*x( n, p );
      }
      
      // Calculate ratio
      double r = help::log_like_cpp( x, branch_loc, node_children_pointer, branch_counts, subtree_counts, loggamma_proposal ) + help::log_phi_cpp( wrap( phi_proposal ), sigma2_phi, wrap( zeta_proposal ) ) + help::log_zeta_G_cpp( a_G, zeta_proposal, b_G, G , branch_j ) - ( help::log_like_cpp( x, branch_loc, node_children_pointer, branch_counts, subtree_counts, loggamma ) + help::log_phi_cpp( wrap( phi ), sigma2_phi, wrap( zeta ) ) + help::log_zeta_G_cpp( a_G, zeta, b_G, G , branch_j ) );
      
      // Calculate acceptance probability
      double a  = log( runif( 1 )[ 0 ] );
      
      // Determine acceptance
      if( a < r ){
        zeta( branch_j, p ) = zeta_proposal( branch_j, p );
        phi( branch_j, p ) = phi_proposal( branch_j, p );
        
        loggamma = loggamma_proposal;
      }
    }
  }  // Add/delete
  
  // Return output
  List between( 3 );
  between[ 0 ] = zeta;
  between[ 1 ] = phi;
  between[ 2 ] = loggamma;
  return between;
}

// Function :: Between Step (jointly update phi and zeta) with add/delete or swap with BB
List between_phi_zeta_swap_update_BB_cpp(
    arma::mat x,
    arma::mat zeta,
    arma::vec alpha,
    arma::mat phi,
    arma::vec branch_location,
    List node_children_pointer,
    arma::mat branch_counts,
    arma::mat subtree_counts,
    double aa,
    double bb,
    double sigma2_phi,
    arma::mat loggamma
){
  
  // Between Model Step
  // Initiate memory space
  int B = zeta.n_rows;
  int P = zeta.n_cols;
  int dim = B*P;
  int N = x.n_rows;
  
  // Set proposal to current
  arma::mat zeta_proposal( B, P );
  arma::mat phi_proposal( B, P );
  zeta_proposal = zeta;
  phi_proposal = phi;
  
  // Add/Delete or swap
  int swap = rbinom( 1, 1, 0.5 )[ 0 ];
  IntegerVector add_vcp;
  IntegerVector delete_vcp;
  
  // Swap
  if( swap == 1 ){
    
    // Get which included and excluded
    IntegerVector included( 0 );
    IntegerVector excluded( 0 );
    for( int bp = 0; bp < dim; ++bp ) {
      if( zeta[ bp ] == 1 ){
        included.push_back( bp );
      }else{
        excluded.push_back( bp );
      }
    }
    
    // Adjust missing indicator
    if( ( included.size() == 0 ) | ( excluded.size() == 0 ) ){
      swap = 0;
    }
    
    // If there is at least one included and excluded - perform a swap
    if( swap == 1 ){
      
      // Select covariate to include
      if( excluded.size() > 1 ){
        add_vcp =  sample_cpp( excluded );
      }else{
        add_vcp = excluded;
      }
      
      // Get a covariate to delete
      if( included.size() > 1 ){
        delete_vcp =  sample_cpp( included );
      }else{
        delete_vcp = included;
      }
      
      // Update proposal
      // For added
      zeta_proposal[ add_vcp[ 0 ] ] = 1;
      double current_phi = phi[ add_vcp[ 0 ] ];
      phi_proposal[ add_vcp[ 0 ] ] = help::phi_prop_cpp( current_phi );
      
      // For deleted
      zeta_proposal[ delete_vcp[ 0 ] ] = 0;
      phi_proposal[ delete_vcp[ 0 ] ] = 0;
      
      // Get branch location of added and deleted to simplify likelihood calculation.
      int add_mod = add_vcp[ 0 ] % B;
      int del_mod = delete_vcp[ 0 ] % B;
      int add_mod_P = floor( add_vcp[ 0 ]/B );
      int del_mod_P = floor( delete_vcp[ 0 ]/B );
      int branch_loc_add = branch_location[ add_mod ] - 1 ;
      int branch_loc_del = branch_location[ del_mod ] - 1 ;
      
      // Adjust loggamma for proposal
      arma::mat loggamma_proposal = loggamma;
      
      for( int n = 0; n < N; ++n ){
        loggamma_proposal( n , add_mod ) = loggamma( n , add_mod ) +  phi_proposal( add_mod, add_mod_P )*x( n, add_mod_P ) - phi( add_mod, add_mod_P )*x( n, add_mod_P ) ;
        loggamma_proposal( n , del_mod ) = loggamma_proposal( n , del_mod ) +  phi_proposal( del_mod, del_mod_P )*x( n, del_mod_P ) - phi( del_mod, del_mod_P )*x( n, del_mod_P ) ;
      }
      
      double r = 0;
      
      // Calculate ratio
      if( add_mod != del_mod){
        r = help::log_like_cpp( x, branch_loc_add, node_children_pointer, branch_counts, subtree_counts, loggamma_proposal ) + help::log_like_cpp( x, branch_loc_del, node_children_pointer, branch_counts, subtree_counts, loggamma_proposal ) + help::log_phi_cpp( wrap( phi_proposal ), sigma2_phi, wrap( zeta_proposal ) ) - ( help::log_like_cpp( x, branch_loc_add, node_children_pointer, branch_counts, subtree_counts, loggamma ) + help::log_like_cpp( x, branch_loc_del, node_children_pointer, branch_counts, subtree_counts, loggamma ) + help::log_phi_cpp( wrap( phi ), sigma2_phi, wrap( zeta ) )  );
      }else{
        r = help::log_like_cpp( x, branch_loc_add, node_children_pointer, branch_counts, subtree_counts, loggamma_proposal ) + help::log_phi_cpp( wrap( phi_proposal ), sigma2_phi, wrap( zeta_proposal ) )  - ( help::log_like_cpp( x, branch_loc_add, node_children_pointer, branch_counts, subtree_counts, loggamma ) + help::log_phi_cpp( wrap( phi ), sigma2_phi, wrap( zeta ) ));
      }
      
      // Calculate acceptance probability
      double a  = log( runif( 1 )[ 0 ] );
      
      
      // Determine acceptance
      if( a < r ){
        zeta[ add_vcp[ 0 ] ] = zeta_proposal[ add_vcp[ 0 ] ];
        phi[ add_vcp[ 0 ] ] = phi_proposal[ add_vcp[ 0 ] ];
        zeta[ delete_vcp[ 0 ] ] = zeta_proposal[ delete_vcp[ 0 ] ];
        phi[ delete_vcp[ 0 ] ] = phi_proposal[ delete_vcp[ 0 ] ];
        
        loggamma = loggamma_proposal;
      }
      
    } // If not missing (still swap)
    
  } // If swap
  
  if( swap == 0 ){
    
    // Choose a branch 
    IntegerVector branches = seq( 0, B-1 );
    int branch_j = sample_cpp( branches );
    
    // Choose a covariate in branch
    IntegerVector cov = seq( 0, P-1 );
    int p = sample_cpp( cov );
    
    // Set proposal to current
    arma::mat zeta_proposal( B, P );
    arma::mat phi_proposal( B, P );
    zeta_proposal = zeta;
    phi_proposal = phi;
    
    // If included propose delete
    if( zeta( branch_j, p ) == 1 ){
      
      // Update proposal
      zeta_proposal( branch_j, p ) = 0;
      phi_proposal( branch_j, p ) = 0;
      
      // Get branch location to simplify likelihood calculation.
      int branch_loc = branch_location[ branch_j ] - 1 ;
      
      // Adjust loggamma for proposal
      arma::mat loggamma_proposal = loggamma;
      
      for( int n = 0; n < N; ++n ){
        loggamma_proposal( n , branch_j ) = loggamma( n , branch_j ) +   phi_proposal( branch_j, p )*x( n, p )  -   phi( branch_j, p )*x( n, p ) ;
      }
      
      // Calculate ratio
      double r = help::log_like_cpp( x, branch_loc, node_children_pointer, branch_counts, subtree_counts, loggamma_proposal ) + help::log_phi_cpp( wrap( phi_proposal ), sigma2_phi, wrap( zeta_proposal ) ) + help::log_zeta_pj_cpp( 0, aa, bb ) - ( help::log_like_cpp( x, branch_loc, node_children_pointer, branch_counts, subtree_counts, loggamma ) + help::log_phi_cpp( wrap( phi ), sigma2_phi, wrap( zeta ) ) + help::log_zeta_pj_cpp( 1, aa, bb ) ) ;
      
      // Calculate acceptance probability
      double a  = log( runif( 1 )[ 0 ] );
      
      // Determine acceptance
      if( a < r ){
        zeta( branch_j, p ) = zeta_proposal( branch_j, p );
        phi( branch_j, p ) = phi_proposal( branch_j, p );
        
        loggamma = loggamma_proposal;
      }
      
    }else{
      // If excluded propose inclusion
      
      // Update proposal
      zeta_proposal( branch_j, p ) = 1;
      double current_phi = phi( branch_j, p );
      phi_proposal( branch_j, p ) = help::phi_prop_cpp( current_phi );
      
      // Get branch location of added to simplify likelihood calculation.
      int branch_loc = branch_location[ branch_j ] - 1 ;   
      
      // Adjust loggamma for proposal
      arma::mat loggamma_proposal = loggamma;
      
      for( int n = 0; n < N; ++n ){
        loggamma_proposal( n , branch_j ) = loggamma( n , branch_j ) +  phi_proposal( branch_j, p )*x( n, p )  - phi( branch_j, p )*x( n, p );
      }
      
      // Calculate ratio
      double r = help::log_like_cpp( x, branch_loc, node_children_pointer, branch_counts, subtree_counts, loggamma_proposal ) + help::log_phi_cpp( wrap( phi_proposal ), sigma2_phi, wrap( zeta_proposal ) ) + help::log_zeta_pj_cpp( 1, aa, bb ) - ( help::log_like_cpp( x, branch_loc, node_children_pointer, branch_counts, subtree_counts, loggamma ) + help::log_phi_cpp( wrap( phi ), sigma2_phi, wrap( zeta ) ) + help::log_zeta_pj_cpp( 0, aa, bb) );
      
      // Calculate acceptance probability
      double a  = log( runif( 1 )[ 0 ] );
      
      // Determine acceptance
      if( a < r ){
        zeta( branch_j, p ) = zeta_proposal( branch_j, p );
        phi( branch_j, p ) = phi_proposal( branch_j, p );
        
        loggamma = loggamma_proposal;
      }
    }
  }  // Add/delete
  
  // Return output
  List between( 3 );
  between[ 0 ] = zeta;
  between[ 1 ] = phi;
  between[ 2 ] = loggamma;
  return between;
}


// Function :: Within Step ( update phi for zeta == 1)
List within_phi_update_cpp(
    arma::mat x,
    arma::mat zeta,
    arma::vec alpha,
    arma::mat phi,
    arma::vec branch_location,
    List node_children_pointer,
    arma::mat branch_counts,
    arma::mat subtree_counts,
    double sigma2_phi,
    arma::mat loggamma
){
  
  // Initiate memory space
  int N = x.n_rows;
  int B = zeta.n_rows;
  int P = zeta.n_cols;
  int dim = B*P;
  arma::mat phi_proposal( B, P );
  
  // For each covariate that is included
  for( int bp = 0; bp < dim; ++bp ) {
    phi_proposal = phi;
    
    if( zeta[ bp ] == 1 ){
      
      // Propose phi
      double current_phi = phi[ bp ];
      phi_proposal[ bp ] = help::phi_prop_cpp( current_phi );
      
      // Get branch location of added and deleted to simplify likelihood calculation.
      int which_mod = bp % B;
      int which_mod_P = floor( bp/B );
      int branch_loc = branch_location[ which_mod ] - 1 ;
      
      // Adjust loggamma for proposal
      arma::mat loggamma_proposal = loggamma;
      
      for( int n = 0; n < N; ++n ){
        loggamma_proposal( n , which_mod ) = loggamma( n , which_mod ) + phi_proposal( which_mod, which_mod_P )*x( n, which_mod_P )  -  phi( which_mod, which_mod_P )*x( n, which_mod_P ) ;
      }
      
      // Calculate ratio
      double r = help::log_like_cpp( x, branch_loc, node_children_pointer, branch_counts, subtree_counts, loggamma_proposal ) + help::log_phi_cpp( wrap( phi_proposal ), sigma2_phi, wrap(zeta) ) - ( help::log_like_cpp( x, branch_loc, node_children_pointer, branch_counts, subtree_counts, loggamma ) + help::log_phi_cpp( wrap( phi ), sigma2_phi, wrap(zeta) )  );
      
      // Calculate acceptance probability
      double a  = log( runif( 1 )[ 0 ] );
      
      // Determine acceptance
      if(a < r){
        phi[ bp ] = phi_proposal[ bp ];
        loggamma = loggamma_proposal;
      }
    }
  }
  
  // Return output
  List phi_return( 2 );
  phi_return[ 0 ] = phi;
  phi_return[ 1 ] = loggamma;
  return phi_return;
}

// Learn graphical structure among covariates
// Simulate MVT normal data
arma::mat mvrnormArma( int n, arma::vec mu, arma::mat sigma ) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn( n, ncols );
  return arma::repmat( mu, 1, n ).t()  + Y * arma::chol( sigma );
}

// Set V matrix based on current G
arma::mat set_V_cpp( arma::mat G, double v0, double v1 ){
  arma::mat G_in = G;
  arma::mat G_out = 1 - G;
  arma::mat V = v1*G_in + v0*G_out;
  int I = G.n_rows;
  
  for( int i = 0; i < I; ++i){
    V( i, i ) = 0;
  }
  
  return V;
}

// Update G matrix of inclusion indicators
arma::mat update_G_cpp(
    arma::mat Omega,
    double v0,
    double v1,
    double pie
){
  // Initiate Memory
  int p = Omega.n_cols;
  arma::mat newG( p, p );
  newG.zeros();
  
  double sd_v1 = sqrt( v1 );
  double sd_v0 = sqrt( v0 );
  
  // Update G
  for( int i = 0; i < p; ++i ){
    for( int j = 0; j < p; ++j ){
      if( i == j ){
        newG( i, j ) = 0;
      }else{
        double N_out = R::dnorm( Omega( i, j ), 0, sd_v0, false );
        double N_in = R::dnorm( Omega( i, j ), 0, sd_v1, false  );
        double prob = ( N_in*pie )/( N_in*pie + N_out*( 1 - pie ) );
        newG( i, j) = rbinom( 1, 1, prob )[ 0 ];
      }
    }
  }
  return newG;
}


// Update precision matrix Omega
arma::mat Omega_update_cpp(
    arma::mat Omega,
    arma::mat S,
    arma::mat V,
    double lambda,
    int n
){
  
  //Initiate Memory
  int p = Omega.n_cols;
  arma::mat Omega_update = Omega;
  
  for( int j = 0; j < p; ++j ){
    Omega_update.shed_col( j );
    Omega_update.shed_row( j );
    
    arma::mat S_temp = S;
    S_temp.shed_row( j );
    arma::vec S_12 = S_temp.col( j );
    
    double S_22 = S( j, j );
    
    arma::mat V_hold = V;
    V_hold.shed_row( j );
    arma::vec V_12 = V_hold.col( j );
    
    arma::mat diag_V_12( p - 1, p - 1 );
    diag_V_12.zeros();
    for( int i = 0; i < ( p - 1 ); ++i ){
      diag_V_12( i, i ) = pow( V_12[ i ], -1 );
    }
    
    arma::mat inv_Omega = inv_sympd( Omega_update );
    
    arma::mat C = inv_sympd( ( S_22 + lambda )*inv_Omega + diag_V_12 );
    
    arma::vec mu = -C*S_12;
    
    arma::vec new_u = help::mvrnormArma( 1, mu, C ).t();
    
    double new_v = rgamma( 1, n/2 + 1, 2/( S_22 + lambda )  )[ 0 ];
    arma::vec new_w22 = new_v + new_u.t()*inv_Omega*new_u;
    
    new_u.insert_rows( j, 1 );
    new_u( j ) = new_w22[ 0 ];
    
    Omega_update.insert_rows( j, 1 );
    Omega_update.insert_cols( j, 1 );
    
    Omega_update.row( j ) = new_u.t();
    Omega_update.col( j ) = new_u;
    
  }
  // Return output
  return Omega_update;
}

// sort of following http://www.johndcook.com/blog/standard_deviation/
double online_mean( int iteration, double last_mean, double curr_obs ){
  
  // don't fall off by one
  int n_vals = iteration++ ;
  
  // assume that iteration > 0 since there's a burnin proposal
  double curr_mean = last_mean + (curr_obs - last_mean)/(double)n_vals;
  
  return( curr_mean );
}

// Calculate online variance for adaptive sampler 
double online_var( int iteration, double last_mean, double last_var, double curr_mean, double curr_obs ){
  
  // note that iteration == n - 1 and that (n - 1) * last_var == last_ss
  // assume that iteration > 0 since there's a burnin proposal
  double curr_ss = (double)iteration * last_var + (curr_obs - last_mean) * (curr_obs - curr_mean);
  
  return( curr_ss/(double)iteration );
}

// simple adaptive MH proposal
// this uses equation (3) from Roberts & Rosenthal (2009) but without the full
// covariance of the target which makes it similar to the component-wise update
// of Haario et al. (2005)
double adap_prop( double curr_var, int n_cats, int n_vars ){
  
  // need dimension (???)
  double dd = (double)n_cats * (double)n_vars;
  
  // ensures bounded convergence
  int safeguard =  rbinom( 1, 1, 0.05 )[ 0 ];
  double usually = pow( 2.38, 2.0 ) * curr_var/dd;
  double unusually = pow( 0.1, 2.0 )/dd;
  
  double prop_var = (1 - safeguard) * rnorm( 1, 0, pow( usually, 0.5 ) )[ 0 ]  + safeguard * rnorm( 1, 0, pow( unusually, 0.5 ) )[ 0 ];
  
  return(prop_var);
}

// log Beta-Bernoulli spike-and-slab prior evaluation
double lprior_bbsas( double betajk, int sjk, double sig_bejk, double mu_bejk, double aa_hp, double bb_hp){
  
  // calculate additional beta factor
  double aa = sjk + aa_hp;
  double bb = 1 - sjk + bb_hp;
  double lbeta_factor = R::lbeta( aa, bb )  - R::lbeta( aa_hp, bb_hp );
  
  // piece-wise function
  if( sjk == 0 ){
    return( log( 1.0 - exp( lbeta_factor) ) );
  } else{
    return( lbeta_factor - 0.5 * log( 2.0 * M_PI * sig_bejk ) - 1.0/( 2.0 * sig_bejk)  * pow( betajk - mu_bejk, 2.0 ) );
  }
}

// calculates the linear predictor for each category of the Dirichlet-Multinomial
double calculate_gamma(arma::mat XX, arma::vec alpha, arma::vec beta, int jj, int ii, int Log, int n_vars){
  
  // loop indices and out
  int hh = 0;
  double out;
  
  out = alpha[ jj ];
  for( int kk = 0 ; kk < n_vars ; ++kk){
    
    hh = kk + jj * n_vars;
    out = out + beta[ hh ] * XX( ii, kk );
  }
  
  if(Log == 1){
    return(out);
  }else{
    return( exp( out ) ) ;
  }
}


// MH update for the regression parameters, Gibbs version
List update_beta_jj( arma::mat XX, arma::mat JJ, arma::mat loggamma,
                     arma::vec beta_temp, arma::vec inclusion_indicator,
                     arma::vec prop_per_beta, arma::mat mu_be, arma::mat sig_be , double aa_hp, double bb_hp, int jj, int n_vars ){
  
  // this function loops through n_vars so it is called for one taxa at a time
  int n_cats = JJ.n_cols; 
  int n_obs = XX.n_rows;
  double sig_prop;
  arma::vec loggamma_p( n_obs );
  loggamma_p.zeros();
  int hh, ii, kk;
  
  // define proposal variable and acceptance ratio variables
  double beta_p;
  double lnu, ln_acp;
  
  double log_full_beta, log_full_beta_p;
  
  for( kk = 0 ; kk < n_vars ; ++kk){
    
    // Stride for full (n_vars * n_cats) vector
    hh = kk + jj * n_vars;
    
    if( inclusion_indicator[hh] == 1 ){
      
      log_full_beta = 0;
      
      for( ii = 0 ; ii < n_obs ; ++ii ){
        log_full_beta = log_full_beta - lgamma( exp( loggamma( ii, jj ) ) );
        log_full_beta = log_full_beta + exp( loggamma( ii, jj ) ) * log( JJ( ii, jj ) );
      }
      
      log_full_beta = log_full_beta + lprior_bbsas( beta_temp[ kk ], inclusion_indicator[ hh ], sig_be( jj, kk ), mu_be( jj, kk ), aa_hp, bb_hp );
      sig_prop = prop_per_beta[ hh ];
      beta_p = beta_temp[ kk ] + help::adap_prop( sig_prop, n_cats, n_vars  );
      
      
      for( ii = 0 ; ii < n_obs ; ++ii ){
        loggamma_p[ ii ] = loggamma( ii, jj ) - beta_temp[ kk ] * XX( ii, kk ) + beta_p * XX( ii, kk );
      }
      
      // calculate proposal probability
      log_full_beta_p = 0;
      for( ii = 0 ; ii < n_obs ; ++ii ){
        log_full_beta_p = log_full_beta_p - lgamma( exp( loggamma_p[ ii ] ) );
        log_full_beta_p = log_full_beta_p + exp( loggamma_p[ ii ] ) * log( JJ( ii, jj ) );
      }
      log_full_beta_p = log_full_beta_p + lprior_bbsas( beta_p, inclusion_indicator[ hh ], sig_be( jj, kk ), mu_be( jj, kk ), aa_hp, bb_hp );
      
      ln_acp = log_full_beta_p - log_full_beta;
      lnu = log( runif( 1 )[ 0 ] );
      
      if(lnu < ln_acp){
        
        beta_temp[ kk ] = beta_p;
        
        for( ii = 0 ; ii < n_obs ; ++ii ){
          
          loggamma( ii, jj ) = loggamma_p[ ii ];
          
        }
        
      } // Close MH accept
    } // Close inclusion_indicator[hh] == 1
  } // Close for kk
  // Return output
  List out( 2 );
  out[ 0 ] = beta_temp;
  out[ 1 ] = loggamma;
  return out ;
  
} // Close function


// MH update for the regression parameters, stochastic search version
List update_beta_jj_ss( arma::mat XX, arma::mat JJ, arma::mat loggamma,
                        arma::vec beta_temp, arma::vec inclusion_indicator,
                        arma::vec prop_per_beta, arma::mat mu_be, arma::mat sig_be, double aa_hp, double bb_hp, int jj, int n_vars ){
  
  // this function loops through n_vars so it is called for one taxa at a time
  int n_cats = JJ.n_cols; 
  int n_obs = XX.n_rows;
  double sig_prop;
  arma::vec loggamma_p( n_obs );
  loggamma_p.zeros();
  int hh, ii, kk;
  
  // define proposal variable and acceptance ratio variables
  double beta_p;
  double lnu, ln_acp;
  
  double log_full_beta, log_full_beta_p;
  
  kk = help::sample_cpp( help::myseq( 0, n_vars - 1) );//for(kk = 0 ; kk < n_vars ; kk++){
  
  // Stride for full (n_vars * n_cats) vector
  hh = kk + jj * n_vars;
  
  if( inclusion_indicator[ hh ] == 1 ){
    
    log_full_beta = 0;
    
    for(ii = 0 ; ii < n_obs ; ++ii){
      log_full_beta = log_full_beta - lgamma( exp( loggamma( ii, jj ) ) ) ;
      log_full_beta = log_full_beta + exp( loggamma( ii, jj ) ) * log( JJ( ii, jj ) );
    }
    
    log_full_beta = log_full_beta + lprior_bbsas( beta_temp[ kk ], inclusion_indicator[ hh ], sig_be( jj, kk ), mu_be( jj, kk ), aa_hp, bb_hp );
    
    sig_prop = prop_per_beta[ hh ];
    beta_p = beta_temp[ kk ] + help::adap_prop( sig_prop, n_cats, n_vars  );
    
    
    for( ii = 0 ; ii < n_obs ; ++ii){
      loggamma_p[ii] = loggamma( ii, jj ) - beta_temp[ kk ] * XX( ii, kk ) + beta_p * XX( ii, kk );
    }
    
    // calculate proposal probability
    log_full_beta_p = 0;
    for( ii = 0 ; ii < n_obs ; ii++ ){
      log_full_beta_p = log_full_beta_p - lgamma( exp( loggamma_p[ ii ] ) );
      log_full_beta_p = log_full_beta_p + exp( loggamma_p[ii] ) * log( JJ( ii, jj ) ) ;
    }
    log_full_beta_p = log_full_beta_p + help::lprior_bbsas( beta_p, inclusion_indicator[ hh ], sig_be( jj, kk ), mu_be( jj, kk ), aa_hp, bb_hp );
    
    ln_acp = log_full_beta_p - log_full_beta;
    lnu = log( runif( 1 )[ 0 ] );
    
    if(lnu < ln_acp){
      
      beta_temp[ kk ] = beta_p;
      
      for( ii = 0 ; ii < n_obs ; ++ii ){
        
        loggamma( ii, jj ) = loggamma_p[ ii ];
        
      }
      
    } // Close MH accept
  } // Close inclusion_indicator[hh] == 1
  
  // Return output
  List out( 2 );
  out[ 0 ] = beta_temp;
  out[ 1 ] = loggamma;
  return out ;
  
  //} // Close for kk
} // Close function

// MH update for the intercept parameters
List update_alpha_jj( arma::mat JJ, arma::mat loggamma, arma::vec alpha, arma::vec prop_per_alpha, arma::vec mu_al, arma::vec sig_al, int jj){
  
  int ii;
  
  double sig_prop = prop_per_alpha[ jj ];
  int n_obs = loggamma.n_rows;
  arma::vec loggamma_p(n_obs);
  loggamma_p.zeros();
  
  double alpha_p;
  double lnu, ln_acp;
  
  // Prepare the current and proposed full conditional values
  double log_full_alpha, log_full_alpha_p;
  
  // Calculate the full conditional for the current value
  log_full_alpha = 0;
  for( ii = 0 ; ii < n_obs ; ++ii ){
    
    log_full_alpha = log_full_alpha - lgamma( exp( loggamma( ii, jj ) ) );
    log_full_alpha = log_full_alpha + exp( loggamma( ii, jj ) ) * log( JJ( ii, jj ) );
    
  }
  
  log_full_alpha = log_full_alpha - 1.0/( 2.0 * sig_al[ jj ] ) * pow( alpha[ jj ] - mu_al[ jj ], 2.0 );
  
  // Propose a new value for alpha[jj] using a random walk proposal centered on
  // the current value of alpha[jj]
  alpha_p = alpha[ jj ] + rnorm( 1, 0, sig_prop )[ 0 ];
  
  // Gamma must be updated too
  for(ii = 0 ; ii < n_obs ; ii++){
    loggamma_p[ ii ] = loggamma( ii, jj ) - alpha[ jj ] + alpha_p;
  }
  
  // Calculate the full conditional for the proposed value
  log_full_alpha_p = 0;
  for( ii = 0 ; ii < n_obs ; ++ii ){
    log_full_alpha_p = log_full_alpha_p - lgamma( exp( loggamma_p[ ii ] ) );
    log_full_alpha_p = log_full_alpha_p + exp( loggamma_p[ ii ] ) * log( JJ( ii, jj ) );
  }
  
  log_full_alpha_p = log_full_alpha_p - 1.0/(2.0 * sig_al[ jj ] ) * pow( alpha_p - mu_al[ jj ], 2.0);
  ln_acp = log_full_alpha_p - log_full_alpha;
  lnu = log( runif( 1 )[ 0 ] );
  
  if(lnu < ln_acp){
    
    // If accepted, update both alpha[jj] and loggamma[ii][jj], and keep
    alpha[ jj ] = alpha_p;
    
    for( ii = 0 ; ii < n_obs ; ++ii ){
      loggamma( ii, jj ) = loggamma_p[ ii ];
    }
  }
  
  // Return output
  List out( 2 );
  out[ 0 ] = alpha;
  out[ 1 ] = loggamma;
  return out ;
} // Close function


// for the Savitsky et al. inclusion proposal step
List between_models_jj( arma::mat XX, arma::mat JJ, arma::mat loggamma,
                        arma::vec beta_temp, arma::vec inclusion_indicator, arma::mat mu_be, 
                        arma::mat sig_be, double aa_hp, double bb_hp, int jj ){
  
  // Metropolis-Hastings proposals
  int n_vars = XX.n_cols;
  double sig_prop;
  double mu_prop;
  
  int hh, kk, ii;
  
  // proposed value vectors
  int inclusion_indicator_p;
  int n_obs = XX.n_rows;
  arma::vec loggamma_p( n_obs );
  loggamma_p.zeros();
  double beta_p;
  
  // acceptance ratio variables
  double lnu, ln_acp;
  
  // must update beta_jk for kk = 1, ..., n_vars
  double log_full_beta, log_full_beta_p;
  
  for( kk = 0 ; kk < n_vars ; ++kk ){
    
    hh = kk + jj * n_vars;
    
    // calculate the current log full conditional
    log_full_beta = 0;
    for( ii = 0 ; ii < n_obs ; ++ii ){
      log_full_beta = log_full_beta - lgamma( exp( loggamma( ii, jj ) ) );
      log_full_beta = log_full_beta + exp( loggamma( ii, jj ) ) * log( JJ( ii, jj ) );
    }
    log_full_beta = log_full_beta + help::lprior_bbsas( beta_temp[ kk ], inclusion_indicator[ hh ], sig_be( jj, kk ), mu_be( jj, kk ), aa_hp, bb_hp );
    
    // proposing a new value for beta[jj][kk] using the prior mean and standard deviation
    
    sig_prop = pow( sig_be( jj, kk ), 0.5);
    mu_prop = mu_be( jj, kk );
    
    // swap and sample beta
    if( inclusion_indicator[ hh ] == 0 ){
      beta_p = mu_prop + rnorm( 1, 0, sig_prop)[ 0 ];
      inclusion_indicator_p = 1;
    }else{
      beta_p = 0.0;
      inclusion_indicator_p = 0;
    }
    
    // update proposal gamma (linear predictor)
    for( ii = 0 ; ii < n_obs ; ++ii ){
      loggamma_p[ ii ] = loggamma( ii, jj ) - beta_temp[ kk ] * XX( ii, kk ) + beta_p * XX( ii, kk );
    }
    
    // calculate the proposed log full conditional
    log_full_beta_p = 0;
    for( ii = 0 ; ii < n_obs ; ++ii ){
      log_full_beta_p = log_full_beta_p - lgamma( exp( loggamma_p[ ii ] ) );
      log_full_beta_p = log_full_beta_p + exp( loggamma_p[ ii ] ) * log( JJ( ii, jj ) );
    }
    log_full_beta_p = log_full_beta_p + lprior_bbsas( beta_p, inclusion_indicator_p, sig_be( jj, kk ), mu_be( jj, kk ), aa_hp, bb_hp );
    
    ln_acp = log_full_beta_p - log_full_beta;
    lnu = log( runif( 1 )[ 0 ] );
    
    // Metropolis-Hastings ratio
    if( lnu < ln_acp ){
      
      // update accepted beta into beta_temp
      beta_temp[ kk ] = beta_p;
      inclusion_indicator[ hh ] = inclusion_indicator_p;
      
      for(ii = 0 ; ii < n_obs ; ii++){
        // also update the gamma (linear predictor)
        loggamma( ii, jj ) = loggamma_p[ ii ];
      }
    } // Close MH if
  } // Close for kk
  // Return output
  List out( 3 );
  out[ 0 ] = beta_temp;
  out[ 1 ] = inclusion_indicator;
  out[ 2 ] = loggamma;
  return out;
} // Close function

} // For namespace 'help'

// Function :: MCMC algorithm
// [[Rcpp::export]]
List DTMbvs(
    int iterations,
    int thin,
    String prior,
    arma::mat x,
    arma::vec branch_location,
    List node_children_pointer,
    arma::mat branch_counts,
    arma::mat subtree_counts,
    arma::mat alpha,
    arma::cube phi,
    arma::cube zeta,
    double sigma2_alpha,
    double sigma2_phi,
    double aa,
    double bb,
    arma::cube Omega,
    arma::cube G,
    arma::cube Var,
    arma::mat S,
    double v0,
    double v1,
    double a_G,
    double b_G,
    double pie,
    double lambda
  
){
  // Initiate memory
  List alpha_return( 2 );
  List between_phi_zeta( 3 );
  List phi_return( 2 );
  int B = alpha.n_rows;
  int N = x.n_rows;
  
  // Set current loggamma
  arma::mat loggamma( N, B );
  loggamma.zeros();
  
  loggamma = x*phi.slice( 0 ).t();
  
  for( int n = 0; n < N; ++n ){
    for( int b = 0; b < B; ++b ){
      loggamma( n, b ) = loggamma( n, b ) + alpha.col( 0 )[ b ] ;
    }
  }
  
  // Set temporary data to enable thinning
  arma::vec temp_alpha = alpha.col( 0 );
  arma::mat temp_zeta = zeta.slice( 0 );
  arma::mat temp_phi = phi.slice( 0 );
  arma::mat temp_Omega = Omega.slice( 0 );
  arma::mat temp_Var = Var.slice( 0 );
  arma::mat temp_G = G.slice( 0 );
  
  // Looping over the number of iterations specified by user
  for( int iter = 0; iter < iterations; ++iter ){
    
    alpha_return =  help::alpha_update_cpp( x, temp_alpha, temp_phi, branch_location, node_children_pointer, branch_counts, subtree_counts, sigma2_alpha, loggamma );
    temp_alpha = as<arma::vec>( alpha_return[ 0 ] );
    loggamma = as<arma::mat>( alpha_return[ 1 ] );
    
    // Swap/Add/Delete or just Add/Delete
    // adjust between based on the prior used
    if( prior == "BB" ){
      between_phi_zeta = help::between_phi_zeta_swap_update_BB_cpp( x, temp_zeta, temp_alpha, temp_phi, branch_location, node_children_pointer, branch_counts, subtree_counts, aa, bb, sigma2_phi, loggamma );
    }
    
    if( prior == "MRF_fixed" ){
      between_phi_zeta = help::between_phi_zeta_swap_update_MRF_cpp( x, temp_zeta, temp_alpha, temp_phi, branch_location, node_children_pointer, branch_counts, subtree_counts, sigma2_phi, a_G, b_G, temp_G, loggamma );
    }
    
    if( prior == "MRF_unknown" ){
      between_phi_zeta = help::between_phi_zeta_swap_update_MRF_cpp( x, temp_zeta, temp_alpha, temp_phi, branch_location, node_children_pointer, branch_counts, subtree_counts, sigma2_phi, a_G, b_G, temp_G, loggamma );
    }
    
    temp_zeta = as<arma::mat>( between_phi_zeta[ 0 ] );
    temp_phi = as<arma::mat>( between_phi_zeta[ 1 ] );
    loggamma =  as<arma::mat>( between_phi_zeta[ 2 ] );
    
    // Update within
    phi_return = help::within_phi_update_cpp(x, temp_zeta, temp_alpha, temp_phi, branch_location, node_children_pointer, branch_counts, subtree_counts, sigma2_phi, loggamma );
    temp_phi = as<arma::mat>( phi_return[ 0 ] );
    loggamma = as<arma::mat>( phi_return[ 1 ] );
    
    // Update Gaussian graphical model if unknown
    if( prior == "MRF_unknown" ){
      temp_Omega = help::Omega_update_cpp( temp_Omega, S, temp_Var, lambda, N );
      temp_G = help::update_G_cpp( temp_Omega, v0, v1, pie );
      temp_Var = help::set_V_cpp( temp_G, v0, v1 );
    }
    
    // Set the starting values for the next iteration
    if( ( iter + 1 ) % thin == 0 ){
      alpha.col( ( iter + 1 )/thin - 1 ) = temp_alpha;
      zeta.slice( ( iter + 1 )/thin - 1 ) = temp_zeta;
      phi.slice( ( iter + 1 )/thin - 1 ) = temp_phi;
      Omega.slice( ( iter + 1 )/thin - 1 ) =  temp_Omega;
      Var.slice( ( iter + 1 )/thin - 1 ) = temp_Var;
      G.slice( ( iter + 1 )/thin - 1 ) = temp_G;
    }
    
    // Print out progress
    double printer = iter % 250;
    
    if( printer == 0 ){
      Rcpp::Rcout << "Iteration = " << iter << std::endl;
    }
  }
  
  // Return output
  List output( 5 );
  output[ 0 ] = alpha;
  output[ 1 ] = zeta;
  output[ 2 ] = phi;
  output[ 3 ] = Omega;
  output[ 4 ] = G;
  return output ;
}

// SSVS Sampler function
// [[Rcpp::export]]
List dmbvs_ss( arma::mat XX, arma::mat YY, arma::vec alpha, arma::vec beta,
               arma::vec mu_al, arma::vec sig_al,
               arma::mat mu_be, arma::mat sig_be,
               double aa_hp, double bb_hp, arma::vec prop_per_alpha,
               arma::vec prop_per_beta, int GG,
               int thin, int Log ){
  
  // loop indices
  int hh, ii, jj, kk, gg;
  
  // adaptive proposal values
  double last_mean;
  double last_var;
  
  ///////// Initialize variables that don't require input or can be obtained
  ///////// from inputs
  
  int n_vars = XX.n_cols;
  int n_cats = YY.n_cols;
  int n_obs = YY.n_rows;
  
  // for adaptive MH need to keep track of target density mean and variance
  arma::vec curr_mean( n_vars * n_cats );
  arma::vec curr_var( n_vars * n_cats );
  
  curr_mean.zeros();
  curr_var.zeros();
  for( hh = 0 ; hh < (n_vars * n_cats) ; ++hh){
    curr_var[ hh ] = 0.5;
  }
  
  // variable inclusion indicator
  arma::vec inclusion_indicator( n_vars * n_cats );
  inclusion_indicator.zeros();
  
  for(hh = 0 ; hh < (n_vars * n_cats) ; ++hh){
    if( beta[ hh ] == 0 ){
      inclusion_indicator[ hh ] = 0;
    }else{
      inclusion_indicator[ hh ] = 1;
    }
  }
  
  // row sums of YY
  arma::vec Ypiu( n_obs );
  Ypiu.zeros();
  
  for( ii = 0 ; ii < n_obs ; ++ii){
    for(jj = 0 ; jj < n_cats ; ++jj){
      Ypiu[ ii ] = Ypiu[ ii ] + YY( ii, jj );
    }
  }
  
  // the beta_temp vector is an n_vars long vector used for updating each of the
  // jj-th rows in Beta corresponding the jj-th category in YY
  arma::vec beta_temp( n_vars);
  beta_temp.zeros(); 
  
  // the linear predictor matrix, represented as a vector, is in the log scale
  arma::mat loggamma( n_obs, n_cats );
  loggamma.zeros();
  for( ii = 0 ; ii < n_obs ; ++ii ){
    for( jj = 0 ; jj < n_cats ; ++jj ){
      loggamma( ii, jj ) = help::calculate_gamma( XX, alpha, beta, jj, ii, Log, n_vars );
    }
  }
  
  // initialize latent variable: JJ ~ gamma()
  arma::mat JJ( n_obs, n_cats );
  JJ.zeros();
  // TT is the vector of normalization constants
  arma::vec TT( n_obs );
  TT.zeros();
  
  for( ii = 0 ; ii < n_obs ; ++ii ){
    for( jj = 0 ; jj < n_cats ; ++jj ){
      JJ( ii, jj ) =  YY( ii, jj );
      // a very small number
      if( JJ( ii, jj )  < pow(10.0, -100.0 ) ){
        JJ( ii, jj ) = pow(10.0, -100.0);
      }
      TT[ ii ] = TT[ ii ] + JJ( ii, jj );
    }
  }
  
  // initialize the other clever latent variable: uu
  // note that this uses the data to initialize
  arma::vec uu( n_obs );
  uu.zeros();
  for(ii = 0 ; ii < n_obs ; ++ii){
    uu[ ii ] = rgamma( 1, Ypiu[ ii ], 1/TT[ ii ])[ 0 ] ;
  }
  
  arma::mat out_alpha( n_cats, floor( GG/thin ) ); 
  arma::mat out_beta( n_cats*n_vars , floor( GG/thin ) ); 
  
  out_alpha.zeros();
  out_beta.zeros();
  
  // the magic starts here
  List between_return( 3 );
  List beta_return( 2 );
  List alpha_return( 2 );
  
  for( gg = 0 ; gg < GG ; ++gg){
    
    // Print out progress
    double printer = gg % 250;
    
    if( printer == 0 ){
      Rcpp::Rcout << "Iteration = " << gg << std::endl;
    }
    
    // first a round of the between-model step for every covariate within every taxa
    jj = help::sample_cpp( help::myseq( 0, n_cats - 1) );
    
    // fill in beta_temp
    for( kk = 0 ; kk < n_vars ; ++kk ){
      hh = kk + jj * n_vars;
      beta_temp[ kk ] = beta[ hh ];
    }
    
    between_return = help::between_models_jj(XX, JJ, loggamma, beta_temp,
                                             inclusion_indicator, mu_be, sig_be, aa_hp, bb_hp, jj);
    
    beta_temp = as<arma::vec>( between_return[ 0 ] );
    inclusion_indicator = as<arma::vec>( between_return[ 1 ] );
    loggamma = as<arma::mat>( between_return[ 2 ] );
    
    // update beta with beta_temp
    for( kk = 0 ; kk < n_vars ; ++kk ){
      hh = kk + jj * n_vars;
      beta[ hh ] = beta_temp[ kk ];
    }
    //}
    
    // now an accelerating within-model step
    for( jj = 0 ; jj < n_cats ; ++jj ){
      
      // fill in beta_temp and update the proposal variance
      for( kk = 0 ; kk < n_vars ; ++kk ){
        hh = kk + jj * n_vars;
        beta_temp[ kk ] = beta[ hh ];
        // wait until each parameter has had a few iterations
        if( gg > 2 * n_vars * n_cats ){
          // calculate the online variance
          last_mean = curr_mean[ hh ];
          last_var = curr_var[ hh ];
          curr_mean[ hh ] = help::online_mean( gg, last_mean, beta[ hh ] );
          curr_var[ hh ] = help::online_var( gg, last_mean, last_var, curr_mean[ hh ], beta[ hh ] );
          // update proposal variance
          prop_per_beta[ hh ] = curr_var[ hh ];
        }
      }
      
      beta_return = help::update_beta_jj_ss(XX, JJ, loggamma, beta_temp, inclusion_indicator,
                                            prop_per_beta, mu_be, sig_be, aa_hp, bb_hp, jj, n_vars );
      
      beta_temp = as<arma::vec>( beta_return[ 0 ] );
      loggamma = as<arma::mat>( beta_return[ 1 ] );
      
      // update beta with beta_temp and write to file
      for( kk = 0 ; kk < n_vars ; ++kk ){
        hh = kk + jj * n_vars;
        beta[ hh ] = beta_temp[ kk ];
      }
    }
    
    
    // update alpha and write to file
    for( jj = 0 ; jj < n_cats ; ++jj ){
      
      alpha_return = help::update_alpha_jj(JJ, loggamma, alpha, prop_per_alpha,
                                           mu_al, sig_al, jj);
      alpha = as<arma::vec>( alpha_return[ 0 ] );
      loggamma = as<arma::mat>( alpha_return[ 1 ] );
      
    }
    
    // update JJ and consequently TT
    for( ii = 0 ; ii < n_obs ; ++ii ){
      TT[ ii ] = 0.0;
      for( jj = 0 ; jj < n_cats ; ++jj ){
        JJ( ii, jj ) = rgamma(1, YY( ii, jj ) + exp(loggamma( ii, jj )), 1/(uu[ii] + 1.0) )[ 0 ] ;
        if(JJ( ii, jj ) < pow(10.0, -100.0)){
          JJ( ii, jj ) = pow(10.0, -100.0);
        }
        TT[ ii ] = TT[ ii ] + JJ( ii, jj );
      }
    }
    
    // update latent variables uu
    for( ii = 0 ; ii < n_obs ; ++ii ){
      uu[ ii ] = rgamma( 1, Ypiu[ii], 1/TT[ ii ])[ 0 ];
    }
    
    
    if( ( gg + 1 ) % thin == 0 ){
      out_alpha.col( ( gg + 1 )/thin - 1  ) = alpha;
      out_beta.col( ( gg + 1 )/thin - 1  ) = beta;
    }
    
    
    
  } // end of iterations
  
  
  // Return output
  List output( 2 );
  output[ 0 ] = out_alpha;
  output[ 1 ] = out_beta;
  
  return output ;
  
} // close SSVS sampler function



// Function :: MCMC algorithm
// Principal function using Gibbs variable selection
// [[Rcpp::export]]
List dmbvs_gibbs(arma::mat XX, arma::mat YY, arma::vec alpha, arma::vec beta,
                 arma::vec mu_al, arma::vec sig_al , arma::mat mu_be, arma::mat sig_be,
                 double aa_hp, double bb_hp, arma::vec prop_per_alpha,
                 arma::vec prop_per_beta, int GG,
                 int thin, int Log ){
  
  // loop indices
  int hh, ii, jj, kk, gg;
  
  // adaptive proposal values
  double last_mean;
  double last_var;
  
  ///////// Initialize variables that don't require input or can be obtained
  ///////// from inputs
  int n_vars = XX.n_cols;
  int n_cats = YY.n_cols;
  int n_obs = YY.n_rows;
  
  // for adaptive MH need to keep track of target density mean and variance
  arma::vec curr_mean( n_vars * n_cats );
  arma::vec curr_var( n_vars * n_cats );
  
  curr_mean.zeros();
  curr_var.zeros();
  for( hh = 0 ; hh < (n_vars * n_cats) ; ++hh){
    curr_var[ hh ] = 0.5;
  }
  
  // variable inclusion indicator
  arma::vec inclusion_indicator( n_vars * n_cats );
  inclusion_indicator.zeros();
  
  for(hh = 0 ; hh < (n_vars * n_cats) ; ++hh){
    if( beta[ hh ] == 0 ){
      inclusion_indicator[ hh ] = 0;
    }else{
      inclusion_indicator[ hh ] = 1;
    }
  }
  
  // row sums of YY
  arma::vec Ypiu( n_obs );
  Ypiu.zeros();
  
  for( ii = 0 ; ii < n_obs ; ++ii){
    for(jj = 0 ; jj < n_cats ; ++jj){
      Ypiu[ ii ] = Ypiu[ ii ] + YY( ii, jj );
    }
  }
  
  // the beta_temp vector is an n_vars long vector used for updating each of the
  // jj-th rows in Beta corresponding the jj-th category in YY
  arma::vec beta_temp( n_vars);
  beta_temp.zeros(); 
  
  // the linear predictor matrix, represented as a vector, is in the log scale
  arma::mat loggamma( n_obs, n_cats );
  loggamma.zeros();
  for( ii = 0 ; ii < n_obs ; ++ii ){
    for( jj = 0 ; jj < n_cats ; ++jj ){
      loggamma( ii, jj ) = help::calculate_gamma( XX, alpha, beta, jj, ii, Log, n_vars );
    }
  }
  
  // initialize latent variable: JJ ~ gamma()
  arma::mat JJ( n_obs, n_cats );
  JJ.zeros();
  
  // TT is the vector of normalization constants
  arma::vec TT( n_obs );
  TT.zeros();
  
  for( ii = 0 ; ii < n_obs ; ++ii ){
    for( jj = 0 ; jj < n_cats ; ++jj ){
      JJ( ii, jj ) =  YY( ii, jj );
      // a very small number
      if( JJ( ii, jj )  < pow(10.0, -100.0 ) ){
        JJ( ii, jj ) = pow(10.0, -100.0);
      }
      TT[ ii ] = TT[ ii ] + JJ( ii, jj );
    }
  }
  
  // initialize the other clever latent variable: uu
  // note that this uses the data to initialize
  arma::vec uu( n_obs );
  uu.zeros();
  for(ii = 0 ; ii < n_obs ; ++ii){
    uu[ ii ] = rgamma( 1, Ypiu[ ii ], 1/TT[ ii ])[ 0 ] ;
  }
  
  
  arma::mat out_alpha( n_cats, floor( GG/thin ) ); 
  arma::mat out_beta( n_cats*n_vars , floor( GG/thin ) ); 
  
  out_alpha.zeros();
  out_beta.zeros();
  
  // the magic starts here
  List between_return( 3 );
  List beta_return( 2 );
  List alpha_return( 3 );
  
  for( gg = 0 ; gg < GG ; ++gg ){
    
    // timing info to the printed output
    // Print out progress
    double printer = gg % 250;
    
    if( printer == 0 ){
      Rcpp::Rcout << "Iteration = " << gg << std::endl;
    }
    
    // first a round of the between-model step for every covariate within every taxa
    for( jj = 0 ; jj < n_cats ; ++jj ){
      
      // fill in beta_temp
      for( kk = 0 ; kk < n_vars ; ++kk ){
        hh = kk + jj * n_vars;
        beta_temp[ kk ] = beta[ hh ];
      }
      
      between_return = help::between_models_jj(XX, JJ, loggamma, beta_temp,
                                               inclusion_indicator, mu_be, sig_be, aa_hp, bb_hp, jj);
      
      beta_temp = as<arma::vec>( between_return[ 0 ] );
      inclusion_indicator = as<arma::vec>( between_return[ 1 ] );
      loggamma = as<arma::mat>( between_return[ 2 ] );
      
      // update beta with beta_temp
      for(kk = 0 ; kk < n_vars ; kk++){
        hh = kk + jj * n_vars;
        beta[ hh ] = beta_temp[ kk ];
      }
      
    }
    
    // now an accelerating within-model step
    for( jj = 0 ; jj < n_cats ; ++jj ){
      
      // fill in beta_temp and update the proposal variance
      for( kk = 0 ; kk < n_vars ; ++kk ){
        hh = kk + jj * n_vars;
        beta_temp[ kk ] = beta[ hh ];
        
        // wait until each parameter has had a few iterations
        if(gg > 2 * n_vars * n_cats){
          // calculate the online variance
          last_mean = curr_mean[ hh ];
          last_var = curr_var[ hh ];
          curr_mean[ hh ] = help::online_mean( gg, last_mean, beta[ hh ] );
          curr_var[ hh ] = help::online_var( gg, last_mean, last_var, curr_mean[ hh ], beta[ hh ] );
          // update proposal variance
          prop_per_beta[ hh ] = curr_var[ hh ];
        }
      }
      
      beta_return = help::update_beta_jj(XX, JJ, loggamma, beta_temp, inclusion_indicator,
                                         prop_per_beta, mu_be, sig_be, aa_hp, bb_hp, jj, n_vars );
      
      beta_temp = as<arma::vec>( beta_return[ 0 ] );
      loggamma = as<arma::mat>( beta_return[ 1 ] );
      
      // update beta with beta_temp and write to file
      for( kk = 0 ; kk < n_vars ; ++kk){
        hh = kk + jj * n_vars;
        beta[ hh ] = beta_temp[ kk ];
      }
    }
    
    // update alpha and write to file
    for( jj = 0 ; jj < n_cats ; ++jj){
      alpha_return = help::update_alpha_jj(JJ, loggamma, alpha, prop_per_alpha,
                                           mu_al, sig_al, jj);
      alpha = as<arma::vec>( alpha_return[ 0 ] );
      loggamma = as<arma::mat>( alpha_return[ 1 ] );
      
    } 
    
    // update JJ and consequently TT
    for( ii = 0 ; ii < n_obs ; ++ii ){
      TT[ ii ] = 0.0;
      for( jj = 0 ; jj < n_cats ; ++jj ){
        JJ( ii, jj ) = rgamma(1, YY( ii, jj ) + exp(loggamma( ii, jj )), 1/(uu[ii] + 1.0) )[ 0 ] ;
        if(JJ( ii, jj ) < pow(10.0, -100.0)){
          JJ( ii, jj ) = pow(10.0, -100.0);
        }
        TT[ ii ] = TT[ ii ] + JJ( ii, jj );
      }
    }
    
    // update latent variables uu
    for( ii = 0 ; ii < n_obs ; ++ii ){
      uu[ ii ] = rgamma( 1, Ypiu[ii], 1/TT[ ii ])[ 0 ];
    }
    
    if( ( gg + 1 ) % thin == 0 ){
      out_alpha.col( ( gg + 1 )/thin - 1  ) = alpha;
      out_beta.col( ( gg + 1 )/thin - 1  ) = beta;
    }
    
  } // end of iterations
  
  
  // Return output
  List output( 2 );
  output[ 0 ] = out_alpha;
  output[ 1 ] = out_beta;
  
  return output ;
}
// close Gibbs (stochastic search) sampler function