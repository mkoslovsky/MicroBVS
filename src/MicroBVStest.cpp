// Code for A joint model for predicting phenotypic responses with human microbiome data
// Author: Matt Koslovsky

// Performs Bayesian variable selection for multivariate count data to predict
// continuous response while simulataneously identifying factors associated with
// multivariate count data


#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
// [[Rcpp::plugins("cpp11")]]

// Define helper functions
namespace help{

// Scale and center a matrix
mat scale_cpp( mat Bal ){
  int J = Bal.n_cols;
  int I = Bal.n_rows;
  vec means( J );
  
  for( int j = 0; j < J; ++j ){
    means[ j ] = sum( Bal.col( j ) )/I;
    double vars = sum(pow( Bal.col( j ) - means[ j ], 2 ) )/( I - 1 );
    Bal.col( j ) = ( Bal.col( j ) - means[ j ] )/sqrt( vars );
  }
  return( Bal );
  
}

// Call sample function from R
double sample_cpp( IntegerVector x ){
  Function f("sample");
  IntegerVector sampled = f(x, Named("size") = 1);
  return sampled[0];
}

// Construct Helmert Contrast Matrix for balances
arma::mat helmert_cpp( int p ){
  arma::mat H( p, p-1 );
  H.zeros();
  
  for( int j = 0; j < (p - 1); ++j ){
    for( int i = 0; i < p; ++i ){
      if( i <= j){
        H( i, j ) = -1;
      }
      if( i == ( j + 1 ) ){
        H( i, j ) = ( j + 1 );
      }
    }
  }
  
  return H;
}

// Get log derivative of matrix
double log_det_cpp( arma::mat x ){
  int p = x.n_cols;
  double sum = 0;
  
  arma::mat chole = chol(x);
  
  for( int j = 0; j < p; ++j ){
    sum += 2*log( chole( j, j ) );
  }
  
  return sum;
}

// Perform isometric log-ratio transformation on percentages (ie construct balances)
arma::mat ilr_fun_cpp( arma::mat x ){
  int n = x.n_rows;
  int p = x.n_cols;
  arma::mat y = log( x );
  arma::vec y_centered( n );
  
  // Center y
  for( int i = 0; i < n; ++i){
    arma::vec rows = y.row( i ).t();
    double means = mean( rows );
    y_centered( i ) = means;
  }
  
  for( int j = 0; j < p; ++j ){
    y.col( j ) = y.col(j) - y_centered;
  }
  
  arma::mat H = help::helmert_cpp( p );
  
  for( int k = 0; k < ( p - 1 ); ++k ){
    H.col( k ) = H.col( k )/sqrt( (k+1)*(k+2) );
  }
  
  arma::mat ilr = y*H;
  
  return ilr;
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
// Function :: Calculate log likelihood of Y ( continuous )
double log_mvtT_cpp( arma::mat sigma, arma::vec mu, double v, arma::vec y){
  int p = sigma.n_cols;
  arma::vec meanC = y - mu;
  arma::vec inside = meanC.t()*inv_sympd( sigma )*meanC;
  
  double log_mvtT =  - 0.5*help::log_det_cpp( sigma ) - ( ( v + p )/2 )*log( 1 + ( 1/v )*inside[0] );
  return log_mvtT;
}

// Function :: Calculate log likelihood of X ( compositional ) across i
double log_multinomial_cpp( arma::mat z, arma::mat psi ){
  int n = psi.n_rows;
  arma::vec log_mult( n );
  
  for( int i = 0; i < n; ++i ){
    log_mult( i ) = lgamma( sum( z.row( i ) ) + 1 ) + sum( - lgamma( z.row( i ) + 1 ) + z.row( i )%log( psi.row( i )) );
  }
  
  double sum_log_mult = sum( log_mult );
  
  return sum_log_mult;
}

// Function :: Calculate log psi contribution for all i
double log_dirichlet_psi_cpp( arma::mat x, arma::vec alpha, arma::mat phi, arma::mat psi ){
  int J = alpha.size();
  int N = x.n_rows;
  arma::mat gamma( N, J );
  arma::mat alpha_expand( N, J );
  arma::vec log_psi( N );
  log_psi.zeros();
  gamma.zeros();
  alpha_expand.zeros();
  
  for( int n = 0; n < N; ++n ){
    alpha_expand.row( n ) = alpha.t();
  }
  
  gamma =  exp( alpha_expand + x*phi.t() ) ;
  
  for( int n = 0; n < N; ++n ){
    log_psi( n ) =  lgamma( sum( gamma.row( n ) ) ) - sum( lgamma( gamma.row( n ) ) ) + sum( ( gamma.row( n ) - 1 )%log( psi.row( n ) ) );
    
  }
  
  
  double sum_log_psi = sum( log_psi );
  return sum_log_psi;
}

// Function :: Calculate log psi contribution for an individual
double log_dirichlet_psi_i_cpp( arma::mat x, arma::vec alpha, arma::mat phi, arma::vec psi, int i ){
  int J = alpha.size();
  arma::vec gamma( J );
  double log_psi = 0;
  
  gamma =  exp( alpha + phi*x.row( i ).t() ) ;
  
  log_psi =  lgamma( sum( gamma ) ) - sum( lgamma( gamma ) ) + sum( ( gamma - 1 )%log( psi ) );
  
  return log_psi;
}

// Function :: Calculate proposal ratio for update of psi
double log_dirichlet_psi_proposal_cpp( double delta, double rho, arma::vec psi, arma::vec psi_proposal ){
  arma::vec inside_proposal = rho*psi_proposal + delta;
  arma::vec inside = rho*psi + delta;
  
  double log_psi_num =  lgamma( sum( inside_proposal ) ) - sum( lgamma( inside_proposal ) ) + sum( ( inside_proposal - 1 )%log( psi ) );
  double log_psi_denom =  lgamma( sum( inside ) ) - sum( lgamma( inside ) ) + sum( ( inside - 1 )%log( psi_proposal ) );
  
  double ratio = log_psi_num - log_psi_denom;
  
  return ratio;
}
// Function :: Sample dirichlet distribution
arma::vec rdirichlet_cpp( arma::vec param, double delta, double rho ){
  int J = param.size();
  arma::vec dir_sample( J );
  dir_sample.zeros();
  double sum = 0;
  arma::vec param_scale_shift(J);
  param_scale_shift = param*rho + delta;
  
  for( int j = 0; j < J; ++j ){
    dir_sample[ j ] = rgamma(1, param_scale_shift[ j ], 1 )[0];
    sum += dir_sample[ j ];
  }
  
  // Protect against dividing by zero
  if(sum == 0){
    sum = 1;
  }
  
  dir_sample = dir_sample/sum;
  
  return dir_sample;
}

// Function :: Calculate log alpha contribution
double log_alpha_cpp( double alpha, double sigma2_alpha){
  // sum the log p(alpha) across all branches
  // initiate log_alpha and set size of vector
  double log_alpha = 0; 
  
  // Sum over all alpha
  log_alpha += -0.50*log( 2*atan(1)*4*sigma2_alpha ) - 1/( 2*sigma2_alpha )*pow( alpha, 2 );
  
  // Return output
  return log_alpha;
}

// Function :: Propose a new alpha
double alpha_prop_cpp( double alpha ){
  // Simulate proposal with jitter around current alpha
  // Make sure that alpha is specified with iteration iterations and branch b
  double alpha_prop = alpha + sqrt(0.5)*rnorm( 1 )[ 0 ];
  
  // Return output
  return alpha_prop ;
}

// Function :: Propose a new phi
double phi_prop_cpp( double phi ){
  // Simulate proposal with jitter around current phi
  // Make sure that phi is specified with iteration iterations and branch b
  double phi_prop = phi + sqrt(0.5)*rnorm( 1 )[ 0 ];
  
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

// Function :: Calculate log phi contribution for individual pj
double log_phi_pj_cpp( double phi, double sigma2_phi){
  
  double log_phi = -0.50*log( 2*atan(1)*4*sigma2_phi ) - 1/( 2*sigma2_phi )*pow(phi,2 );
  
  // Return output
  return log_phi;
}

// Function :: Calculate individual zeta
double log_zeta_pj_cpp( double t_pj, double a, double b
){
  double post_a = t_pj + a;
  double post_b = 1 - t_pj + b;
  double log_zeta_pj = lgamma( post_a ) + lgamma( post_b ) - lgamma( post_a + post_b ) - ( lgamma( a ) + lgamma( b ) - lgamma( a + b ) );
  
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


// Function :: Calculate individual xi
double log_xi_m_cpp( double t_m, double a_m, double b_m
){
  double post_a = t_m + a_m;
  double post_b = 1 - t_m + b_m;
  double log_xi_m = lgamma( post_a ) + lgamma( post_b ) - lgamma( post_a + post_b ) - ( lgamma( a_m ) + lgamma( b_m ) - lgamma( a_m + b_m ) );
  
  // Return output
  return log_xi_m ;
}


// Function :: Update alpha augment
// Make sure arguments are indexed by iteration
List alpha_update_augment_cpp(
    arma::mat loggamma,
    arma::mat cc, 
    arma::vec alpha,
    arma::mat phi, 
    double sigma2_alpha
){
  // Initiate memory space
  int J = phi.n_rows;
  int N = loggamma.n_rows;
  
  // Update each alpha
  for( int j = 0; j < J; ++j ){
    
    double log_like = 0;
    // Calculate log_like
    for( int n = 0; n < N; ++n ){
      log_like = log_like - lgamma( exp( loggamma( n , j ) ) ) + exp( loggamma( n , j ) )*log( cc( n, j ) );
    }
    
    // Adjust alpha and calculate log_like_prop
    // Proposal alpha
    double current_alpha = alpha[ j ];
    double alpha_proposal  = help::alpha_prop_cpp( current_alpha );
    
    
    arma::vec loggamma_prop( N );
    loggamma_prop.zeros();
    
    for( int n = 0 ; n < N ; ++n ){
      loggamma_prop[ n ] = loggamma( n, j ) - current_alpha + alpha_proposal;
    }
    
    double log_like_prop = 0; 
    for( int n = 0; n < N; ++n ){
      log_like_prop = log_like_prop - lgamma( exp( loggamma_prop[ n ] ) ) + exp( loggamma_prop[ n ] )*log( cc( n, j ) );
    }
    
    // Calculate ratio
    double r = log_like_prop + help::log_alpha_cpp( alpha_proposal, sigma2_alpha ) - ( log_like + help::log_alpha_cpp(  current_alpha, sigma2_alpha ) );
    
    // Calculate acceptance probability
    double a  = log( runif( 1 )[ 0 ] );
    
    // Determine acceptance
    if(a < r){
      alpha[ j ] = alpha_proposal;
      loggamma.col( j ) = loggamma_prop;
    }
  }
  // Return output
  List alpha_out( 2 );
  alpha_out[ 0 ] = alpha; 
  alpha_out[ 1 ] = loggamma;
  return alpha_out;
}





// Function :: Between Step (jointly update phi and zeta) with add/delete using BB prior
List between_phi_zeta_augment_update_BB_cpp(
    arma::mat x,
    arma::mat zeta,
    arma::vec alpha,
    arma::mat phi, 
    arma::mat loggamma, 
    arma::mat cc,
    double sigma2_phi,
    double a,
    double b
){
  
  // Between Model Step
  // Initiate memory space
  int J = zeta.n_rows;
  int P = zeta.n_cols; 
  int N = x.n_rows;
  
  // Get a covariate p to update in branch j 
  
  int j = help::sample_cpp( help::myseq( 0, J - 1) ); 
  for( int p = 0; p < P; ++p ){
    
    // Calculate current log likelihood 
    double log_like = 0; 
    
    for( int n = 0; n < N; ++n ){
      log_like = log_like - lgamma( exp( loggamma( n , j ) ) ) + exp( loggamma( n , j ) )*log( cc( n, j ) );
    }
    
    // Calculate proposed likelihood 
    double phi_prop;
    int zeta_prop; 
    
    if( zeta( j, p ) == 1 ){
      phi_prop = 0.0;
      zeta_prop = 0;
    }else{
      phi_prop = help::phi_prop_cpp( phi( j, p ) ); 
      zeta_prop = 1;
    }
    
    arma::vec loggamma_prop( N );
    loggamma_prop.zeros();
    
    for( int n = 0 ; n < N ; ++n ){
      loggamma_prop[ n ] = loggamma( n, j ) - phi( j, p ) * x( n, p ) + phi_prop * x( n, p );
    }
    
    double log_like_prop = 0; 
    for( int n = 0; n < N; ++n ){
      log_like_prop = log_like_prop - lgamma( exp( loggamma_prop[ n ] ) ) + exp( loggamma_prop[ n ] )*log( cc( n, j ) );
    }
    
    // Calculate ratio for single zeta
    double r;
    
    // If BB
      if( zeta( j, p ) == 1 ){
       r = log_like_prop + help::log_zeta_pj_cpp( 0, a, b) - ( log_like + help::log_phi_pj_cpp( phi( j, p ), sigma2_phi ) + help::log_zeta_pj_cpp( 1, a, b ) );
      }else{
       r = log_like_prop + help::log_phi_pj_cpp( phi_prop, sigma2_phi ) + help::log_zeta_pj_cpp( 1, a, b) - ( log_like + help::log_zeta_pj_cpp( 0, a, b ) );
      }
 
    
    // Calculate acceptance probability
    double aa  = log( runif( 1 )[ 0 ] );
    
    // Determine acceptance
    if( aa < r ){
      zeta( j, p ) = zeta_prop;
      phi( j, p ) = phi_prop;
      loggamma.col( j ) = loggamma_prop;
    }
    
  }
  
  // Return output
  List between( 3 );
  between[ 0 ] = zeta;
  between[ 1 ] = phi;
  between[ 2 ] = loggamma;
  return between;
}

// Function :: Between Step (jointly update phi and zeta) with add/delete using MRF prior
List between_phi_zeta_augment_update_MRF_cpp(
    arma::mat x,
    arma::mat zeta,
    arma::vec alpha,
    arma::mat phi, 
    arma::mat loggamma, 
    arma::mat cc,
    double sigma2_phi,
    double a_G,
    double b_G,
    arma::mat temp_G
){
  
  // Between Model Step
  // Initiate memory space
  int J = zeta.n_rows;
  int P = zeta.n_cols; 
  int N = x.n_rows;
  
  // Get a covariate p to update in branch j 
  
  int j = help::sample_cpp( help::myseq( 0, J - 1) ); 
  for( int p = 0; p < P; ++p ){
    
    // Calculate current log likelihood 
    double log_like = 0; 
    
    for( int n = 0; n < N; ++n ){
      log_like = log_like - lgamma( exp( loggamma( n , j ) ) ) + exp( loggamma( n , j ) )*log( cc( n, j ) );
    }
    
    // Calculate proposed likelihood 
    double phi_prop;
    int zeta_prop; 
    
    if( zeta( j, p ) == 1 ){
      phi_prop = 0.0;
      zeta_prop = 0;
    }else{
      phi_prop = help::phi_prop_cpp( phi( j, p ) ); 
      zeta_prop = 1;
    }
    
    arma::vec loggamma_prop( N );
    loggamma_prop.zeros();
    
    for( int n = 0 ; n < N ; ++n ){
      loggamma_prop[ n ] = loggamma( n, j ) - phi( j, p ) * x( n, p ) + phi_prop * x( n, p );
    }
    
    double log_like_prop = 0; 
    for( int n = 0; n < N; ++n ){
      log_like_prop = log_like_prop - lgamma( exp( loggamma_prop[ n ] ) ) + exp( loggamma_prop[ n ] )*log( cc( n, j ) );
    }
    
    // Calculate ratio for single zeta
    double r;
    
    // If MRF
    arma::mat zeta_prop_mat = zeta;
    zeta_prop_mat( j, p ) = zeta_prop;
    
    if( zeta( j, p ) == 1 ){
      r = log_like_prop + help::log_zeta_G_cpp( a_G, zeta_prop_mat, b_G, temp_G, j ) - ( log_like + help::log_phi_pj_cpp( phi( j, p ), sigma2_phi ) + help::log_zeta_G_cpp( a_G, zeta, b_G, temp_G, j ) );
    }else{
      r = log_like_prop + help::log_phi_pj_cpp( phi_prop, sigma2_phi ) + help::log_zeta_G_cpp( a_G, zeta_prop_mat, b_G, temp_G, j ) - ( log_like + help::log_zeta_G_cpp( a_G, zeta, b_G, temp_G, j ) );
    }
    
    // Calculate acceptance probability
    double aa  = log( runif( 1 )[ 0 ] );
    
    // Determine acceptance
    if( aa < r ){
      zeta( j, p ) = zeta_prop;
      phi( j, p ) = phi_prop;
      loggamma.col( j ) = loggamma_prop;
    }
    
  }
  
  // Return output
  List between( 3 );
  between[ 0 ] = zeta;
  between[ 1 ] = phi;
  between[ 2 ] = loggamma;
  return between;
}


// Function :: Between Step (jointly update phi and zeta) with add/delete or swap
List between_phi_zeta_swap_update_cpp(
    arma::mat x,
    arma::mat zeta,
    arma::vec alpha,
    arma::mat phi,
    arma::mat psi,
    double sigma2_phi,
    double a,
    double b
){
  
  // Between Model Step
  // Initiate memory space
  int J = zeta.n_rows;
  int P = zeta.n_cols;
  int dim = J*P;
  IntegerVector included( 0 );
  IntegerVector excluded( 0 );
  
  // Set proposal to current
  arma::mat zeta_proposal( J, P );
  arma::mat phi_proposal( J, P );
  zeta_proposal = zeta;
  phi_proposal = phi;
  
  
  // Get which included and excluded
  for( int jp = 0; jp < dim; ++jp ) {
    if( zeta[ jp ] == 1 ){
      included.push_back( jp );
    }else{
      excluded.push_back( jp );
    }
  }
  
  // Add/Delete or swap
  double swap = rbinom( 1, 1, 0.5 )[ 0 ];
  IntegerVector add_vcp;
  IntegerVector delete_vcp;
  
  // Protect against missingness
  if( ( included.size() < 1 ) | ( excluded.size() < 1 ) ){
    swap = 0;
  }
  
  if( swap == 1 ){
    // Swap
    // Select covariate to include
    if( excluded.size() == 1 ){
      add_vcp = excluded;
    }else{
      add_vcp =  help::sample_cpp( excluded );
    }
    
    // Get a covariate to delete
    if( included.size() == 1 ){
      delete_vcp = included;
    }else{
      delete_vcp =  help::sample_cpp( included );
    }
    
    // Update proposal
    // For added
    zeta_proposal[ add_vcp[ 0 ] ] = 1;
    double current_phi = phi[ add_vcp[ 0 ] ];
    double proposed_phi = help::phi_prop_cpp( current_phi );
    phi_proposal[ add_vcp[ 0 ] ] = proposed_phi;
    
    // For deleted
    zeta_proposal[ delete_vcp[ 0 ] ] = 0;
    double phi_removed =  phi_proposal[ delete_vcp[ 0 ] ];
    phi_proposal[ delete_vcp[ 0 ] ] = 0;
    
    // Calculate ratio for single zeta
    double r = help::log_dirichlet_psi_cpp( x, alpha,  phi_proposal, psi ) + help::log_phi_pj_cpp( proposed_phi, sigma2_phi )  - ( help::log_dirichlet_psi_cpp( x, alpha,  phi, psi ) + help::log_phi_pj_cpp( phi_removed, sigma2_phi ) );
    
    // Calculate acceptance probability
    double a  = log( runif( 1 )[ 0 ] );
    
    // Determine acceptance
    if( a < r ){
      zeta[ add_vcp[ 0 ] ] = zeta_proposal[ add_vcp[ 0 ] ];
      phi[ add_vcp[ 0 ] ] = phi_proposal[ add_vcp[ 0 ] ];
      zeta[ delete_vcp[ 0 ] ] = zeta_proposal[ delete_vcp[ 0 ] ];
      phi[ delete_vcp[ 0 ] ] = phi_proposal[ delete_vcp[ 0 ] ];
    }
    
  }else{
    
    // Add or delete
    double add = rbinom( 1, 1, 0.5 )[ 0 ];
    IntegerVector add_vcp;
    IntegerVector delete_vcp;
    
    // Protect against missingness
    if( ( included.size() < 1 ) & ( add == 0 ) ){
      add = 1;
    }
    if( ( excluded.size() < 1) & (add == 1)){
      add = 0;
    }
    
    if( add == 1 ){
      // Select covariate to include
      if( excluded.size() == 1 ){
        add_vcp = excluded;
      }else{
        add_vcp =  help::sample_cpp( excluded );
      }
      
      // Update proposal
      zeta_proposal[ add_vcp[ 0 ] ] = 1;
      double current_phi = phi[ add_vcp[ 0 ] ];
      double phi_proposed = help::phi_prop_cpp( current_phi );
      phi_proposal[ add_vcp[ 0 ] ] = phi_proposed;
      
      // Calculate ratio for single zeta
      double r = help::log_dirichlet_psi_cpp( x, alpha,  phi_proposal, psi ) + help::log_phi_pj_cpp( phi_proposed, sigma2_phi ) + help::log_zeta_pj_cpp( 1, a, b) - ( help::log_dirichlet_psi_cpp( x, alpha,  phi, psi ) + help::log_zeta_pj_cpp( 0, a, b ) );
      
      // Calculate acceptance probability
      double a  = log( runif( 1 )[ 0 ] );
      
      // Determine acceptance
      if( a < r ){
        zeta[ add_vcp[ 0 ] ] = zeta_proposal[ add_vcp[ 0 ] ];
        phi[ add_vcp[ 0 ] ] = phi_proposal[ add_vcp[ 0 ] ];
      }
    }else{
      // Get a covariate to delete
      if( included.size() == 1 ){
        delete_vcp = included;
      }else{
        delete_vcp =  sample_cpp( included );
      }
      
      // Update proposal
      zeta_proposal[ delete_vcp[ 0 ] ] = 0;
      double phi_removed = phi_proposal[ delete_vcp[ 0 ] ];
      phi_proposal[ delete_vcp[ 0 ] ] = 0;
      
      // Calculate ratio for zeta single
      double r = help::log_dirichlet_psi_cpp( x, alpha,  phi_proposal, psi ) + help::log_zeta_pj_cpp( 0, a, b ) - ( help::log_dirichlet_psi_cpp( x, alpha,  phi, psi ) + help::log_phi_pj_cpp( phi_removed, sigma2_phi ) + help::log_zeta_pj_cpp( 1, a, b ) );
      
      // Calculate acceptance probability
      double a  = log( runif( 1 )[ 0 ] );
      
      // Determine acceptance
      if( a < r ){
        zeta[ delete_vcp[ 0 ] ] = zeta_proposal[ delete_vcp[ 0 ] ];
        phi[ delete_vcp[ 0 ] ] = phi_proposal[ delete_vcp[ 0 ] ];
      }
      
    }
  }
  
  // Return output
  List between( 2 );
  between[ 0 ] = zeta;
  between[ 1 ] = phi;
  return between;
}

// Function :: Within Step ( update phi )
arma::mat within_phi_update_cpp(
    arma::mat x,
    arma::mat zeta,
    arma::vec alpha,
    arma::mat phi, 
    arma::mat psi,
    double sigma2_phi
){
  
  // Initiate memory space
  int J = zeta.n_rows;
  int P = zeta.n_cols;
  int dim = J*P;
  arma::mat phi_proposal( J, P );
  
  // Get which covariates are included
  for( int jp = 0; jp < dim; ++jp ) {
    phi_proposal = phi;
    
    if( zeta[ jp ] == 1 ){
      
      // Propose phi
      double current_phi = phi[ jp ];
      double phi_proposed = help::phi_prop_cpp( current_phi );
      phi_proposal[ jp ] = phi_proposed;
      
      // Calculate ratio
      double r = help::log_dirichlet_psi_cpp( x, alpha,  phi_proposal, psi ) + help::log_phi_pj_cpp( phi_proposed, sigma2_phi) - ( help::log_dirichlet_psi_cpp( x, alpha,  phi, psi ) + help::log_phi_pj_cpp( current_phi, sigma2_phi )  );
      
      // Calculate acceptance probability
      double a  = log( runif( 1 )[ 0 ] );
      
      // Determine acceptance
      if(a < r){
        phi[ jp ] = phi_proposal[ jp ];
      }
    }
  }
  
  // Return output
  return phi;
}

// Function :: Within Step ( update phi )
List within_phi_augment_cpp(
    arma::mat x,
    arma::mat zeta,
    arma::vec alpha,
    arma::mat phi, 
    arma::mat loggamma,
    arma::mat cc,
    double sigma2_phi
){
  
  // Initiate memory space
  int J = zeta.n_rows;
  int P = zeta.n_cols;
  int N = x.n_rows; 
  arma::mat phi_proposal( J, P );
  
  // Get which covariates are included
  for( int j = 0; j < J; ++j ){
    for( int p = 0; p < P; ++p ){
      if( zeta( j , p ) == 1 ){
        
        // Calculate current log likelihood 
        double log_like = 0; 
        
        for( int n = 0; n < N; ++n ){
          log_like = log_like - lgamma( exp( loggamma( n , j ) ) ) + exp( loggamma( n , j ) )*log( cc( n, j ) );
        }
        
        // Calculate proposed likelihood 
        double phi_prop = help::phi_prop_cpp( phi( j, p ) );
        
        arma::vec loggamma_prop( N );
        loggamma_prop.zeros();
        
        for( int n = 0 ; n < N ; ++n ){
          loggamma_prop[ n ] = loggamma( n, j ) - phi( j, p ) * x( n, p ) + phi_prop * x( n, p );
        }
        
        double log_like_prop = 0; 
        for( int n = 0; n < N; ++n ){
          log_like_prop = log_like_prop - lgamma( exp( loggamma_prop[ n ] ) ) + exp( loggamma_prop[ n ] )*log( cc( n, j ) );
        }
        
        // Calculate ratio
        double r = log_like_prop + help::log_phi_pj_cpp( phi_prop, sigma2_phi) - ( log_like + help::log_phi_pj_cpp( phi( j, p ), sigma2_phi )  );
        
        // Calculate acceptance probability
        double a  = log( runif( 1 )[ 0 ] );
        
        // Determine acceptance
        if(a < r){
          phi( j, p ) = phi_prop;
          loggamma.col( j ) = loggamma_prop;
        }
      }
    }
  }
  
  // Return output
  List phi_out( 2 );
  phi_out[ 0 ] = phi;
  phi_out[ 1 ] = loggamma;
  return phi_out;
}

// Function :: Update psi
arma::mat psi_update_cpp(
    arma::vec y,
    arma::mat Z,
    arma::mat psi,
    arma::vec alpha,
    arma::mat x,
    arma::mat phi,
    arma::vec xi,
    double delta,
    double rho,
    double a_0,
    double b_0,
    double h_alpha,
    double h_beta
){
  // Initiate memory space
  
  
  int N = psi.n_rows;
  int J = psi.n_cols;
  int M = xi.size();
  arma::mat I( N, N );
  I.zeros();
  arma::vec I_vec( N );
  I_vec.ones();
  arma::vec mu( N );
  mu.zeros();
  
  arma::mat psi_proposal( N, J );
  
  for( int i = 0; i < N; ++i ){
    I( i, i ) = 1;
  }
  
  // Update every psi
  for( int i = 0; i < N; ++i ){
    
    psi_proposal = psi;
    
    // Propose a psi
    arma::vec current_psi = psi.row( i ).t();
    arma::vec psi_proposal_get =  help::rdirichlet_cpp( current_psi, delta, rho );
    
    // Make sure that psi_proposal_get does not have any zeros
    for( int j = 0; j < J; ++j ){
      if( psi_proposal_get[ j ] == 0 ){
        psi_proposal_get[ j ] = 0.5/7500;
      }
    }
    // Make proposal sum to 1
    psi_proposal_get = psi_proposal_get/sum( psi_proposal_get );
    psi_proposal.row( i ) = psi_proposal_get.t();
    
    // Transform psi into balances
    arma::mat current_balance = help::ilr_fun_cpp( psi );
    arma::mat proposal_balance = help::ilr_fun_cpp( psi_proposal );
    
    // Get covariates that are selected
    IntegerVector which_balance( 0 );
    for( int m = 0; m < M; ++m ){
      if( xi[ m ] == 1 ){
        which_balance.push_back( m );
      }
    }
    
    arma::mat current_sigma( N , N);
    arma::mat proposal_sigma( N , N);
    
    if( sum( xi ) != 0){
      // Remove balances that are currently not selected
      uvec which = as<uvec>( which_balance );
      arma::mat current_balance_reduced = help::scale_cpp( current_balance.cols( which ) ) ;
      arma::mat proposal_balance_reduced = help::scale_cpp( proposal_balance.cols( which ) ) ;
      current_sigma = ( b_0/a_0 )*( I + h_alpha*I_vec*I_vec.t() + h_beta*current_balance_reduced*current_balance_reduced.t() );
      proposal_sigma = ( b_0/a_0 )*( I + h_alpha*I_vec*I_vec.t() + h_beta*proposal_balance_reduced*proposal_balance_reduced.t() );
      
    }else{
      current_sigma = ( b_0/a_0 )*( I + h_alpha*I_vec*I_vec.t() );
      proposal_sigma = ( b_0/a_0 )*( I + h_alpha*I_vec*I_vec.t() );
    }
    
    //Calculate ratio
    double r = help::log_mvtT_cpp( proposal_sigma, mu, 2*a_0, y) + help::log_multinomial_cpp( Z, psi_proposal ) + help::log_dirichlet_psi_i_cpp( x, alpha, phi, psi_proposal_get, i ) - ( help::log_mvtT_cpp( current_sigma, mu, 2*a_0, y) + help::log_multinomial_cpp( Z, psi ) + help::log_dirichlet_psi_i_cpp( x, alpha, phi, current_psi, i ) ) + help::log_dirichlet_psi_proposal_cpp( delta, rho, current_psi, psi_proposal_get );
    
    // Calculate acceptance probability
    double a  = log( runif( 1 )[ 0 ] );
    
    // Determine acceptance
    if(a < r){
      psi.row( i ) = psi_proposal.row( i );
    }
  }
  
  // Return output
  return psi;
}

// Function :: Update psi augment
arma::mat psi_augment_cpp(
    arma::mat cc
){
  // Initiate memory space
  int N = cc.n_rows;
  int J = cc.n_cols;
  arma::mat psi( N, J );
  psi.zeros();
  
  arma::vec TT( N );
  TT.zeros();
  arma::vec TT_norm(N);
  TT_norm.zeros();
  
  for( int n = 0; n < N; ++n ){
    for( int j = 0; j < J; ++j ){
      TT[ n ] += cc( n, j );
    }
    psi.row( n ) = cc.row( n )/TT[ n ];
    
    // Control for really small psi
    for( int j = 0; j < J; ++j ){
      if( psi( n, j ) < 0.5/7500){
        psi( n, j ) = 0.5/7500;
      }
      TT_norm[ n ] += psi( n, j );
    }    
    psi.row( n ) = psi.row( n )/TT_norm[ n ];
  }
  
  // Return output
  return psi;
}

// Update latent variable cc
arma::mat cc_update_cpp( arma::mat z, arma::mat loggamma, arma::vec uu ){
  int N = z.n_rows;
  int J = z.n_cols;
  
  arma::mat cc( N, J );
  cc.zeros();
  
  for( int n = 0; n < N; ++n ){
    for( int j = 0; j < J; ++j ){
      cc( n, j ) = rgamma( 1, z( n, j ) + exp(loggamma( n, j ) ), 1/( uu[ n ] + 1 )  )[ 0 ];
      if( cc( n, j ) < pow(10.0, -100.0)){
        cc( n, j ) = pow(10.0, -100.0);
      }
    }
  } 
  
  return cc;
}

// Update latent variable uu 
arma::vec uu_update_cpp( arma::mat z, arma::mat cc ){
  int N = z.n_rows; 
  
  arma::vec uu( N );
  uu.zeros();
  
  for( int n = 0; n < N; ++n ){
    double z_dot_n = sum( z.row( n ) );
    double T_n = sum( cc.row( n ) );
    
    uu[ n ] = rgamma( 1, z_dot_n, 1/T_n )[ 0 ]; 
  }
  
  return uu; 
}

// Update xi for selection of balances
arma::vec xi_update_cpp(
    arma::vec y,
    arma::mat psi,
    arma::vec xi,
    double a_m,
    double b_m,
    double a_0,
    double b_0,
    double h_alpha,
    double h_beta
){
  // Initiate memory space
  int N = psi.n_rows;
  int M = xi.size();
  arma::mat I( N, N );
  I.zeros();
  arma::vec I_vec( N );
  I_vec.ones();
  arma::vec mu( N );
  mu.zeros();
  for( int i = 0; i < N; ++i ){
    I( i, i ) = 1;
  }
  
  // Set proposal to current
  arma::vec xi_proposal( M );
  xi_proposal = xi;
  
  // Select a covariate
  int m = help::sample_cpp( help::myseq( 0, M - 1) ); 
  
  
  if( xi[ m ] == 1 ){
    xi_proposal[ m ] = 0;
  }else{
    xi_proposal[ m ] = 1;
  }
  
  // Transform psi into balances
  arma::mat current_balance = help::ilr_fun_cpp( psi );
  
  // Get covariates that are selected for proposal
  IntegerVector which_balance_proposal( 0 );
  for( int m = 0; m < M; ++m ){
    if( xi_proposal[ m ] == 1 ){
      which_balance_proposal.push_back( m );
    }
  }
  
  // Get covariates that are selected for current
  IntegerVector which_balance_current( 0 );
  for( int m = 0; m < M; ++m ){
    if( xi[ m ] == 1 ){
      which_balance_current.push_back( m );
    }
  }
  
  // Set sigma ( Make sure current is not empty )
  arma::mat current_sigma( N , N);
  if( sum( xi ) != 0){
    // Remove balances that are currently not selected
    uvec which_current = as<uvec>( which_balance_current );
    arma::mat current_balance_reduced = help::scale_cpp( current_balance.cols( which_current ) );
    current_sigma = ( b_0/a_0 )*( I + h_alpha*I_vec*I_vec.t() + h_beta*current_balance_reduced*current_balance_reduced.t() );
  }else{
    current_sigma = ( b_0/a_0 )*( I + h_alpha*I_vec*I_vec.t() );
  }
  
  // Set sigma ( Make sure proposal is not empty )
  arma::mat proposal_sigma( N , N);
  if( sum( xi_proposal ) != 0){
    // Remove balances that are currently not selected
    uvec which_proposal = as<uvec>(which_balance_proposal);
    arma::mat proposal_balance_reduced = help::scale_cpp( current_balance.cols( which_proposal ) );
    proposal_sigma = ( b_0/a_0 )*( I + h_alpha*I_vec*I_vec.t() + h_beta*proposal_balance_reduced*proposal_balance_reduced.t() );
  }else{
    proposal_sigma = ( b_0/a_0 )*( I + h_alpha*I_vec*I_vec.t() );
  }
  
  double r = 0;
  
  if( xi[ m ] == 1 ){
    r = help::log_mvtT_cpp( proposal_sigma, mu, 2*a_0, y) + help::log_xi_m_cpp( 0 , a_m, b_m)  - ( help::log_mvtT_cpp( current_sigma, mu, 2*a_0, y) + help::log_xi_m_cpp( 1 , a_m, b_m) );
  }else{
    r = help::log_mvtT_cpp( proposal_sigma, mu, 2*a_0, y) + help::log_xi_m_cpp( 1 , a_m, b_m)  - ( help::log_mvtT_cpp( current_sigma, mu, 2*a_0, y) + help::log_xi_m_cpp( 0 , a_m, b_m) );
  }
  
  // Calculate acceptance probability
  double a  = log( runif( 1 )[ 0 ] );
  
  // Determine acceptance
  if( a < r ){
    xi[ m ] = xi_proposal[ m ];
  }
  
  
  // Return output
  return xi;
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
    for( int j = i; j < p; ++j ){
      if( i == j ){
        newG( i, j ) = 0;
      }else{
        double N_out = R::dnorm( Omega( i, j ), 0, sd_v0, false );
        double N_in = R::dnorm( Omega( i, j ), 0, sd_v1, false  );
        double prob = ( N_in*pie )/( N_in*pie + N_out*( 1 - pie ) );
        newG( i, j) = rbinom( 1, 1, prob )[ 0 ];
        newG( j, i) = newG( i, j);
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
    double r = help::log_like_cpp( x, branch_loc, node_children_pointer,  branch_counts, subtree_counts, loggamma_proposal ) + help::log_alpha_cpp( alpha_proposal[ b ], sigma2_alpha ) - ( help::log_like_cpp( x, branch_loc, node_children_pointer, branch_counts, subtree_counts, loggamma ) + help::log_alpha_cpp( current_alpha, sigma2_alpha ) );
    
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

// Function :: Between Step (jointly update phi and zeta) with add/delete or swap with BB DTM
List between_phi_zeta_swap_update_BB_tree_cpp(
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
        add_vcp =  help::sample_cpp( excluded );
      }else{
        add_vcp = excluded;
      }
      
      // Get a covariate to delete
      if( included.size() > 1 ){
        delete_vcp =  help::sample_cpp( included );
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

// Function :: Between Step (jointly update phi and zeta) with add/delete or swap with MRF DTM
List between_phi_zeta_swap_update_MRF_tree_cpp(
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



// Function :: Within Step ( update phi for zeta == 1)
List within_phi_update_tree_cpp(
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

} // For namespace 'help'

// Function :: MCMC algorithm
// Make sure arguments are indexed by iteration
// [[Rcpp::export]]
List dm_lm_bvs(
    int iterations,
    int thin,
    arma::mat alpha,
    arma::vec y,
    arma::mat z,
    arma::mat x,
    arma::cube phi,
    arma::cube psi,
    arma::mat temp_cc,
    arma::vec temp_uu, 
    double sigma2_alpha,
    arma::cube zeta,
    arma::mat xi,
    double sigma2_phi,
    double a, // Beta-Bin
    double b, // Beta-Bin
    double a_0, // Inverse-Gamma
    double b_0, // Inverse-Gamma
    double h_alpha,
    double h_beta, 
    double a_m,
    double b_m
){
  
  // Initiate memory
  List between_phi_zeta( 3 );
  List alpha_out( 2 );
  List phi_out( 2 );
  int B = alpha.n_rows;
  int P = zeta.n_cols;
  int N = y.n_rows;  
  arma::vec temp_alpha( B );
  arma::mat temp_zeta( B, P );
  arma::mat temp_phi( B, P );
  arma::mat temp_psi( N, B );
  arma::vec temp_xi( B - 1 );
  
  temp_alpha = alpha.col( 0 );
  temp_zeta = zeta.slice( 0 );
  temp_phi = phi.slice( 0 );
  temp_psi = psi.slice( 0 );
  temp_xi = xi.col( 0 );
  
  arma::mat temp_loggamma( N, B );
  temp_loggamma.zeros();
  for( int n = 0; n < N; ++n ){
    for( int j = 0; j < B; ++j ){
      temp_loggamma( n, j ) = temp_alpha[ j ] + ( x.row( n )*temp_phi.row( j ).t() ).eval()( 0, 0 );
    }
  }
  
  // Looping over the number of iterations specified by user
  for( int iter = 0; iter < iterations; ++iter ){
    
    // Add/Delete
    between_phi_zeta = help::between_phi_zeta_augment_update_BB_cpp( x, temp_zeta, temp_alpha, temp_phi, temp_loggamma, temp_cc, sigma2_phi, a, b );
    temp_zeta = as<arma::mat>( between_phi_zeta[ 0 ] );
    temp_phi = as<arma::mat>( between_phi_zeta[ 1 ] );
    temp_loggamma = as<arma::mat>( between_phi_zeta[ 2 ] );
    
    // Update within
    phi_out = help::within_phi_augment_cpp( x, temp_zeta, temp_alpha, temp_phi, temp_loggamma, temp_cc, sigma2_phi );
    temp_phi =  as<arma::mat>( phi_out[ 0 ] );
    temp_loggamma = as<arma::mat>( phi_out[ 1 ] );
    
    // Update alphas
    alpha_out = help::alpha_update_augment_cpp( temp_loggamma, temp_cc, temp_alpha, temp_phi, sigma2_alpha );
    temp_alpha = as<arma::vec>( alpha_out[ 0 ] );
    temp_loggamma = as<arma::mat>( alpha_out[ 1 ] );
    
    // Update cc 
    temp_cc = help::cc_update_cpp( z, temp_loggamma, temp_uu );
    
    // Update uu
    temp_uu = help::uu_update_cpp( z, temp_cc );
    
    // Update psi
    temp_psi = help::psi_augment_cpp( temp_cc );

    // Update xi
    temp_xi = help::xi_update_cpp( y, temp_psi, temp_xi, a_m, b_m, a_0, b_0, h_alpha, h_beta);
    

    
    // Set the starting values for the next iteration
    if( (iter + 1) % thin == 0 ){
      alpha.col( (iter + 1)/thin - 1 ) = temp_alpha;
      zeta.slice( (iter + 1)/thin - 1 ) = temp_zeta;
      phi.slice( (iter + 1)/thin - 1 ) = temp_phi;
      psi.slice( (iter + 1)/thin - 1 ) = temp_psi;
      xi.col( (iter + 1)/thin - 1 ) = temp_xi;
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
  output[ 3 ] = psi;
  output[ 4 ] = xi;
  
  return output ;
}

// Function :: MCMC algorithm
// Make sure arguments are indexed by iteration
// [[Rcpp::export]]
List dm_bvs(
    int iterations,
    int thin,
    arma::mat alpha, 
    arma::mat z,
    arma::mat x,
    arma::cube phi, 
    arma::mat temp_cc,
    arma::vec temp_uu, 
    double sigma2_alpha,
    arma::cube zeta, 
    double sigma2_phi,
    String prior,
    double a, // Beta-Bin
    double b, // Beta-Bin 
    arma::cube Omega, // MRF
    arma::cube G, // MRF
    arma::cube Var, // MRF
    arma::mat S, // MRF
    double v0, // MRF
    double v1, // MRF
    double a_G, // MRF
    double b_G, // MRF
    double pie, // MRF
    double lambda // MRF
){
  
  // Initiate memory
  List between_phi_zeta( 3 );
  List alpha_out( 2 );
  List phi_out( 2 );
  int B = alpha.n_rows;
  int P = zeta.n_cols;
  int N = x.n_rows;  
  arma::vec temp_alpha( B );
  arma::mat temp_zeta( B, P );
  arma::mat temp_phi( B, P ); 
  
  temp_alpha = alpha.col( 0 );
  temp_zeta = zeta.slice( 0 );
  temp_phi = phi.slice( 0 ); 
  
  arma::mat temp_Omega = Omega.slice( 0 );
  arma::mat temp_Var = Var.slice( 0 );
  arma::mat temp_G = G.slice( 0 );
  
  arma::mat temp_loggamma( N, B );
  temp_loggamma.zeros();
  for( int n = 0; n < N; ++n ){
    for( int j = 0; j < B; ++j ){
      temp_loggamma( n, j ) = temp_alpha[ j ] + ( x.row( n )*temp_phi.row( j ).t() ).eval()( 0, 0 );
    }
  }
  
  // Looping over the number of iterations specified by user
  for( int iter = 0; iter < iterations; ++iter ){
    
    // Add/Delete
    if( prior == "BB" ){
      between_phi_zeta = help::between_phi_zeta_augment_update_BB_cpp( x, temp_zeta, temp_alpha, temp_phi, temp_loggamma, temp_cc, sigma2_phi, a, b );
    }

    if( prior == "MRF_fixed" | prior == "MRF_unknown" ){
      between_phi_zeta = help::between_phi_zeta_augment_update_MRF_cpp( x, temp_zeta, temp_alpha, temp_phi, temp_loggamma, temp_cc, sigma2_phi, a_G, b_G, temp_G );
    }

    temp_zeta = as<arma::mat>( between_phi_zeta[ 0 ] );
    temp_phi = as<arma::mat>( between_phi_zeta[ 1 ] );
    temp_loggamma = as<arma::mat>( between_phi_zeta[ 2 ] );

    // Update within
    phi_out = help::within_phi_augment_cpp( x, temp_zeta, temp_alpha, temp_phi, temp_loggamma, temp_cc, sigma2_phi );
    temp_phi =  as<arma::mat>( phi_out[ 0 ] );
    temp_loggamma = as<arma::mat>( phi_out[ 1 ] );

    // Update alphas
    alpha_out = help::alpha_update_augment_cpp( temp_loggamma, temp_cc, temp_alpha, temp_phi, sigma2_alpha );
    temp_alpha = as<arma::vec>( alpha_out[ 0 ] );
    temp_loggamma = as<arma::mat>( alpha_out[ 1 ] );

    // Update cc
    temp_cc = help::cc_update_cpp( z, temp_loggamma, temp_uu );

    // Update uu
    temp_uu = help::uu_update_cpp( z, temp_cc );
    
    // Update Gaussian graphical model if unknown
    if( prior == "MRF_unknown" ){
      temp_Omega = help::Omega_update_cpp( temp_Omega, S, temp_Var, lambda, N );
      temp_G = help::update_G_cpp( temp_Omega, v0, v1, pie );
      temp_Var = help::set_V_cpp( temp_G, v0, v1 );
    }
    
    int jk = help::sample_cpp( help::myseq( 0, 10) );
    for( int m = 0; m < jk; ++m){
      double ok = rbinom( 1, 1, 0.5 )[ 0 ];
      ok = ok + 1;
    }
    
    // Set the starting values for the next iteration
    if( (iter + 1) % thin == 0 ){
      alpha.col( (iter + 1)/thin - 1 ) = temp_alpha;
      zeta.slice( (iter + 1)/thin - 1 ) = temp_zeta;
      phi.slice( (iter + 1)/thin - 1 ) = temp_phi; 
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

// Function :: MCMC algorithm
// Make sure arguments are indexed by iteration
// [[Rcpp::export]]
List wang(
    int iterations,
    int thin,
    arma::mat x,
    arma::cube Omega, // MRF
    arma::cube G, // MRF
    arma::cube Var, // MRF
    arma::mat S, // MRF
    double v0, // MRF
    double v1, // MRF
    double pie, // MRF
    double lambda // MRF
){
  
  // Initiate memory
  int N = x.n_rows;  
  
  arma::mat temp_Omega = Omega.slice( 0 );
  arma::mat temp_Var = Var.slice( 0 );
  arma::mat temp_G = G.slice( 0 );

  // Looping over the number of iterations specified by user
  for( int iter = 0; iter < iterations; ++iter ){
    
      temp_Omega = help::Omega_update_cpp( temp_Omega, S, temp_Var, lambda, N );
      temp_G = help::update_G_cpp( temp_Omega, v0, v1, pie );
      temp_Var = help::set_V_cpp( temp_G, v0, v1 );
  
    
    // Set the starting values for the next iteration
    if( (iter + 1) % thin == 0 ){
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
  List output( 2 );
  output[ 0 ] = Omega;
  output[ 1 ] = G;
  
  return output ;
}

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
      between_phi_zeta = help::between_phi_zeta_swap_update_BB_tree_cpp( x, temp_zeta, temp_alpha, temp_phi, branch_location, node_children_pointer, branch_counts, subtree_counts, aa, bb, sigma2_phi, loggamma );
    }
    
    if( prior == "MRF_fixed" ){
      between_phi_zeta = help::between_phi_zeta_swap_update_MRF_tree_cpp( x, temp_zeta, temp_alpha, temp_phi, branch_location, node_children_pointer, branch_counts, subtree_counts, sigma2_phi, a_G, b_G, temp_G, loggamma );
    }
    
    if( prior == "MRF_unknown" ){
      between_phi_zeta = help::between_phi_zeta_swap_update_MRF_tree_cpp( x, temp_zeta, temp_alpha, temp_phi, branch_location, node_children_pointer, branch_counts, subtree_counts, sigma2_phi, a_G, b_G, temp_G, loggamma );
    }
    
    temp_zeta = as<arma::mat>( between_phi_zeta[ 0 ] );
    temp_phi = as<arma::mat>( between_phi_zeta[ 1 ] );
    loggamma =  as<arma::mat>( between_phi_zeta[ 2 ] );
    
    // Update within
    phi_return = help::within_phi_update_tree_cpp(x, temp_zeta, temp_alpha, temp_phi, branch_location, node_children_pointer, branch_counts, subtree_counts, sigma2_phi, loggamma );
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

////////////////////////////////////////////////////////////////////////////////////////////////////////
// Perform isometric log-ratio transformation on percentages (ie construct balances)
// [[Rcpp::export]]
arma::mat ilr_fun_cpp( arma::mat x ){
  int n = x.n_rows;
  int p = x.n_cols;
  arma::mat y = log( x );
  arma::vec y_centered( n );
  
  // Center y
  for( int i = 0; i < n; ++i){
    arma::vec rows = y.row( i ).t();
    double means = mean( rows );
    y_centered( i ) = means;
  }
  
  for( int j = 0; j < p; ++j ){
    y.col( j ) = y.col(j) - y_centered;
  }
  
  arma::mat H = help::helmert_cpp( p );
  
  for( int k = 0; k < ( p - 1 ); ++k ){
    H.col( k ) = H.col( k )/sqrt( (k+1)*(k+2) );
  }
  
  arma::mat ilr = y*H;
  
  return ilr;
}
// Calculate log posterior distribution
// [[Rcpp::export]]
arma::vec log_post(
    int iterations,
    arma::mat alpha,
    arma::vec y,
    arma::mat z,
    arma::mat x,
    arma::cube phi,
    arma::cube psi,
    double sigma2_alpha,
    arma::cube zeta,
    arma::mat xi,
    double sigma2_phi,
    double a, // Beta-Bin
    double b, // Beta-Bin
    double a_0, // Inverse-Gamma
    double b_0, // Inverse-Gamma
    double h_alpha,
    double h_beta, 
    double a_m,
    double b_m
){
  
  int N = x.n_rows;
  int P = x.n_cols;
  int J = psi.n_cols;
  int M = xi.n_rows;
  arma::vec log_post( iterations );
  log_post.zeros();
  arma::mat I( N, N );
  I.zeros();
  arma::vec I_vec( N );
  I_vec.ones();
  
  for( int n = 0; n < N; ++n ){
    I( n, n ) = 1;
  }
  
  arma::vec mu( N );
  mu.zeros();
  
  // Get log_post for each iteration
  for( int i = 0; i < iterations; ++i){
    arma::mat balance = ilr_fun_cpp( psi.slice(i) );
    
    IntegerVector which_balance_current( 0 );
    vec curr_xi = xi.col( i );
    
    for( int m = 0; m < M; ++m ){
      if( curr_xi[ m ] == 1 ){
        which_balance_current.push_back( m );
      }
    }
    
    uvec which_current = as<uvec>( which_balance_current );
    arma::mat current_balance_reduced = help::scale_cpp( balance.cols( which_current ) );
    arma::mat current_sigma = ( b_0/a_0 )*( I + h_alpha*I_vec*I_vec.t() + h_beta*current_balance_reduced*current_balance_reduced.t() );
    
    log_post[ i ] +=  help::log_mvtT_cpp( current_sigma, mu, a_0*2, y) + help::log_multinomial_cpp( z, psi.slice(i) ) + help::log_dirichlet_psi_cpp( x, alpha.col( i ), phi.slice( i ), psi.slice( i ) );
    
    arma::mat cur_phi = phi.slice( i );
    arma::mat cur_zeta = zeta.slice( i );
    
    for( int j = 0; j < J; ++j){
      
      log_post[ i ] += help::log_alpha_cpp( alpha( j, i ), sigma2_alpha );
      
      for( int p = 0; p < P; ++p){
        log_post[ i ] += help::log_phi_pj_cpp( cur_phi( j, p ), sigma2_phi) + help::log_zeta_pj_cpp( cur_zeta( j, p), a, b);
        
      }
    }
    for( int m = 0; m < M; ++m){
      
      log_post[ i ] += help::log_xi_m_cpp( curr_xi[ m ] , a_m, b_m);
      
    }
  }
  return( log_post );
}

// Preps input for 'loo' package in R to obtain approximate leave-one-out cross-validation
// [[Rcpp::export]]
arma::mat loo_prep(
    arma::vec y,
    arma::cube psi,
    arma::mat xi,
    double a_0,
    double b_0,
    double h_alpha,
    double h_beta
){
  int M = xi.n_rows;
  int N = y.size();
  int iterations = xi.n_cols;
  
  arma::mat loo_mat( iterations, N );
  
  for( int iter = 0; iter < iterations; ++iter){
    arma::vec xi_col = xi.col( iter );
    
    arma::mat balance = ilr_fun_cpp( psi.slice( iter ) );
    
    IntegerVector which_balance_current( 0 );
    vec curr_xi = xi.col( iter );
    
    for( int m = 0; m < M; ++m ){
      if( curr_xi[ m ] == 1 ){
        which_balance_current.push_back( m );
      }
    }
    
    uvec which_current = as<uvec>( which_balance_current );
    arma::mat current_balance_reduced = help::scale_cpp( balance.cols( which_current ) );
    
    for( int n = 0; n < N; ++n){
      
      arma::mat current_balance_n = current_balance_reduced.row( n );
      arma::vec current_sigma = ( b_0/a_0 )*( 1 + h_alpha + h_beta*current_balance_n*current_balance_n.t() );
      double meanC = y[ n ];
      
      double log_T = lgamma( ( 2*a_0 + 1 )/2 ) - lgamma( a_0 ) - 0.5*log( 2*a_0*atan(1)*4 ) - log( sqrt( current_sigma[0] ) ) - ( ( 2*a_0 + 1 )/2 )*log( 1 + ( 1/( 2*a_0 ) )*pow( meanC, 2 )/current_sigma[ 0 ] );
      
      loo_mat( iter, n ) = log_T;
      
    }
  }
  return( loo_mat );
}

// [[Rcpp::export]]
double log_zeta_pj_cpp( double t_pj, double a, double b
){
  double post_a = t_pj + a;
  double post_b = 1 - t_pj + b;
  double log_zeta_pj = lgamma( post_a ) + lgamma( post_b ) - lgamma( post_a + post_b ) - ( lgamma( a ) + lgamma( b ) - lgamma( a + b ) );
  
  // Return output
  return log_zeta_pj ;
}

