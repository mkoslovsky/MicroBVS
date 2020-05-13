

out <- DTMbvs::DMbvs_R(Y = data$Y, X = data$X[,-1], MCMC = "SSVS", warmstart = F, iterations = 20000)
xit <- apply(out$beta[,1000:2000]!=0,c(1),mean)
View(matrix(xit,nrow = 40,ncol = 100, byrow = TRUE))
View(matrix((xit > 0.5)*1,nrow = 40,ncol = 100, byrow = TRUE))
View( 1*(data$betas[,-1] != 0) ) 

G <- diag(20)
G[1:5,1:5] <- 1
ok <- simulate_DM( n_vars = 20,  rho = 0.7, G = G )
out2 <- DMbvs_R(iterations = 10000, thin = 10, z = ok$Y, x = ok$X[,-1], alpha = NULL, phi = NULL,
                sigma2_alpha = sqrt(10), sigma2_phi = sqrt(10), prior = "MRF_unknown", G = matrix(0, ncol( ok$X[,-1] ), ncol( ok$X[,-1] )), a = 1, b = 4, seed = 1212 )
xit <- apply( out2$zeta[,,500:1000]!=0,c(1,2),mean )
View( ( xit > 0.5)*1 )
zeta_means <- apply( out2[[2]][, ,500:1000], c(1, 2), mean )
selected_zeta <- which( zeta_means >= 0.5, arr.ind = T )
this <- apply( out2$G[,,500:1000], c(1,2), mean)
View( (this > 0.5)*1)


# Test DM BB (GOOD TO GO)
sourceCpp("MicroBVStest.cpp")
data <- simulate_DM( rho = 0.3 )
set.seed(1)
out2 <- DMbvs_R(iterations = 1000, thin = 10, z = data$Y, x = data$X[,-1], alpha = NULL, phi = NULL,
                sigma2_alpha = sqrt(10), sigma2_phi = sqrt(10), prior = "BB", a = 1, b = 9, seed = 1212 ) # 1 minute to run 

sele <- selected_DM(out2, burnin = 500) # From selectedDM.R 

# Test DM MRF fixed (GOOD TO GO) Look at sensitivity to b_G
sourceCpp("MicroBVStest.cpp")
Sigma <- matrix(0,30,30)
Sigma[1:15,1:15] <- 0.7
diag(Sigma) <- 1
G <- (Sigma != 0 )*1
set.seed(1)
data <- simulate_DM( Sigma = Sigma )
out2 <- DMbvs_R(iterations = 10000, thin = 10, z = data$Y, x = data$X[,-1], alpha = NULL, phi = NULL,
                sigma2_alpha = sqrt(10), sigma2_phi = sqrt(10), prior = "MRF_fixed", a_G = log(0.1/0.9), b_G = 0.2, seed = 1212, G = G ) # 2 minute to run 

sele_2 <- selected_DM(out2, burnin = 500) # From selectedDM.R 

# Test DM MRF unknown (GOOD TO GO) Look at sensitivity to b_G
sourceCpp("MicroBVStest.cpp")
Sigma <- matrix(0,30,30)
Sigma[1:15,1:15] <- 0.7
diag(Sigma) <- 1
G <- (Sigma != 0 )*1
set.seed(1)
data <- simulate_DM( n_obs= 250, Sigma = Sigma )
start_time <- Sys.time()
out2 <- DMbvs_R(iterations = 10000, thin = 10, z = data$Y, x = data$X[,-1], alpha = NULL, phi = NULL,
                sigma2_alpha = sqrt(10), sigma2_phi = sqrt(10), prior = "MRF_unknown", a_G = log(0.1/0.9), b_G = 0.1, seed = 1212, v0 = 0.0001, v1 = 10 ) # 3.87 minute to run 
end_time <- Sys.time()
end_time - start_time
sele_3 <- selected_DM(out2, burnin = 500, G = T) # From selectedDM.R 
View( 1*( sele_3$estimated_G > 0.5 ) )
sele_3$selected_zeta

# Test DMLM BB (GOOD TO GO)
sourceCpp("MicroBVStest.cpp")
data <- simulate_DMLM(rho = 0.3)
set.seed(1)
out2 <- dm_lm_bvs_R( iterations = 10000, thin = 10, y = data$Y, z = data$Z, x = data$X, sigma2_alpha = sqrt( 10 ), sigma2_phi = sqrt( 10 ),
                     h_alpha = 1, h_beta = 1, a_m = 1, b_m = 9, a = 1, b = 9,
                     a_0 = 2, b_0 = 2, seed = 1, prior = "BB" )
sele <- selected_DMLM(out2, burnin = 500) # From selectedDM.R 

# Test DMLM MRF fixed (GOOD TO GO)
sourceCpp("MicroBVStest.cpp")
Sigma <- matrix(0,30,30)
Sigma[1:15,1:15] <- 0.7
diag(Sigma) <- 1
G <- (Sigma != 0 )*1
set.seed(1)
data <- simulate_DMLM( n_vars= 30,  Sigma = Sigma )
out2 <- dm_lm_bvs_R( iterations = 10000, thin = 10, y = data$Y, z = data$Z, x = data$X, sigma2_alpha = sqrt( 10 ), sigma2_phi = sqrt( 10 ),
                     h_alpha = 1, h_beta = 1, a_m = 1, b_m = 9, a = 1, b = 9,
                     a_0 = 2, b_0 = 2, seed = 1, prior = "MRF_fixed", G= G, a_G = log(0.1/0.9), b_G = 0.2 )
sele <- selected_DMLM(out2, burnin = 500) # From selectedDMLM.R 

# Test DMLM MRF unknown (GOOD TO GO)
sourceCpp("MicroBVStest.cpp")
Sigma <- matrix(0,30,30)
Sigma[1:15,1:15] <- 0.7
diag(Sigma) <- 1
G <- (Sigma != 0 )*1
set.seed(1)
data <- simulate_DMLM( n_vars= 30,  Sigma = Sigma )
out2 <- dm_lm_bvs_R( iterations = 10000, thin = 10, y = data$Y, z = data$Z, x = data$X, sigma2_alpha = sqrt( 10 ), sigma2_phi = sqrt( 10 ),
                     h_alpha = 1, h_beta = 1, a_m = 1, b_m = 9, a = 1, b = 9,
                     a_0 = 2, b_0 = 2, seed = 1,  prior = "MRF_unknown", a_G = log(0.1/0.9), b_G = 0.2, v0 = 0.0001, v1 = 10 )
sele <- selected_DMLM(out2, burnin = 500, G = T) # From selectedDMLM.R 

# Test DTM BB 





# Test Graphical Model 
n_obs <- 100
covariates_sim <- 30
iterations <- 10000
thin <- 10
samples <- floor( iterations/thin )
Sigma <- matrix( 0, covariates_sim, covariates_sim)
Sigma[ 1:15, 1:15] <- 0.7
diag(Sigma) <- 1 
v0 <- 0.0001
v1 <- 10
lambda <- 1
set.seed(211) 
X <- scale(MASS::mvrnorm(n = n_obs, mu = rep(0, covariates_sim), Sigma = Sigma))
Omega. <- array( 0, dim = c( covariates_sim, covariates_sim, samples ) )
Var. <- array( 0, dim = c( covariates_sim, covariates_sim, samples ) )
G. <- array( 0, dim = c( covariates_sim, covariates_sim, samples ) )

# Adjust inital values for graphical parameters if they are still NULL
pie <- .5 #2/( covariates_sim - 1) 

set_V <- function( G. = G, v0. = v0, v1. = v1 ){
  G_in <- G.
  G_out <- 1 - G.
  
  V <- v1.*G_in + v0.*G_out
  diag( V ) <- 0
  return(V)
}

Omega.[ , , 1 ] <- diag( covariates_sim )  
G.[ , , 1 ] <- 0  

Var.[ , , 1] <- set_V( G.[ , , 1 ], v0, v1 )
S. <- t( X ) %*% X

graph <- wang( iterations, thin, data$X[,-1], Omega., G., Var., S., v0, v1, pie, lambda )
View( 1*( 0.5 < apply( graph[[ 2 ]][ , , 500:1000],c( 1, 2), mean ) ))





