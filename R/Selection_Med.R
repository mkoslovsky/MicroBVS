Selection_Med1 <- function( model = NULL, burnin = 250, threshold = c(0.5, 0.5), quantile = c(0.025,0.975), 
                           covar_select = FALSE, interest_covar = 1){
  library(ggplot2)
  
  if(burnin%%1 != 0){
    stop("Bad input: burn-in should be an integer")
  }
  
  if( length( threshold ) != 2 ){
    stop("Bad input: threshold must be a two-dimensional vector")
  }
  
  if( threshold[ 1 ] > 1 | threshold[ 1 ] < 0 | threshold[ 2 ] > 1 | threshold[ 2 ] < 0){
    stop("Bad input: threshold should be between 0 and 1")
  }

  
  samples<- dim(model[[1]])[2]
  B_sim <- dim(model[[1]])[1]
  tot <- dim(model[[5]])[1]
  balance_structure <- matrix(0,nrow=B_sim,ncol=(B_sim-1))
  for(i in 1:B_sim){
    for(j in 1:(B_sim-1)){
      if(j < i){
        balance_structure[i,j] = -1/(B_sim-j)*sqrt((B_sim-i)/(B_sim+1-i))
      }
      if(j == i){
        balance_structure[i,j] = sqrt((B_sim-j)/(B_sim+1-j))
      }
    }
  }
  
  #test if the taxa in the first play of the roll is selected
  selected_taxa = mean(model[[2]][1, 1 , (burnin + 1):samples]) > threshold[1] & mean(model[[5]][(tot-B_sim+2), ( burnin + 1 ):samples]) > threshold[2]
  
  covar_DM = NULL
  covar_LM = NULL
  if(covar_select == TRUE & dim(model[[2]])[2] >= (interest_covar + 1)){
    covar_DM = mean(model[[2]][1, (interest_covar + 1), (burnin + 1):samples]) > threshold[1]
    covar_LM = mean(model[[5]][2+interest_covar, ( burnin + 1 ):samples]) > threshold[2]
  }
  if(covar_select == TRUE & dim(model[[2]])[2] < (interest_covar + 1)){
    stop("Bad input: the interest covariate is not included in the data")
  }
  
  #Direct Effect:  
  mean_direct <- mean(model[[6]][2, ( burnin + 1 ):samples])
  quantile_direct <- quantile(model[[6]][2, ( burnin + 1 ):samples],quantile)
  
  ##calculate the \gamma and digamma function
  temp_psi1 <- exp(model[[3]][,1,( burnin + 1 ):samples] + model[[1]][,( burnin + 1 ):samples])
  temp_psi0 <- exp(model[[1]][,( burnin + 1 ):samples])
  tot_psi1 <- digamma(apply(temp_psi1,2,sum))
  tot_psi0 <- digamma(apply(temp_psi0,2,sum))
  
  E_log <- matrix(ncol = B_sim, nrow = (samples - burnin))
  temp_balance <- matrix(ncol = (B_sim - 1), nrow = (samples - burnin))
  for(l in 1:B_sim){
    E_log[,l] = digamma(temp_psi1[l,]) - tot_psi1 - (digamma(temp_psi0[l,]) - tot_psi0)
  }
  
  #calculate E[balance_1(T=1)] - E[balance_1(T=0)]
  marginal_mediation <- array(0,dim=c(B_sim,(samples - burnin)))
  for(i in 1:(B_sim)){
    for(j in 1:(samples - burnin)){
      marginal_mediation[i,j] <- sum(balance_structure[i,]*model[[6]][(tot-B_sim+2):tot, (j+burnin)])*E_log[j,i]
    } 
  }
  
    global_mediation <- apply(marginal_mediation,2,sum)
    global_mean <- mean(global_mediation)
    global_quantile <- quantile(global_mediation,c(quantile[1],quantile[2]))
    marginal_mean <- mean(marginal_mediation[1,])
    marginal_quantile <- quantile(marginal_mediation[1,],c(quantile[1],quantile[2]))

  return(list(selected = selected_taxa, global_mean = global_mean, global_quantile = global_quantile, 
              marginal_mean = marginal_mean, marginal_quantile = marginal_quantile,
              direct_mean = mean_direct, direct_quantile = quantile_direct,
              covar_result = c(covar_DM, covar_LM)))
}

Selection_Med2 <- function( model = NULL, burnin = 250, threshold = c(0.5, 0.5), quantile = c(0.025,0.975), 
                            covar_select = FALSE, interest_covar = 1){
  library(ggplot2)
  
  if(burnin%%1 != 0){
    stop("Bad input: burn-in should be an integer")
  }
  
  if( length( threshold ) != 2 ){
    stop("Bad input: threshold must be a two-dimensional vector")
  }
  
  if( threshold[ 1 ] > 1 | threshold[ 1 ] < 0 | threshold[ 2 ] > 1 | threshold[ 2 ] < 0){
    stop("Bad input: threshold should be between 0 and 1")
  }
  
  
  samples<- dim(model[[1]])[2]
  B_sim <- dim(model[[1]])[1]
  tot <- dim(model[[5]])[1]
  
  balance_structure <- matrix(0,nrow=B_sim,ncol=(B_sim-1))
  for(i in 1:B_sim){
    for(j in 1:(B_sim-1)){
      if(j < i){
        balance_structure[i,j] = -1/(B_sim-j)*sqrt((B_sim-i)/(B_sim+1-i))
      }
      if(j == i){
        balance_structure[i,j] = sqrt((B_sim-j)/(B_sim+1-j))
      }
    }
  }
  
  for(j in 1:(B_sim-1)){
    balance_structure[B_sim,j] = -1/(B_sim-j)*sqrt((B_sim-j)/(B_sim+1-j))
  }
  
  covar_DM = NULL
  covar_LM = NULL
  if(covar_select == TRUE & dim(model[[2]])[2] >= (interest_covar + 1)){
    covar_DM = mean(model[[2]][1, (interest_covar + 1), (burnin + 1):samples]) > threshold[1]
    covar_LM = mean(model[[5]][2+interest_covar, ( burnin + 1 ):samples]) > threshold[2]
  }
  if(covar_select == TRUE & dim(model[[2]])[2] < (interest_covar + 1)){
    stop("Bad input: the interest covariate is not included in the data")
  }
  
  #Direct Effect:  
  mean_direct <- mean(model[[6]][2, ( burnin + 1 ):samples])
  quantile_direct <- quantile(model[[6]][2, ( burnin + 1 ):samples],quantile)
  
  ##calculate the \gamma and digamma function
  temp_psi1 <- exp(model[[3]][,1,( burnin + 1 ):samples] + model[[1]][,( burnin + 1 ):samples])
  temp_psi0 <- exp(model[[1]][,( burnin + 1 ):samples])
  tot_psi1 <- digamma(apply(temp_psi1,2,sum))
  tot_psi0 <- digamma(apply(temp_psi0,2,sum))
  
  E_log <- matrix(ncol = B_sim, nrow = (samples - burnin))
  temp_balance <- matrix(ncol = (B_sim - 1), nrow = (samples - burnin))
  for(l in 1:B_sim){
    E_log[,l] = digamma(temp_psi1[l,]) - tot_psi1 - (digamma(temp_psi0[l,]) - tot_psi0)
  }
  
  #calculate E[balance_1(T=1)] - E[balance_1(T=0)]
  marginal_mediation <- array(0,dim=c(B_sim,(samples - burnin)))
  for(i in 1:(B_sim)){
    for(j in 1:(samples - burnin)){
      marginal_mediation[i,j] <- sum(balance_structure[i,]*model[[6]][(tot-B_sim+2):tot, (j+burnin)])*E_log[j,i]
    } 
  }
  
  global_mediation <- apply(marginal_mediation,2,sum)
  global_mean <- mean(global_mediation)
  global_quantile <- quantile(global_mediation,c(quantile[1],quantile[2]))
  marginal_mean <- apply(marginal_mediation,1,mean)
  marginal_lowerquantile <- apply(marginal_mediation,1,quantile,c(quantile[1]))
  marginal_upperquantile <- apply(marginal_mediation,1,quantile,c(quantile[2]))
  selected_taxa <- marginal_lowerquantile < 0  & marginal_upperquantile > 0
  
  return(list(selected = !selected_taxa, global_mean = global_mean, global_quantile = global_quantile, 
              marginal_mean = marginal_mean, marginal_lowerquantile = marginal_lowerquantile, marginal_upperquantile = marginal_upperquantile,
              direct_mean = mean_direct, direct_quantile = quantile_direct,
              covar_result = c(covar_DM, covar_LM)))
}
Selection_Med3 <- function( model = NULL, burnin = 250, threshold = c(0.5, 0.5), quantile = c(0.025,0.975), 
                            covar_select = FALSE, interest_covar = 1){
  library(ggplot2)
  
  if(burnin%%1 != 0){
    stop("Bad input: burn-in should be an integer")
  }
  
  if( length( threshold ) != 2 ){
    stop("Bad input: threshold must be a two-dimensional vector")
  }
  
  if( threshold[ 1 ] > 1 | threshold[ 1 ] < 0 | threshold[ 2 ] > 1 | threshold[ 2 ] < 0){
    stop("Bad input: threshold should be between 0 and 1")
  }
  
  
  samples<- dim(model[[1]])[2]
  B_sim <- dim(model[[1]])[1]
  tot <- dim(model[[5]])[1]

  mediators <- 1:B_sim
  est_zeta_temp <- apply( model[[2]][,,251:500], c(1), mean )
  potential_mediator <- mediators[est_zeta_temp >=0.5]
  
  covar_DM = NULL
  covar_LM = NULL
  if(covar_select == TRUE & dim(model[[2]])[2] >= (interest_covar + 1)){
    covar_DM = mean(model[[2]][1, (interest_covar + 1), (burnin + 1):samples]) > threshold[1]
    covar_LM = mean(model[[5]][2+interest_covar, ( burnin + 1 ):samples]) > threshold[2]
  }
  if(covar_select == TRUE & dim(model[[2]])[2] < (interest_covar + 1)){
    stop("Bad input: the interest covariate is not included in the data")
  }
  
  mean_direct <- mean(model[[6]][2, ( burnin + 1 ):samples])
  quantile_direct <- quantile(model[[6]][2, ( burnin + 1 ):samples],quantile)
  selected_taxa <- est_zeta_temp >=0.5
  
  if(length(potential_mediator) == 0){
    return(list(direct_mean = mean_direct, direct_quantile = quantile_direct,
                covar_result = c(covar_DM, covar_LM)), message = "No mediator found")
  }
  
  if(length(potential_mediator) == 1){
    
    ##calculate the \gamma and digamma function
    temp_psi1 <- exp(model[[3]][,1,( burnin + 1 ):samples] + model[[1]][,( burnin + 1 ):samples])
    temp_psi0 <- exp(model[[1]][,( burnin + 1 ):samples])
    tot_psi1 <- digamma(apply(temp_psi1,2,sum))
    tot_psi0 <- digamma(apply(temp_psi0,2,sum))
    
    E_log <- matrix(ncol = B_sim, nrow = (samples - burnin))
    temp_balance <- matrix(ncol = (B_sim - 1), nrow = (samples - burnin))
    for(l in 1:B_sim){
      E_log[,l] = digamma(temp_psi1[l,]) - tot_psi1 - (digamma(temp_psi0[l,]) - tot_psi0)
    }
    
    #calculate E[balance_1(T=1)] - E[balance_1(T=0)]
    marginal_mediation <- array(0,dim=c(B_sim,(samples - burnin)))
    for(i in 1:(B_sim)){
      for(j in 1:(samples - burnin)){
        marginal_mediation[i,j] <- sum(balance_structure[i,]*model[[6]][3:(B_sim+1), (j+burnin)])*E_log[j,i]
      } 
    }
    
    global_mediation <- apply(marginal_mediation,2,sum)
    global_mean <- mean(global_mediation)
    global_quantile <- quantile(global_mediation,c(quantile[1],quantile[2]))
    marginal_mean <- mean(marginal_mediation[potential_mediator,])
    marginal_quantile <- quantile(marginal_mediation[potential_mediator,],c(quantile[1],quantile[2]))
    
    return(list(selected = selected_taxa, global_mean = global_mean, global_quantile = global_quantile, 
                marginal_mean = marginal_mean, marginal_quantile = marginal_quantile,
                direct_mean = mean_direct, direct_quantile = quantile_direct,
                covar_result = c(covar_DM, covar_LM)),message = "1 mediator found")
  }
  
  if(length(potential_mediator) > 1){
    return(list(selected = selected_taxa, message = "more than 1 mediators founded, process the input data"))
  }
  
}
