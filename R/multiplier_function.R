##################################################
#Script: this function return the value of h and K for vectors of multipliers


# Juliette Champagnat
# started 2/07/21
# 17/08/22 upload and adapt for github 
##################################################


"mult_fun" <- function(lambda_DI,lambda_DD,lambda_surf,h,K,Mlarvae,deltaT,Wbar){
  
  # retrocalculation of alpha and M_DI from h
  retro_alpha <- as.numeric(4*h/(exp(-Mlarvae)*Wbar*(1-h)))
  M_DI_0 <- -log(retro_alpha)/deltaT
  
  # retrocalculation of M_DD from K and retro_M_DI
  M_DD_0 <- M_DI_0/(K*(exp(M_DI_0*deltaT)-1))
  
  #creation of a vector of M_DI based on retrocalculated value & multiplicator vec
  vec_M_DI <- ifelse(lambda_DI==0,M_DI_0,lambda_DI*M_DI_0)#(M_DI_0+lambda_quali*M_DI_0) #
  
  #idem for M_DD
  vec_M_DD <- ifelse(lambda_DD==0,M_DD_0,lambda_DD*M_DD_0)#(M_DI_0+lambda_quali*M_DI_0) #
  
  #computation of vector of h and K associated with the multiplicators
  vec_h <- exp(-vec_M_DI*deltaT)*exp(-Mlarvae)*Wbar/(4+(exp(-vec_M_DI*deltaT)*exp(-Mlarvae)*Wbar))
  vec_K_temp <- vec_M_DI/(vec_M_DD*(exp(vec_M_DI*deltaT)-1))
  
  #adding lambda_surf effects
  vec_K <-vec_K_temp/lambda_surf
  
  #defining sim type according to multiplier
  sim_type <- NULL
  for (i in seq_along(lambda_DI)){
    if(lambda_DI[i]!=1){
      if(lambda_surf[i]==1){sim_type_temp="quality"}else{sim_type_temp="both"}
    }else{
      if(lambda_surf[i]==1){sim_type_temp="ref"}else{sim_type_temp="surface"}}
    
    sim_type <- c(sim_type,sim_type_temp)}
  
  return(data.frame(type=sim_type,lambda_DI=lambda_DI,lambda_DD=lambda_DD,lambda_surf=lambda_surf,h=vec_h,K=vec_K))
}


