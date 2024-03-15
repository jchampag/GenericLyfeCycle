##################################################
#Function: return the value of h and K (e.g., Beverton-Holt stock-recruitement parameters) for a vector of quality multiplier

# Juliette Champagnat
##################################################

#arguments
#' @lambda_DI multiplier applied to density-independent mortality rate
#' @lambda_DD multiplier applied to density-dependent mortality rate
#' @h  initial vqlue of steepness
#' @K initial value of carrying capacity
#' @deltaT duration of the recruitement phase
#' @Wbar the average number of eggs produced by a recruit during its lifetime



"mult_quality" <- function(lambda_DI,lambda_DD,h,K,Mlarvae,deltaT=1,Wbar){
  
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
  vec_K <- vec_M_DI/(vec_M_DD*(exp(vec_M_DI*deltaT)-1))
  
  
  return(data.frame(type="quality",lambda_DI=lambda_DI,lambda_DD=lambda_DD,lambda_surf=1,h=vec_h,K=vec_K))
}
