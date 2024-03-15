##################################################
#Function: for a list of parameters for given species compute
# the average number of eggs produced by an individual of age a

# Juliette Champagnat
##########################################

#arguments
#' @a  age for which W is computed
#' @data dataset used for grabbing parameters

"Wbar_computing" <- function(a,data){

  source("R/survival_fun.R",local=T)

  #compute survival until age a
  S_a <- survival_computingSIMPLE(a=a,Fval=0,data=data,first_age=1)

  
  #the average number of eggs produced by an individual of age a during it's lifetime
  W_a = S_a*data$Pf[a-min(data$ages)+1]*data$mat[a-min(data$ages)+1]*data$fec[a-min(data$ages)+1] ##a+1 because those vectors start at age 0

return(W_a)
}