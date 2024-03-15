##################################################
#Funtion: for a given F combination (3 values), this script:
# simulates population dynamic for a given number of year 
# then computes # mean Catch and/or SSB for the last half of the simulation period

# Juliette Champagnat
##########################################

#arguments
#' @F_fixed_value a value of fishing mortality for which the simulation is run
#' @ntraj numnber of trajectories to run
#' @years year duration of the simulation
#' @age.plus is the last age an age plus group?
#' @comments something to pass out in the simulation name #maybe useless
#' @data dataset used to grab population parameters


"MSY_wraper" <- function(F_fixed_value,
                             n_traj=1,
                             years=c(1:500),
                             age.plus=F,
                             data
                             
){

  #create Fval from the wanted F combination and selectivity at age  
  Fval <- matrix(data=NA,nrow=length(data$ages),ncol=length(data$regions))

  for (a in 1:length(data$ages)){
    #for (z in 1:length(regions)){ 
    Fval[a] <- as.numeric(data$Selec[a] * F_fixed_value)
  }
  
  #format and arrange data 
  newdata <- format_data(data=data,yr=years,Fvalue=Fval)
  
  #### load simulation function #####
  source("R/sim_function_new.R",local=T)
  
  #### run simulation ######
  S1 <-  life_cycle_sim(data=newdata,traj=n_traj,param_R = "h",age.plus = age.plus)
  
  ## return output
  df <- data.frame(Fval=F_fixed_value,Catch=S1$C_tot[max(years)-1],SSB=S1$SSB[max(years)-1])
  
  
  return(df)
}
