##################################################
#Function: from a list of input data, format them according to simulation function needs and return a new list

# Juliette Champagnat
##########################################

#arguments
#' @data dataset used to grab population parameters
#' @yr year duration of the simulation
#' @Fvalue input fishing mortality value for the simulation


"format_data" <- function(data, yr=1:1000,Fvalue){
  
  # generate a new list
  data_formatted <- data 
  
  ########## Duration of the simulation ###########
  data_formatted$years <- yr
  
  ##### Nursery parameters #######
  data_formatted$Surf_nurs <- array(rep(data$Surf,each=length(yr)), dim=c(length(yr),length(data$nurseries)),
                                    dimnames=list(yr,data$nurseries_names))
  
  data_formatted$K <- array(rep(data$K,each=length(yr)), dim=c(length(yr),length(data$nurseries)),
                            dimnames=list(yr,data$nurseries_names))
  
  data_formatted$h <- array(rep(data$h,each=length(yr)), dim=c(length(yr),length(data$nurseries)),
                            dimnames=list(yr,data$nurseries_names))
  
  ####### Fishing mortality #########
  data_formatted$F_nurseries <- array(0,dim=c(length(data$ages),length(yr),length(data$nurseries)),
                                      dimnames = list(data$ages,yr,data$nurseries_names))#
  data_formatted$F_nurseries[1:data$age_sa-min(data$ages)+1,,] <- Fvalue[1:data$age_sa-min(data$ages)+1]
  
  data_formatted$F_regions <- array(0,dim=c(length(data$ages),length(yr),length(data$regions)),
                                    dimnames = list(data$ages,yr,data$regions_names))#
  
  #taking Fval for ages present in adult regions
  #/!\ Fval is now considered as a vector, should be modified when simulating with multiple region/spawning ground
  for(r in seq_along(data$regions)){
    data_formatted$F_regions[(data$age_sa-min(data$ages)+1):(max(data$ages)-min(data$ages)+1),,r] <- 
      Fvalue[(data$age_sa-min(data$ages)+1):(max(data$ages)-min(data$ages)+1)]
  }
  
  
  data_formatted$F_spwn_grd <- array(0,dim=c(length(data$ages),length(yr),length(data$spwn_grd)),
                                     dimnames = list(data$ages,yr,data$spwn_grd_names))#
  
  for(s in seq_along(data$spwn_grd)){
    data_formatted$F_spwn_grd[(data$age_mat-min(data$ages)+1):(max(data$ages)-min(data$ages)+1),,s] <- 
      Fvalue[(data$age_mat-min(data$ages)+1):(max(data$ages)-min(data$ages)+1)]
  }
  
  ########### Natural Mortality ##########
  ##on nurseries
  data_formatted$M_nurseries <- array(0, dim=c(length(data$ages),length(yr),length(data$nurseries)),
                                      dimnames = list(data$ages,yr,data$nurseries))#
  
  data_formatted$M_nurseries[1:(data$age_sa-min(data$ages)+1),,] <- data$M_at_age[1:(data$age_sa-min(data$ages)+1)] 
  
  
  
  ##on regions
  data_formatted$M_regions <- array(0, dim=c(length(data$ages),length(yr),length(data$regions)),
                                    dimnames = list(data$ages,yr,data$regions))#
  
  for (a in (data$age_sa-min(data$ages)+1):(max(data$ages)-min(data$ages)+1)){
    data_formatted$M_regions[a,,] <- data$M_at_age[a]}
  
  
  #on spawning ground
  data_formatted$M_spwn_grd <- data_formatted$M_regions#*M_spwgrd_mult
  data_formatted$M_spwn_grd[1:(data$age_mat-min(data$ages)+1-1),,] <- 0 #for fish unmature fish => no mortality on spw ground
  data_formatted$M_spwgrd_surv <- 1 #multiplier
  
  ###### initialization values #########
  data_formatted$Na_year1 <- data$NatA_0
  data_formatted$Egg_0 <- data$Egg_0
  if(!is.null(data$Mlarvae)){data_formatted$Larvae_year1 <- data$Egg_0*exp(-data$Mlarvae)}
  
  ######## standard deviation for process error #########
  data_formatted$sd_Sad <- 0#sd_Sad
  data_formatted$sd_Sjuv <- 0#sd_Sjuv
  data_formatted$sd_R <- 0#sd_R
  
  return(data_formatted)
  
}
