##################################################
#Function: for a given parameter, make a MSY search: simulates population dynamic for a given Fseq

# Juliette Champagnat
##########################################

#arguments
#' @ligne indexing number for subseting parameters matrix
#' @par a matrix of parameters for a giver species
#' @Fseq  a range of fishing mortality values covered to search MSY
#' @ntraj numnber of trajectories to run
#' @years year duration of the simulation
#' @comments something to pass out in the simulation name #maybe useless
#' @data dataset used to grab population parameters


"MSYanalysis_wraper" <- function(ligne,par,
                                     F_seq,
                                     n_traj=1,
                                     years=c(1:500),
                                     age.plus=F,
                                     data
){
  
  par <- par[ligne,]
  
  data[[as.numeric(par[14])]]$h <- as.numeric(par[5])
  data[[as.numeric(par[14])]]$K <- as.numeric(par[6])
  
  # run the analysis
  source("R/MSY_wraper.R",local=T)
  
  dfs <- purrr::map_df(F_seq,MSY_wraper,age.plus=age.plus,data= data[[as.numeric(par[14])]],
                       n_traj=n_traj,years=years)
  
  dfs <- dfs %>% 
    filter(Catch==max(Catch)) %>% 
    mutate(type=as.character(par[1]),lambda_DI=as.numeric(par[2]),lambda_DD=as.numeric(par[3]), lambda_surf=as.numeric(par[4]),
           h=as.numeric(par[5]),K=as.numeric(par[6]),
           Mlarvae=data[[as.numeric(par[14])]]$Mlarvae,draw=as.numeric(par[14])) %>% 
    rename(F_MSY=Fval,C_MSY=Catch,SSB_MSY=SSB)
  
  return(dfs)
  
}


