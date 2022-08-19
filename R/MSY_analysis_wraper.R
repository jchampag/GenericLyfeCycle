##################################################
#Script: for a given parameter, make a MSY search: simulates population dynamic for a given Fseq


# Juliette Champagnat
# 17/08/22 upload and adapt for github 
##########################################

"MSYanalysis_wraper" <- function(par, which_par ,
                                 hval=NULL,Kval=NULL,
                                 F_seq,
                                 simulation.name = 'test1',
                                 what_save = 'stat',   #'all'
                                 value_of_interest = "both",#'C',# "SSB" ,  
                                 path_simulation=path_simulation,
                                 n_traj=1,
                                 years=c(1:100),
                                 comments=NULL,
                                 age.plus=F
){
  

  #set few things 
  file_path=paste0(path_simulation,which_par,"_MSY/",genus_name,"_",sp_name) #where the output will be save
  if(!dir.exists(file_path)){dir.create(file_path,recursive = T)}
 
  sc_name<- paste0("DI_",par[2],"_DD_",par[3],"_Surf_",par[4])
  
  if(which_par == "both"){h <- as.numeric(par[5]); K <- as.numeric(par[6])}
  
  # run the analysis
  source("R/MSY_wraper.R",local=T)
  
  purrr::map_dbl(F_seq,MSY_wraper,value_of_interest = value_of_interest ,age.plus=age.plus,
                 n_traj=n_traj,years=years,path_simulation=file_path,simulation.name=sc_name,comments=paste0("K_",K,"_h_",h))
  
  return(paste0(which_par,"_",par))
  
}
