##########################################################
#Script : MAIN SCRIPT
# for an archetype perform analysis of responses to quality and quantity multiplier

# Juliette Champagnat
# started 28/09/2021
# 17/08/22 upload and adapt for github 
###########################################################

rm(list=ls())
library(tidyverse);library(reshape2) 


#####################################
####### Load and format data ########
#####################################

#specify species studied
genus_name <- "Leucoraja"#"Merluccius"#"Engraulis"#"Pleuronectes"#
sp_name <-  "naevus"#"merluccius" #"encrasicolus"#"platessa"#

#load data
fec_study <- ifelse(sp_name=="merluccius",3,0)
source("R/get_data.R")


## set and format data for simulation
years <- c(1:1000) #specify simulation duration
Fval <- Selec*0.2;CWatAval <- CWatA;SWatAval <- SWatA #initilize some data frame
source("R/make_data.R",local=T) #format data
source("R/Wbar_function.R",local=T) #compute Wbar according to data


##########################################################
####### Modify SR parameters according to habitat ########
##########################################################

## setting habitat multipliers
mult_values_DI <- mult_values_DD <- c(1,0.5,0.75,0.85,0.95,1.05,1.15,1.25,1.5)
mult_values_surf <- c(0.5,0.75,0.85,0.95,1.05,1.15,1.25,1.5)

# create a table with couples of modified (h,K) values in relation to multipliers (DD,DI and surf) entered
source("R/multiplier_function.R",local=T)

## marginal scenarios 
simu_par<- mult_fun(lambda_DI=c(mult_values_DI,rep(1,length(mult_values_surf))),
                            lambda_DD=c(mult_values_DD,rep(1,length(mult_values_surf))),
                            lambda_surf = c(rep(1,length(mult_values_DD)),mult_values_surf),
                            deltaT=1,Mlarvae=Mlarvae,Wbar= mydata$Wbar,h=h,K=K)

## or all combinations 
all_crossover <- expand_grid(DI=seq(0.5,1.5,0.05),surf=seq(0.5,1.5,0.05)) 
simu_par<- mult_fun(lambda_DI=all_croisement$DI,lambda_DD=all_croisement$DI,
                           lambda_surf = all_croisement$surf,
                           deltaT=1,Mlarvae=Mlarvae,Wbar= mydata$Wbar,h=h,K=K)


#########################################################################
####### run a series of MSY analysis with modified SR parameters ########
#########################################################################


source('R/MSY_analysis_wraper.R')


apply(simu_par[c(1,2),],1,MSYanalysis_wraper,which_par = "both",  F_seq = seq(0,1.5,0.01),
      path_simulation=paste0("Simulation/"), years=c(1:1000),n_traj=1,value_of_interest="both")

######################################
####### basic plot of results ########
######################################

#grab results
which_par ="both"
path_simulation <- paste0("Simulation/",which_par,"_MSY/",genus_name,"_",sp_name)

## grab all MSY run from txt in a df
MSYfiles <- list.files(path_simulation,pattern = paste0("MSY_"))
df_MSY <- do.call("rbind",purrr::map(MSYfiles,function(x){read.table(paste0(path_simulation,"/",x) ,header=T)%>% mutate(simu=x)}))


##plot productivity curves
df_MSY %>% 
  filter(Zone =="ALL") %>% #,grepl(sim,pattern="0.9")) %>% 
  pivot_wider(id_cols = c(Zone,Fvalue,sim),names_from=type,values_from=c(mean))%>%
  filter(Zone=="ALL") %>% ggplot()+geom_line(aes(x=Fvalue,y=C,col=sim))


