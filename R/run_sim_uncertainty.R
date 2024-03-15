##########################################################
#Script : for a species perform analysis of population responses to nursery habitat 
# degradation/restoration together with considering population parameters uncertainty

# Juliette Champagnat
###########################################################

rm(list=ls())
library(tidyverse);library(rfishbase);library(FishLife)

#########################################################
############ 1. Load data ###############################
#########################################################

#specify species of interest name
genus_name <- "Pleuronectes"#"Leucoraja"# "Engraulis"#
sp_name <- "platessa"#"naevus"#"encrasicolus"#


#grab from rfishbase et Fishlife
df_estimates <- estimate(paste0(genus_name," ",sp_name))
Predict = Plot_taxa( Search_species(Genus=genus_name,Species=sp_name)$match_taxonomy, mfrow=c(2,3),Database=FishLife::FishBase_and_RAM)

## uncertainty reduction
source("R/varcovar_reduction.R",local=T)

#use Fishlife prediciton to generate an array of parameters using the varcovar matrix just created 
nb_tirage <- 10
varcovar_sp <- varcovar_ple # /!\don't forget to change input matrix according to species
predictions <- MASS::mvrnorm(n=nb_tirage, mu=Predict[[1]]$Mean_pred[c(1,2,6,13)],Sigma =varcovar_sp ) %>% as.data.frame()

## add parameters considered as fixed in the table 
predictions <- predictions %>% mutate(tmax=Predict[[1]]$Mean_pred[4],tm=Predict[[1]]$Mean_pred[5],
                                      Winf=Predict[[1]]$Mean_pred[3])

#check parameters values
predictions %>%  mutate(across(.cols=!c(h),.fns =exp )) %>%  summary()

# plot parameter uncertainty
predictions %>% #names()
  # select("Loo","K","tmax","tm","M","h") %>% 
  mutate(across(.cols=!h,.fns =exp )) %>% 
  # mutate(across(.cols=c(tmax,tm),.fns = round)) %>% 
  mutate(tirage=c(1:nb_tirage)) %>%
  pivot_longer(-tirage) %>% 
  ggplot()+geom_boxplot(aes(x=1,y=value))+facet_wrap(~name,scale="free_y")

######################################################
####### 2. Generate data list ########################
######################################################


## create a list with all inputs for each draw
source("R/make_data_fun.R",local=T)
list_params <- apply(predictions,1,make_data)


## generate Wbar values
source("R/Wbar_computing_fun.R",local=T)
for(i in 1:nb_tirage){
  list_params[[i]]$Wbar <- sum(map_dbl(list_params[[i]]$ages,Wbar_computing,data=list_params[[i]]),na.rm=T)} #,M_spwgrd_surv=M_spwgrd_surv,

# list_params[[1]]$Wbar;list_params[[2]]$Wbar


###############################################################
######## 3. Apply multipliers #################################
###############################################################


##compute new values of (h,K) with multiplier + check viability
source("R/quality_multiplier_fun.R",local=T)

simu_par_all <- NULL

for(i in 1:nb_tirage){
  #for 1 draw, computes new values of h and K according to mutipliers & others parameters values
  
  ## uses mult_quality() for habitat quality multipliers
  simu_par_qual <- mult_quality(lambda_DI=c(0.75,0.85,0.95,1,1.05,1.15,1.25),
                                lambda_DD=c(0.75,0.85,0.95,1,1.05,1.15,1.25),#seq(0.75,1.25,0.05)
                                deltaT=1,Mlarvae=list_params[[i]]$Mlarvae,Wbar= list_params[[i]]$Wbar,
                                h=list_params[[i]]$h,K=list_params[[i]]$K)
  
  ## makes a simple division for habitat surface area multiplier
  mult_values_surf <- seq(0.75,1.25,0.1)
  simu_par_surf <- data.frame(type="surface",lambda_DI=1,lambda_DD=1,
                              lambda_surf=mult_values_surf,
                              h=rep(list_params[[i]]$h,length(mult_values_surf)),
                              K=list_params[[i]]$K/mult_values_surf)
  
  ## bind the tables + check population viability with parameters modification
  simu_par <- rbind(simu_par_qual,simu_par_surf) %>% 
    mutate(deltaT=1,Mlarvae=list_params[[i]]$Mlarvae,Wbar=list_params[[i]]$Wbar) %>% 
    mutate(surv_max = 4*h/(Wbar*(1-h))) %>% 
    mutate(alpha = 4*h/(exp(-Mlarvae)*Wbar*(1-h))) %>% 
    mutate(M_DI = -log(alpha)/deltaT) %>% 
    mutate(M_DD = M_DI/(K*(exp(M_DI*deltaT)-1))) %>% #%>% summary()
    mutate(draw=i)
  
  if(any(simu_par %>% pull(alpha)<0|simu_par %>% pull(alpha)>1)){
    print(paste0("Draw ",i, " generates inapropriate values of alpha"))}
  if(any(simu_par %>% pull(M_DI)<0)){
    print(paste0("Draw ",i, " generates inapropriate values of M_DI"))}
  if(any(simu_par %>% pull(M_DD)>1)){
    print(paste0("Draw ",i, " generates inapropriate values of M_DD"))}
  
  simu_par_all <- bind_rows(simu_par_all,simu_par) 
}

# can check population parameter sets
simu_par_all %>% summary()

###############################################
####### 4. Simulations ########################
###############################################

# main analysis function + underlying function
source('R/MSYanalysis_wraper.R') #main
source('R/MSY_wraper.R') # for a set of F value, simulate population
source("R/data_formating_fun.R") # data formating 
source('R/sim_function_new.R') # simulates population dynamics


## WITH NO PARALLELIZATION

df_results <- data.table::rbindlist(lapply(seq(1,nrow(simu_par_all)),FUN=MSYanalysis_wraper,
                                           par=simu_par_all,F_seq = seq(0,1.5,0.01),
                             years=c(1:1000),data=list_params))

## WITH PARALLELIZATION (Advised if using multiple parameters draw for 1 species)
library(parallel)
cl <- makeCluster(2)

clusterEvalQ(cl, {library(tidyverse)})

clusterExport(cl,
              c("MSYanalysis_wraper","MSY_wraper","format_data","life_cycle_sim", #functions
                   "list_params","simu_par_all" #objects
))


results <- parLapply(cl,seq(1,nrow(simu_par_all)),fun=MSYanalysis_wraper,par=simu_par_all,F_seq = seq(0,1.5,0.01),
                     years=c(1:1000),
                     data=list_params)
stopCluster(cl)

df_results <- data.table::rbindlist(results)


############################################
######## 5. Save parameters & results ######
############################################

write.table(simu_par_all,
            file = paste0("Simulation/",genus_name,"_",sp_name,"_simu_parameters.txt"),row.names = F)

save(list_params,file= paste0("Simulation/",genus_name,"_",sp_name,"_simu_data.Rdata"))

write.table(data.table::rbindlist(results),
            file = paste0("Simulation/",genus_name,"_",sp_name,"_res.txt"),row.names = F)

