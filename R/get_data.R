###########################################################
# Script: for a given species name grab related data needed for simulation

# Juliette C
# started 28/09/21
# 17/08/22 upload and adapt for github 
###########################################################

## first on Fishbase
library(rfishbase);library(FishLife)

DYN=F
#########
# vignette("tutorial","FishLife")

#grab from rfishbase
df_estimates <- estimate(paste0(genus_name," ",sp_name))
df_fecundity <- fecundity(paste0(genus_name," ",sp_name))

#grab from Fishlife
Predict = Plot_taxa( Search_species(Genus=genus_name,Species=sp_name)$match_taxonomy, mfrow=c(2,3),Database=FishLife::FishBase_and_RAM)

#####Age related#####
age_max <- as.numeric(round(exp(Predict[[1]]$Mean_pred[4])))
age_rec <- 1 #age de recrutement ie 1er qui n'est pas englobé dans le h
ages <-  c(age_rec:age_max) # is age_max the age we want to use here? 
age_sa <- age_rec # how to find it?
age_mat <- as.numeric(round(exp(Predict[[1]]$Mean_pred[5]))) #
age_fullselec <- age_mat #how to proceed ? same for all fish? 

#######Spawning related##########

d_sp <- 100 #data of spawning peak 
ts_sp <- 30 #time span of spawning season

######Nurseries related##########
d_sa <- 1 #date (in day) of nurseries leaving for subadult fish 

######Zones###########
spwn_grd <- 1
spwn_grd_names <- 'spw'

nurseries <-1
nurseries_names <- 'nurs'

regions <-1
regions_names <- "region"

#########Female proportion######
# if(is.na(spawning(paste0(genus_name," ",sp_name)) %>% pull(SexRatiomid))#,SexRmodRef))){
#   
# }else{
Pf <- rep(0.5,length(ages)) #generic value}

########## Spawning ground related ###########
# controle de la densite dependance de F_spw_grd #avec 0 pas d'effet DD
beta_DD <- 0
#surfaces of spawning grounds #with 1 nothing happens
Surf_spwn_grd <- rep(1,length(spwn_grd))

##########Process errors########
# toutes les mettre a 0 pour le moment 
sd_Sad <- 0
sd_Sjuv <- 0
sd_R <- 0


###### Lenght & Weight at age######
# parameter of Von B relationship 

# for length L(a)=L_inf * (1-exp(-K_growth*a))
K_growth = as.numeric(exp(Predict[[1]]$Mean_pred[2]))
L_inf = as.numeric(exp(Predict[[1]]$Mean_pred[1]))
LatA <- L_inf * (1-exp(-K_growth*(ages-(-0.1))))

# for weight W(a) = w_cst * L(a) ^ w_pow
if(!is.na(df_estimates$a)){w_cst <- df_estimates$a}
if(!is.na(df_estimates$b)){w_pow <- df_estimates$b}
if(!any(!is.na(df_estimates$a))&!any(!is.na(df_estimates$b))){print("You need some length-weight parameters")}

## Sometimes the weight constant has to be rescale by 10^3 sometime no
SWatA <- CWatA <- (w_cst/1000) * LatA ^ w_pow

#I'm using a twisted way to determine is the order of magnitude of Winfinity (from Thorson)
# and the weight of max age fish match

if(round(log10(as.numeric(exp(Predict[[1]]$Mean_pred[3]))))!=round(log10(CWatA[age_max]))){
  SWatA <- CWatA <- (w_cst) * LatA ^ w_pow
}
# if not i'm writing it for user to check
if(round(log10(as.numeric(exp(Predict[[1]]$Mean_pred[3]))))!=round(log10(CWatA[age_max]))){
  print("Seems to have a weight issue")
}


######## Maturity & Fecondity#############

## maturity 

# from a maturation probability function (Kinsvater) mat(a=) = 1/(1+e(-q*(a-a_mat)))

# from a knife edge function
mat <- ifelse(ages<age_mat,0,1)

## fecundity

# from a fonction of length fec(a) = fec_cst * L(a) ^ fec_pow (estimates from rfishbase)

if(any(!is.na(df_fecundity$a))&any(!is.na(df_fecundity$b))){ # there are fec-L estimates in fish base
  
  if(sum(!is.na(df_fecundity$a))>1|sum(!is.na(df_fecundity$b))>1) #there are more than one L-fec estimate
    
    if(is.null(fec_study)){ #need to choose which fec study will be used
      print("Decide which fec study you want to use and enter fec_study number (line of the study in the df_fecundity)")
      print(df_fecundity %>% select(a,b,r2,Number))
      
    }else{ #take a and b estimates from this study
      fec_cst <- df_fecundity[fec_study,]$a
      fec_pow <-  df_fecundity[fec_study,]$b
    }
}
# if(any(!is.na(df_fecundity$a))){fec_cst <- mean(df_fecundity$a,na.rm = T)}
# if(any(!is.na(df_fecundity$b))){fec_pow <- mean(df_fecundity$b,na.rm = T)}
# #if(!any(!is.na(df_fecundity$a))&!any(!is.na(df_fecundity$b))){print("You need some fec-length parameters")}

#for some species i'm gathering fec information from litterature
if(genus_name=="Solea"&sp_name=="solea"){
  fec <- ifelse(ages<age_mat,0,exp(5.619+1.17*log(SWatA*1000))) #using Rochette equation
  
}else if(genus_name=="Cetorhinus"&&sp_name=="maximus"){
  fec <- ifelse(ages<age_mat,0,1.5) # according to litterature mean is 1.5 max is 2
  
}else if (genus_name=="Limanda"&&sp_name=="limanda"){
  fec<-ifelse(ages<age_mat,0, 7.2783 * LatA[ages] ^ 3.4525) # from Bohl (1057) estimate for North sea 
  
}else if(genus_name=="Coryphaenoides"&&sp_name=="rupestris"){
  fec <- ifelse(ages<age_mat,0,28.54 * SWatA[ages] -2187.3) #from Allain, V. (2001)
  
}else if(genus_name=="Scophthalmus"&&sp_name=="maximus"){
  fec <- ifelse(ages<age_mat,0,exp(5.820)+SWatA[ages]^1.298) #from Nissling et al 2013
  
}else if(genus_name=="Pleuronectes"&&sp_name=="platessa"){
  plaice_fecdata <- data.frame(Age =c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30),
                               Length=c(0.0,0.8,19.4,26.2,31.5,35.7,39.0,41.6,43.7,45.3,46.5,47.5,48.3,49.0,49.4,49.8,50.1,50.4,50.6,50.7,50.8,50.9,51.0,51.0,51.1,51.1,51.1,51.2,51.2,51.2,51.2),
                               Weight=c(0.000,0.013,0.073,0.180,0.313,0.456,0.594,0.721,0.833,0.928,1.009,1.075,1.129,1.173,1.209,1.237,1.260,1.278,1.292,1.303,1.312,1.320,1.325,1.330,1.333,1.336,1.338,1.340,1.341,1.342,1.343),
                               mat=c(0,0,0.01,0.06,0.29,0.80,0.99,rep(1,24)),
                               fec=c(0,0,0,0.1,3.5,38.6,77.2,95.5,110.3,123.0,133.6,142.4,149.6,155.5,160.1,163.9,166.9,169.3,171.2,172.7,173.9,174.8,175.6,176.2,176.6,177.0,177.3,177.5,177.7,177.8,177.9))
  
  fec <- (plaice_fecdata %>% filter(Age%in%ages) %>% pull(fec))*10^3
  
}else if(genus_name=="Amblyraja"&&sp_name=="radiata"){
  fec <- ifelse(ages<age_mat,0,32)
  
}else if(genus_name=="Leucoraja"&&sp_name=="naevus"){
  fec <- ifelse(ages<age_mat,0,80)
  
}else if(genus_name=="Engraulis"&sp_name=="encrasicolus"){
  # fec <- ifelse(ages == 1,11188*SWatA,8750*SWatA) #Motos 1996
  # fec <- ifelse(LatA<14,11188*SWatA,8750*SWatA) #Motos 1996
  #fec <- ifelse(ages == 1,110000,350000) #Motos 1996
  # fec <- 478.9*SWatA*20 #Gatti et al 2016
  fec <- 550*SWatA*30 #Pethybrigde et al 2013
}else{
  fec <- ifelse(ages<age_mat,0,fec_cst*(LatA[ages]^fec_pow))
}
if(sp_name=="merluccius"){fec <- fec*1000}

######Eggs survival#######
p_eggs_surv <- rep(1, length(spwn_grd))


#####Selectivity###### 
# for now using a knife edge starting at maturity (c'est aussi ce que fait Kindsvater)
Selec <- ifelse(ages<age_mat,0,1)

########Natural mortality############
M_at_age <- rep(exp(Predict[[1]]$Mean_pred[6]),length(ages))

if(genus_name %in% c("Amblyraja","Leucoraja") &&sp_name%in%c("radiata","naevus")){
  Mlarvae <- 1
}else {Mlarvae <- 3*log(10)}

###### Nurseries parameters######
alpha <- NULL
Surf <- 1
K <- 10^6
h <- Predict[[1]]$Mean_pred[13]

#######Dispersion matrices############
##useless for now because no multiples zones
D_larvae <- diag(length(spwn_grd)) #Larval key
D_subadult <- diag(length(nurseries)) # Sub adult matrix
D_spawner <- diag(length(spwn_grd)) # Spawner matrix 
D_return <- diag(length(spwn_grd)) # Post spawning matrix 
D_adult <- diag(length(regions)) # Adult connectivity 

####### Inits ##############
# initial numbers as a decreasing function of K and natural mortality
#/!\ Eggs estimated from adult 
NatA_0 <- purrr::map_dbl(ages,function(x){K*exp(-sum(M_at_age[1:x]))})
Egg_0 <- sum(NatA_0*mat*fec*Pf,na.rm = T)

if(genus_name %in% c("Amblyraja","Leucoraja") &&sp_name%in%c("radiata","naevus")){
 NatA_0=NatA_0/10;Egg_0=Egg_0/10}

#######################################################################################################################


