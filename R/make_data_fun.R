##################################################
# Function: for a vector of parameters for a given species generate related data needed for simulation

# Juliette Champagnat
##########################################

#arguments
#' @pred predicted population parameters based on FishLife estimation


"make_data" <- function(preds){
  
  mydata_uncertain <- list()
  
  #####Age related#####
  mydata_uncertain$age_max <- as.numeric(round(exp(preds[5])))
  mydata_uncertain$age_rec <- 1 # recruitment age e.g., 1st mot include in h
  mydata_uncertain$ages <-  c(mydata_uncertain$age_rec:mydata_uncertain$age_max)
  mydata_uncertain$age_sa <- mydata_uncertain$age_rec # 
  mydata_uncertain$age_mat <-  as.numeric(round(exp(preds[6])))
  mydata_uncertain$age_fullselec <- mydata_uncertain$age_mat 
  
  #######Spawning related##########
  
  mydata_uncertain$d_sp <- 100 #data of spawning peak 
  mydata_uncertain$ts_sp <- 30 #time span of spawning season
  
  ######Nurseries related##########
  mydata_uncertain$d_sa <- 1 #date (in day) of nurseries leaving for subadult fish 
  
  ######Zones###########
  mydata_uncertain$spwn_grd <- 1
  mydata_uncertain$spwn_grd_names <- 'spw'
  
  mydata_uncertain$nurseries <-1
  mydata_uncertain$nurseries_names <- 'nurs'
  
  mydata_uncertain$regions <-1
  mydata_uncertain$regions_names <- "region"
  
  #########Female proportion######
  mydata_uncertain$Pf <- rep(0.5,length(mydata_uncertain$ages)) #generic value
  
  ########## Spawning ground related ###########
  # control of density dependence on F_spw_grd #no DD with 0 
  mydata_uncertain$beta_DD <- 0
  #surfaces of spawning grounds #with 1 nothing happens
  mydata_uncertain$Surf_spwn_grd <- rep(1,length(mydata_uncertain$spwn_grd))
  
  ##########Process errors########
  # toutes les mettre a 0 pour le moment 
  mydata_uncertain$sd_Sad <- 0
  mydata_uncertain$sd_Sjuv <- 0
  mydata_uncertain$sd_R <- 0
  
  
  ###### Lenght & Weight at age######
  # parameter of Von B relationship 
  mydata_uncertain$K_growth = as.numeric(exp(preds[2]))
  mydata_uncertain$L_inf = as.numeric(exp(preds[1]))
  mydata_uncertain$LatA <- mydata_uncertain$L_inf * (1-exp(-mydata_uncertain$K_growth*(mydata_uncertain$ages-(-0.1))))
  
  # for weight W(a) = w_cst * L(a) ^ w_pow
  if(!is.na(df_estimates$a)){mydata_uncertain$w_cst <- df_estimates$a}
  if(!is.na(df_estimates$b)){mydata_uncertain$w_pow <- df_estimates$b}
  if(!any(!is.na(df_estimates$a))&!any(!is.na(df_estimates$b))){print("You need some length-weight parameters")}
  
  ## Sometimes the weight constant has to be rescale by 10^3 sometime no
  mydata_uncertain$SwatA <- mydata_uncertain$CwatA <- (mydata_uncertain$w_cst/1000) * mydata_uncertain$LatA ^ mydata_uncertain$w_pow
  
  #I'm using a twisted way to determine is the order of magnitude of Winfinity (from Thorson)
  # and the weight of max age fish match
  
  if(round(log10(as.numeric(exp(preds[7]))))!=round(log10(mydata_uncertain$CwatA[mydata_uncertain$age_max]))){
    mydata_uncertain$SwatA <- mydata_uncertain$CwatA <- (mydata_uncertain$w_cst) * mydata_uncertain$LatA ^ mydata_uncertain$w_pow
  }
  # if not i'm writing it for user to check
  if(round(log10(as.numeric(exp(preds[7]))))!=round(log10(mydata_uncertain$CwatA[mydata_uncertain$age_max]))){
    print("Seems to have a weight issue")
  }
  
  ######## Maturity & Fecondity#############
  
  ## maturity from a knife edge function
  mydata_uncertain$mat <- ifelse(mydata_uncertain$ages<mydata_uncertain$age_mat,0,1)
  
  ## fecundity (values from litterature or rfishbase package)
  
  if(genus_name=="Pleuronectes"&&sp_name=="platessa"){
    
    mydata_uncertain$fec_cst <- 2.33;mydata_uncertain$fec_pow<- 3.10
    
    mydata_uncertain$fec <- ifelse(mydata_uncertain$ages<mydata_uncertain$age_mat,0,
                                   mydata_uncertain$fec_cst*(mydata_uncertain$LatA[mydata_uncertain$ages]^mydata_uncertain$fec_pow))
    
    
  }else if(genus_name=="Leucoraja"&&sp_name=="naevus"){
    mydata_uncertain$fec <- ifelse(mydata_uncertain$ages<mydata_uncertain$age_mat,0,80)
    
  }else if(genus_name=="Engraulis"&sp_name=="encrasicolus"){
    mydata_uncertain$fec  <- 478.9*mydata_uncertain$SwatA*20 #Gatti et al 2016
  }else{
    mydata_uncertain$fec <- ifelse(mydata_uncertain$ages<mydata_uncertain$age_mat,0,
                                   mydata_uncertain$fec_cst*(mydata_uncertain$LatA[mydata_uncertain$ages]^mydata_uncertain$fec_pow))
  }

  ######Eggs survival#######
  mydata_uncertain$p_eggs_surv <- rep(1, length(mydata_uncertain$spwn_grd))
  
  
  #####Selectivity###### 
  # for now using a knife edge starting at maturity (c'est aussi ce que fait Kindsvater)
  mydata_uncertain$Selec <- ifelse(mydata_uncertain$ages<mydata_uncertain$age_mat,0,1)
  
  ########Natural mortality############
  mydata_uncertain$M_at_age <- rep(as.numeric(exp(preds[3])),length(mydata_uncertain$ages))
  
  if(genus_name %in% c("Amblyraja","Leucoraja") &&sp_name%in%c("radiata","naevus")){
    mydata_uncertain$Mlarvae <- 0
  }else {mydata_uncertain$Mlarvae <- 3*log(10)}
  
  ###### Nurseries parameters######
  mydata_uncertain$alpha <- NULL
  mydata_uncertain$Surf <- 1
  mydata_uncertain$K <- 10^6
  mydata_uncertain$h <- as.numeric(preds[4])
  
  #######Dispersion matrices############
  ##useless for now because no multiples zones
  mydata_uncertain$D_larvae <- diag(length(mydata_uncertain$spwn_grd)) #Larval key
  mydata_uncertain$D_subadult <- diag(length(mydata_uncertain$nurseries)) # Sub adult matrix
  mydata_uncertain$D_spawner <- diag(length(mydata_uncertain$spwn_grd)) # Spawner matrix 
  mydata_uncertain$D_return <- diag(length(mydata_uncertain$spwn_grd)) # Post spawning matrix 
  mydata_uncertain$D_adult <- diag(length(mydata_uncertain$regions)) # Adult connectivity 
  
  ####### Inits ##############
  # initial numbers as a decreasing function of K and natural mortality
  #/!\ Eggs estimated from adult 
  mydata_uncertain$NatA_0 <- purrr::map_dbl(mydata_uncertain$ages,
                                            function(x){mydata_uncertain$K*exp(-sum(mydata_uncertain$M_at_age[1:x]))})
  
  mydata_uncertain$Egg_0 <- sum(mydata_uncertain$NatA_0*mydata_uncertain$mat*mydata_uncertain$fec*mydata_uncertain$Pf,na.rm = T)
  
  
  return(mydata_uncertain)
}

