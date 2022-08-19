##################################################
#Script: life cycle simulation with inputs from Make_data.R

# Juliette Champagnat
# 17/08/22 upload and adapt for github 
#################################################

"life_cycle_sim" <- function(data,traj= 1,param_R="h",age.plus=FALSE){
  
  #------------------------------
  ##df creation to save outputs##
  #------------------------------
  
  ## total catches
  C_tot <- array(NA,dim=c(length(data$years), max(traj)), dimnames = list(data$years,1:max(traj)))
  
  # related to spawning grounds
  ##eggs and SSB
  SSB <- Eggs <- array(NA,dim=c(length(data$years), length(data$spwn_grd), max(traj)),
                       dimnames = list(data$years,data$spwn_grd,1:max(traj)))
  
  ##abundance 
  N_spwn_grd_bf <- N_spwn_grd_end<-  C_spwn_grd <- array(data=NA, dim=c(length(data$ages),length(data$years),length(data$spwn_grd),max(traj)),
                                                         dimnames = list(data$ages,data$years,data$spwn_grd,1:max(traj)))
  
  
  # related to nurseries 
  dL <- L <- array(NA,dim=c(length(data$years), length(data$nurseries), max(traj)),
                   dimnames = list(data$years,data$nurseries,1:max(traj))) #larvaes
  
  d_nurseries <- C_nurseries <- array(data=NA, dim=c(length(data$ages),length(data$years),length(data$nurseries),max(traj)),
                                      dimnames = list(data$ages,data$years,data$nurseries,1:max(traj)))
  
  N_nurseries <- array(data=NA, dim=c(length(data$ages),(length(data$years)+1),length(data$nurseries),max(traj)),
                       dimnames = list(data$ages,c(data$years,(max(data$years)+1)),data$nurseries,1:max(traj)))
  
  #related to stock regions
  C_regions <- C_regions_bf<- array(data=NA, dim=c(length(data$ages),length(data$years),length(data$regions),max(traj)),
                                    dimnames = list(data$ages,data$years,data$regions,1:max(traj)))
  
  N_regions <- N_regions_af_mvt <- array(data=NA, dim=c(length(data$ages),(length(data$years)+1),length(data$regions),max(traj)),
                                         dimnames = list(data$ages,c(data$years,max(data$years)+1),data$regions,1:max(traj)))
  
  for(s in 1:traj){
    
    # initilization
    N_regions[(data$age_sa-min(data$ages)+2):(max(data$ages)-min(data$ages)+1),1,,s] <- data$Na_year1[(data$age_sa-min(data$ages)+2):(max(data$ages)-min(data$ages)+1)]
    #N_nurseries[1,1,,s] <-
    N_nurseries[1:(data$age_sa-min(data$ages)+1),1,,s] <-  data$Na_year1[1:(data$age_sa-min(data$ages)+1)]
    #N_nurseries[3,1,,s] <- data$Na_year1[age_sa] #ça devrait pouvoir être enlevé 30/04/21
    
    L[1,,s] <- sum(data$Larvae_year1)
    
    for (y in 1:max(data$years)){
      
      for (a in 1:(max(data$ages)-min(data$ages)+1)){
        
        #############
        ##nursery phase
        #############
        
        if(a <data$age_sa-min(data$ages)+1){ #for fish of age < age_subadult 
          if (param_R=="BH"){
            
            for(i in seq_along(data$nurseries)){
              
              
              #number of larvae arriving at each nursery
              
              L[y,i,s] <- ifelse(y==1, L[1,,s] ,sum(data$D_larvae[,i]*Eggs[y-1,,s]*exp(-data$Mlarvae)))
              
              
              # number of age 0 fish produced by larval input
              # Beverton Holt formulation
              N_nurseries[1,y,i,s] <- data$alpha[y,i]*L[y,i,s]/(1+L[y,i,s]*(data$alpha[y,i]/(data$K[y,i]*data$Surf_nurs[y,i])))
              
            }#closing loop on nurseries  
            
            
            
            if(a==1){
              #computing number at age in each nursery
              # for age 1 the mortality is divide by three because it happens between the end of the summer and December
              N_nurseries[a+1,y+1,,s] <- N_nurseries[a,y,,s]*
                exp(-(data$M_nurseries[a,y,]+data$F_nurseries[a,y,])/1)#*
              # exp(pe_Sjuv[a,y,i,s]-data$sd_Sjuv/2)
              
              #compute catches
              C_nurseries[a,y,,s] <- N_nurseries[a,y,,s]*(data$F_nurseries[a,y,]/(data$F_nurseries[a,y,]+data$M_nurseries[a,y,]))*
                (1-exp(-(data$F_nurseries[a,y,]+data$M_nurseries[a,y,])/1))*(data$CwatA[a])
              
            }else{ #Closing if on age 1  
              
              
              #computing number at age in each nursery
              N_nurseries[a+1,y+1,,s] <- N_nurseries[a,y,,s]*
                exp(-(data$M_nurseries[a,y,]+data$F_nurseries[a,y,]))#*
              #exp(pe_Sjuv[a,y,i,s]-data$sd_Sjuv/2)
              
              #compute catches
              C_nurseries[a,y,,s] <- N_nurseries[a,y,,s]*(data$F_nurseries[a,y,]/(data$F_nurseries[a,y,]+data$M_nurseries[a,y,]))*
                (1-exp(-(data$F_nurseries[a,y,]+data$M_nurseries[a,y,])))*(data$CwatA[a])
              
            }
            
            #}#closing if on juvenile age
            
          } else { #closing if on recrutement parametrization
            # else = param_h
            
            if(a==(data$age_rec-min(data$ages)+1)){ #pour l'age du recrutement je calcule avec la SR en h 
              
              for(i in seq_along(data$nurseries)){
                # /!\ en spatialisant je n'ai aucune idée de comment je vais faire
                
                # /!\ il me faut un Egg_y0 pour initiliser !
                numerateur<- 4*data$h[y,i]*data$p_eggs_surv*ifelse(y==1,data$Egg_0,sum(data$D_larvae[,i]*Eggs[y-1,,s]))/(data$Wbar*(1-data$h[y,i]))
                
                denominateur <- 1+(4*data$h[y,i]*data$p_eggs_surv*ifelse(y==1,data$Egg_0,sum(data$D_larvae[,i]*Eggs[y-1,,s]))/(data$Wbar*(1-data$h[y,i])*data$K[y,i]*data$Surf_nurs[y,i]))
                
                N_nurseries[(data$age_rec-min(data$ages)+1),y,i,s] <- numerateur/denominateur
                
              } #close loop on nurseries
              
            } #close if on age_rec
            
            if(a>=(data$age_rec-min(data$ages)+1) & a<(data$age_sa-min(data$ages)+1)){
              #computing number at age in each nursery
              N_nurseries[a+1,y+1,,s] <- N_nurseries[a,y,,s]*
                exp(-(data$M_nurseries[a,y,]+data$F_nurseries[a,y,]))#*
              #exp(pe_Sjuv[a,y,i,s]-data$sd_Sjuv/2)
              
              #compute catches
              C_nurseries[a,y,,s] <- N_nurseries[a,y,,s]*(data$F_nurseries[a,y,]/(data$F_nurseries[a,y,]+data$M_nurseries[a,y,]))*
                (1-exp(-(data$F_nurseries[a,y,]+data$M_nurseries[a,y,])))*(data$CwatA[a])
              
            } #close if on age<age_sa
          } #closing else on param_R
          
        }else if(a==(data$age_sa-min(data$ages)+1)& data$age_sa!=data$age_mat){#closing if on juvenile age + if cas simple si age_sa < age_mat 
          
          ###########
          ##adults pop
          ###########
          
          
          #####################
          ##arrivée des juv depuis les nour => on les sommes ds chaque region
          #####################
          
          if(data$age_sa==data$age_rec){ #need to compute a number on nursery ground, because we did not enter the case before
            for(i in seq_along(data$nurseries)){
              
              numerateur<- 4*data$h[y,i]*data$p_eggs_surv*ifelse(y%in%c(1:data$age_rec),data$Egg_0,sum(data$D_larvae[,i]*Eggs[y-data$age_rec,,s]))/(data$Wbar*(1-data$h[y,i]))
              
              denominateur <- 1+(4*data$h[y,i]*data$p_eggs_surv*ifelse(y%in%c(1:data$age_rec),data$Egg_0,sum(data$D_larvae[,i]*Eggs[y-data$age_rec,,s]))/(data$Wbar*(1-data$h[y,i])*data$K[y,i]*data$Surf_nurs[y,i]))
              
              N_nurseries[(data$age_rec-min(data$ages)+1),y,i,s] <- numerateur/denominateur
              
            } #close loop on nurseries
          } #close if age_sa=age_rec
          
          
          for (r in seq_along(data$regions)){
            
            N_juv_in <- sum(data$D_subadult[,r]*N_nurseries[a,y,,s]*
                              exp(-data$d_sa*(data$M_nurseries[a,y,]+data$F_nurseries[a,y,])/365))#*
            #exp(pe_Sjuv[a,y,,s]-data$sd_Sjuv/2))
            
            
            
            N_regions[a+1,y+1,r,s] <- N_juv_in*
              exp(-(365-data$d_sa)*(data$M_regions[a,y,r]+data$F_regions[a,y,r])/365)#*
            #exp(pe_Sad[a,y,r,s]-data$sd_Sad/2)
            
            
            #compute catches on nurseries and on regions
            C_nurseries[a,y,,s] <- N_nurseries[a,y,,s]*(data$F_nurseries[a,y,]/(data$F_nurseries[a,y,]+data$M_nurseries[a,y,]))*
              (1-exp(-data$d_sa*(data$F_nurseries[a,y,]+data$M_nurseries[a,y,])/365))*(data$CwatA[a])
            
            C_regions[a,y,r,s] <- N_juv_in* (data$F_regions[a,y,r]/(data$F_regions[a,y,r]+data$M_regions[a,y,r]))*
              (1-exp(-((365-data$d_sa)/365)*(data$F_regions[a,y,r]+data$M_regions[a,y,r])))*(data$CwatA[a])
            
            
          }#closing loop region
          
        }else{ #end of if a=age_sa & age_sa!=age_mat 
          
          ###############
          ##Dynamic of adult pop
          ###############
          
          
          #first making all fish move through a matrix from CMR connectivity model (Lecompte et al 2020)
          N_regions_af_mvt[a,y,,s] <-  apply(data$D_adult*N_regions[a,y,,s],2,sum,na.rm=T) #N_regions[a,y,,s]#
          
          
          if(a<(data$age_mat-min(data$ages)+1)){  #for non mature fish, no movement on spawning grounds
            
            N_regions[a+1,y+1,,s] <- N_regions_af_mvt[a,y,,s]*
              exp(-(data$M_regions[a,y,]+data$F_regions[a,y,]))#*
            # exp(pe_Sad[a,y,r,s]-data$sd_Sad/2)
            
            #compute catches
            C_regions[a,y,,s] <- N_regions_af_mvt[a,y,,s]* (data$F_regions[a,y,]/(data$F_regions[a,y,]+data$M_regions[a,y,]))*
              (1-exp(-(data$F_regions[a,y,]+data$M_regions[a,y,])))*(data$CwatA[a])
            
          }else{ #mature adults spend time in adult region and spawning ground
            
            #compute catches for 1st part of the year (before going to spawning grounds)
            C_regions_bf[a,y,,s] <- N_regions_af_mvt[a,y,,s]* (data$F_regions[a,y,]/(data$F_regions[a,y,]+data$M_regions[a,y,]))*
              (1-exp(-(data$d_sp-data$ts_sp/2)*(data$F_regions[a,y,]+data$M_regions[a,y,])/365))*(data$CwatA[a])
            
            
            if(a==(data$age_sa-min(data$ages)+1)&data$age_sa==data$age_mat){ #special case for those fish
              ##arrivée des juv depuis les nour => on les sommes ds chaque region
              
              if(data$age_sa==data$age_rec){ #need to compute a number on nursery ground, because we did not enter the case before
                for(i in seq_along(data$nurseries)){
                  
                  numerateur<- 4*data$h[y,i]*data$p_eggs_surv*ifelse(y%in%c(1:data$age_rec),data$Egg_0,sum(data$D_larvae[,i]*Eggs[y-data$age_rec,,s]))/(data$Wbar*(1-data$h[y,i]))
                  
                  denominateur <- 1+(4*data$h[y,i]*data$p_eggs_surv*ifelse(y%in%c(1:data$age_rec),data$Egg_0,sum(data$D_larvae[,i]*Eggs[y-data$age_rec,,s]))/(data$Wbar*(1-data$h[y,i])*data$K[y,i]*data$Surf_nurs[y,i]))
                  
                  N_nurseries[(data$age_rec-min(data$ages)+1),y,i,s] <- numerateur/denominateur
                  
                }
              }
              
              
              for (r in seq_along(data$regions)){
                
                N_juv_in <- sum(data$D_subadult[,r]*N_nurseries[a,y,,s]*
                                  exp(-data$d_sa*(data$M_nurseries[a,y,]+data$F_nurseries[a,y,])/365))#*
                #exp(pe_Sjuv[a,y,,s]-data$sd_Sjuv/2))
                
                
                
                # N_regions[a+1,y+1,r,s] <- N_juv_in*
                #   exp(-(365-data$d_sa)*(data$M_regions[a,y,r]+data$F_regions[a,y,r])/365)#*
                #exp(pe_Sad[a,y,r,s]-data$sd_Sad/2)
                
                
                #compute catches on nurseries and on regions
                C_nurseries[a,y,,s] <- N_nurseries[a,y,,s]*(data$F_nurseries[a,y,]/(data$F_nurseries[a,y,]+data$M_nurseries[a,y,]))*
                  (1-exp(-data$d_sa*(data$F_nurseries[a,y,]+data$M_nurseries[a,y,])/365))*(data$CwatA[a])
                
                # C_regions[a,y,r,s] <- N_juv_in* (data$F_regions[a,y,r]/(data$F_regions[a,y,r]+data$M_regions[a,y,r]))*
                #   (1-exp(-((365-data$d_sa)/365)*(data$F_regions[a,y,r]+data$M_regions[a,y,r])))*(data$CwatA[a])
                
                # on ajoute ces poissons subadult qui arrivent à N_region_af_mvt  => ils bougeront donc sur les spw ground
                N_regions_af_mvt[a,y,r,s] <-  N_juv_in
                
              }#closing loop region
              
            } #closing case age_sa=age_mat
            
            
            for(j in seq_along(data$spwn_grd)){
              
              #Compute the number of fish arriving on a spawning ground (summing accross dispersion matrix)
              # has suffered Znsp before moving to spw grounds
              
              if(j>1){ #pas hyper propre mais je ne sais pas comment m'en sortir quand on a pas de multiple j et r
                N_spwn_grd_bf[(data$age_mat-min(data$ages)+1):(max(data$ages)-min(data$ages)+1),y,j,s]  <- apply(apply(N_regions_af_mvt[(data$age_mat-min(data$ages)+1):(max(data$ages)-min(data$ages)+1),y,,s]*
                                                                                                                         exp(-(data$d_sp-data$ts_sp/2)*(data$M_regions[(data$age_mat-min(data$ages)+1):(max(data$ages)-min(data$ages)+1),y,]+data$F_regions[(data$age_mat-min(data$ages)+1):(max(data$ages)-min(data$ages)+1),y,])/365),1,
                                                                                                                       function(x){x*data$D_spawner[,j]}),2,sum)
                
              } else {
                N_spwn_grd_bf[(data$age_mat-min(data$ages)+1):(max(data$ages)-min(data$ages)+1),y,j,s]  <- N_regions_af_mvt[(data$age_mat-min(data$ages)+1):(max(data$ages)-min(data$ages)+1),y,,s]*
                  exp(-(data$d_sp-data$ts_sp/2)*(data$M_regions[(data$age_mat-min(data$ages)+1):(max(data$ages)-min(data$ages)+1),y,]+data$F_regions[(data$age_mat-min(data$ages)+1):(max(data$ages)-min(data$ages)+1),y,])/365)
              }
              
              
              
              if(a==(max(data$ages)-min(data$ages)+1)){
                #compute nb of eggs produces - hypo: computation at mid time ie suffered half Zsp
                Eggs[y,j,s] <- sum(N_spwn_grd_bf[,y,j,s]*
                                     exp(-data$ts_sp/(2*365)*(data$F_spwn_grd[,y,j]+data$M_spwn_grd[,y,j]))*
                                     data$fec*data$Pf*data$mat,na.rm = T)*data$p_eggs_surv[j]
                
                
                #compute SSB - hypo: computation at mid time ie suffered half Zsp
                SSB[y,j,s] <- sum(N_spwn_grd_bf[,y,j,s]*data$SwatA*mat*
                                    exp(-data$ts_sp/(2*365)*(data$F_spwn_grd[,y,j]+data$M_spwn_grd[,y,j])),na.rm = T)
              }
              
              #Compute the number of fish at the end of spawning season ie ready to move back to regions
              # has suffered all Zsp
              N_spwn_grd_end[a,y,j,s] <- N_spwn_grd_bf[a,y,j,s]*
                exp(-data$ts_sp/365*(data$F_spwn_grd[a,y,j]+data$M_spwn_grd[a,y,j])) 
              
              #compute catches
              C_spwn_grd[a,y,j,s] <- N_spwn_grd_bf[a,y,j,s]*  
                ((data$F_spwn_grd[a,y,j])/(data$F_spwn_grd[a,y,j]+data$M_spwn_grd[a,y,j]))*
                (1-exp(-(data$ts_sp/365)*(data$F_spwn_grd[a,y,j]+data$M_spwn_grd[a,y,j]))) *(data$CwatA[a])
              
              
            } #end of loop on spawning grounds
            
            for (r in seq_along(data$regions)){ #mature fish finish the year on adult feeding region
              
              if(a<(max(data$ages)-min(data$ages)+1)){ #for adult fish < age max
                #compute numbers of fish coming back to region after spawning
                N_regions[a+1,y+1,r,s] <- sum(data$D_return[,r]*N_spwn_grd_end[a,y,,s],na.rm=T)*
                  exp(-(data$M_regions[a,y,r]+data$F_regions[a,y,r])*((365-(data$d_sp+data$ts_sp/2))/365))#*
                #exp(pe_Sad[a,y,r,s]-data$sd_Sad/2)
                
                #compute catches
                C_regions[a,y,r,s] <- C_regions_bf[a,y,r,s]+ #some catches has already happened before spawning
                  sum(data$D_return[,r]*N_spwn_grd_end[a,y,r,s],na.rm=T)* #take fish who came back
                  (data$F_regions[a,y,r]/(data$F_regions[a,y,r]+data$M_regions[a,y,r]))*
                  (1-exp(-((365-(data$d_sp+data$ts_sp/2))/365)*(data$F_regions[a,y,r]+data$M_regions[a,y,r])))*(data$CwatA[a])
                
              }else{#close if in adult age < age max
                
                if(age.plus){# & a==(max(data$ages)-min(data$ages)+1)){  # if an age + should be computed 
                  #you just add the plus grup abundance to the last age
                  
                  N_plusgrp <-   sum(data$D_return[,r]*N_spwn_grd_end[a,y,,s],na.rm=T)*
                    exp(-(data$M_regions[a,y,r]+data$F_regions[a,y,r])*((365-(data$d_sp+data$ts_sp/2))/365))#*
                  #   #exp(pe_Sad[a,y,r,s]-data$sd_Sad/2)
                  
                  N_regions[a,y+1,r,s] <-  N_regions[a,y+1,r,s]+ N_plusgrp
                } # closing if age.plus 
                
                C_regions[a,y,r,s] <- C_regions_bf[a,y,r,s]+ #some catches has already happened before spawning
                  sum(data$D_return[,r]*N_spwn_grd_end[a,y,r,s],na.rm=T)* #take fish who came back
                  (data$F_regions[a,y,r]/(data$F_regions[a,y,r]+data$M_regions[a,y,r]))*
                  (1-exp(-((365-(data$d_sp+data$ts_sp/2))/365)*(data$F_regions[a,y,r]+data$M_regions[a,y,r])))*(data$CwatA[a])
              } #close else of a=age max
              
              
            }#closing loop on regions
          } #closing else on mature age 
          # } #closing else on subadult age 
        } #closing else on adult age
      } #closing loop on age
      
      #computing the sum of catches
      C_tot[y,s] <- sum(C_regions[,y,,s],na.rm=T)+sum(C_nurseries[,y,,s],na.rm=T)+sum(C_spwn_grd[,y,,s],na.rm=T)
      
    } #closing loop on year
  } #closing loop on traj
  
  
  #compute biomass
  # B_nurseries <- N_nurseries*data$SwatA
  # B_regions <- N_regions*data$SwatA
  # B_spwn_grd_bf <- N_spwn_grd_bf*data$SwatA
  # B_spwn_grd_end <- N_spwn_grd_end*data$SwatA
  
  
  res <- list(N_nurseries=N_nurseries,
              N_regions=N_regions,
              N_regions_af_mvt=N_regions_af_mvt,
              N_spwn_grd_bf=N_spwn_grd_bf,
              N_spwn_grd_end=N_spwn_grd_end,
              # B_nurseries=B_nurseries,
              # B_regions=B_regions,
              # B_spwn_grd_bf=B_spwn_grd_bf,
              # B_spwn_grd_end=B_spwn_grd_end,
              C_nurseries=C_nurseries,
              C_regions=C_regions,
              C_spwn_grd=C_spwn_grd,
              Larvae=L,
              Eggs=Eggs,
              SSB=SSB,
              C_tot =C_tot)#,
  # num=numerateur,
  # deno=denominateur)
  
  
  return(res)
  
  
  
} 
