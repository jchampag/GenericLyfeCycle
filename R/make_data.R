##################################################
#Script: from loaded data, format them according to simulation needs


# Juliette Champagnat
# started 01/07/20
# 17/08/22 upload and adapt for github 
##########################################

#creation of a data list
mydata <- list()


###### Number and names of spatial blocs #######
mydata$regions <- regions 
mydata$regions_names <- regions_names

mydata$spwn_grd <- spwn_grd
mydata$spwn_grd_names <- spwn_grd_names

mydata$nurseries <- nurseries
mydata$nurseries_names <- nurseries_names


###### ages ########
mydata$ages <- ages
mydata$age_sa <- age_sa #age when subadults leave nurseries
mydata$age_rec <- age_rec
mydata$age_mat <- age_mat #age of sexual maturity
mydata$age_fullselec <- age_fullselec #age of full selectivity

##### years of simulation ########
mydata$years <- years


######### Population parameters #############
mydata$CwatA <- CWatAval #weight at age
mydata$SwatA <- SWatAval #weight at age
mydata$fec <- fec #fecundity at age
mydata$mat <- mat #fecundity at age
mydata$Pf <- Pf # female proportion at age
mydata$p_eggs_surv <- p_eggs_surv #spwn grd specific egg survival


##### nursery parameters #######
if(!is.null(Mlarvae)){mydata$Mlarvae <- Mlarvae}
mydata$d_sa <- d_sa #mean date of subadult leaving nurseries for adult regions

mydata$Surf_nurs <- array(rep(Surf,each=length(years)), dim=c(length(years),length(nurseries)),
                          dimnames=list(years,nurseries_names))

if(!is.null(alpha)){ mydata$alpha <- array(rep(alpha,each=length(years)), dim=c(length(years),length(nurseries)),
                                           dimnames=list(years,nurseries_names))}

mydata$K <- array(rep(K,each=length(years)), dim=c(length(years),length(nurseries)),
                  dimnames=list(years,nurseries_names))

mydata$h <- array(rep(h,each=length(years)), dim=c(length(years),length(nurseries)),
                  dimnames=list(years,nurseries_names))

##### spawning ground parameters ######
mydata$d_sp <- d_sp #mean reproduction date 
mydata$ts_sp <- ts_sp #mean reproduction duration


####### Fishing mortality #########
mydata$F_nurseries <- array(0,dim=c(length(ages),length(years),length(nurseries)),dimnames = list(ages,years,nurseries_names))#
mydata$F_nurseries[1:age_sa-min(ages)+1,,] <- Fval[1:age_sa-min(ages)+1]

if(DYN){
  for(a in 1:age_sa-min(ages)+1){
    if(age_rec==0){
      mydata$F_nurseries[a,,] <- FatA %>% filter(Age==0) %>% pull(F)
    }
    mydata$F_nurseries[a,,] <- FatA %>% filter(Age==a) %>% pull(F)
  }
}

mydata$F_regions <- array(0,dim=c(length(ages),length(years),length(regions)),dimnames = list(ages,years,regions_names))#

#taking Fval for ages present in adult regions
#/!\ Fval is now considered as a vector, should be modified when simulating with multiple region/spawning ground
for(r in seq_along(regions)){
  mydata$F_regions[(age_sa-min(ages)+1):(max(ages)-min(ages)+1),,r] <- Fval[(age_sa-min(ages)+1):(max(ages)-min(ages)+1)]
}


mydata$F_spwn_grd <- array(0,dim=c(length(ages),length(years),length(spwn_grd)),dimnames = list(ages,years,spwn_grd_names))#

for(s in seq_along(spwn_grd)){
  mydata$F_spwn_grd[(age_mat-min(ages)+1):(max(ages)-min(ages)+1),,s] <- Fval[(age_mat-min(ages)+1):(max(ages)-min(ages)+1)]
}


if(DYN){
  mydata$F_regions[,,1]<- as.matrix(FatA %>% pivot_wider(names_from=Year,values_from=F) %>% select(-Age))
  mydata$F_regions[1,,] <-0
  
  mydata$F_spwn_grd <-  mydata$F_regions
  mydata$F_spwn_grd[1:(age_mat-min(ages)),,] <- 0
}

########### Natural Mortality ##########
##on nurseries
mydata$M_nurseries <- array(0, dim=c(length(ages),length(years),length(nurseries)), dimnames = list(ages,years,nurseries))#

mydata$M_nurseries[1:(age_sa-min(ages)+1),,] <- M_at_age[1:(age_sa-min(ages)+1)] 



##on regions
mydata$M_regions <- array(0, dim=c(length(ages),length(years),length(regions)), dimnames = list(ages,years,regions))#

for (a in (age_sa-min(ages)+1):(max(ages)-min(ages)+1)){
  mydata$M_regions[a,,] <- M_at_age[a]}


#on spawning ground
mydata$M_spwn_grd <- mydata$M_regions


###### initialization values #########
mydata$Na_year1 <- NatA_0
mydata$Egg_0 <- Egg_0
if(!is.null(Mlarvae)){mydata$Larvae_year1 <- Egg_0*exp(-Mlarvae)}

####### dispersion matrices #########
# these matrices allocate individual from x departure areas to y arrival areas
mydata$D_larvae <- D_larvae # from spawning grounds to nurseries
mydata$D_subadult <- D_subadult # from nurseries to adult regions
mydata$D_spawner <- D_spawner # from adult regions to spawning grounds
mydata$D_return <- D_return # from spawning grounds to adult region
mydata$D_adult <- D_adult # connectivity between adult zones

######## standard deviation for process error #########
mydata$sd_Sad <- 0#sd_Sad
mydata$sd_Sjuv <- 0#sd_Sjuv
mydata$sd_R <- 0#sd_R


