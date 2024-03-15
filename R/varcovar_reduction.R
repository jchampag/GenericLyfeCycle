##########################################################
#Script : manipulate parameters varcovar matrices from FishLife to reduce uncertainty
# uses the matrix of 1 species which has smaller parameters uncertainty (here plaice) to 
# rescale other species parameter uncertainty (here amchovy and cuckoo ray) so that for a
# a given parameter, each species have the same relative level of uncertainty (e.g., same CV)

# Juliette Champagnat
###########################################################

Predict_ray = Plot_taxa( Search_species(Genus="Leucoraja",Species= "naevus")$match_taxonomy,
                         mfrow=c(2,3),Database=FishLife::FishBase_and_RAM)

Predict_anc = Plot_taxa( Search_species(Genus="Engraulis",Species= "encrasicolus")$match_taxonomy,
                         mfrow=c(2,3),Database=FishLife::FishBase_and_RAM)

Predict_ple = Plot_taxa( Search_species(Genus="Pleuronectes",Species= "platessa")$match_taxonomy,
                         mfrow=c(2,3),Database=FishLife::FishBase_and_RAM)


## reduce the matrix to parameters of interest only
varcovar_ray <-Predict_ray[[1]]$Cov_pred[rownames(Predict_ray[[1]]$Cov_pred)%in%c("Loo","K","M","h"),#"tmax","tm"
                                         colnames(Predict_ray[[1]]$Cov_pred)%in%c("Loo","K","M","h") ]

varcovar_ple <-Predict_ple[[1]]$Cov_pred[rownames(Predict_ple[[1]]$Cov_pred)%in%c("Loo","K","M","h"),
                                         colnames(Predict_ple[[1]]$Cov_pred)%in%c("Loo","K","M","h") ]


varcovar_anc <-Predict_anc[[1]]$Cov_pred[rownames(Predict_anc[[1]]$Cov_pred)%in%c("Loo","K","M","h"),
                                         colnames(Predict_anc[[1]]$Cov_pred)%in%c("Loo","K","M","h") ]

#### REDUCTION /!\ differs for parameters in log space (M,Loo et k)

# compute initial sigma value
sigma_ple <- sqrt(diag(varcovar_ple))
sigma_ray <- sqrt(diag(varcovar_ray))
sigma_anc <- sqrt(diag(varcovar_anc))

# compute plaice CV
CV_ple <- sigma_ple/abs(Predict_ple[[1]]$Mean_pred[c(1,2,6,13)]) #4,5,
CV_ray <- sigma_ray/abs(Predict_ray[[1]]$Mean_pred[c(1,2,6,13)])
CV_anc <- sigma_anc/abs(Predict_anc[[1]]$Mean_pred[c(1,2,6,13)])


# plot CVs
data.frame(bind_rows(CV_ray,CV_ple,CV_anc)) %>%
  mutate(sp=c("ray","plaice","anchovy")) %>%
  pivot_longer(cols=-sp) %>%
  # filter(name=="Loo") %>% 
  ggplot()+geom_point(aes(name,value,col=sp),size=5,shape=4)+
  geom_hline(yintercept=0)+
  labs(x="Parameters",y="CV")

# compute the lambda/reduction factor
red_diag_ray <- red_diag_anc <- NULL
## for h
red_diag_ray[4] <- CV_ple[4]/CV_ray[4]#Predict_ray[[1]]$Mean_pred[c(13)]*CV_ple[4]/sigma_ray[4]
red_diag_anc[4] <-  CV_ple[4]/CV_anc[4]#Predict_anc[[1]]$Mean_pred[c(13)]*CV_anc[4]/sigma_anc[4]

#for LogM logMoo log k
red_diag_ray[c(1:3)] <- sigma_ple[c(1:3)]/sigma_ray[c(1:3)]
red_diag_anc[c(1:3)] <- sigma_ple[c(1:3)]/sigma_anc[c(1:3)]

# compute new sigma to have 
new_sigma_ray <- red_diag_ray*sigma_ray
new_sigma_anc <- red_diag_anc*sigma_anc
# reset Loo sigma to original value
# new_sigma_anc[1:2] <- sigma_anc[1:2]


# built new varcovar matrix
new_varcovar_anc <- new_varcovar_ray <- temp_ray <- temp_anc <- matrix(0,4,4,dimnames=dimnames(varcovar_ray))#varcovar_ray

# compute the reduction factor used on the upper triangle
red_up_tri_ray <- as.numeric(c(red_diag_ray[1]*red_diag_ray[2],red_diag_ray[1]*red_diag_ray[3],red_diag_ray[2]*red_diag_ray[3],
                               red_diag_ray[1]*red_diag_ray[4],red_diag_ray[2]*red_diag_ray[4],
                               red_diag_ray[3]*red_diag_ray[4]))

red_up_tri_anc <- as.numeric(c(red_diag_anc[1]*red_diag_anc[2],red_diag_anc[1]*red_diag_anc[3],red_diag_anc[2]*red_diag_anc[3],
                               red_diag_anc[1]*red_diag_anc[4],red_diag_anc[2]*red_diag_anc[4],
                               red_diag_anc[3]*red_diag_anc[4]))

# fill a symetrical matrice
temp_ray[upper.tri(temp_ray)] <- varcovar_ray[upper.tri(varcovar_ray)]*red_up_tri_ray
temp_anc[upper.tri(temp_anc)] <- varcovar_anc[upper.tri(varcovar_anc)]*red_up_tri_anc

# build the new matrix
new_varcovar_ray <- temp_ray+t(temp_ray)+diag(new_sigma_ray*new_sigma_ray, 4) #varcovar_ray;new_varcovar_ray
new_varcovar_anc <- temp_anc+t(temp_anc)+diag(new_sigma_anc*new_sigma_anc, 4) #varcovar_anc;new_varcovar_anc


### plot old and new parameters distribution
# nb_tirage <- 1000
# 
# predictions_all <- MASS::mvrnorm(n=nb_tirage, mu=Predict_ray[[1]]$Mean_pred[c(1,2,6,13)],Sigma = varcovar_ray) %>% as.data.frame() %>%
#   mutate(type="ray_raw") %>%
#   bind_rows(MASS::mvrnorm(n=nb_tirage, mu=Predict_ray[[1]]$Mean_pred[c(1,2,6,13)],Sigma = new_varcovar_ray) %>% as.data.frame() %>%
#               mutate(type="ray_reduced")) %>%
#   bind_rows(MASS::mvrnorm(n=nb_tirage, mu=Predict_ple[[1]]$Mean_pred[c(1,2,6,13)],Sigma = varcovar_ple) %>% as.data.frame() %>%
#               mutate(type="plaice_raw")) %>%
#   bind_rows(MASS::mvrnorm(n=nb_tirage, mu=Predict_anc[[1]]$Mean_pred[c(1,2,6,13)],Sigma = varcovar_anc) %>% as.data.frame() %>%
#               mutate(type="anc_raw")) %>%
#   bind_rows(MASS::mvrnorm(n=nb_tirage, mu=Predict_anc[[1]]$Mean_pred[c(1,2,6,13)],Sigma = new_varcovar_anc) %>% as.data.frame() %>%
#               mutate(type="anc_reduced"))
# 
# predictions_all %>%
#   mutate(across(.cols=!c(h,type),.fns =exp )) %>%
#   pivot_longer(cols=-type) %>%
#   # filter(type!="anc_reduced") %>%
#   # filter(name=="M") %>%
#   ggplot()+geom_density(aes(y=value,fill=type),alpha=0.5)+facet_wrap(~name,scale="free")#+
#   # coord_cartesian(xlim=c(0,20))
# 