##################################################
#Script: function which for a given a computes survival until this age


# Juliette Champagnat
# started 09/9/21
# 17/08/22 upload and adapt for github 
##################################################

#arguments
#' @a  age until with survival is computed
#' @F fishing mortality value (will be multiplied with a selectivity vector)
#' @data dataset used for grabbing parameters
#' @first_age first_sge considered for survival computation, if 1 the computation starts at an age post recruitment, of 0 the survival density independant from recruitment phase is also considered


"survival_computing" <- function(a,F,data=mydata,first_age=1){
  
  
  S_a <- exp(ifelse(a<(age_rec-min(ages)+1),0, #si on est < age rec on ne subit que la mortalité de la phase de recrutement (qui est codée en dessous)
                    
                    ifelse(a<(age_sa-min(ages)+1), #si on est inf à l'age de sortie de nourricerie, on ne subit que la mortalité nat des nourriceries
                           -sum(data$M_nurseries[1:(a-min(ages)+1),1,1]-F*Selec[1:(a-min(ages)+1)]), 
                           
                           ifelse(a==(age_sa-min(ages)+1),#si on est = à l'age de sortie des nour on subit les 2 mortalités pondérée par la durée
                                  -sum(data$M_nurseries[1:(age_sa-min(ages)),1,1]+F*Selec[1:(age_sa-min(ages))])-
                                    (d_sa*(data$M_nurseries[(age_sa-min(ages)+1),1,1]+F*Selec[(age_sa-min(ages)+1)])/365)-
                                    ((365-d_sa)/365*(data$M_regions[(age_sa-min(ages)+1),1,1]+F*Selec[(age_sa-min(ages)+1)])),
                                  
                                  #ifelse((a-min(ages)+1)<max(ages), # si on est > à l'age de sortie des nour on subit la juv, celle de transition pondéreé, et la suite
                                  -sum(data$M_nurseries[1:(age_sa-min(ages)),1,1]+F*Selec[1:(age_sa-min(ages))])-
                                    (d_sa*(data$M_nurseries[(age_sa-min(ages)+1),1,1]+F*Selec[(age_sa-min(ages)+1)])/365)-
                                    ((365-d_sa)/365*(data$M_regions[(age_sa-min(ages)+1),1,1]+F*Selec[(age_sa-min(ages)+1)]))-
                                    sum(data$M_regions[(age_sa-min(ages)+2):(a-min(ages)+1),1,1]+F*Selec[(age_sa-min(ages)+2):(a-min(ages)+1)])))))#,
  

  if(first_age==0){S_a <- S_a* (4*h/(mydata$Wbar*((1-h))))} #si on veut prendre en compte la mortalité DI de la phase de recrutement on l'ajoute a la survie
  
  
  return(S_a)
}
