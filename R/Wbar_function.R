##################################################
#Script: function which computes Wbar and adds it to mydata() list


# Juliette Champagnat
# started sometimes in early 2021
# 17/08/22 upload and adapt for github 
##################################################

"Wbar_computing" <- function(a,data=mydata,output = "W_a"){
  
  if(grepl(pattern="/R",getwd())){source("../survival_function.R",local=T)
  }else{source("R/survival_function.R",local=T)}
  
 S_a <- survival_computing(a=a,data=mydata,F=0,first_age = 1)
  
  #puis le nombre d'oeuf moyen
  W_a = S_a*data$Pf[a-min(ages)+1]*data$mat[a-min(ages)+1]*data$fec[a-min(ages)+1] ##a+1 car ces vecteurs commencent à l'age 0
  
  if(output == "S_a") {return(S_a)}
  else{return(W_a)}
}

mydata$Wbar <- sum(map_dbl(mydata$ages,Wbar_computing),na.rm=T)

