##################################################
#Function: for a given age computes survival until this age

# Juliette Champagnat
##################################################

#arguments
#' @a  age until which survival is computed
#' @Fval fishing mortality value (will be multiplied with a selectivity vector)
#' @data dataset used for grabbing parameters
#' @first_age first_sge considered for survival computation, if 1 the computation starts at an age post recruitment, of 0 the survival density independant from recruitment phase is also considered
                             
"survival_computingSIMPLE" <- function(a,Fval=0,data=mydata,first_age=1){
  if(is.numeric(data$M_at_age)){
    allM <- data$M_at_age
  }else{ #if M_at_age is a matrix age*year, we average over years
    allM <- as.numeric(apply(data$M_at_age,1,mean,na.rm=T))}
  
  if(min(data$ages)==0){
    S_a = ifelse(a==0,1,exp(-sum(allM[1:a])))
  } else{
    S_a = ifelse(a==1,1,exp(-sum(allM[1:a-1])))
  }
  return(S_a)
}
  