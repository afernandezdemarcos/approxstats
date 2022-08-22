# Circular statistic computation using (n, alpha)-stabilization

library(sphunif)
library(Rcpp)
sourceCpp("./src/cir_stat_mod_newton.cpp")

cir_stat_mod <- function(type, statistic, n, alpha, Stephens = FALSE){
  
  if(type=="KS"){mod.stat <- cir_stat_Kuiper_mod(statistic, n, alpha, 
                                                 KS = TRUE, Stephens = Stephens)}
  
  if(type=="Kuiper"){mod.stat <- cir_stat_Kuiper_mod(statistic, n, alpha, 
                                                     KS = FALSE, Stephens = Stephens)}
  
  if(type=="CvM"){mod.stat <- cir_stat_Watson_mod(statistic, n, alpha, 
                                                  CvM = TRUE, Stephens = Stephens)}
  
  if(type=="Watson"){mod.stat <- cir_stat_Watson_mod(statistic, n, alpha, 
                                                     CvM = FALSE, Stephens = Stephens)}
  
  if(type=="AD"){mod.stat <- cir_stat_AD_mod(statistic, n, alpha, 
                                                          Stephens = Stephens)}
  
  return(mod.stat)
  
}
