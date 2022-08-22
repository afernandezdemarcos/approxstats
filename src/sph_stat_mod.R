library(sphunif)
library(Rcpp)
sourceCpp("./src/sph_stat_mod_newton.cpp")

sph_stat_mod <- function(type, statistic, n, alpha, p){
  
  if(type=="PCvM"){mod.stat <- sph_stat_PCvM_mod(statistic, n, alpha, p)}
  
  if(type=="PAD"){mod.stat <- sph_stat_PAD_mod(statistic, n, alpha, p)}
  
  if(type=="PAD_2"){mod.stat <- sph_stat_PAD_mod(statistic, n, alpha, p, PAD_2 = TRUE)}
  
  if(type=="Bakshaev"){mod.stat <- sph_stat_Bakshaev_mod(statistic, n, alpha, p)}
  
  return(mod.stat)
  
}