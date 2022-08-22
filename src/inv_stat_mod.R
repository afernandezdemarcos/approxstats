inv_cir_stat_Kuiper_mod <- function(asymp_statistic, n, alpha, KS = FALSE, Stephens = FALSE){
  
  if (Stephens){
    if (KS){
      f <- (1 + 0.12/sqrt(n) + 0.11/n)
    }else{
      f <- (1 + 0.155/sqrt(n) + 0.24/n)
    }
  }else{
    if (KS){
      f <- (1 + 0.1575/sqrt(n) + 0.0192/(n*sqrt(alpha)) - 0.0051/sqrt(alpha*n))
    }else{
      f <- (1 + 0.2330/sqrt(n) + 0.0276/(n*sqrt(alpha)) - 0.0068/sqrt(alpha*n))
    }  
  }
  
  return(asymp_statistic / f)
  
}

inv_cir_stat_Watson_mod <- function(asymp_statistic, n, alpha, CvM = FALSE, Stephens = FALSE){
  
  if (Stephens){
    if (CvM){
      g <- -0.4/n +0.6/(n*n)
      f <- (1 + 1/n)
    }else{
      g <- -0.1/n +0.1/(n*n)
      f <- (1 + 0.8/n)
    }
  }else{
    if (CvM){
      g <- 0
      f <- (1 - 0.1651/n + 0.0749/(n*sqrt(alpha)) - 0.0014/(n*alpha))
    }else{
      g <- 0
      f <- (1 - 0.1505/n + 0.0917/(n*sqrt(alpha)) - 0.0018/(n*alpha))
    }
  }
  
  return((asymp_statistic / f) - g)
  
}

inv_cir_stat_AD_mod <- function(asymp_statistic, n, alpha, Stephens = FALSE){
  
  if (Stephens){
    f <- 1
  }else{
    f <- (1 + 0.0356/n - 0.0234/(n*sqrt(alpha)) + 0.0006/(n*alpha))
  }
  
  return(asymp_statistic / f)
  
}

inv_cir_stat_PAD_2_mod <- function(asymp_statistic, n, alpha){
  
  f <- (1 - 0.0683/n + 0.0692/(n*sqrt(alpha)) - 0.0014/(n*alpha))
  
  return(asymp_statistic / f)
  
}

inv_sph_stat_PCvM_mod <- function(asymp_statistic, n, alpha, p){
  
  f <- (1 + (+0.113/sqrt(p) -0.5415/p) / n
              + (+0.1438/sqrt(p)) / (n*sqrt(alpha))
              + (-0.0031/sqrt(p)) / (n*alpha))
  
  return(asymp_statistic / f)
}

inv_sph_stat_PAD_mod <- function(asymp_statistic, n, alpha, p){
  
  f <- (1 + (+0.0978/sqrt(p) -0.3596/p) / n
              + (0 +0.1126/sqrt(p)) / (n*sqrt(alpha))
              + (0 -0.0025/sqrt(p)) / (n*alpha))
  
  return(asymp_statistic / f)
}

inv_sph_stat_Bakshaev_mod <- function(asymp_statistic, n, alpha, p){
  
  f <- (1 + (+0.1189/sqrt(p) -0.5838/p) / n
              + (0 +0.1210/sqrt(p) +0.0385/p) / (n*sqrt(alpha))
              + (0 -0.0030/sqrt(p)) / (n*alpha)
  )
  
  return(asymp_statistic / f)
}

inv_cir_stat_mod <- function(type, asymp_statistic, n, alpha, Stephens = FALSE){
  
  if(type=="KS"){mod.stat <- inv_cir_stat_Kuiper_mod(asymp_statistic, n, alpha, 
                                                 KS = TRUE, Stephens = Stephens)}
  
  if(type=="Kuiper"){mod.stat <- inv_cir_stat_Kuiper_mod(asymp_statistic, n, alpha, 
                                                     KS = FALSE, Stephens = Stephens)}
  
  if(type=="CvM"){mod.stat <- inv_cir_stat_Watson_mod(asymp_statistic, n, alpha, 
                                                  CvM = TRUE, Stephens = Stephens)}
  
  if(type=="Watson"){mod.stat <- inv_cir_stat_Watson_mod(asymp_statistic, n, alpha, 
                                                     CvM = FALSE, Stephens = Stephens)}
  
  if(type=="AD"){mod.stat <- inv_cir_stat_AD_mod(asymp_statistic, n, alpha, 
                                             Stephens = Stephens)}
  
  if(type=="PAD"){mod.stat <- inv_cir_stat_PAD_2_mod(asymp_statistic, n, alpha)}
  
  return(mod.stat)
  
}

inv_sph_stat_mod <- function(type, asymp_statistic, n, alpha, p){
  
  if(type=="PCvM"){mod.stat <- inv_sph_stat_PCvM_mod(asymp_statistic, n, alpha, p)}
  
  if(type=="PAD"){mod.stat <- inv_sph_stat_PAD_mod(asymp_statistic, n, alpha, p)}
  
  if(type=="Bakshaev"){mod.stat <- inv_sph_stat_Bakshaev_mod(asymp_statistic, n, alpha, p)}
  
  return(mod.stat)
  
}

get_p_real <- function(q, n, crit_val_n, statistic){
  
  # Search p - values
  n_filter <- (crit_val_n$n == n)
  MC_q <- crit_val_n[n_filter, statistic]
  names(MC_q) <- crit_val_n[n_filter, "alpha"]
  MC_q <- sort(MC_q)
  alpha <- as.numeric(names(MC_q))
  idx <- findInterval(q, MC_q)
  
  # Interpolate p - value
  x <- c(MC_q[idx], MC_q[idx + 1])
  y <- c(alpha[idx], alpha[idx + 1])
  alpha_real <- approx(x, y, xout = q)$y
  
  return(alpha_real)
  
}