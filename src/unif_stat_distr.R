# (Algorithm 1) p-value approximation using (n, alpha)-stabilization
# Wrapper for arbitrary length of statistic values

library(Rcpp)
sourceCpp("src/cir_stat_mod_newton.cpp")
sourceCpp("src/sph_stat_mod_newton.cpp")
source("src/approxstats_utils.R")

# Asymptotic quantiles of hyperspherical statistics 
# (Interpolate non available quantiles for other dimensions p)
sph_asymp_q <- list("PCvM" = asymptotic_quantiles_table("PCvM"),
                    "PAD" = asymptotic_quantiles_table("PAD"),
                    "Bakshaev" = asymptotic_quantiles_table("Bakshaev"))

p_value_grid <- function(T_n, statistic, n, p = 2){
  
  # Set arguments for distribution and arguments name
  
  KS <- FALSE
  CvM <- FALSE
  PAD_2 <- FALSE
  
  if (statistic == "KS") {
    statistic <- "Kuiper" 
    KS <- TRUE
  }
  if (statistic == "CvM") {
    statistic <- "Watson" 
    CvM <- TRUE
  }
  
  if (statistic == "PAD_2") {
    statistic <- "PAD" 
    PAD_2 <- TRUE
  }
  
  if (statistic %in% avail_cir_tests){
    prefix <- "p_cir_stat_"
  } else {
    prefix <- "p_sph_stat_"
  }
  
  args <- list("KS" = KS, "CvM" = CvM, "PAD_2" = PAD_2, "n" = n, "p" = p)

  m <- if (length(T_n) > 1) {q <- length(T_n)} else {q <- 1}
  p.val <- numeric(m)
  
  # Hyperspherical statistics: Pre-compute \mathcal{T}_{\infinity} set
  if (statistic %in% avail_sph_tests){
    
    alpha.grid <- seq(0.001, 0.25, 0.001)
    
    stat_sph_asymp_q <- sph_asymp_q[[statistic]]
    avail.p <- as.numeric(colnames(stat_sph_asymp_q))
    
    # If p is not available, then interpolate from adjacent p's
    if(p %in% avail.p){
      
      T_inf <- stat_sph_asymp_q[1:length(alpha.grid), as.character(p)]
  
    }else{
      
      avail.alpha <- rownames(stat_sph_asymp_q)
      T_inf <- numeric(length(avail.alpha))
      names(T_inf) <- avail.alpha
      idx.upp <- which.min(p > avail.p)
      p.low <- avail.p[idx.upp - 1]
      p.upp <- avail.p[idx.upp]
      
      for(alpha in avail.alpha){
        T_inf[alpha] <- approx(c(p.low, p.upp), 
                               c(stat_sph_asymp_q[alpha, as.character(p.low)], stat_sph_asymp_q[alpha, as.character(p.upp)]), p)$y
      }
      T_inf <- T_inf[as.character(alpha.grid)]
      
    }
    
    args <- c(list("T_inf" = T_inf), args)
    
  }
  
  name_distr <- paste0(prefix, statistic)
  names_args <- names(args)
  
  distr_args <- args[names_args %in% names(formals(name_distr))]

  # Compute p-value grid approximation for each statistic value
  for(i in 1:q){
    p.val[i] <- do.call(what = name_distr, 
                        args = c(list("T_n" = T_n[i]), distr_args))
  }

  return(p.val)
  
}
