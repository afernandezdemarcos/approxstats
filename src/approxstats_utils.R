# Utils

library(stringr)
library(bootstrap)

avail_cir_tests <- c("KS", "Kuiper", "CvM", "Watson", "AD")
avail_sph_tests <- c("PAD", "PCvM", "Bakshaev", "PAD_2")

get_stepAIC_penalty <- function(k, data){
  # Penalty measure
  penalties <- list("AIC" = 2, "BIC" = log(nrow(data)))
  split.k <- strsplit(k, "\\*")[[1]]
  if(length(split.k) > 1){
    k <- as.numeric(split.k[1])*penalties[[split.k[2]]]
  }else{
    k <- penalties[[k]]
  }
  return(k)
}

get_stepAIC_direction <- function(fit.direction){
  # Fit direction
  fit.directions <- c("forward", "both")
  if(length(fit.direction) == 1){
    if(fit.direction == "all"){
      fit.direction <- fit.directions
    }
  }
  return(fit.direction)
}

# Compute weights for (n, alpha, T_inf, T_n, f_Tn, p) applying weight.fun
get_weights_from_function_n <- function(weights.fun, n_vector, alpha_vector, 
                                        T_inf_vector, T_n_vector, f_Tn_vector,
                                        p_vector = NULL){
  if(!is.null(weights.fun)){
    
    if(typeof(weights.fun) == "closure"){
      
      weights <- if(is.null(p_vector)) {
        
        weights.fun(n = n_vector, 
                    alpha = alpha_vector, 
                    T_inf = T_inf_vector, 
                    T_n = T_n_vector, 
                    f_Tn = f_Tn_vector)
        
      }else{
        
        weights.fun(n = n_vector, 
                    alpha = alpha_vector, 
                    p_vector)
        
      }
      
    }else{
      
      stop("weights.fun must be a function of n")
      
    }
    
  }else{
    
    weights <- weights.fun
    
  }
  
  return(weights)
  
}

# Compute pdf from sample quantiles
compute_pdf <- function(T_n){
  
  alpha_vector <- as.numeric(colnames(T_n))
  n_vector <- as.numeric(row.names(T_n))
  
  pdf_Tn <- matrix(nrow = length(n_vector), ncol = length(alpha_vector), 
                   dimnames = list(n_vector, alpha_vector))
  
  for (i in 1:length(n_vector)){

    cdf_spline <- splinefun(x = T_n[i,], y = 1 - alpha_vector)
    # Interpolate cdf with splines
    cdf <- function(x) {
      
      pmin(pmax(cdf_spline(x), 0), 1)
      
    }
    
    # Numerical differentiation
    pdf <- function(x) {
      
      pmax(numDeriv::grad(func = cdf, x = x), 0)
      
    }
    pdf_Tn[i, ] <- pdf(T_n[i,])

  }
  
  return(pdf_Tn)
  
}

get_N_index <- function(N.list, N.start, N.end = NULL){
  
  N.list.num <- as.numeric(N.list) # Numeric list
  if(is.null(N.end)){N.end <- max(N.list.num)} # Check NULL of N
  
  # Sample sizes greater than selected minimum
  N.mask <- (N.list.num >= N.start) & (N.list.num <= N.end)
  N.index <- N.list[N.mask] # Index from ratios
  N.num <- as.numeric(N.index) # x axis
  
  return(list(index = N.index, num = N.num, end = N.end))
}

get_formula_lm <- function(model, intercept = TRUE, precision = 4){
  if(intercept){
    form <- paste0(round(coefficients(model)[1],2), " + ", 
                   paste(sprintf(paste0("%.", precision,"f * %s"), 
                                 coefficients(model)[-1],  
                                 names(coefficients(model)[-1])), 
                         collapse=" + "))
  }else{
    form <- paste0(round(coefficients(model)[1],precision), " + ", 
                   paste(sprintf(paste0("%.", precision,"f * %s"), 
                                 coefficients(model)[-1],  
                                 names(coefficients(model)[-1])), 
                         collapse=" + "))
  }
  
  return(form)
}

subsetList <- function(myList, elementNames) {
  lapply(elementNames, FUN=function(x) myList[[x]])
}

process_MC_results <- function(results, B = 100, ci.limits = c(0.025, 0.975)){
  
  results.dim <- dim(results)
  results.dimnames <- dimnames(results)
  
  cir <- length(results.dim) == 3
  
  alpha.list <- results.dimnames[[3]]
  n.list <- results.dimnames[[2]]

  rejection.results <- array(dim = results.dim[2:length(results.dim)], 
                             dimnames = subsetList(results.dimnames, 
                                                   2:length(results.dim)))
  ci.lower <- array(dim = results.dim[2:length(results.dim)], 
                    dimnames = subsetList(results.dimnames, 
                                          2:length(results.dim)))
  ci.upper <- array(dim = results.dim[2:length(results.dim)], 
                    dimnames = subsetList(results.dimnames, 
                                          2:length(results.dim)))
  
  if(cir){
    
    for (alpha in alpha.list){
      
      rejection.results[, alpha] <- apply(results[, , alpha], 2, mean)
      
      for(n in n.list){
        
        rej.star <- bootstrap(results[, n, alpha], B, mean)$thetastar
        ci <- quantile(rej.star, ci.limits)
        ci.lower[n, alpha] <- ci[1]
        ci.upper[n, alpha] <- ci[2]
        
      }
      
    }
    
  } else {
    
    p.list <- results.dimnames[[4]]
    
    for (p in p.list){
      
      for (alpha in alpha.list){
        rejection.results[, alpha, p] <- apply(results[, , alpha, p], 2, mean)

      }
      
    }
    
  }
  
  return(list(rejection.results = rejection.results, 
              ci.lower = ci.lower, ci.upper = ci.upper))
}

asymptotic_quantiles_table <- function(statistic){
  
  if(statistic == "PAD_2"){statistic <- "PAD"}
  
  e <- new.env()
  load(paste0("distributions/asymptotic/", 
              if (statistic %in% avail_sph_tests) "n500_" else "",
              "asymp_qua_", statistic, ".RData"), 
       envir = e)
  return(e$T_inf)
  
}
