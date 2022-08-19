# Quantile ratios T_{n;\alpha}/T_{n;\alpha_0} and T_{\infty;\alpha}/T_{n;\alpha}

source("src/asymptotic_distribution.R")

# Returns the quantile ratios T_{n;\alpha}/T_{n;\alpha_ref}
get_quantile_ratios <- function(quantiles, alpha_ref = 0.1){
  
  N <- nrow(quantiles)
  N.list <- rownames(quantiles)# List of n sample sizes
  alpha <- as.numeric(colnames(quantiles))# alpha-quantiles levels
  which.alphas <- !(alpha %in% alpha_ref)# indicator of alpha other than reference
  
  # Obtain ratios between alpha-quantiles and alpha_ref-quantiles
  quantile_ratios <- matrix(nrow = N, ncol = sum(which.alphas),
                            dimnames = list(N.list, alpha[which.alphas])
  )
  
  for (n in rownames(quantile_ratios)){
    
    for (q in colnames(quantile_ratios)){
      
      quantile_ratios[n, q] <- quantiles[n, q]/quantiles[n, as.character(alpha_ref)]
    
    }
    
  }
  
  return(quantile_ratios)
  
}

# Returns the ratios T_{\infty;\alpha}/T_{n;\alpha} for statistic.
# p indicates the dimension of the ambient space R^p.
# compute.asymp indicates to compute the asymptotic cdf if TRUE,
#   otherwise it contains the asymptotic quantiles
# return.T.Tn indicates whether to return T and T_n values (TRUE)
#   or T/T_n ratios (FALSE)
get_T_Tn <- function(quantiles, statistic, p, compute.asymp = TRUE, 
                     T_infty = seq(0, 9, 0.0001), tol = 1e-6,
                     return.T.Tn = FALSE){
  
  N <- nrow(quantiles)
  N.list <- rownames(quantiles) # List of n sample sizes
  alpha <- as.numeric(colnames(quantiles)) # alpha-quantiles levels
  
  if (isTRUE(compute.asymp)){
    
    # Asymptotic cdf computation
    P_T_infty <- asymptotic_distribution(statistic, p, T_infty, tol)
    
  }else{
    
    asymp.quantiles <- compute.asymp
    
  }
  
  # Obtain ratios T_infty/T_n as a function of n
  T_Tn <- matrix(nrow=N, ncol=length(alpha),
                 dimnames = list(N.list, alpha[1:length(alpha)])
  )
  if (isTRUE(compute.asymp)){
    
    asymp.quantiles <- matrix(nrow = N, ncol = length(alpha),
                              dimnames = list(N.list, alpha[1:length(alpha)]))
    
    for (q in 1:length(alpha)){
      
      asymp.quantiles[,q] <- T_infty[which.max(P_T_infty <= alpha[q])]
      T_Tn[,q] <- asymp.quantiles[,q] / quantiles[,q]
      
    }
    
  } else {
    
    for (q in 1:length(alpha)){
      
      T_Tn[,q] <- asymp.quantiles[q] / quantiles[,q]
      
    }
    
  }
  
  # Return the ratios or the quantiles (both asymptotic and exact)
  if(return.T.Tn){
    
    return(list("T_infty" = asymp.quantiles, "T_n" = quantiles))
    
  }else{
    
    return(T_Tn)
    
  }
}
