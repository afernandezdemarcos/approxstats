# Weight functions

# Asymptotic Monte Carlo sampling variance
AVar <- function(alpha, T_inf, T_n, f_Tn, B = as.numeric(M), ...){
  ((T_inf^2*alpha*(1-alpha))/(B*f_Tn^2*T_n^4))^(-1/2)
}

# Upper-tail quantiles
weight.fun.1 <- function(n, alpha, ...){as.numeric(alpha<=0.25)}

# Upper-tail focused on small sample sizes
weight.fun.2 <- function(n, alpha, ...){n^(-1/2)*(alpha<=0.25)}

# Upper-tail balanced with Monte Carlo AVar
weight.fun.3 <- function(n, alpha, T_inf, T_n, f_Tn, B = as.numeric(M)){
  (alpha<=0.25) * AVar(alpha = alpha, T_inf = T_inf, 
                       T_n = T_n, f_Tn = f_Tn, B = B)
}

# Upper-tail focused on small sample sizes balanced with Monte Carlo AVar
weight.fun.4 <- function(n, alpha, T_inf, T_n, f_Tn, B = as.numeric(M)){
  (n)^(-1/2) * (alpha<=0.25) * AVar(alpha = alpha, T_inf = T_inf, 
                                    T_n = T_n, f_Tn = f_Tn, B = B)
}

# All quantiles balanced with Monte Carlo AVar
weight.fun.5 <- function(alpha, T_inf, T_n, f_Tn, B = as.numeric(M), ...){
  AVar(alpha = alpha, T_inf = T_inf, T_n = T_n, f_Tn = f_Tn, B = B)
}

# Asymptotic variance focused on small sample sizes
weight.fun.6 <- function(n, alpha, T_inf, T_n, f_Tn, B = as.numeric(M)){
  (n)^(-1/2) * AVar(alpha = alpha, T_inf = T_inf, 
                    T_n = T_n, f_Tn = f_Tn, B = B)
}

# All quantiles weighted by the inverse of alpha focused on small sample sizes
weight.fun.7 <- function(n, alpha, ...){(n*alpha)^(-1/2)}
