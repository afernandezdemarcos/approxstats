source("src/inv_stat_mod.R")
library(sphunif)

statistic_list <- c("PCvM", "Bakshaev", "PAD")
M <- "1e7"

alpha_list <- c(0.25, 0.2, 0.15, 0.1, 0.05, 0.02, 0.01)
n_list <- c(5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 100, 200, 300)
p_list <- c(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 21, 51, 101, 151)

replicates <- 1e3
M.mc <- 1e4 * replicates

crit_val_MC <- array(dim = c(length(n_list), length(alpha_list), length(p_list), replicates),
                     dimnames = list(n_list, alpha_list, p_list, NULL))

mod.results <- array(dim = c(length(n_list), length(alpha_list), length(p_list), replicates, length(statistic_list)),
                     dimnames = list(n_list, alpha_list, p_list, NULL, statistic_list))

for (p in p_list){
  
  print(paste("p = ", p))
  
  load(paste0("distributions/qua_p", p, "_M", M, ".RData"))
  
  for (n in n_list){
    
    print(paste("n = ", n))
    
    unif_stat_MC_n <- unif_stat_MC(n, type = statistic_list, p = p, return_stats = TRUE,
                                   cores = 3, M = M.mc, chunks = 10)
    
    for (statistic in statistic_list){
      print(statistic)
      
      apply_get_p_real <- function(x){
        get_p_real(q = x, n = n, 
                   crit_val_n = crit_val_n, statistic = statistic)
      }
      
      stat_MC_n <- matrix(unif_stat_MC_n$stats_MC[, statistic], ncol = replicates)
      crit_val <- rbind(apply(stat_MC_n, 2, quantile, probs = 1 - alpha_list,
                              na.rm = TRUE))
      rownames(crit_val) <- alpha_list
      
      crit_val_MC[as.character(n), ,as.character(p),] <- crit_val
      
      for (alpha in alpha_list){
        
        mod_crit_val <- crit_val_MC[as.character(n), as.character(alpha), as.character(p),]
        
        # Real MC quantile
        alpha_real <- sapply(mod_crit_val, apply_get_p_real)
        
        # Save alpha_real
        mod.results[as.character(n), as.character(alpha), as.character(p), , statistic] <- alpha_real
        
      }
    }
  }
}

for (statistic in statistic_list){
  # Save results into a table
  save.res <- mod.results[,,,,statistic]
  save(save.res, file = paste0("results/crit_val_MC_", M, "_",
                               statistic, "_MC.RData"))
}
