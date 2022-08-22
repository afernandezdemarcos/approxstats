source("src/cir_stat_mod.R")
source("src/asymptotic_distribution.R")

library(Hmisc)
library(goftest)

M <- 1e6

statistic.list <- c("KS", "Kuiper", "CvM", "Watson", "AD")
alpha.list <- seq(0.01, 0.25, 0.01)
colors <- c("#5BBEDA", "#3C92BF", "#ABABAB", "#C0D52B", "#D5AC2B", "#D5802B", "#D55C2B")
n.list <- c(5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 100, 200, 300)

samples <- array(dim = c(M, length(statistic.list), length(n.list)),
                 dimnames = list(1:M, statistic.list, n.list))

for (n in n.list){
  print(paste0("n = ", n))
  
  cir <- unif_stat_MC(n = n, M = M, type = statistic.list, p = 2, 
                      chunks = 100, seeds = 1:100, cores = 3)$stats_MC
  
  for (statistic in statistic.list){
    samples[, statistic, as.character(n)] <- cir[, statistic]
  }
}

for (statistic in statistic.list){
  
  print(paste0("<=========== ", statistic, " ===========>"))
  
  if (statistic %in% c("KS", "CvM", "AD")){
    Stephens.list <- c("TRUE", "FALSE", "particular_method")
  }else{
    Stephens.list <- c("TRUE", "FALSE")
  }
  
  for (Stephens in Stephens.list){
    
    print(paste0("<=========== ", if (Stephens == "particular_method") "Particular" else if (Stephens == "TRUE") "Stephens" else "Ours", " ===========>"))
    
    mod.results <- array(dim = c(M, length(n.list), length(alpha.list)),
                         dimnames = list(1:M, n.list, alpha.list))
    
    for (alpha in alpha.list){
    
      print(paste0("alpha = ", alpha))
      
      # Asymptotic statistic critical value
      asymp.crit.val <- asymptotic_quantiles(statistic = statistic, alpha = alpha)
      
      for (n in n.list){
        
        print(paste0("n = ", n))
        
        cir <- samples[, statistic, as.character(n)]
        
        if(Stephens == "particular_method"){
          
          # approximated alpha-quantile by the corresponding method
          p_val <- if (statistic == "KS") {
            function(x) 1 - .Call(stats:::C_pKolmogorov2x, x/sqrt(n), n) - alpha
          } else if (statistic == "CvM"){
            function(x) goftest::pCvM(q = x, n = n, lower.tail = F) - alpha
          } else if (statistic == "AD"){
            function(x) goftest::pAD(q = x, n = n, lower.tail = F) - alpha
          }
          
          q <- uniroot(p_val, c(0,10))$root
          
          # Actual proportion of rejected statistics (alpha.mod)
          alpha.mod <- cir > q
          
        }else{
          
          # Modification of the statistic values obtained
          cir.mod <- unlist(lapply(cir,
                           function(x) return(cir_stat_mod(type = statistic, x, n, alpha, Stephens = as.logical(Stephens)))
          ))
          
          # Actual proportion of rejected statistics (alpha.mod)
          alpha.mod <- cir.mod > asymp.crit.val
          
        }

        # Store results in matrix
        mod.results[, as.character(n), as.character(alpha)] <- alpha.mod
        
      }
    }
    
    # Save results into a table
    save(mod.results, file = paste0("results/crit_val_samples_MC_", M,"_",
                                    statistic, "_", 
                                    if (Stephens == "particular_method") Stephens else if (Stephens == "TRUE") "stephens" else "ours", 
                                    ".RData"))
    
  }
}
