library(goftest)
library(sphunif)

M <- 1e5
Stephens <- FALSE

alpha.list <- seq(0.01, 0.25, 0.01)
colors <- c("#5BBEDA", "#3C92BF", "#ABABAB", "#C0D52B", "#D5AC2B", "#D5802B", "#D55C2B")
n.list <- c(5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 100, 200, 300)

for (statistic in c("KS", "CvM", "AD")){
  
  print(paste0("<=========== ", statistic, " ===========>"))
  
  mod.results <- array(dim = c(M, length(n.list), length(alpha.list)),
                       dimnames = list(1:M, n.list, alpha.list))
  
  for (alpha in alpha.list){
    
    print(paste0("alpha = ", alpha))
    
    for (n in n.list){
      
      print(paste0("n = ", n))
      
      # Monte Carlo statistic simulation
      cir <- unif_stat_MC(n = n, M = M, type = statistic, p = 2, alpha = alpha, 
                            chunks = 100, seeds = 1:100)$stats_MC
      
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
      
      # Store results in matrix
      mod.results[, as.character(n), as.character(alpha)] <- alpha.mod
      
    }
  }
  
  # Save results into a table
  save(mod.results, file = paste0("results/crit_val_samples_MC_", M,"_",
                                  statistic, "_particular_method.RData"))
  
}
