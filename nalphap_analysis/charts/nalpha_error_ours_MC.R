# Error analysis comparison between (n, alpha)-stabilization and Monte Carlo

library(Hmisc)

M <- 1e7
M.name <- "1e7"

alpha.plot <- c("0.25", "0.2", "0.15", "0.1", "0.05", "0.02", "0.01")

# Overall errors of Stephens and ours stabilizations

for (statistic in c("KS", "Kuiper", "CvM", "Watson", "AD")){
  
  # Load Stephens MC results
  load(file = paste0("results/crit_val_MC_", M.name,"_", statistic, "_stephens.RData"))
  mod.results.steph <- mod.results[, alpha.plot]
  
  # Load Ours MC results
  load(file = paste0("results/crit_val_MC_", M.name,"_", statistic, "_ours.RData"))
  mod.results.ours <- mod.results[, alpha.plot]
  
  alpha <- as.numeric(colnames(mod.results.steph))
  alpha.matrix <- alpha[col(mod.results.ours)]
  
  # Absolute errors
  error.steph <- abs(mod.results.steph - alpha.matrix)
  error.ours <- abs(mod.results.ours - alpha.matrix)

  # Relative errors
  rel.error.steph <- error.steph / alpha.matrix * 100
  rel.error.ours <- error.ours / alpha.matrix * 100

  # Marginals
  margin.error.steph <- apply(error.steph, 2, mean)
  margin.error.ours <- apply(error.ours, 2, mean)

  alpha.margin.rel.error.steph <- matrix(apply(rel.error.steph, 2, mean),
                                   nrow = 1,
                                   dimnames = list(NULL, colnames(error.steph))
                                   )
  n.margin.rel.error.steph <- matrix(apply(rel.error.steph, 1, mean),
                                   ncol = 1,
                                   dimnames = list(rownames(error.steph), NULL)
                                   )
  alpha.margin.rel.error.ours <- matrix(apply(rel.error.ours, 2, mean),
                                  nrow = 1,
                                  dimnames = list(NULL, colnames(error.ours))
                                   )
  n.margin.rel.error.ours <- matrix(apply(rel.error.ours, 1, mean),
                                  ncol = 1,
                                  dimnames = list(rownames(error.ours), NULL)
                                   )
  print(statistic)
  print(round(mean(rel.error.steph), 2))
  print(round(mean(rel.error.ours), 2))

}

################################################
#####    Charts alpha Stephens vs. ours    #####
################################################

legend_x.list <- list("Kuiper" = 0.05, "Watson" = 0.05)

for (statistic in c("Kuiper", "Watson")){
  
  # Load Stephens MC results
  load(file = paste0("results/crit_val_MC_", M.name,"_", statistic, "_stephens.RData"))
  mod.results.steph <- mod.results
  
  # Load Ours MC results
  load(file = paste0("results/crit_val_MC_", M.name,"_", statistic, "_ours.RData"))
  mod.results.ours <- mod.results
  
  alpha <- as.numeric(colnames(mod.results.steph))
  alpha.matrix <- alpha[col(mod.results.ours)]
  
  # Load MC 1e4 MC results
  load(file = paste0("results/crit_val_MC_", M.name,"_", statistic, "_MC.RData"))
  mod.results.MC <- save.res
  error.MC.m <- array(0, dim = dim(mod.results.MC), dimnames = dimnames(mod.results.MC))
  for(m in 1:dim(mod.results.MC)[3]){
    error.MC.m[,,m] <- abs(mod.results.MC[,,m] - alpha.matrix)
  }
  
  # MC confidence intervals
  ci.alpha <- binconf(alpha*M, M, alpha = 0.05)
  ci.rel.low <- abs(ci.alpha[,2]-ci.alpha[,1])/ci.alpha[,1]*100
  ci.rel.upp <- abs(ci.alpha[,3]-ci.alpha[,1])/ci.alpha[,1]*100
  ci.rel <- rowMeans(cbind(ci.rel.low, ci.rel.upp), na.rm=TRUE)
  
  # Absolute errors
  error.steph <- abs(mod.results.steph - alpha.matrix)
  error.ours <- abs(mod.results.ours - alpha.matrix)
  error.MC <- apply(error.MC.m, c(1,2), mean)
  
  # Relative errors
  rel.error.steph <- error.steph / alpha.matrix * 100
  rel.error.ours <- error.ours / alpha.matrix * 100
  rel.error.MC <- error.MC / alpha.matrix * 100
  
  bins <- c("n: [5, 10)", "n: [10, 100)", "n: [100, 300]")
  bins.limits <- c(10, 100, 500)
  rel.error.steph.summ <- matrix(nrow = length(bins), ncol = length(alpha), 
                                dimnames = list(bins, alpha))
  rel.error.ours.summ <- matrix(nrow = length(bins), ncol = length(alpha), 
                                dimnames = list(bins, alpha))
  rel.error.MC.summ <- matrix(nrow = length(bins), ncol = length(alpha), 
                                dimnames = list(bins, alpha))
  
  for (i in 1:length(bins)){
    index <- as.numeric(rownames(rel.error.steph)) < bins.limits[i]
    rel.error.steph.summ[bins[i],] <- apply(rel.error.steph[index, ], 2, mean)
    rel.error.ours.summ[bins[i],] <- apply(rel.error.ours[index, ], 2, mean)
    rel.error.MC.summ[bins[i],] <- apply(rel.error.MC[index, ], 2, mean)
  }
  
  colors <- c("#f0a500", "#8ac4d0", "#28527a")
  
  # Plot setup
  pdf(file = paste0("results/rel_error_alpha_", M.name, "_", statistic, ".pdf"),
      width = 6,
      height = 6)
  
  par(mfrow = c(1,1),
      xpd = FALSE)
  
  plot(range(as.numeric(colnames(rel.error.steph.summ))), 
       c(0, 8), type="n",
       # main = statistic
       ylab = "Relative Error (%)",
       xlab = expression(alpha),
       yaxt = "n"
  )
  ytick <- 0:8
  axis(side = 2, at = ytick, labels = FALSE)
  text(par("usr")[1], ytick+0.1, labels = ytick, adj = 0.5, pos = 2, offset = 1, srt = 90, xpd = TRUE)
  
  x.poly <- c(alpha[1], alpha, alpha[length(alpha)])
  y.poly <- c(range(rel.error.ours.summ, rel.error.steph.summ, rel.error.MC.summ)[1], ci.rel, range(rel.error.ours.summ, rel.error.steph.summ, rel.error.MC.summ)[1])
  polygon(x.poly, y.poly, col=gray(0.9), border=NA)
  
  i <- 0
  for (bin in bins){
    i <- i + 1
    
    errors.steph <- rel.error.steph.summ[bin, ]
    errors.ours <- rel.error.ours.summ[bin, ]
    errors.MC <- rel.error.MC.summ[bin, ]
    lines(as.numeric(names(errors.ours)), errors.ours, pch = 16, col = colors[i], cex = 1, lty = 1, type = "b")
    lines(as.numeric(names(errors.steph)), errors.steph, pch = 2, col = colors[i], cex = 1, lty = 2, type = "b")
    lines(as.numeric(names(errors.MC)), errors.MC, pch = 0, col = colors[i], cex = 1, lty = 2, type = "b")
  }
    
  dev.off()

}
