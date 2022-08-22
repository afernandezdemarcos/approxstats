# Error analysis comparison between (n, alpha)-stabilization, Monte Carlo, 
# and particular methods

library(Hmisc)

M <- 1e7
M.name <- "1e7"

statistic.list <- c("CvM", "KS", "AD")
method_filename.list <- list("CvM"= "Csorgo_and_Faraway_1996", 
                        "KS" = "Marsaglia_et_al_2003", 
                        "AD" = "Marsaglia_and_Marsaglia_2004")
method_name.list <- list("CvM" = "Csörgo & Faraway (1996)",
                    "KS" = "Marsaglia et al. (2003)",
                    "AD" = "Marsaglia & Marsaglia (2004)")
legend_x.list <- list("CvM" = 0.025, "KS" = 0.05, "AD" = 0.1)

for (statistic in statistic.list){
  
  method_filename <- method_filename.list[[statistic]]
  method_name <- method_name.list[[statistic]]
  
  # Load Stephens MC results
  load(file = paste0("results/crit_val_MC_", M.name,"_", statistic, "_stephens.RData"))
  mod.results.steph <- mod.results
  
  # Load particular MC results
  load(file = paste0("results/crit_val_MC_", M.name,"_", statistic, "_particular_method.RData"))
  mod.results.part <- mod.results
  
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
  error.part <- abs(mod.results.part - alpha.matrix)
  error.ours <- abs(mod.results.ours - alpha.matrix)
  error.MC <- apply(error.MC.m, c(1,2), mean)
  
  # Relative errors
  rel.error.steph <- error.steph / alpha.matrix * 100
  rel.error.part <- error.part / alpha.matrix * 100
  rel.error.ours <- error.ours / alpha.matrix * 100
  rel.error.MC <- error.MC / alpha.matrix * 100
  
  bins <- c("n: [5, 10)", "n: [10, 100)", "n: [100, 300]")
  bins.limits <- c(10, 100, 500)
  rel.error.steph.summ <- matrix(nrow = length(bins), ncol = length(alpha), 
                                 dimnames = list(bins, alpha))
  rel.error.part.summ <- matrix(nrow = length(bins), ncol = length(alpha), 
                                dimnames = list(bins, alpha))
  rel.error.ours.summ <- matrix(nrow = length(bins), ncol = length(alpha), 
                                dimnames = list(bins, alpha))
  rel.error.MC.summ <- matrix(nrow = length(bins), ncol = length(alpha), 
                                dimnames = list(bins, alpha))
  
  for (i in 1:length(bins)){
    index <- as.numeric(rownames(rel.error.steph)) < bins.limits[i]
    rel.error.steph.summ[bins[i],] <- apply(rel.error.steph[index, ], 2, mean)
    rel.error.part.summ[bins[i],] <- apply(rel.error.part[index, ], 2, mean)
    rel.error.ours.summ[bins[i],] <- apply(rel.error.ours[index, ], 2, mean)
    rel.error.MC.summ[bins[i],] <- apply(rel.error.MC[index, ], 2, mean)
  }
  
  colors <- c("#f0a500", "#8ac4d0", "#28527a")
  
  # Plot setup
  pdf(file = paste0("results/rel_error_alpha_", M.name, "_", statistic, ".pdf"),
      width = 6,#7
      height = 6)#5
  
  par(mfrow = c(1,1),
      xpd = FALSE)
  
  plot(range(as.numeric(colnames(rel.error.steph.summ))), 
       c(0, 8), type="n",
       ylab = "Relative Error (%)",
       xlab = expression(alpha),
       yaxt = "n"
  )
  ytick <- 0:8
  axis(side = 2, at = ytick, labels = FALSE)
  text(par("usr")[1], ytick+0.1, labels = ytick, adj = 0.5, pos = 2, offset = 1, srt = 90, xpd = TRUE)
  
  x.poly <- c(alpha[1], alpha, alpha[length(alpha)])
  y.poly <- c(0, ci.rel, 0)
  polygon(x.poly, y.poly, col=gray(0.9), border=NA)
  i <- 0
  for (bin in bins){
    i <- i + 1
    
    errors.steph <- rel.error.steph.summ[bin, ]
    errors.part <- rel.error.part.summ[bin, ]
    errors.ours <- rel.error.ours.summ[bin, ]
    errors.MC <- rel.error.MC.summ[bin, ]
    lines(as.numeric(names(errors.ours)), errors.ours, pch = 16, col = colors[i], cex = 1, lty = 1, type = "b")
    lines(as.numeric(names(errors.steph)), errors.steph, pch = 2, col = colors[i], cex = 1, lty = 2, type = "b")
    lines(as.numeric(names(errors.part)), errors.part, pch = 4, col = colors[i], cex = 1, lty = 3, type = "b")
    lines(as.numeric(names(errors.MC)), errors.MC, pch = 0, col = colors[i], cex = 1, lty = 3, type = "b")
  }
  
  if(statistic == "KS"){
    legend("topright",
           legend = c("Stephens", method_name, "Ours", expression(paste("Monte Carlo")), bins), 
           col = c(rep("black", 4), colors),
           pch = c(2, 4, 16, 0, rep(16, length(bins))),
           lty = c(2, 3, 1, 3, rep(NA, length(bins))))
  } else {
    legend("topright", 
           legend = c(method_name), 
           bg = "white",
           col = c("black"),
           pch = c(4),
           lty = c(3))
  }
  
  dev.off()
  
}
