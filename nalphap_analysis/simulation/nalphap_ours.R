source("src/sph_stat_mod.R")

M <- 1e5

alpha.list <- seq(0.01, 0.25, 0.01)
colors <- c("#5BBEDA", "#3C92BF", "#ABABAB", "#C0D52B", "#D5AC2B", "#D5802B", "#D55C2B")

for (statistic in c("PCvM", "Bakshaev", "PAD")) {
  print(paste0("<========== ", statistic, " ==========>"))

  for (range_p_n in c("gt_p21_n50", "leq_p21_n300")) {
    if (range_p_n == "gt_p21_n50") {
      n.list <- c(5, 6, 7, 8, 9, 10, 25, 50)
      p.list <- c(51, 101, 151, 201, 301)
    }
    if (range_p_n == "leq_p21_n300") {
      n.list <- c(5, 6, 7, 8, 9, 10, 25, 50, 100, 200, 300)
      p.list <- c(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 21)
    }

    # Asymptotic critical values
    load(paste0("distributions/n500_asymp_qua_", statistic, ".RData"))
    asymp.quantiles <- quantiles

    mod.results <- array(
      dim = c(M, length(n.list), length(alpha.list), length(p.list)),
      dimnames = list(1:M, n.list, alpha.list, p.list)
    )
    for (p in p.list) {
      print(paste0("p = ", p))

      for (alpha in alpha.list) {
        print(paste0("alpha = ", alpha))

        # Asymptotic statistic critical value
        asymp.crit.val <- asymp.quantiles[as.character(alpha), as.character(p)] # n = 500

        for (n in n.list) {
          print(paste0("n = ", n))

          # Monte Carlo statistic simulation
          sph <- unif_stat_MC(
            n = n, M = M, type = statistic, p = p, alpha = alpha,
            chunks = 100, seeds = 1:100
          )

          # Modification of the statistic values obtained
          sph.mod <- apply(
            sph$stats_MC,
            c(1, 2),
            function(x) {
              return(sph_stat_mod(type = statistic, x, n, alpha, p))
            }
          )

          # Actual proportion of rejected statistics (alpha.mod)
          alpha.mod <- (sph.mod > asymp.crit.val)

          # Store results in matrix
          mod.results[, as.character(n), as.character(alpha), as.character(p)] <- alpha.mod
        }
      }
    }

    # Save results into a table
    save(mod.results, file = paste0(
      "results/crit_val_samples_MC_", M, "_",
      statistic, "_", range_p_n, ".RData"
    ))
  }
}
