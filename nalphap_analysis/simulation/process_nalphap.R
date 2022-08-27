source("./src/approxstats_utils.R")

M <- 1e5

for (statistic in c("PCvM", "Bakshaev", "PAD")) {
  print(paste0("<========= ", statistic, " =========>"))

  for (range_p_n in c("gt_p21_n50", "leq_p21_n300")) {

    # Load Stephens MC results
    load(file = paste0("./results/crit_val_samples_MC_", M, "_", statistic, "_", range_p_n, ".RData"))

    results <- process_MC_results(mod.results)
    mod.results <- results$rejection.results
    ci.lower <- results$ci.lower
    ci.upper <- results$ci.upper

    # Save summary results
    save(mod.results, ci.lower, ci.upper,
      file = paste0("results/crit_val_MC_", M, "_", statistic, "_", range_p_n, ".RData")
    )
  }
}
