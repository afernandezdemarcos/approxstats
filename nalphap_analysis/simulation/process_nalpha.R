source("./src/approxstats_utils.R")

M <- 1e5

for (statistic in c("KS", "Kuiper", "CvM", "Watson", "AD")) {
  print(paste0("<========= ", statistic, " =========>"))

  for (method in c("stephens", "ours")) {
    print(paste0("<===== ", method, " =====>"))

    # Load Stephens MC results
    load(file = paste0("./results/crit_val_samples_MC_", M, "_", statistic, "_", method, ".RData"))

    results <- process_MC_results(mod.results)
    mod.results <- results$rejection.results
    ci.lower <- results$ci.lower
    ci.upper <- results$ci.upper

    # Save summary results
    save(mod.results, ci.lower, ci.upper,
      file = paste0("results/crit_val_MC_", M, "_", statistic, "_", method, ".RData")
    )
  }
}
