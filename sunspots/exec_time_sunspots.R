# Comparison of Monte Carlo and (n, alpha)-stabilization execution times

cycles <- c(21, 22, 23)
tot_ours <- 0
tot_MC <- 0

for (cycle in cycles) {
  load(paste0(
    "./sunspots/results/",
    cycle, "/sunspots_cir_test_all_10_exec_times.RData"
  ))
  ours <- execution.times

  load(paste0(
    "./sunspots/results/",
    cycle, "/sunspots_cir_test_all_10_MC_exec_times.RData"
  ))
  MC <- execution.times

  tot_MC <- tot_MC + sum(MC)
  tot_ours <- tot_ours + sum(ours)
}

tot_MC / length(cycles)
tot_ours / length(cycles)
