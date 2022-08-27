# Sunspots PAD uniformity test
#
# Sunspot births uniformity are analyzed with the PAD test over a 1-Carrington day
# rolling window of prefixed intervals (rotation.k, expressed in Carrington rotations)
# Solar cycles analyzed are specified in variable cycles.
# The test is tagged with the date in which the observation was made.

library(rotasym)

source("./src/unif_test.R")

dir.create(file.path("./sunspots/results/"), showWarnings = FALSE)

# Carrington rotation
carrington.rotation <- 27.2753 * 3600 * 24 # seconds

# Analysis variables
cycles <- c(21, 22, 23)
rotation.k <- c(10)
obs.window.list <- rotation.k * carrington.rotation

# Choose the uniformity tests to be performed
statistic.list <- c("PAD_2")
# Available statistics: ("PAD", "PCvM", "Bakshaev", "Watson", "Kuiper")

for (cycle in cycles) {
  print(paste0("<============== Cycle ", cycle, " ==============>"))

  for (hemisphere in c("N", "S", "all")) {
    print(paste0("<==== ", hemisphere, " hemisphere ====>"))

    data <- rotasym::sunspots_births
    cycle.filter <- data$cycle == cycle
    hemisphere.filter <- if (hemisphere == "N") {
      data$phi >= 0
    } else if (hemisphere == "S") {
      data$phi < 0
    } else {
      TRUE
    }
    data <- data[cycle.filter & hemisphere.filter, c("date", "theta")]

    first.date <- min(data$date)
    last.date <- max(data$date)

    for (obs.window in obs.window.list) {
      print(paste0(
        obs.window / carrington.rotation,
        " Carrington rotation periods"
      ))

      total.cycles <- floor(as.numeric(difftime(last.date,
        first.date + obs.window,
        units = "secs"
      )) / carrington.rotation)

      test.results <- array(
        dim = c(total.cycles, length(statistic.list) + 1),
        dimnames = list(
          1:total.cycles,
          c("date", statistic.list)
        )
      )

      execution.times <- array(
        dim = c(total.cycles, length(statistic.list)),
        dimnames = list(1:total.cycles, statistic.list)
      )

      for (statistic in statistic.list) {
        asymp.quantiles <- asymptotic_quantiles_table(statistic)

        ini <- first.date

        for (i in 1:total.cycles) {
          end <- ini + obs.window

          Theta <- as.matrix(data[
            (data$date >= ini) & (data$date < end),
            "theta"
          ])

          # Clock start
          start_time <- Sys.time()
          if (statistic %in% avail_cir_tests) {
            test <- unif_test_mod(Theta, statistic)
          } else if (statistic %in% avail_sph_tests) {
            Theta <- array(c(apply(Theta, 2, cos), apply(Theta, 2, sin)),
              dim = c(length(Theta), 2, 1)
            )
            test <- unif_test_mod(Theta, statistic)
          } else {
            stop(paste0(statistic, " statistic is not available."))
          }

          # Clock end
          end_time <- Sys.time()
          exec_time <- end_time - start_time
          execution.times[i, statistic] <- exec_time
          test.results[i, statistic] <- if (is.na(test$p_value)) {
            0.26
          } else {
            as.numeric(test$p_value)
          }

          ini <- ini + carrington.rotation
          test.results[i, "date"] <- as.numeric(
            as.POSIXct(strptime(end, "%Y-%m-%d %H:%M:%S"))
          )
        }
      }

      # Save results
      test.results <- data.frame(test.results)
      test.results$date <- as.POSIXct(test.results$date, origin = "1970-01-01")
      dir.create(file.path(paste0("./sunspots/results/", cycle)),
        showWarnings = FALSE
      )

      save(test.results, file = paste0(
        "./sunspots/results/",
        cycle, "/sunspots_cir_test_",
        hemisphere, "_",
        obs.window / carrington.rotation,
        "_results.RData"
      ))

      save(execution.times, file = paste0(
        "./sunspots/results/",
        cycle, "/sunspots_cir_test_",
        hemisphere, "_",
        obs.window / carrington.rotation,
        "_exec_times.RData"
      ))
    }
  }
}
