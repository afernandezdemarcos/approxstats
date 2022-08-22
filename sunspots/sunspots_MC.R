# Sunspots PAD uniformity test with MC (1e4)
#       
# Only execution times are saved to be compared with Algorithm 1 approximation

library(rotasym)
library(sphunif)

dir.create(file.path("./sunspots/results/"), showWarnings = FALSE)

# Carrington rotation
carrington.rotation <- 27.2753 * 3600 * 24# seconds

# Analysis variables
cycles <- c(21, 22, 23)
rotation.k <- c(10)
obs.window.list <-rotation.k * carrington.rotation

statistic.list <- c("PAD")
# Available statistics: ("PCvM", "Bakshaev", "Watson", "Kuiper")

for (cycle in cycles){
  
  print(paste0("<============== Cycle ", cycle, " ==============>"))
  
  for(hemisphere in c("N", "S", "all")){
    
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
    
    for (obs.window in obs.window.list){
      
      print(paste0(obs.window/carrington.rotation, 
                   " Carrington rotation periods"))
      
      total.cycles <- floor(as.numeric(difftime(last.date, 
                                                first.date + obs.window, 
                                                units = "secs")
                                       )/carrington.rotation)
      
      test.results <- array(dim = c(total.cycles, length(statistic.list) + 1),
                            dimnames = list(1:total.cycles, 
                                            c("date", statistic.list)))
      
      execution.times <- array(dim = c(total.cycles, length(statistic.list)),
                               dimnames = list(1:total.cycles, statistic.list))
      
      for (statistic in statistic.list){
        
        ini <- first.date
        
        for (i in 1:total.cycles){
          
          end <- ini + obs.window
          
          Theta <- as.matrix(data[(data$date >= ini) & (data$date < end), 
                                  "theta"])
          
          Theta <- array(c(apply(Theta, 2, cos), apply(Theta, 2, sin)), 
                         dim = c(length(Theta), 2, 1))
          # Clock start
          start_time <- Sys.time()
          test <- unif_test(Theta, type = statistic, p_value = "MC", M = 5e3,
                            chunks = 50)
          # Clock end
          end_time <- Sys.time()
          exec_time <- end_time - start_time
          execution.times[i, statistic] <- exec_time
          test.results[i, statistic] <- as.numeric(test$p.value)
          
          ini <- ini + carrington.rotation
          test.results[i, "date"] <- as.numeric(
            as.POSIXct(strptime(end, "%Y-%m-%d %H:%M:%S"))
          )
          
        }
      }
      
      save(execution.times, file = paste0("./sunspots/results/", 
                                          cycle,"/sunspots_cir_test_", 
                                          hemisphere, "_", 
                                          obs.window/carrington.rotation, 
                                          "_MC_exec_times.RData"))
      
    }
    
  }
  
}
