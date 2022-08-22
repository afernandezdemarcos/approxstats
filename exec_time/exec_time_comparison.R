source("exec_time/exact_from_asymp_stats.R")
source("src/approxstats_utils.R")
library(microbenchmark)
library(Rcpp)
library(sphunif)

sourceCpp("src/cir_stat_mod_newton.cpp")

dir.create(file.path("./results/"), showWarnings = FALSE)

# Compute for pairs of (n, alpha) p-values, and save execution statistics.
n.list <- c(5, 6, 7, 8, 9, 10, 20, 30, 40, 50)
alpha.list <- c(0.01, 0.02, 0.05, 0.1, 0.15, 0.25)

# Methods to be compared for each statistic
D.exec_time <- array(dim = c(length(n.list), length(alpha.list), 2),
                   dimnames = list(n.list, alpha.list, NULL))
W2.exec_time <- array(dim = c(length(n.list), length(alpha.list), 3),
                   dimnames = list(n.list, alpha.list, NULL))
A2.exec_time <- array(dim = c(length(n.list), length(alpha.list), 2),
                   dimnames = list(n.list, alpha.list, NULL))

D_inf <- asymptotic_quantiles_table("KS")
W2_inf <- asymptotic_quantiles_table("CvM")
A2_inf <- asymptotic_quantiles_table("AD")

newton_raphson_quadratic_R <- function(T_n, n, 
                               beta_n, beta_n_alpha_1_2, beta_nalpha, 
                               T_inf2, T_inf1, alpha2, alpha1, thr = 1e-7){
  
  fobj <- function(alpha){
    m <- (T_inf2 - T_inf1)/(alpha2 - alpha1);
    h <- (T_n - T_inf2 + m*alpha2 + beta_n*T_n/n)*alpha + 
      (beta_n_alpha_1_2*T_n/n)*sqrt(alpha) + beta_nalpha*T_n/n - m*alpha*alpha;
    return(h)
  }
  derobj <- function(alpha){
    m <- (T_inf2 - T_inf1)/(alpha2 - alpha1);
    der <- (T_n - T_inf2 + m*alpha2 + beta_n*T_n/sqrt(n)) + 
      1/2*(beta_n_alpha_1_2*T_n/n)/sqrt(alpha) - 2*m*alpha;
    return(der)
  }
  
  a <- alpha1 + (alpha2-alpha1)/2
  z <- a
  e <- 0
  while((abs(fobj(z))>thr) && (e < 7)){
    e <- e + 1
    z <- a - (fobj(a)/derobj(a))
    a <- z
  }
  
  return(a)
  
}

p_cir_stat_CvM <- function(T_n, n){
  alpha_grid <- seq(0.001, 0.25, 0.001)
  m <- length(alpha_grid)
  T_inf <- c(1.167869, 1.03846, 0.9633351, 0.9103154, 0.8693864, 0.8360821, 
             0.8080174, 0.7837828, 0.7624706, 0.7434891, 0.7262938, 0.7106819, 
             0.6963499, 0.6831086, 0.6708072, 0.6593525, 0.6485588, 0.6384363, 
             0.6288909, 0.6197956, 0.6112114, 0.6030277, 0.5952203, 0.5877707, 
             0.5806214, 0.5737654, 0.5671797, 0.5608443, 0.5547412, 0.5488543, 
             0.5431421, 0.5376563, 0.5323439, 0.5271954, 0.5222018, 0.5173548, 
             0.5126467, 0.5080701, 0.5036184, 0.4992852, 0.4950577, 0.4909746, 
             0.4869476, 0.4830278, 0.4792077, 0.4754568, 0.471827, 0.4682607, 
             0.464771, 0.4613538, 0.4580049, 0.4547201, 0.4515563, 0.4484017, 
             0.445322, 0.4423021, 0.4393398, 0.4364336, 0.4335817, 0.4307826, 
             0.4280349, 0.4252765, 0.4226228, 0.4200009, 0.4174821, 0.4149427, 
             0.4124424, 0.4099876, 0.4075772, 0.4052044, 0.4028679, 0.4005669,
             0.3983004, 0.3960669, 0.3938653, 0.3916953, 0.3895561, 0.3874467, 
             0.3853665, 0.3833147, 0.3812906, 0.3792935, 0.3773227, 0.3753777, 
             0.3734577, 0.3715622, 0.3696907, 0.3678424, 0.3660171, 0.364214, 
             0.3624327, 0.3606728, 0.3589337, 0.357215, 0.3555163, 0.3538371, 
             0.352177, 0.3505357, 0.3489127, 0.3473077, 0.3457203, 0.3441502, 
             0.3425969, 0.3410603, 0.3395399, 0.3380354, 0.3365466, 0.3350731, 
             0.3336147, 0.3321711, 0.3307419, 0.329324, 0.3279232, 0.3265362, 
             0.3251626, 0.3238023, 0.3224549, 0.3211204, 0.3197984, 0.3184887, 
             0.3171912, 0.3159057, 0.3146319, 0.3133696, 0.3121186, 0.3108789, 
             0.3096501, 0.3084321, 0.3072248, 0.306028, 0.3048414, 0.303665, 
             0.3024986, 0.3013421, 0.3001952, 0.2990579, 0.2979299, 0.2968113, 
             0.2957017, 0.2945745, 0.2934823, 0.2923989, 0.2913241, 0.2902577, 
             0.2891997, 0.28815, 0.2871083, 0.2861357, 0.28511, 0.2840921, 
             0.2830818, 0.2820791, 0.2810839, 0.2800829, 0.2791023, 0.278129, 
             0.2771627, 0.2762035, 0.2752512, 0.2743058, 0.2733671, 0.272435, 
             0.2715095, 0.2705906, 0.269678, 0.2687718, 0.2678718, 0.266978, 
             0.2660903, 0.2652087, 0.2643329, 0.2634631, 0.2625991, 0.2617408,
             0.2608882, 0.2600413, 0.2591998, 0.2583639, 0.2575333, 0.2567081,
             0.2558882, 0.2550735, 0.254264, 0.2534596, 0.2526602, 0.2518659, 
             0.2510765, 0.250292, 0.2495123, 0.2487374, 0.2479672, 0.2472017, 
             0.2464408, 0.2456845, 0.2449327, 0.2441854, 0.2434426, 0.2427041,
             0.24197, 0.2412401, 0.2405145, 0.2397931, 0.2390759, 0.2383628,
             0.2376537, 0.2369487, 0.2362477, 0.2355507, 0.2348575, 0.2341683, 
             0.2334828, 0.2328012, 0.2321252, 0.2314533, 0.2307851, 0.2301204, 
             0.2294593, 0.2288017, 0.2281476, 0.227497, 0.2268497, 0.2262058,
             0.2255653, 0.2249281, 0.2242943, 0.2236636, 0.2230362, 0.2224121,
             0.22173, 0.2211122, 0.2204974, 0.2198858, 0.2192772, 0.2186717, 
             0.2180692, 0.2174697, 0.2168731, 0.2162795, 0.2156888, 0.2151009, 
             0.214516, 0.2139339, 0.2133546, 0.2127781, 0.2122043, 0.2116333, 
             0.2110651, 0.2104995, 0.2099366, 0.2093764)
  
  for(i in 1:m){
    alpha <- alpha_grid[i]
    T_mod <- cir_stat_Watson_mod(T_n, n, alpha, CvM = FALSE)
    
    if(T_mod > T_inf[i]){
      if(i > 1){
        alpha <- newton_raphson_quadratic_R(T_n, n, -0.1651, 0.0749, -0.0014, 
                                       T_inf[i-1], T_inf[i], alpha_grid[i-1], alpha)
      }
      return(alpha)
    }
  }
  return(NA);
}

for (n in n.list){
  
  print(paste0("n = ", n))
  
  for (alpha in alpha.list){
    
    D_n <- exact_quantile_from_asymp(T_inf = D_inf[[as.character(alpha)]], 
                                     statistic = "KS", n = n, alpha = alpha)
    W2_n <- exact_quantile_from_asymp(T_inf = W2_inf[[as.character(alpha)]], 
                                      statistic = "CvM", n = n, alpha = alpha)
    A2_n <- exact_quantile_from_asymp(T_inf = A2_inf[[as.character(alpha)]], 
                                      statistic = "AD", n = n, alpha = alpha)

    KS <- summary(
      microbenchmark(
        "Ours" = {
          p_cir_stat_Kuiper(T_n = D_n, n = n, KS = TRUE)
        },
        "Marsaglia et al. (2003)" = {
          # the statistic for this cdf is Dn/sqrt(n)
          1 - .Call(stats:::C_pKolmogorov2x, D_n/sqrt(n), n)
        },
        unit = "us",
        control = list(warmup = 10),
        times = 1e3
      )
    )
    
    CvM <- summary(
      microbenchmark(
        "Ours_R" = {
          p_cir_stat_CvM(T_n = W2_n, n = n)
        },
        "Ours_C++" = {
          p_cir_stat_Watson(T_n = W2_n, n = n, CvM = T)
        },
        "Csorgo & Faraway (1996)" = {
          goftest::pCvM(q = W2_n, n = n, lower.tail = F)
        },
        unit = "us",
        control = list(warmup = 10),
        times = 1e3
      )
    )
    
    AD <- summary(
      microbenchmark(
        "Ours" = {
          p_cir_stat_AD(T_n = A2_n, n = n)
        },
        "Marsaglia & Marsaglia (2004)" = {
          goftest::pAD(q = A2_n, n = n, lower.tail = F)
        },
        unit = "us",
        control = list(warmup = 10),
        times = 1e3
      )
    )
    
    D.exec_time[as.character(n), as.character(alpha), ] <- KS$median
    W2.exec_time[as.character(n), as.character(alpha), ] <- CvM$median
    A2.exec_time[as.character(n), as.character(alpha), ] <- AD$median

  }
}

dimnames(D.exec_time)[[3]] <- levels(KS$expr)
dimnames(W2.exec_time)[[3]] <- levels(CvM$expr)
dimnames(A2.exec_time)[[3]] <- levels(AD$expr)

save(D.exec_time, W2.exec_time, A2.exec_time, 
     file = "./results/execution_time_comparison.RData")

####################################
# MC = 1e4
alpha <- 0.05

D.exec_time.MC <- numeric(length(n.list))
W2.exec_time.MC <- numeric(length(n.list))
A2.exec_time.MC <- numeric(length(n.list))

names(D.exec_time.MC) <- n.list
names(W2.exec_time.MC) <- n.list
names(A2.exec_time.MC) <- n.list

for (n in n.list){
  
  D_n <- exact_quantile_from_asymp(T_inf = D_inf[[as.character(alpha)]], 
                                   statistic = "KS", n = n, alpha = alpha)
  W2_n <- exact_quantile_from_asymp(T_inf = W2_inf[[as.character(alpha)]], 
                                    statistic = "CvM", n = n, alpha = alpha)
  A2_n <- exact_quantile_from_asymp(T_inf = A2_inf[[as.character(alpha)]], 
                                    statistic = "AD", n = n, alpha = alpha)
  
  KS_MC <- summary(
    microbenchmark(
      "MC" = {
        MC_Tn <- unif_stat_MC(n = n, statistic = "KS", p = 2, M = 1e4,
                              cores = 1, return_stats = TRUE)$stats_MC
        mean(D_n > MC_Tn)
      },
      unit = "us",
      control = list(warmup = 10),
      times = 1e3
    )
  )
  
  CvM_MC <- summary(
    microbenchmark(
      "MC" = {
        MC_Tn <- unif_stat_MC(n = n, statistic = "CvM", p = 2, M = 1e4,
                              cores = 1, return_stats = TRUE)$stats_MC
        mean(W2_n > MC_Tn)
      },
      unit = "us",
      control = list(warmup = 10),
      times = 1e3
    )
  )
  
  AD_MC <- summary(
    microbenchmark(
      "MC" = {
        MC_Tn <- unif_stat_MC(n = n, statistic = "AD", p = 2, M = 1e4,
                              cores = 1, return_stats = TRUE)$stats_MC
        mean(A2_n > MC_Tn)
      },
      unit = "us",
      control = list(warmup = 10),
      times = 1e3
    )
  )
  
  D.exec_time.MC[as.character(n)] <- KS_MC$median
  W2.exec_time.MC[as.character(n)] <- CvM_MC$median
  A2.exec_time.MC[as.character(n)] <- AD_MC$median
  
}

save(D.exec_time.MC, W2.exec_time.MC, A2.exec_time.MC, 
     file = "./results/execution_time_MC_comparison.RData")

### Time Analysis
load("./results/execution_time_comparison.RData")
mean(D.exec_time[,,"Marsaglia et al. (2003)"]/D.exec_time[,,"Ours"])
mean(W2.exec_time[,,"Csorgo & Faraway (1996)"]/W2.exec_time[,,"Ours_R"])
mean(A2.exec_time[,,"Marsaglia & Marsaglia (2004)"]/A2.exec_time[,,"Ours"])

load("./results/execution_time_MC_comparison.RData")
mean(D.exec_time.MC/D.exec_time[,"0.05","Ours"])
mean(W2.exec_time.MC/W2.exec_time[,"0.05","Ours_C++"])
mean(A2.exec_time.MC/A2.exec_time[,"0.05","Ours"])
