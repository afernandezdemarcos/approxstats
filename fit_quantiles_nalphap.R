# Load src functions and packages

source("src/get_quantile_ratios.R")# Quantile ratios
source("src/fit_T_Tn_stable_ratios.R")# Fit model functions
source("src/weights.R")# Weight functions
source("src/drop_terms.R")# Further simplification from BIC-optimal model
library(tidyselect)
library(dplyr)
library(car)

# Parameters
# Dimension of ambient space R^p
p.list <- c(seq(2, 11, 1), seq(21, 101, 10), seq(151, 301, 50))
# Sample sizes to load
N.start <- 1
N.end <- 500
# Quantiles alpha values to load
q.start <- 0.001
q.end <- 1

# Choose a statistic in avail_sph_tests
statistic <- "Bakshaev"

# Prepare data
data <- data.frame(n=integer(),
                   alpha=double(),
                   T_Tn=double(),
                   inv_n_5_2=double(),
                   inv_n_4_2=double(),
                   inv_n_3_2=double(),
                   inv_n_2_2=double(),
                   inv_n_1_2=double(),
                   alpha_1_2=double(),
                   inv_alpha_2_2=double(),
                   inv_alpha_1_2=double(),
                   p = integer(),
                   inv_p_1_2 = double(), 
                   inv_p_2_2 = double(),
                   inv_p_3_2 = double(),
                   inv_p_4_2 = double(),
                   stringsAsFactors=FALSE)

for (p in p.list){
  
  # Load asymptotic distributions previously stored
  load(paste0("distributions/n500_asymp_qua_", statistic, ".RData"))
  asymp.quantiles <- quantiles[(as.numeric(rownames(quantiles) >= q.start)) 
                               & (as.numeric(rownames(quantiles) <= q.end)), 
                               as.character(p)]
  
  # Load Monte Carlo null quantiles for up to N sample sizes and quantiles
  load(paste0("distributions/qua_p", p, "_M1e7.RData"))
  
  alpha_mask <- (crit_val_n["alpha"] >= q.start) & (crit_val_n["alpha"] <= q.end)
  n_mask <- (crit_val_n["n"] >= N.start) & (crit_val_n["n"] <= N.end)
  quantiles_long <- crit_val_n[alpha_mask & n_mask, 
                               c("n", "alpha", all_of(statistic))]
  
  quantiles_wide <- pivot_wider(quantiles_long, 
                                names_from = alpha, 
                                values_from = all_of(statistic))
  
  quantiles <- as.matrix(select(quantiles_wide, !c("n")))
  rownames(quantiles) <- pull(quantiles_wide, 'n')
  
  # Compute ratios and T_Tn
  T_and_Tn <- get_T_Tn(quantiles = quantiles,
                       statistic = statistic,
                       compute.asymp = asymp.quantiles,
                       p = p,
                       return.T.Tn = TRUE)
  T_Tn <- get_T_Tn(quantiles = quantiles,
                   statistic = statistic,
                   compute.asymp = asymp.quantiles,
                   p = p)
  
  T_inf <- T_and_Tn$T_infty
  T_n <- T_and_Tn$T_n
  pdf_T_n <- compute_pdf(T_n)
  
  # Fit for different p values
  data.p <- prepare_data(T_Tn = T_Tn, T_inf = T_inf, T_n = T_n, 
                         f_T_n = pdf_T_n,
                         N.start = 5,
                         N.end = N.end,
                         fit.alpha = TRUE,
                         lambda = 5, mu = 2)
  data.p$p <- p
  data.p$inv_p_1_2 <- p^(-1/2)
  data.p$inv_p_2_2 <- p^(-2/2)
  data.p$inv_p_3_2 <- p^(-3/2)
  data.p$inv_p_4_2 <- p^(-4/2)
  
  data <- rbind(data, data.p)
  
}

# Fit T_Tn
fit.p <- TRUE

T_Tn_manual <- fit_T_Tn(T_Tn = NULL, T_inf = NULL, T_n = NULL, f_T_n = NULL, 
                        fit.direction = "both",
                        N.start = 5,
                        fit.p = fit.p,
                        data = data,
                        k = "BIC",
                        model.reference = "T_Tn ~ 0 + (inv_p_2_2 + inv_p_1_2):(
                        inv_n_2_2 + inv_n_2_2:inv_alpha_2_2 + inv_n_2_2:inv_alpha_1_2)",
                        weights.fun = weight.fun.2)

summary(T_Tn_manual$both)
model <- T_Tn_manual$both

# Drop terms based on R^2
R2_original <- summary(T_Tn_manual$both)$r.squared
new.formula <- drop_term_r_squared(model, delta.max = 1e-3, 
                                   R2_0 = R2_original, scope = NULL)
i <- 0
while (!is.null(new.formula)) {

  T_Tn_particular <- fit_T_Tn(T_Tn = T_Tn, T_n = T_n,
                              T_inf = T_inf, f_T_n = pdf_T_n,
                              fit.direction = c("both"),
                              N.start = 5,
                              fit.p = fit.p,
                              data = data,
                              model.formula = new.formula,
                              weights.fun = weight.fun.2
  )

  new.formula <- drop_term_r_squared(T_Tn_particular$particular,
                                     delta.max = 1.5e-3, R2_0 = R2_original)

  i <- i + 1
  print(paste0("Nb. of variables dropped: ", i))

}

summary(T_Tn_particular$particular)
