# Load src functions and packages

source("src/get_quantile_ratios.R")# Quantile ratios
source("src/fit_T_Tn_stable_ratios.R")# Fit model functions ########################## TODO: Ajustar automáticamente potencias de n, alpha
source("src/weights.R")# Weight functions
source("src/drop_terms.R")# Further simplification from BIC-optimal model
library(tidyselect)
library(dplyr)
library(car)

# Parameters
# Dimension of ambient space R^p (In this case, circular p = 2)
p <- 2
# Sample sizes to load
N.start <- 1
N.end <- 500
# Quantiles alpha values to load
q.start <- 0.001
q.end <- 1

# Choose a statistic in avail_cir_tests
statistic <- "KS"

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
                     p = p,
                     return.T.Tn = TRUE
                     )
T_Tn <- get_T_Tn(quantiles = quantiles,
                 statistic = statistic,
                 p = p)

T_inf <- T_and_Tn$T_infty
T_n <- T_and_Tn$T_n
pdf_T_n <- compute_pdf(T_n)

rm(quantiles_long, quantiles_wide, n_mask, alpha_mask, T_and_Tn, crit_val_n)

# Fit T_Tn
T_Tn_fitted <- fit_T_Tn(T_Tn = T_Tn, T_n = T_n,
                        T_inf = T_inf, f_T_n = pdf_T_n,
                        fit.direction = "both",
                        N.start = 5,
                        k = "BIC",
                        weights.fun = weight.fun.2,
                        lambda = 2, mu = 2)

summary(T_Tn_fitted$both)
BIC(T_Tn_fitted$both)
vif(T_Tn_fitted$both)

# Drop terms one-by-one based on R^2 decrease
R2_original <- summary(T_Tn_fitted$both)$r.squared
new.formula <- drop_term_r_squared(T_Tn_fitted$both, 
                                   delta.max = 5e-3, R2_0 = R2_original)
i <- 0
while (!is.null(new.formula)) {
  
  T_Tn_particular <- fit_T_Tn(T_Tn = T_Tn, T_n = T_n,
                              T_inf = T_inf, f_T_n = pdf_T_n,
                              fit.direction = c("both"),
                              N.start = 5,
                              model.formula = new.formula,
                              weights.fun = weight.fun.1
  )
  
  new.formula <- drop_term_r_squared(T_Tn_particular$particular, 
                                     delta.max = 1.5e-3, R2_0 = R2_original)
  
  i <- i + 1
  print(paste0("Nb. of variables dropped: ", i))
  
}

summary(T_Tn_particular$particular)
