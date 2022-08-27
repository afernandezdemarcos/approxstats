# Asymptotic null cdf and quantile distributions

library(sphunif)
library(goftest)

avail_statistics <- c(
  "KS", "Kuiper", "CvM", # Classical
  "Watson", "AD", # Circular
  "PAD", "PCvM", "Bakshaev"
) # Hyperspherical

# Returns the asymptotic null distribution of statistic for T_infty values of
# the statistic. p indicates the dimension of the ambient space R^p
asymptotic_distribution <- function(statistic, p, T_infty = seq(0, 9, 0.0001),
                                    tol = 1e-6) {
  if (!(statistic %in% avail_statistics)) {
    stop(paste0(
      statistic,
      " statistic is not available among avail_statistics."
    ))
  }

  # Asymptotic null distribution
  if (statistic == "KS") {
    P_T_infty <- 1 - .Call(stats:::C_pKS2, T_infty, tol)
  }

  if (statistic == "Kuiper") {
    P_T_infty <- as.numeric(1 - sphunif::p_cir_stat_Kuiper(T_infty, n = 1e10))
  }

  if (statistic == "CvM") {
    P_T_infty <- 1 - goftest::pCvM(T_infty)
  }

  if (statistic == "Watson") {
    P_T_infty <- as.numeric(1 - sphunif::p_cir_stat_Watson(T_infty))
  }

  if (statistic == "AD") {
    P_T_infty <- 1 - goftest::pAD(T_infty)
  }

  if (statistic == "PAD") {
    P_T_infty <- 1 - p_sph_stat_PAD(T_infty, p = p)
  }

  if (statistic == "PCvM") {
    P_T_infty <- 1 - p_sph_stat_PCvM(T_infty, p = p)
  }

  return(P_T_infty)
}

# Returns the asymptotic null quantiles of statistic for alpha values.
# p indicates the dimension of the ambient space R^p
asymptotic_quantiles <- function(statistic, alpha, p,
                                 T_infty = seq(0, 9, 0.0001),
                                 tol = 1e-6) {
  if (!(statistic %in% avail_statistics)) {
    stop(paste0(
      statistic,
      " statistic is not available among avail_statistics."
    ))
  }

  if (statistic %in% c("KS", "Kuiper", "Watson")) {

    # Asymptotic probability distribution
    P_T_infty <- asymptotic_distribution(statistic, alpha, T_infty, tol)

    asymp.quantile <- numeric(length(alpha))
    names(asymp.quantile) <- alpha

    for (q in alpha) {
      asymp.quantile[as.character(q)] <- T_infty[which.max(P_T_infty <= q)]
    }
  }

  if (statistic == "CvM") {
    asymp.quantile <- goftest::qCvM(p = 1 - alpha)
  }

  if (statistic == "AD") {
    asymp.quantile <- goftest::qAD(p = 1 - alpha)
  }

  if (statistic == "PAD") {
    # asymp.quantile <- 1 - p_sph_stat_PAD(T_infty, p = p)
  }

  if (statistic == "PCvM") {
    asymp.quantile <- q_Sobolev(
      u = 1 - alpha, p = p, type = "PCvM",
      K_max = 1e5, thre = 0
    )
  }

  return(asymp.quantile)
}
