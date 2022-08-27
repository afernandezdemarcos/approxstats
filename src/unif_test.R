# Uniformity test with p-value approximation computed by Algorithm 1

library(sphunif)
library(Rcpp)
source("src/unif_stat_distr.R")

unif_test_mod <- function(Theta, statistic, sorted = FALSE, method = "grid") {

  # Set arguments for distribution and arguments name

  KS <- FALSE
  CvM <- FALSE
  AD <- FALSE
  name_statistic <- statistic

  if (statistic == "KS") {
    name_statistic <- "Kuiper"
    KS <- TRUE
  }

  if (statistic == "CvM") {
    name_statistic <- "Watson"
    CvM <- TRUE
  }

  if (statistic == "AD") {
    name_statistic <- "PAD"
    AD <- TRUE
  }

  if (statistic == "PAD_2") {
    name_statistic <- "PAD"
  }

  if (statistic %in% avail_cir_tests) {
    prefix <- "cir_stat_"
  } else {
    prefix <- "sph_stat_"
  }

  name_distr <- paste0(prefix, name_statistic)

  args <- list(
    "X" = Theta, "Theta" = Theta, "sorted" = sorted,
    "KS" = KS, "CvM" = CvM, "AD" = AD
  )
  names_args <- names(args)

  distr_args <- args[names_args %in% names(formals(name_distr))]

  n <- nrow(Theta) # sample size

  if (statistic %in% avail_sph_tests) {
    p <- ncol(Theta) # Hypersphere dimension
  } else if (statistic %in% avail_cir_tests) {
    p <- 2
  } else {
    stop(paste0(
      statistic,
      " statistic is not available among avail_statistics."
    ))
  }

  # Compute (n,alpha)-stabilized statistic
  T_n <- do.call(what = name_distr, args = c(distr_args))

  # (Algorithm 1) p-value approximation using (n, alpha)-stabilization
  if (method == "grid") {
    p.val <- p_value_grid(T_n = T_n, statistic = statistic, n = n, p = p)
    gt0.25 <- if (is.na(p.val)) TRUE else FALSE
  }

  return(list("statistic" = T_n, "p_value" = p.val, "greater_0.25" = gt0.25))
}
