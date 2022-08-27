exact_quantile_from_asymp <- function(T_inf, statistic, n, alpha, Stephens = FALSE) {
  if (statistic == "KS") {
    exact.stat <- T_inf / Kuiper_mod(T_inf, n, alpha,
      KS = TRUE, Stephens = Stephens
    )
  }
  if (statistic == "Kuiper") {
    exact.stat <- T_inf / Kuiper_mod(T_inf, n, alpha,
      KS = FALSE, Stephens = Stephens
    )
  }
  if (statistic == "CvM") {
    exact.stat <- T_inf / Watson_mod(T_inf, n, alpha,
      CvM = TRUE, Stephens = Stephens
    )
  }
  if (statistic == "Watson") {
    exact.stat <- T_inf / Watson_mod(T_inf, n, alpha,
      CvM = FALSE, Stephens = Stephens
    )
  }
  if (statistic == "AD") {
    exact.stat <- T_inf / AndersonDarling_mod(T_inf, n, alpha,
      Stephens = Stephens
    )
  }
  return(exact.stat)
}

Kuiper_mod <- function(statistic, n, alpha, KS = FALSE, Stephens = FALSE) {
  if (Stephens) {
    if (KS) {
      f <- (1 + 0.12 * n^(-1 / 2) + 0.11 * n^(-1))
    } else {
      f <- (1 + 0.155 * n^(-1 / 2) + 0.24 * n^(-1))
    }
  } else {
    if (KS) {
      f <- (1 + 0.1570 * n^(-1 / 2) + 0.0195 * (n)^(-1) * (alpha)^(-1 / 2) - 0.0052 * (alpha * n)^(-1 / 2))
    } else {
      f <- 1 + 0.2342 * n^(-1 / 2) + 0.0268 * (alpha)^(-1 / 2) * (n)^(-1) - 0.0066 * (alpha * n)^(-1 / 2)
    }
  }

  return(f)
}

Watson_mod <- function(statistic, n, alpha, CvM = FALSE, Stephens = FALSE) {
  if (Stephens) {
    if (CvM) {
      statistic <- statistic - 0.4 * n^(-1) + 0.6 * n^(-2)
      f <- (1 + 1 * n^(-1))
    } else {
      statistic <- statistic - 0.1 * n^(-1) + 0.1 * n^(-2)
      f <- (1 + 0.8 * n^(-1))
    }
  } else {
    if (CvM) {
      f <- (1 - 0.1654 * n^(-1) + 0.0749 * n^(-1) * (alpha)^(-1 / 2) - 0.0014 * n^(-1) * alpha^(-1))
    } else {
      f <- (1 - 0.1517 * n^(-1) + 0.0917 * n^(-1) * alpha^(-1 / 2) - 0.0018 * (n * alpha)^-1)
    }
  }

  return(f)
}

AndersonDarling_mod <- function(statistic, n, alpha, Stephens = FALSE) {
  if (Stephens) {
    f <- 1
  } else {
    f <- (1 + 0.0356 * n^(-1) - 0.0234 * n^(-1) * alpha^(-1 / 2) + 0.0006 * (n * alpha)^-1)
  }

  return(f)
}
