# Kde for circular-linear data
kde_cir_lin <- function(th, z, data_cir, data_lin, h, g) {
  n <- length(data_cir)
  stopifnot(n == length(data_lin))
  log_c <- -log((2 * pi) * besselI(x = 1 / h^2, nu = 0, expon.scaled = TRUE))
  kde_cir <- exp(log_c - (1 - cos(outer(th, data_cir, "-"))) / h^2)
  kde_lin <- dnorm(outer(z, data_lin, "-"), sd = g)
  return(rowSums(kde_cir * kde_lin) / n)
}

# Estimation of level sets
get_c_alpha <- function(data_cir, data_lin, h, g, alpha) {
  kde <- kde_cir_lin(
    th = data_cir, z = data_lin,
    data_cir = data_cir, data_lin = data_lin,
    h = h, g = g
  )
  return(unname(quantile(kde, probs = 1 - alpha)))
}
