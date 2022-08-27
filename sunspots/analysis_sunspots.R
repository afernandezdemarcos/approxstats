library(rotasym)
library(np)
library(NPCirc)
library(sdetorus)
library(DirStats)
library(colorspace)

source("sunspots/kde_level_set.R")

save.plots <- TRUE

carrington.rotation <- 27.2753 * 3600 * 24 # seconds
cycles <- c(21, 22, 23)

rotation.k <- c(10) # window in carrington rotations
step <- 1 # rolling window step in carrington rotations
obs.window.list <- rotation.k * carrington.rotation # window in seconds

# Bandwidth for N-W nonparametric regression
adjust.bw <- TRUE
# if TRUE then optimal bw is computed by CV, otherwise bw0 is set as bandwidth
# (following previously computed bandwidth:)
bw0 <- rotation.k
bw0.list <- list(
  "21" = list("N" = 10.26768, "S" = 10.82924, "all" = 5.483238),
  "22" = list("N" = 2.687911, "S" = 5.273556, "all" = 2.022129),
  "23" = list("N" = 2.854751, "S" = 6.777612, "all" = 8.820482)
)

scale.factor <- 1e7 # Scale factor used for density ridges
alpha_level_set <- 0.2 # alpha level to compute highest density level sets

statistic <- "PAD_2"

colors <- list(
  "N" = rgb(
    red = 18, green = 93, blue = 152,
    alpha = 80, maxColorValue = 255
  ),
  "S" = rgb(
    red = 177, green = 99, blue = 110,
    alpha = 80, maxColorValue = 255
  ),
  "all" = rgb(
    red = 0, green = 0, blue = 0,
    alpha = 80, maxColorValue = 255
  )
)

line.colors <- list(
  "N" = rgb(
    red = 18, green = 93, blue = 152,
    alpha = 255, maxColorValue = 255
  ),
  "S" = rgb(
    red = 177, green = 99, blue = 110,
    alpha = 255, maxColorValue = 255
  ),
  "all" = "black"
)

i <- 0

for (obs.window in obs.window.list) {
  print(paste0(
    "<======== Obs. window: ",
    obs.window / carrington.rotation, " ========>"
  ))

  i <- i + 1
  if (save.plots) {
    pdf(paste0(
      "sunspots/results/sunspots_cir_", statistic,
      "_test_scatter_", obs.window / carrington.rotation,
      "_bw_", if (adjust.bw) "automated" else bw0,
      "_level_set_", alpha_level_set * 100, ".pdf"
    ),
    width = 25, height = 13
    )
  }

  # Plot overall layout
  par(mfcol = c(4, length(cycles)))
  par(
    oma = c(4, 4, .5, .5),
    mgp = c(2, .6, 0)
  )

  for (cycle in cycles) {
    print(paste0("<====== Cycle: ", cycle, " ======>"))

    data <- rotasym::sunspots_births
    cycle.filter <- data$cycle == cycle
    data <- data[cycle.filter, ]
    scaler_date <- scale(data$date)
    scaled_date <- c(scaler_date)

    first.date <- min(data$date)
    last.date <- max(data$date)

    first.date.z <- min(scaled_date)
    last.date.z <- max(scaled_date)


    ############################
    # Non-uniformity test plot #
    ############################

    par(mar = c(
      -0.1, if (cycle == 21) 4 else -0.1, -0.1,
      if (cycle == 23) 2 else -0.1
    ) + 0.1)

    plot(c(first.date, last.date), c(0, 1),
      type = "n",
      ylab = "Adjusted p-values",
      axes = FALSE
    )
    if (cycle == 21) axis(2L)
    box()

    for (hemisphere in c("N", "S", "all")) {

      # Non-uniformity test p-values
      load(paste0(
        "sunspots/results/",
        cycle, "/sunspots_cir_test_",
        hemisphere, "_",
        obs.window / carrington.rotation, "_results", ".RData"
      ))

      # Adjust p-values using Benjamini and Yekutielli FDR
      adj.results <- p.adjust(test.results[, statistic],
        method = "BY"
      )

      lines(test.results$date - obs.window / (2 * carrington.rotation),
        adj.results,
        col = line.colors[[hemisphere]],
        pch = if (hemisphere == "N") 19 else if (hemisphere == "S") 17 else 4,
        cex = 1,
        type = "b",
        lty = if (hemisphere != "all") 1 else 2
      )
    }

    abline(h = 0.1, col = "black", lty = 2)
    abline(h = 0.05, col = "darkgray", lty = 2)

    if (cycle == 21) {
      legend(
        legend = c("North", "South", "Both"),
        bg = "white",
        pch = c(19, 17, 4),
        lty = c(1, 1, 2),
        col = unlist(line.colors),
        "topleft",
        inset = 0.012,
        title = "Hemisphere"
      )
    }

    #############################
    #   Circ-lin density plot   #
    #############################

    dens.bw <- list("N" = NULL, "S" = NULL, "all" = NULL)

    for (hemisphere in c("N", "S", "all")) {
      print(hemisphere)

      # Filter data by hemisphere
      if (hemisphere == "N") {
        mask <- data$phi > 0
      } else if (hemisphere == "S") {
        mask <- data$phi < 0
      } else {
        mask <- TRUE
      }
      data_hemis <- data[mask, ]
      data_hemis_scaler <- (as.double(data_hemis$date)
      - attr(scaler_date, "scaled:center")
      ) / attr(scaler_date, "scaled:scale")
      data_hemis$z <- data_hemis_scaler

      # Compute joint level sets using rule of thumb bandwidth selectors
      data_cir <- sdetorus::toPiInt(data_hemis$theta)
      data_lin <- data_hemis$z
      data_date <- data_hemis$date
      h <- DirStats::bw_dir_rot(to_cir(data_cir))
      g <- ks::hns(data_lin)
      alpha <- seq(0, 1, by = 0.1)
      c_alpha <- get_c_alpha(
        data_cir = data_cir, data_lin = data_lin,
        h = h, g = g,
        alpha = alpha
      )
      c_alpha <- c_alpha + c(0.1, rep(0, length(alpha) - 1))

      # Plot
      L <- 100
      th <- seq(-1.25 * pi, 1.25 * pi, l = 1.25 * L)
      z <- seq(first.date.z, last.date.z, l = L)
      dates <- seq(first.date, last.date, l = L)
      thz <- as.matrix(expand.grid(th, z))
      dens.bw[[hemisphere]] <- (g * attr(scaler_date, "scaled:scale")
      ) / carrington.rotation

      if (hemisphere != "all") {
        kde_vec <- kde_cir_lin(
          th = thz[, 1], z = thz[, 2],
          data_cir = data_cir, data_lin = data_lin,
          h = h, g = g
        )

        kde_mat <- matrix(kde_vec, nrow = length(th), ncol = length(z))

        par(mar = c(
          -0.1, if (cycle == 21) 4 else -0.1, -0.1,
          if (cycle == 23) 2 else -0.1
        ) + 0.1)

        plot(c(first.date.z, last.date.z), c(-pi, pi),
          type = "n",
          ylab = if (hemisphere == "N") {
            "Northern longitudes (º)"
          } else {
            "Southern longitudes(º)"
          },
          axes = FALSE
        )
        if (cycle == 21) {
          axis(2L,
            at = c(-pi, -pi / 2, 0, pi / 2, pi),
            c("-180", "-90", "0", "90", "180")
          )
        }

        col <- c(
          sequential_hcl(
            length(c_alpha) - 1,
            if (hemisphere == "N") "Blues 3" else "Reds 3"
          ),
          "#FFFFFF"
        )
        .filled.contour(z, th, t(kde_mat),
          levels = rev(c_alpha),
          col = rev(col)
        )
        contour(
          x = z, y = th, z = t(kde_mat),
          levels = c_alpha[2:(length(c_alpha) - 1)],
          ylim = c(-pi, pi),
          col = "black",
          labels = paste(alpha[2:(length(c_alpha) - 1)] * 100, "%"),
          axes = FALSE,
          add = TRUE
        )
        box()

        points(data_lin, data_cir,
          col = colors[[hemisphere]],
          pch = if (hemisphere == "N") 19 else if (hemisphere == "S") 17 else 4,
          cex = 0.8
        )
      }
    }

    ############################
    # Scatter and NW regr plot #
    ############################

    par(mar = c(
      1, if (cycle == 21) 4 else -0.1, -0.1,
      if (cycle == 23) 2 else -0.1
    ) + 0.1)

    plot(c(first.date, last.date), c(-pi, pi),
      type = "n",
      ylab = "Longitudes (º)",
      xlab = "Date",
      yaxt = "n"
    )
    if (cycle == 21) {
      axis(2L,
        at = c(-pi, -pi / 2, 0, pi / 2, pi),
        c("-180", "-90", "0", "90", "180")
      )
    }

    for (hemisphere in c("N", "S", "all")) {
      print(hemisphere)

      # Filter data by hemisphere
      if (hemisphere == "N") {
        mask <- data$phi > 0
      } else if (hemisphere == "S") {
        mask <- data$phi < 0
      } else if (hemisphere == "all") {
        mask <- TRUE
      }

      theta <- toPiInt(circular(data$theta[mask], units = "radians"))
      date <- (data$date[mask])

      if (hemisphere != "all") {
        points(date, theta,
          col = colors[[hemisphere]],
          pch = if (hemisphere == "N") 19 else if (hemisphere == "S") 17 else 4,
          cex = 0.8
        )
      }

      # Circular non-parametric regression
      bw0 <- if (adjust.bw) dens.bw[[hemisphere]] else bw0
      kre0 <- kern.reg.lin.circ(as.double(date) / carrington.rotation, theta,
        method = "NW", bw = bw0, len = 5000
      )
      print(bw0)

      y <- kre0$y
      neg <- y < 0
      y[neg] <- y[neg] + 2 * pi

      linesCirc(as.POSIXct(kre0$x * carrington.rotation + bw0 * carrington.rotation / 2,
        origin = "1970-01-01"
      ),
      y,
      col = line.colors[[hemisphere]],
      lwd = 5,
      lty = if (hemisphere != "all") {
        1
      } else {
        2
      }
      )
    }
  }

  if (save.plots) dev.off()
}
