M <- 1e7
M.name <- "1e7"

p.list <- c(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 21, 51, 101, 151)
alpha.list <- c(0.25, 0.2, 0.15, 0.1, 0.05, 0.02, 0.01)

n.list <- c(5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 100, 200, 300)
n.list.MC <- c(5, 10, 20)
plot.errors <- TRUE

colors <- c(
  "#E74C3C", "#9B59B6", "#3498DB", "#27AE60",
  "#F1C40F", "#E67E22", "#5D6D7E"
)

for (statistic in c("PCvM", "Bakshaev", "PAD")) {
  print(paste0("<======== ", statistic, " ========>"))
  errors.summary <- matrix(
    nrow = length(p.list), ncol = length(alpha.list),
    dimnames = list(p.list, alpha.list)
  )
  errors.summary.MC <- matrix(
    nrow = length(p.list), ncol = length(alpha.list),
    dimnames = list(p.list, alpha.list)
  )

  for (p in p.list) {
    load(file = paste0("results/crit_val_MC_", M, "_", statistic, ".RData"))
    mod.results.ours <- mod.results[as.character(n.list), , as.character(p)]

    alpha <- as.numeric(colnames(mod.results.ours))
    alpha.matrix <- alpha[col(mod.results.ours)]

    load(file = paste0("results/crit_val_MC_", M.name, "_", statistic, "_MC.RData"))
    mod.results.MC <- save.res[as.character(n.list.MC), , as.character(p), ]
    alpha.matrix.MC <- alpha.list[col(mod.results.MC[, , 1])]
    error.MC.m <- array(0,
      dim = dim(mod.results.MC),
      dimnames = dimnames(mod.results.MC)
    )
    for (m in 1:dim(mod.results.MC)[3]) {
      error.MC.m[, , m] <- abs(mod.results.MC[, , m] - alpha.matrix.MC)
    }

    # Absolute errors
    error <- abs(mod.results.ours - alpha.matrix)
    error.MC <- apply(error.MC.m, c(1, 2), function(x) mean(x, na.rm = TRUE))

    # Relative errors
    rel.error <- error / alpha.matrix * 100
    rel.error.MC <- error.MC / alpha.matrix.MC * 100

    errors.summary[as.character(p), ] <- apply(
      rel.error, 2,
      mean
    )[as.character(alpha.list)]
    errors.summary.MC[as.character(p), ] <- apply(
      rel.error.MC, 2,
      function(x) {
        mean(x,
          na.rm = TRUE
        )
      }
    )[as.character(alpha.list)]
  }

  print(errors.summary)
  print(errors.summary.MC)

  if (plot.errors) {

    # Plot setup
    pdf(
      file = paste0("results/rel_error_p_", M.name, "_", statistic, ".pdf"),
      width = 6, # 8
      height = 6
    ) # 6.7

    par(
      mfrow = c(1, 1),
      xpd = FALSE
    )

    plot(range(as.numeric(rownames(errors.summary))),
      c(0, 10.2),
      type = "n",
      ylab = "Relative Error (%)",
      xlab = expression(p),
      log = "x",
      yaxt = "n"
    )
    ytick <- 0:10
    axis(side = 2, at = ytick, labels = FALSE)
    text(2, ytick + 0.2,
      labels = ytick,
      adj = 0.5, pos = 2, offset = 2,
      srt = 90, xpd = TRUE
    )

    i <- 0
    for (alpha in alpha.list) {
      i <- i + 1

      errors <- errors.summary[, as.character(alpha)]
      errors.MC <- errors.summary.MC[, as.character(alpha)]
      lines(as.numeric(names(errors)), errors,
        pch = 16,
        col = colors[i],
        cex = 1,
        lty = 2,
        type = "b"
      )

      # MC regression lines (constant lines)
      err <- errors.summary.MC[, as.character(alpha)]
      p.err <- as.numeric(names(err))
      MC.lin <- lm(err ~ p.err)
      errors.line.MC <- predict(MC.lin, data.frame(list("n.err" = p.list)))

      # empirical points
      points(as.numeric(names(errors.MC)), errors.MC,
        pch = 0,
        col = colors[i],
        cex = 1
      )
      # fitted line
      lines(as.numeric(p.list), errors.line.MC,
        col = colors[i],
        lty = 3,
        cex = 1
      )
    }

    if (statistic == "PCvM") {
      leg.1 <- paste(
        "expression(paste(alpha, \": \", ",
        format(alpha.list, nsmall = 2), "))"
      )
      leg.2 <- paste(leg.1, collapse = ", ")
      leg <- paste("c(", leg.2, ")", sep = "")

      legend("top",
        legend = c(
          "", "", "",
          "Ours", expression(paste("Monte Carlo")),
          c(
            expression(paste(alpha, ": ", "0.25")),
            expression(paste(alpha, ": ", "0.20")),
            expression(paste(alpha, ": ", "0.15")),
            expression(paste(alpha, ": ", "0.10")),
            expression(paste(alpha, ": ", "0.05")),
            expression(paste(alpha, ": ", "0.02")),
            expression(paste(alpha, ": ", "0.01"))
          ), "", "", ""
        ),
        col = c(rep("white", 3), rep("black", 2), colors, rep("white", 3)),
        bg = "white",
        pch = c(rep(1, 3), 16, 0, rep(16, length(alpha.list)), rep(1, 3)),
        lty = c(rep(1, 3), 2, 3, rep(0, length(alpha.list)), rep(1, 3)),
        ncol = 5
      )
    }

    dev.off()
  }
}
