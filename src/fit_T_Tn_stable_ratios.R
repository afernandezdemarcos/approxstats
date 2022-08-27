library(tidyr)
library(glmnet)
library(MASS)
library(stringr)

source("src/approxstats_utils.R")

# Data preparation. Returns a data frame with the following columns
# n | alpha | T_Tn(T_{\infty;\alpha}/T_{\n;\alpha})
# | T_inf (T_{\infty;\alpha}) | T_n (T_{\n;\alpha})
# | f_T_n (pdf[T_{\n}, alpha]) | n^{-1/2} | n^{-1} | ... | n^{-lambda/2}
# | alpha^{-1/2} | alpha^{-1} | ... | alpha^{-mu/2}
prepare_data <- function(T_Tn, T_inf, f_T_n, T_n,
                         N.start, N.end,
                         fit.alpha,
                         lambda = 2, mu = 2) {
  N <- get_N_index(rownames(T_Tn), N.start, N.end) # Sample sizes
  N.index <- N$index
  N.num <- N$num

  data <- data.frame(n = N.num, T_Tn = T_Tn[N.index, ])
  alpha <- as.numeric(colnames(T_Tn)) # alpha-quantiles
  colnames(data) <- c("n", alpha) # Set column names of data
  data$n <- factor(data$n)
  data <- gather(data, alpha, T_Tn,
    as.character(alpha[1]):as.character(alpha[length(alpha)]),
    factor_key = TRUE
  )

  data_T_inf <- data.frame(n = N.num, T_inf = T_inf[N.index, ])
  data_T_n <- data.frame(n = N.num, T_n = T_n[N.index, ])
  data_f_T_n <- data.frame(n = N.num, f_T_n = f_T_n[N.index, ])
  colnames(data_T_inf) <- c("n", alpha) # Set column names of data
  colnames(data_T_n) <- c("n", alpha) # Set column names of data
  colnames(data_f_T_n) <- c("n", alpha) # Set column names of data
  data_T_inf$n <- factor(data_T_inf$n)
  data_T_n$n <- factor(data_T_n$n)
  data_f_T_n$n <- factor(data_f_T_n$n)
  data_T_inf <- gather(data_T_inf, alpha, T_inf,
    as.character(alpha[1]):as.character(alpha[length(alpha)]),
    factor_key = TRUE
  )
  data_T_n <- gather(data_T_n, alpha, T_n,
    as.character(alpha[1]):as.character(alpha[length(alpha)]),
    factor_key = TRUE
  )
  data_f_T_n <- gather(data_f_T_n, alpha, f_T_n,
    as.character(alpha[1]):as.character(alpha[length(alpha)]),
    factor_key = TRUE
  )

  data <- (data %>% inner_join(data_T_inf, by = c("n", "alpha"))
    %>% inner_join(data_T_n, by = c("n", "alpha"))
    %>% inner_join(data_f_T_n, by = c("n", "alpha")))

  # Polynomial forms of n
  data$n <- as.numeric(levels(data$n))[data$n]
  for (l in 1:lambda) {
    eval(parse(text = paste0("data$inv_n_", l, "_2 <- (data$n)^(-", l, "/2)")))
  }

  data$alpha <- as.numeric(levels(data$alpha))[data$alpha]
  if (fit.alpha) {
    # Polynomial forms of alpha
    for (m in 1:mu) {
      eval(parse(text = paste0("data$inv_alpha_", m, "_2 <- (data$alpha)^(-", m, "/2)")))
    }
  }

  return(data)
}

# Weighted least squares regression of T_Tn (T_{\infty;\alpha}/T_{\n;\alpha})
# using model selection based on optimal measure k
# with search direction specified by fit.direction
# with initial searching model the one at model.reference
# and weight function specified by weights.fun
#
# Additional tuning parameters:
#   - fit.direction: c("all", "forward", "both")
#   - (N.start, N.end): Specify the sample sizes to fit the model
#   - fit.alpha: Use alpha and its powers as predictors (TRUE) or not (FALSE).
#       If FALSE, then a particular alpha value must be specified in which.alpha
#   - fit.p: Whether the dimension of ambient space R^p and
#       its powers are used as predictors. If TRUE, data must be given
#   - data: Only required for (n, alpha, p)-model (i.e., if fit.p is TRUE)
#   - model.reference: If "Stephens", the one at Stephens(1970) (n^{-1/2}, n^{-1});
#       if "Lasso", the initial model is one fitted with Lasso predictors;
#       if an object lm() is provided, is used as starting point;
#       if a string is provided, it is interpreted as a formula for lm().
#   - which.alpha: If fit.alpha is FALSE, then a particular alpha value
#       must be given to fit this particular model
#   - k: Optimality criteria for model search, "BIC" or "AIC".
#   - weights.fun: Function for computing weights.
#       If NULL, no weights are applied.
#   - trace: Indicates whether to print the model search procedure (T) or not (F)
#   - model.formula: Specifies a particular formula for the model
#       (Used in drop one-by-one terms to simplify model)
#   - lambda: Maximum power (-lambda/2) to include as predictor of n
#   - mu: Maximum power (-mu/2) to include as predictor of alpha
#
# Output: List consisting on the following elements
#   - data: Data used to fit the model
#   - Stephens: Estimated Stephens (1970) model if used as reference.model
#   - both: Estimated optimal model when searching direction is "both"
#   - forward: Estimated optimal model when searching direction is "forward"
#   - Lasso: Estimated Lasso initial model (fit, lambda.1se, lambda.min)
#   - particular: Estimated particular model when model.formula is specified.
#
fit_T_Tn <- function(T_Tn, T_n, T_inf, f_T_n,
                     fit.direction = "all",
                     N.start = 3, N.end = nrow(T_Tn),
                     fit.alpha = TRUE,
                     fit.p = FALSE,
                     data = FALSE,
                     model.reference = "Stephens",
                     which.alpha = 0.05,
                     k = "BIC",
                     weights.fun = NULL,
                     trace = FALSE,
                     model.formula = NULL,
                     lambda = 2,
                     mu = 2) {

  # Output list
  results <- list(
    data = NULL,
    Stephens = NULL,
    both = NULL,
    forward = NULL,
    Lasso = NULL,
    particular = NULL
  )

  # Prepare data from quantiles ratios T_Tn, quantiles T_inf, T_n,
  # pdf f_T_n, where n ranges from N.start to N.end sample sizes.
  if (!fit.p) {
    data <- prepare_data(
      T_Tn = T_Tn,
      T_inf = T_inf,
      f_T_n = f_T_n,
      T_n = T_n,
      N.start = N.start,
      N.end = N.end,
      fit.alpha = fit.alpha,
      lambda = lambda,
      mu = mu
    )
    results[["data"]] <- data
  } else {
    if (is.null(data)) stop("data is required to fit p")
  }

  # Prepare weights
  weights.vec <- (
    if (!fit.p) {

      # (n, alpha)-model
      get_weights_from_function_n(weights.fun,
        n_vector = data$n,
        alpha_vector = data$alpha,
        T_inf_vector = data$T_inf,
        T_n_vector = data$T_n,
        f_Tn_vector = data$f_T_n
      )
    } else {

      # (n, alpha, p)-model
      get_weights_from_function_n(weights.fun,
        n_vector = data$n,
        alpha_vector = data$alpha,
        p_vector = data$p
      )
    })

  # Specify formula for the saturated model

  # Drop n, alpha, and arguments used for weight computation from predictors
  n.alpha.ind <- c("n", "alpha", "T_inf", "T_n", "f_T_n")
  vars.wo.n <- names(data)[!(names(data) %in% n.alpha.ind)]

  # Interactions for n and alpha blocks
  vars <- str_split(vars.wo.n, "_")
  n.var.set <- c()
  alpha.var.set <- c()
  for (i in 1:length(vars)) {
    var <- vars[[i]]
    var.name <- vars.wo.n[i]
    if ("n" %in% var) {
      n.var.set <- c(n.var.set, var.name)
    }
    if ("alpha" %in% var) {
      alpha.var.set <- c(alpha.var.set, var.name)
    }
  }

  form.n <- paste(n.var.set, collapse = " + ")
  form.alpha <- paste(alpha.var.set, collapse = " + ")

  # Select data from which.alpha and adjust saturated model
  if (!fit.alpha) {

    # Model only with sample size predictors (n)
    data <- data[data$alpha == as.character(which.alpha), ]
    saturated.form <- formula(paste("T_Tn ~", form.n))
  } else {

    # (n, alpha)-model
    saturated.form <- formula(paste0(
      "T_Tn ~ ",
      "0 + ", form.n,
      " + (", form.n, "):(", form.alpha, ")"
    ))

    if (fit.p) {

      # (n, alpha, p)-model
      if (is.null(model.formula)) {
        saturated.form <- formula(model.reference)
      } else {
        saturated.form <- formula("T_Tn ~ 1")
      }
    }
  }

  # Set intercept to 1
  data$T_Tn <- data$T_Tn - 1

  # Specify reference model

  # Dummy model
  T_Tn_Zero <- lm(formula = T_Tn ~ 0, data, weights = weights.vec)

  # Saturated model
  T_Tn_All <- lm(formula = saturated.form, data, weights = weights.vec)

  # Search initial model
  if (model.reference == "Stephens") {

    # Stephens (1970) model (n^{-1/2}, n*{-1})
    T_Tn_initial <- lm(
      formula = formula("T_Tn ~ 0 + inv_n_2_2 + inv_n_1_2"),
      data,
      weights = weights.vec
    )

    results[["Stephens"]] <- T_Tn_initial
  } else if (model.reference == "Lasso") {

    # Lasso predictors as initial model

    # Fit Lasso model by CV (Leave-One-Out)
    lasso.fit <- fit_lasso_T_Tn(saturated.form,
      data = data,
      weights = weights.vec
    )

    # Get formula from list of predictors selected by Lasso
    lasso.preds <- get_lasso_predictors(lasso.fit = lasso.fit)
    T_Tn_initial.form <- formula(paste(
      "T_Tn ~ 0 + ",
      paste(lasso.preds,
        collapse = " + "
      )
    ))

    T_Tn_initial <- lm(T_Tn_initial.form,
      data = data,
      weights = weights.vec
    )

    results[["Lasso"]] <- list(
      kcv = lasso.fit$fit,
      lambda.1se = lasso.fit$lambda.1se,
      lambda.min = lasso.fit$lambda.min
    )
  } else {
    if (typeof(model.reference) == "list") {
      if (is.null(model.reference$call[1])) {
        stop("Reference model must be an lm() fit")
      }

      if (model.reference$call[1] != "lm()") {
        stop("Reference model must be an lm() fit")
      }

      # Particular linear model lm() as a starting model
      T_Tn_initial <- model.reference
    } else if (typeof(model.reference) == "character") {

      # Particular formula to build a linear model lm() as a starting model
      T_Tn_initial <- lm(formula(model.reference),
        data = data,
        weights = weights.vec
      )
    } else {
      stop("Reference model must be an lm() fit")
    }
  }

  if (!is.null(model.formula)) {

    # In case a particular formula is entered
    T_Tn_fit <- lm(formula(model.formula),
      data = data,
      weights = weights.vec
    )
    results[["particular"]] <- T_Tn_fit
  } else {

    # Fitting direction for stepwise regression
    fit.direction <- get_stepAIC_direction(fit.direction)
    k <- get_stepAIC_penalty(k, data) # Penalty measure

    for (direction in fit.direction) {
      T_Tn_fit <- MASS::stepAIC(T_Tn_initial,
        direction = direction,
        scope = list(
          lower = T_Tn_Zero,
          upper = T_Tn_All
        ),
        k = k,
        trace = trace
      )

      results[[direction]] <- T_Tn_fit
    }
  }

  return(results)
}

# Automated process of searching minimum lambda penalty to avoid edge problems
extend_lambda_grid <- function(kcvLasso) {
  upper.lambda <- kcvLasso$lambda[1]
  lower.lambda <- kcvLasso$lambda[length(kcvLasso$lambda)]

  min.lambda.pos <- which(kcvLasso$lambda.min == kcvLasso$lambda)

  # Check if minimum is at a extreme
  minlambda.extreme.lower <- min.lambda.pos == length(kcvLasso$lambda)
  minlambda.extreme.upper <- min.lambda.pos == 1
  minlambda.extreme <- minlambda.extreme.lower + minlambda.extreme.upper

  if (minlambda.extreme) {
    if (minlambda.extreme.lower) {
      print("Minimum lambda is on the lower extreme of the grid")

      # Extend lower extreme
      new.lower.loglambda <- 2 * log10(lower.lambda) - log10(upper.lambda)

      lambdaGrid <- 10^seq(log10(kcvLasso$lambda[1]), new.lower.loglambda,
        length.out = 150
      )
    }

    if (minlambda.extreme.upper) {
      print("Minimum lambda is on the upper extreme of the grid")

      # Extend upper extreme
      new.upper.loglambda <- 2 * log10(upper.lambda) - log10(lower.lambda)

      lambdaGrid <- 10^seq(new.upper.loglambda,
        log10(kcvLasso$lambda[length(kcvLasso$lambda)]),
        length.out = 150
      )
    }
    return(lambdaGrid)
  }
}

# Lasso fit using k.folds-CV
fit_lasso_T_Tn <- function(data, saturated.formula, alpha = 1,
                           weights = NULL,
                           k.folds = nrow(data),
                           max.iter.lambda = 3) {

  # Model matrix
  x <- model.matrix(saturated.formula, data = data)
  y <- data$T_Tn

  # k-CV lasso model
  kcvLasso <- cv.glmnet(
    x = x,
    y = y,
    weights = weights,
    alpha = alpha,
    nfolds = k.folds
  )

  # Automated process of searching minimum lambda penalty,
  # until it's not on an extreme
  lambdaGrid <- extend_lambda_grid(kcvLasso)
  iter <- 1
  while ((!is.null(lambdaGrid)) | (iter < max.iter.lambda)) {
    kcvLasso <- cv.glmnet(
      x = x, y = y,
      alpha = 1,
      nfolds = 10,
      lambda = lambdaGrid
    )

    lambdaGrid <- extend_lambda_grid(kcvLasso)

    iter <- iter + 1
  }

  plot(kcvLasso)
  indMin <- which.min(kcvLasso$cvm)
  abline(h = kcvLasso$cvm[indMin] + c(0, kcvLasso$cvsd[indMin]))

  modLassoCV <- kcvLasso$glmnet.fit # Fitted models

  return(list(
    fit = modLassoCV,
    lambda.1se = kcvLasso$lambda.1se,
    lambda.min = kcvLasso$lambda.min
  ))
}

# Returns lasso.fit predictors for an optimal lambda value in
# opt.lambda (1se or min)
get_lasso_predictors <- function(lasso.fit, opt.lambda = "1se") {
  sel.preds <- predict(lasso.fit$fit,
    type = "coefficients",
    s = c(lasso.fit$lambda.min, lasso.fit$lambda.1se)
  )[-1, ] != 0

  # Get model with optimal lambda. If not 1se then take min.lambda
  name.sel.preds <- names(which(
    sel.preds[, ifelse(opt.lambda == "1se", 2, 1)]
  ))

  return(name.sel.preds)
}
