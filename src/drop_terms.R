# Further simplification of the model if R^2_{Adj} 
# does not decrease significantly

library(MASS)

# Iteratively drop one-by-one the predictors of model that contribute 
# the least to the R^2_{Adj}. The maximum allowed decrease to drop a 
# predictor is set by delta.max, relative to the original model with R^2: R2_0.
drop_term_r_squared <- function(model, delta.max = 1e-3, R2_0 = 1,
                                scope = formula(model), trace = FALSE){
  
  # Set delta.max relative to original R^2
  delta.max <- delta.max * R2_0 
  
  # Step AIC drop term
  if(!is.null(scope)){
    drop.terms <- dropterm(model, scope = scope)
  }else{
    drop.terms <- dropterm(model)
  }
  mod.summary <- summary(model)
  
  if(trace){print(drop.terms)}
  
  # Compute SSE and SST for R squared
  sse <- drop.terms$RSS
  sst <- sse[1]/(1 - mod.summary$r.squared)
  
  r.sq <- 1 - sse/sst
  dfr <- model$df
  p <- length(attr(model$terms, "term.labels")); n <- dfr + p + 1
  # -1 since we drop one variable from original model
  r.sq.adj <- 1 - (1 - r.sq) * ((n - 1)/(dfr - 1))
  delta.r.sq <- abs(r.sq.adj - r.sq.adj[1])[2:length(r.sq.adj)]
  
  # Get variable with minimum change in R^2 always less than delta.max
  mask <- delta.r.sq < delta.max
  drop.variables <- rownames(drop.terms)[2:nrow(drop.terms)]
  drop.variables <- drop.variables[mask]
  drop.var <- drop.variables[which.min(delta.r.sq[mask])]
  
  if (length(drop.var) > 0){
    
    pos.drop <- grep(paste0("^", drop.var, "$"), 
                     attr(model$terms, "term.labels"))
    new.formula <- drop.terms(model$terms, pos.drop, keep.response = T)
  
  }else{
    
    new.formula <- NULL
    
  }
  
  return(new.formula)
  
}
