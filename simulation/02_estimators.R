# =============================================================================
# FILE 2: IPW Estimators
# Implements Propositions 1-5 from Huber (2013)
# =============================================================================
#
# WHAT THESE ESTIMATORS DO:
#
# All estimators take propensity scores as inputs and compute weighted averages
# of the outcome Y to recover direct and indirect effects.
#
# Two propensity scores are needed:
#   p_X    = Pr(D=1 | X)         treatment propensity given covariates only
#   p_MX   = Pr(D=1 | M, X)      treatment propensity given mediator + covariates
#   p_WX   = Pr(D=1 | W, X)      treatment propensity given W + covariates
#   p_MWX  = Pr(D=1 | M, W, X)   treatment propensity given M + W + covariates
#
# The "normalized" versions (used throughout) rescale weights so they sum to 1
# within each treatment group. This improves finite-sample performance.
# =============================================================================

library(glmnet)



# -----------------------------------------------------------------------------
# PROPENSITY SCORE ESTIMATION
# Supports: probit, lasso, random forest, gradient boosting
# -----------------------------------------------------------------------------

estimate_pscore <- function(D, X, method = "probit") {
  X <- as.matrix(X)
  
  # Handle edge cases
  if (ncol(X) == 0) {
    return(rep(mean(D), length(D)))
  }
  
  # LASSO needs at least 2 columns
  if (method == "lasso" && ncol(X) < 2) {
    method <- "logit" 
  }
  
  if (method == "probit") {
    df  <- data.frame(D = D, X)
    fit <- glm(D ~ ., data = df, family = binomial(link = "probit"))
    phat <- predict(fit, type = "response")
    
  } else if (method == "logit") {
    df  <- data.frame(D = D, X)
    fit <- glm(D ~ ., data = df, family = binomial(link = "logit"))
    phat <- predict(fit, type = "response")
    
  } else if (method == "lasso") {
    D_fac <- as.factor(D)
    fit  <- cv.glmnet(X, D_fac, family = "binomial", alpha = 1,
                      nfolds = min(5, nrow(X)), type.measure = "deviance")
    phat <- as.vector(predict(fit, newx = X, s = "lambda.min",
                              type = "response"))
  } else {
    stop("Unknown method: choose from probit, logit, lasso")
  }
  
  # Clip to avoid extreme weights
  phat <- pmax(pmin(phat, 0.999), 0.001)
  return(phat)
}



# -----------------------------------------------------------------------------
# FRAMEWORK 1: Assumptions 1 & 2 (gamma = 0, pre-treatment confounders only)
# Uses propensity scores: p_X and p_MX
# -----------------------------------------------------------------------------

# PROPOSITION 1: Average Direct Effect theta(d)
# theta(d) = E[ (Y*D/p(M,X) - Y*(1-D)/(1-p(M,X))) * p(D=d|M,X)/p(D=d|X) ]
ipw_direct_1 <- function(Y, D, p_X, p_MX) {
  # theta(1): keep mediator at M(1), vary treatment
  w1_num <- D / p_MX                              # treated weight
  w0_num <- (1 - D) * p_MX / ((1 - p_MX) * p_X) # control weight, tilted

  theta1 <- sum(Y * D / p_X)    / sum(D / p_X) -
            sum(Y * w0_num)     / sum(w0_num)

  # theta(0): keep mediator at M(0), vary treatment
  w1_num2 <- D * (1 - p_MX) / (p_MX * (1 - p_X))
  w0_num2 <- (1 - D) / (1 - p_X)

  theta0 <- sum(Y * w1_num2) / sum(w1_num2) -
            sum(Y * w0_num2) / sum(w0_num2)

  c(theta_1 = theta1, theta_0 = theta0)
}

# PROPOSITION 2: Average Indirect Effect delta(d)
# delta(d) = E[ Y*I{D=d}/p(D=d|M,X) * (p(M,X)/p(X) - (1-p(M,X))/(1-p(X))) ]
ipw_indirect_1 <- function(Y, D, p_X, p_MX) {
  # delta(1)
  w1 <- D / p_MX
  delta1 <- sum(Y * w1 * (p_MX / p_X - (1 - p_MX) / (1 - p_X))) /
            sum(w1)

  # delta(0)
  w0 <- (1 - D) / (1 - p_MX)
  delta0 <- sum(Y * w0 * (p_MX / p_X - (1 - p_MX) / (1 - p_X))) /
            sum(w0)

  c(delta_1 = delta1, delta_0 = delta0)
}


# -----------------------------------------------------------------------------
# FRAMEWORK 2: Assumptions 3 & 4 (gamma != 0, post-treatment confounders W)
# Uses propensity scores: p_X, p_WX, p_MWX
# -----------------------------------------------------------------------------

# PROPOSITION 3: Average Direct Effect theta(d) with W
ipw_direct_2 <- function(Y, D, p_X, p_MWX) {
  # Same structure as Proposition 1 but conditioning on M, W, X
  w1_num <- D / p_MWX
  w0_num <- (1 - D) * p_MWX / ((1 - p_MWX) * p_X)

  theta1 <- sum(Y * D / p_X)  / sum(D / p_X) -
            sum(Y * w0_num)   / sum(w0_num)

  w1_num2 <- D * (1 - p_MWX) / (p_MWX * (1 - p_X))
  w0_num2 <- (1 - D) / (1 - p_X)

  theta0 <- sum(Y * w1_num2) / sum(w1_num2) -
            sum(Y * w0_num2) / sum(w0_num2)

  c(theta_1 = theta1, theta_0 = theta0)
}

# PROPOSITION 4: Average Partial Indirect Effect delta_p(d)
# Holds W fixed, only captures the direct D -> M -> Y path
ipw_partial_indirect <- function(Y, D, p_X, p_WX, p_MWX) {
  # delta_p(1)
  w1 <- D / p_MWX * (p_WX / p_X)
  delta_p1 <- sum(Y * w1 * (p_MWX / p_WX - (1 - p_MWX) / (1 - p_WX))) /
              sum(w1)

  # delta_p(0)
  w0 <- (1 - D) / (1 - p_MWX) * ((1 - p_WX) / (1 - p_X))
  delta_p0 <- sum(Y * w0 * (p_MWX / p_WX - (1 - p_MWX) / (1 - p_WX))) /
              sum(w0)

  c(delta_p1 = delta_p1, delta_p0 = delta_p0)
}

# PROPOSITION 5: Total Indirect Effect delta_t(d)
# Requires Assumption 5 (linearity of outcome in M)
# Uses OLS regression within each treatment arm + IPW reweighting
ipw_total_indirect <- function(Y, D, M, W, X, p_X) {
  # Step 1: Fit linear model mu(d, M, W, X) separately for D=1 and D=0
  df <- data.frame(Y = Y, M = M, W = W, X = X)

  mu1 <- lm(Y ~ ., data = df[D == 1, ])  # model for treated
  mu0 <- lm(Y ~ ., data = df[D == 0, ])  # model for control

  # Step 2: Estimate E[M(1-d, W(1-d))] using IPW
  # E[M(1)] = weighted mean of M among treated
  EM1 <- sum(M * D / p_X) / sum(D / p_X)
  # E[M(0)] = weighted mean of M among control
  EM0 <- sum(M * (1 - D) / (1 - p_X)) / sum((1 - D) / (1 - p_X))

  # Step 3: Predict counterfactual outcomes
  # For delta_t(1): plug E[M(0,W(0))] into mu1
  newdat1      <- df[D == 1, ]
  newdat1$M    <- EM0
  mu1_pred     <- predict(mu1, newdata = newdat1)

  # delta_t(1) = E[(Y - mu1(E[M0], W, X)) * D/p(X)]
  resid1  <- Y[D == 1] - mu1_pred
  delta_t1 <- sum(resid1 * (D / p_X)[D == 1]) / sum((D / p_X)[D == 1])

  # For delta_t(0): plug E[M(1,W(1))] into mu0
  newdat0   <- df[D == 0, ]
  newdat0$M <- EM1
  mu0_pred  <- predict(mu0, newdata = newdat0)

  resid0    <- mu0_pred - Y[D == 0]
  delta_t0  <- sum(resid0 * ((1 - D) / (1 - p_X))[D == 0]) /
               sum(((1 - D) / (1 - p_X))[D == 0])

  c(delta_t1 = delta_t1, delta_t0 = delta_t0)
}


# -----------------------------------------------------------------------------
# OLS BENCHMARKS (for comparison with IPW)
# -----------------------------------------------------------------------------

ols_effects <- function(Y, D, M, V, U, include_interaction = FALSE) {
  if (include_interaction) {
    # OLS with D-M interaction (OLS.ia in paper)
    fit_Y <- lm(Y ~ D + M + D:M + V + U)
    fit_M <- lm(M ~ D + V + U)

    b_D   <- coef(fit_Y)["D"]
    b_M   <- coef(fit_Y)["M"]
    b_DM  <- coef(fit_Y)["D:M"]
    b_D_M <- coef(fit_M)["D"]

    # Direct effect: coefficient on D + interaction * E[M(d)]
    EM1 <- mean(M[D == 1])
    EM0 <- mean(M[D == 0])

    theta1 <- b_D + b_DM * EM1
    theta0 <- b_D + b_DM * EM0

    # Indirect effect: first stage * (M coefficient + interaction * d)
    delta1 <- b_D_M * (b_M + b_DM)
    delta0 <- b_D_M * b_M

  } else {
    # OLS without interaction (OLS in paper)
    fit_Y <- lm(Y ~ D + M + V + U)
    fit_M <- lm(M ~ D + V + U)

    b_D   <- coef(fit_Y)["D"]
    b_M   <- coef(fit_Y)["M"]
    b_D_M <- coef(fit_M)["D"]

    # Homogeneous direct and indirect effects
    theta1 <- theta0 <- b_D
    delta1 <- delta0 <- b_M * b_D_M
  }

  c(theta_1 = unname(theta1), theta_0 = unname(theta0),
    delta_1 = unname(delta1), delta_0 = unname(delta0))
}

ols_naive <- function(Y, D, M) {
  # Naive OLS: ignores all confounders V and U
  fit_Y <- lm(Y ~ D + M)
  fit_M <- lm(M ~ D)

  b_D   <- coef(fit_Y)["D"]
  b_M   <- coef(fit_Y)["M"]
  b_D_M <- coef(fit_M)["D"]

  c(theta_1 = unname(b_D),
    theta_0 = unname(b_D),
    delta_1 = unname(b_M * b_D_M),
    delta_0 = unname(b_M * b_D_M))
}
