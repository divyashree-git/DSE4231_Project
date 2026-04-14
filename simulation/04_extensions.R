# =============================================================================
# FILE 7: GUARANTEED EXTENSION RUNNER (Serial & Safe)
# Targets: LASSO vs Probit, Gamma Sensitivity, Sample Size Analysis
# =============================================================================

source("simulation/01_dgp.R")
source("simulation/02_estimators.R")

library(dplyr)
library(tidyr)
library(glmnet)

#1.PROPENSITY SCORE ESTIMATOR

safe_estimate_pscore <- function(D, X, method = "probit") {
  X_mat <- as.matrix(X)
  
  # Fallback to logit if LASSO has only 1 variable
  if (method == "lasso" && ncol(X_mat) < 2) {
    method <- "logit" 
  }
  
  if (method == "probit") {
    df  <- data.frame(D = D, X_mat)
    fit <- glm(D ~ ., data = df, family = binomial(link = "probit"))
    phat <- predict(fit, type = "response")
    
  } else if (method == "logit") {
    df  <- data.frame(D = D, X_mat)
    fit <- glm(D ~ ., data = df, family = binomial(link = "logit"))
    phat <- predict(fit, type = "response")
    
  } else if (method == "lasso") {
    # CRITICAL FIX: glmnet requires D to be a factor for binomial classification
    fit  <- cv.glmnet(X_mat, as.factor(D), family = "binomial", alpha = 1, nfolds = 5)
    phat <- as.numeric(predict(fit, newx = X_mat, s = "lambda.min", type = "response"))
  }
  
  # Clip to prevent extreme weights dividing by zero
  return(pmax(pmin(phat, 0.999), 0.001))
}

# ── 2. SIMULATION LOGIC ────────────────────────────────────────────────────
run_one_sim_flexible <- function(gamma, beta, n_obs, method = "probit") {
  dat <- generate_data(n_obs, gamma = gamma, beta = beta)
  Y <- dat$Y; D <- dat$D; M <- dat$M; V <- dat$V; U <- dat$U
  
  if (gamma == 0) {
    p_X  <- safe_estimate_pscore(D, cbind(V, U), method)
    p_MX <- safe_estimate_pscore(D, cbind(M, V, U), method)
    
    ipw_dir <- ipw_direct_1(Y, D, p_X, p_MX)
    ipw_ind <- ipw_indirect_1(Y, D, p_X, p_MX)
    
    return(c(
      theta_1 = unname(ipw_dir["theta_1"]), theta_0 = unname(ipw_dir["theta_0"]),
      delta_t1 = unname(ipw_ind["delta_1"]), delta_t0 = unname(ipw_ind["delta_0"]),
      delta_p1 = unname(ipw_ind["delta_1"]), delta_p0 = unname(ipw_ind["delta_0"])
    ))
  } else {
    W <- V
    p_X   <- safe_estimate_pscore(D, matrix(U, ncol=1), method)
    p_WX  <- safe_estimate_pscore(D, cbind(W, U), method)
    p_MWX <- safe_estimate_pscore(D, cbind(M, W, U), method)
    
    ipw_dir  <- ipw_direct_2(Y, D, p_X, p_MWX)
    ipw_part <- ipw_partial_indirect(Y, D, p_X, p_WX, p_MWX)
    ipw_tot  <- ipw_total_indirect(Y, D, M, W, U, p_X) # Passed U directly for safety
    
    return(c(
      theta_1 = unname(ipw_dir["theta_1"]), theta_0 = unname(ipw_dir["theta_0"]),
      delta_t1 = unname(ipw_tot["delta_t1"]), delta_t0 = unname(ipw_tot["delta_t0"]),
      delta_p1 = unname(ipw_part["delta_p1"]), delta_p0 = unname(ipw_part["delta_p0"])
    ))
  }
}

# Analytical True Values (Appendix A.3)
get_true_theta1  <- function(g) { 0.5 + 0.5 * (0.5 + 0.5 * g) }
get_true_theta0  <- function(g) { 0.5 }
get_true_deltat1 <- function(g) { 1.5 * (0.5 + 0.5 * g) }
get_true_deltap1 <- function(g) { 1.5 * 0.5 }

# ── 3. EXPERIMENT EXECUTION ────────────────────────────────────────────────
sample_sizes <- c(500, 1000, 2000)
gamma_values <- c(0, 0.2, 0.4)
methods      <- c("probit", "lasso")

ext_grid <- expand.grid(N = sample_sizes, gamma = gamma_values, method = methods, stringsAsFactors = FALSE)
ext_grid$theta_1_bias <- NA; ext_grid$theta_1_MSE <- NA
ext_grid$delta_t1_bias <- NA; ext_grid$delta_t1_MSE <- NA

N_REPS <- 100 # Safe, fast number to guarantee completion

cat("\n=================================================================\n")
cat("STARTING SERIAL EXTENSION RUN (This will take ~60 seconds)\n")
cat("=================================================================\n")

for (i in 1:nrow(ext_grid)) {
  curr <- ext_grid[i, ]
  cat(sprintf("[%d/%d] N=%d | Gamma=%.1f | Method=%s ... ", 
              i, nrow(ext_grid), curr$N, curr$gamma, curr$method))
  
  set.seed(42 + i)
  
  # Because there is no tryCatch, if this breaks, you will see the exact error immediately
  sim_results <- lapply(1:N_REPS, function(s) {
    run_one_sim_flexible(gamma = curr$gamma, beta = 0.5, n_obs = curr$N, method = curr$method)
  })
  
  # Convert list to dataframe safely
  res_df <- do.call(rbind, sim_results)
  
  # Calculate Metrics (using na.rm=TRUE to be absolutely safe)
  t1_ests <- res_df[, "theta_1"]
  dt1_ests <- res_df[, "delta_t1"]
  
  ext_grid$theta_1_bias[i] <- mean(t1_ests, na.rm=TRUE) - get_true_theta1(curr$gamma)
  ext_grid$theta_1_MSE[i]  <- mean((t1_ests - get_true_theta1(curr$gamma))^2, na.rm=TRUE)
  
  ext_grid$delta_t1_bias[i] <- mean(dt1_ests, na.rm=TRUE) - get_true_deltat1(curr$gamma)
  ext_grid$delta_t1_MSE[i]  <- mean((dt1_ests - get_true_deltat1(curr$gamma))^2, na.rm=TRUE)
  
  cat("Done.\n")
}

# ── 4. TABLE GENERATION ────────────────────────────────────────────────────
cat("\n=================================================================\n")

# TABLE 1: ML Comparison (LASSO vs Probit)
table1_ml <- ext_grid %>%
  filter(N == 2000) %>%
  group_by(method) %>%
  summarise(
    Avg_Theta1_Bias = round(mean(theta_1_bias), 4),
    Avg_Theta1_MSE  = round(mean(theta_1_MSE), 4),
    Avg_DeltaT1_Bias = round(mean(delta_t1_bias), 4),
    Avg_DeltaT1_MSE  = round(mean(delta_t1_MSE), 4)
  )

# TABLE 2: Gamma Sensitivity
table2_gamma <- ext_grid %>%
  group_by(gamma, method) %>%
  summarise(
    Theta1_Bias = round(mean(theta_1_bias), 4),
    DeltaT1_Bias = round(mean(delta_t1_bias), 4),
    .groups = "drop"
  ) %>%
  arrange(method, gamma)

# TABLE 3: Sample Size Analysis
table3_N <- ext_grid %>%
  group_by(N, method) %>%
  summarise(
    Theta1_MSE = round(mean(theta_1_MSE), 4),
    DeltaT1_MSE = round(mean(delta_t1_MSE), 4),
    .groups = "drop"
  ) %>%
  arrange(method, N)

cat("\n--- TABLE 1: PROBIT vs LASSO (at N=2000) ---\n")
print(as.data.frame(table1_ml))

cat("\n--- TABLE 2: GAMMA SENSITIVITY (Bias by Gamma) ---\n")
print(as.data.frame(table2_gamma))

cat("\n--- TABLE 3: SAMPLE SIZE ANALYSIS (MSE by N) ---\n")
print(as.data.frame(table3_N))

# Save to output
dir.create("output/tables", showWarnings = FALSE)
write.csv(table1_ml, "output/tables/1_ml_comparison.csv", row.names=FALSE)
write.csv(table2_gamma, "output/tables/2_gamma_sensitivity.csv", row.names=FALSE)
write.csv(table3_N, "output/tables/3_sample_size.csv", row.names=FALSE)
cat("\n✓ Tables saved successfully to output/tables/.\n")