
source("simulation/01_dgp.R")
source("simulation/02_estimators.R")

library(dplyr)

N_SIM <- 5000
N_OBS <- 2000

run_one_sim <- function(gamma, beta, method = "probit") {
  dat <- generate_data(N_OBS, gamma = gamma, beta = beta)
  Y <- dat$Y; D <- dat$D; M <- dat$M; V <- dat$V; U <- dat$U

  if (gamma == 0) {
    X_base <- cbind(V, U)
    X_MX   <- cbind(M, V, U)
    p_X  <- estimate_pscore(D, X_base, method = method)
    p_MX <- estimate_pscore(D, X_MX,   method = method)
    ipw_dir <- ipw_direct_1(Y, D, p_X, p_MX)
    ipw_ind <- ipw_indirect_1(Y, D, p_X, p_MX)
    ols_ia  <- ols_effects(Y, D, M, V, U, include_interaction = TRUE)
    ols_    <- ols_effects(Y, D, M, V, U, include_interaction = FALSE)
    ols_nav <- ols_naive(Y, D, M)
    return(list(
      theta_1_ipw=ipw_dir["theta_1"], theta_0_ipw=ipw_dir["theta_0"],
      delta_1_ipw=ipw_ind["delta_1"], delta_0_ipw=ipw_ind["delta_0"],
      theta_1_ols_ia=ols_ia["theta_1"], theta_0_ols_ia=ols_ia["theta_0"],
      delta_1_ols_ia=ols_ia["delta_1"], delta_0_ols_ia=ols_ia["delta_0"],
      theta_1_ols=ols_["theta_1"], theta_0_ols=ols_["theta_0"],
      delta_1_ols=ols_["delta_1"], delta_0_ols=ols_["delta_0"],
      theta_1_naive=ols_nav["theta_1"], theta_0_naive=ols_nav["theta_0"],
      delta_1_naive=ols_nav["delta_1"], delta_0_naive=ols_nav["delta_0"]
    ))
  } else {
    W <- V
    X_base <- matrix(U, ncol=1)
    X_WX   <- cbind(W, U)
    X_MWX  <- cbind(M, W, U)
    p_X   <- estimate_pscore(D, X_base, method=method)
    p_WX  <- estimate_pscore(D, X_WX,   method=method)
    p_MWX <- estimate_pscore(D, X_MWX,  method=method)
    ipw_dir  <- ipw_direct_2(Y, D, p_X, p_MWX)
    ipw_part <- ipw_partial_indirect(Y, D, p_X, p_WX, p_MWX)
    ipw_tot  <- ipw_total_indirect(Y, D, M, W, matrix(U, ncol=1), p_X)
    p_MX_wrong <- estimate_pscore(D, cbind(M, U), method=method)
    ipw_wrong  <- ipw_indirect_1(Y, D, p_X, p_MX_wrong)
    ols_ia  <- ols_effects(Y, D, M, V, U, include_interaction=TRUE)
    ols_    <- ols_effects(Y, D, M, V, U, include_interaction=FALSE)
    ols_nav <- ols_naive(Y, D, M)
    return(list(
      theta_1_ipw=ipw_dir["theta_1"], theta_0_ipw=ipw_dir["theta_0"],
      theta_1_ols_ia=ols_ia["theta_1"], theta_0_ols_ia=ols_ia["theta_0"],
      theta_1_ols=ols_["theta_1"], theta_0_ols=ols_["theta_0"],
      theta_1_naive=ols_nav["theta_1"], theta_0_naive=ols_nav["theta_0"],
      delta_t1_correct=ipw_tot["delta_t1"], delta_t0_correct=ipw_tot["delta_t0"],
      delta_t1_wrong=ipw_wrong["delta_1"], delta_t0_wrong=ipw_wrong["delta_0"],
      delta_t1_ols_ia=ols_ia["delta_1"], delta_t0_ols_ia=ols_ia["delta_0"],
      delta_t1_ols=ols_["delta_1"], delta_t0_ols=ols_["delta_0"],
      delta_t1_naive=ols_nav["delta_1"], delta_t0_naive=ols_nav["delta_0"],
      delta_p1_correct=ipw_part["delta_p1"], delta_p0_correct=ipw_part["delta_p0"],
      delta_p1_wrong=ipw_wrong["delta_1"], delta_p0_wrong=ipw_wrong["delta_0"],
      delta_p1_ols_ia=ols_ia["delta_1"], delta_p0_ols_ia=ols_ia["delta_0"],
      delta_p1_ols=ols_["delta_1"], delta_p0_ols=ols_["delta_0"],
      delta_p1_naive=ols_nav["delta_1"], delta_p0_naive=ols_nav["delta_0"]
    ))
  }
}

summarise_estimator <- function(estimates, true_val) {
  bias <- mean(estimates, na.rm=TRUE) - true_val
  vari <- var(estimates,  na.rm=TRUE)
  mse  <- bias^2 + vari
  round(c(bias=bias, var=vari, MSE=mse), 3)
}

run_framework1 <- function(beta, method="probit") {
  cat(sprintf("
--- Framework 1 | beta=%.1f | method=%s ---
", beta, method))
  results <- vector("list", N_SIM)
  set.seed(42)
  for (i in seq_len(N_SIM)) {
    results[[i]] <- tryCatch(run_one_sim(gamma=0, beta=beta, method=method), error=function(e) NULL)
  }
  results <- Filter(Negate(is.null), results)
  cat(sprintf("  Completed: %d / %d
", length(results), N_SIM))
  tv <- if (beta==0) true_values$g0_b0 else true_values$g0_b05
  extract <- function(name) sapply(results, function(r) r[[name]])
  list(
    theta_1=rbind(IPW=summarise_estimator(extract("theta_1_ipw"),tv$theta_1),
                  OLS.ia=summarise_estimator(extract("theta_1_ols_ia"),tv$theta_1),
                  OLS=summarise_estimator(extract("theta_1_ols"),tv$theta_1),
                  Naive=summarise_estimator(extract("theta_1_naive"),tv$theta_1)),
    theta_0=rbind(IPW=summarise_estimator(extract("theta_0_ipw"),tv$theta_0),
                  OLS.ia=summarise_estimator(extract("theta_0_ols_ia"),tv$theta_0),
                  OLS=summarise_estimator(extract("theta_0_ols"),tv$theta_0),
                  Naive=summarise_estimator(extract("theta_0_naive"),tv$theta_0)),
    delta_1=rbind(IPW=summarise_estimator(extract("delta_1_ipw"),tv$delta_1),
                  OLS.ia=summarise_estimator(extract("delta_1_ols_ia"),tv$delta_1),
                  OLS=summarise_estimator(extract("delta_1_ols"),tv$delta_1),
                  Naive=summarise_estimator(extract("delta_1_naive"),tv$delta_1)),
    delta_0=rbind(IPW=summarise_estimator(extract("delta_0_ipw"),tv$delta_0),
                  OLS.ia=summarise_estimator(extract("delta_0_ols_ia"),tv$delta_0),
                  OLS=summarise_estimator(extract("delta_0_ols"),tv$delta_0),
                  Naive=summarise_estimator(extract("delta_0_naive"),tv$delta_0))
  )
}

run_framework2 <- function(beta, method="probit") {
  cat(sprintf("
--- Framework 2 | beta=%.1f | method=%s ---
", beta, method))
  results <- vector("list", N_SIM)
  set.seed(42)
  for (i in seq_len(N_SIM)) {
    results[[i]] <- tryCatch(run_one_sim(gamma=0.2, beta=beta, method=method), error=function(e) NULL)
  }
  results <- Filter(Negate(is.null), results)
  cat(sprintf("  Completed: %d / %d
", length(results), N_SIM))
  tv <- if (beta==0) true_values$g02_b0 else true_values$g02_b05
  extract <- function(name) sapply(results, function(r) r[[name]])
  list(
    theta=rbind(IPW=summarise_estimator(extract("theta_1_ipw"),tv$theta_1),
                OLS.ia=summarise_estimator(extract("theta_1_ols_ia"),tv$theta_1),
                OLS=summarise_estimator(extract("theta_1_ols"),tv$theta_1),
                Naive=summarise_estimator(extract("theta_1_naive"),tv$theta_1)),
    delta_t=rbind(Correct_IPW=summarise_estimator(extract("delta_t1_correct"),tv$delta_t1),
                  Wrong_IPW=summarise_estimator(extract("delta_t1_wrong"),tv$delta_t1),
                  OLS.ia=summarise_estimator(extract("delta_t1_ols_ia"),tv$delta_t1),
                  OLS=summarise_estimator(extract("delta_t1_ols"),tv$delta_t1),
                  Naive=summarise_estimator(extract("delta_t1_naive"),tv$delta_t1)),
    delta_p=rbind(Correct_IPW=summarise_estimator(extract("delta_p1_correct"),tv$delta_p1),
                  Wrong_IPW=summarise_estimator(extract("delta_p1_wrong"),tv$delta_p1),
                  OLS.ia=summarise_estimator(extract("delta_p1_ols_ia"),tv$delta_p1),
                  OLS=summarise_estimator(extract("delta_p1_ols"),tv$delta_p1),
                  Naive=summarise_estimator(extract("delta_p1_naive"),tv$delta_p1))
  )
}

dir.create("output", showWarnings=FALSE)

cat("
===== FRAMEWORK 1 (Table 2) =====
")
table2_b0  <- run_framework1(beta=0,   method="probit")
table2_b05 <- run_framework1(beta=0.5, method="probit")

cat("
===== FRAMEWORK 2 (Table 3) =====
")
table3_b0  <- run_framework2(beta=0,   method="probit")
table3_b05 <- run_framework2(beta=0.5, method="probit")

save(table2_b0, table2_b05, table3_b0, table3_b05, file="output/baseline_results.RData")

cat("
--- TABLE 2 beta=0 (theta_1) ---
"); print(table2_b0$theta_1)
cat("
--- TABLE 2 beta=0 (delta_1) ---
"); print(table2_b0$delta_1)
cat("
--- TABLE 2 beta=0.5 (theta_1) ---
"); print(table2_b05$theta_1)
cat("
--- TABLE 2 beta=0.5 (delta_1) ---
"); print(table2_b05$delta_1)
cat("
--- TABLE 3 beta=0 (theta) ---
"); print(table3_b0$theta)
cat("
--- TABLE 3 beta=0 (delta_t) ---
"); print(table3_b0$delta_t)
cat("
--- TABLE 3 beta=0 (delta_p) ---
"); print(table3_b0$delta_p)
cat("
--- TABLE 3 beta=0.5 (theta) ---
"); print(table3_b05$theta)
cat("
--- TABLE 3 beta=0.5 (delta_t) ---
"); print(table3_b05$delta_t)
cat("
--- TABLE 3 beta=0.5 (delta_p) ---
"); print(table3_b05$delta_p)

write.csv(table2_b0$theta_1,  "output/table2_b0_theta1.csv")
write.csv(table2_b0$delta_1,  "output/table2_b0_delta1.csv")
write.csv(table2_b05$theta_1, "output/table2_b05_theta1.csv")
write.csv(table2_b05$delta_1, "output/table2_b05_delta1.csv")
write.csv(table3_b0$theta,    "output/table3_b0_theta.csv")
write.csv(table3_b0$delta_t,  "output/table3_b0_deltat.csv")
write.csv(table3_b0$delta_p,  "output/table3_b0_deltap.csv")
write.csv(table3_b05$theta,   "output/table3_b05_theta.csv")
write.csv(table3_b05$delta_t, "output/table3_b05_deltat.csv")
write.csv(table3_b05$delta_p, "output/table3_b05_deltap.csv")

cat("
✓ All tables saved to output/ folder
")

