# =============================================================================
# FILE 4: ML Extension
# Replaces probit propensity scores with ML alternatives
# This is your main EXTENSION beyond the paper's original results
# =============================================================================
#
# MOTIVATION (from the paper itself, Section 2.2):
#   "IPW is less robust to propensity score misspecification than other
#    classes of estimators" (Kang & Schafer, 2007; Waernbaum, 2012)
#
# RESEARCH QUESTION for the extension:
#   Does replacing the parametric probit model with flexible ML estimators
#   improve or destabilize the IPW estimator's performance?
#   Under what conditions does ML help?
#
# METHODS COMPARED:
#   1. Probit (original paper)
#   2. LASSO logistic regression (regularized, feature selection)
#   3. Random Forest (nonparametric, captures interactions)
#   4. Gradient Boosted Trees (flexible, often best in practice)
# =============================================================================

source("simulation/01_dgp.R")
source("simulation/02_estimators.R")
source("simulation/03_run_simulations.R")

library(parallel)

# Methods to compare
METHODS <- c("probit", "lasso", "rf", "gbm")

cat("\n=================================================================\n")
cat("ML EXTENSION: Comparing propensity score estimators\n")
cat("=================================================================\n")


# =============================================================================
# RUN ALL METHODS FOR FRAMEWORK 1
# =============================================================================

library(parallel)

# 1. Setup Parallelism (detect how many "brains" your Mac has)
num_cores <- detectCores() - 1 

# 2. Define Adaptive Iterations 
# We do 5000 for fast ones, 500 for slow ones to save time
get_n_sim <- function(m) {
  if (m %in% c("rf", "gbm")) return(500) else return(5000)
}

cat("\n--- Framework 1 (gamma=0, beta=0.5): Parallel ML Extension ---\n")

ml_results_f1 <- lapply(METHODS, function(m) {
  n_val <- get_n_sim(m)
  cat(sprintf("\nStarting %s: Running %d iterations on %d cores...\n", m, n_val, num_cores))
  
  # 3. Parallelize the internal loop using mclapply
  # This replaces the 'for' loop inside run_framework1
  set.seed(42)
  results <- mclapply(1:n_val, function(i) {
    tryCatch(run_one_sim(gamma = 0, beta = 0.5, method = m), 
             error = function(e) NULL)
  }, mc.cores = num_cores)
  
  # 4. Clean up results
  results <- Filter(Negate(is.null), results)
  
  # 5. Summarize (using your existing extract logic)
  # Ensure your summarise_estimator function handles the correct true values
  tv <- true_values$g0_b05 
  extract <- function(name) unlist(lapply(results, function(r) r[[name]]))
  
  list(
    theta_1 = rbind(IPW = summarise_estimator(extract("theta_1_ipw"), tv$theta_1),
                    OLS = summarise_estimator(extract("theta_1_ols"), tv$theta_1)),
    delta_1 = rbind(IPW = summarise_estimator(extract("delta_1_ipw"), tv$delta_1),
                    OLS = summarise_estimator(extract("delta_1_ols"), tv$delta_1))
  )
})

names(ml_results_f1) <- METHODS

# =============================================================================
# RUN ALL METHODS FOR FRAMEWORK 2
# =============================================================================

cat("\n--- Framework 2 (gamma=0.2, beta=0.5): All Methods ---\n")

ml_results_f2 <- lapply(METHODS, function(m) {
  cat(sprintf("\nRunning: %s\n", m))
  run_framework2(beta = 0.5, method = m)
})
names(ml_results_f2) <- METHODS


# # =============================================================================
# # SAMPLE SIZE ANALYSIS
# # How does performance change with smaller/larger samples?
# # =============================================================================
# 
cat("\n--- Sample Size Analysis ---\n")

sample_sizes <- c(500, 1000, 2000, 5000)

run_sample_size_analysis <- function(n_obs, method = "probit") {
  # Temporarily override N_OBS
  old_N <- N_OBS
  N_OBS <<- n_obs

  # Run a reduced simulation (500 reps to save time)
  results <- vector("list", 500)
  set.seed(42)
  for (i in 1:500) {
    results[[i]] <- tryCatch(
      run_one_sim(gamma = 0, beta = 0.5, method = method),
      error = function(e) NULL
    )
  }
  results    <- Filter(Negate(is.null), results)
  N_OBS     <<- old_N

  tv <- true_values$g0_b05
  extract <- function(name) sapply(results, function(r) r[[name]])

  data.frame(
    n      = n_obs,
    method = method,
    effect = "theta(1)",
    bias   = mean(extract("theta_1_ipw"), na.rm=TRUE) - tv$theta_1,
    MSE    = mean((extract("theta_1_ipw") - tv$theta_1)^2, na.rm=TRUE)
  )
}

size_results <- do.call(rbind, lapply(sample_sizes, function(n) {
  do.call(rbind, lapply(METHODS, function(m) {
    run_sample_size_analysis(n, m)
  }))
}))


# # =============================================================================
# # MISSPECIFICATION ANALYSIS
# # What happens when the outcome model is nonlinear?
# # This extends Table 4 in the paper
# # =============================================================================
# 
cat("\n--- Misspecification Analysis (nonlinear DGP) ---\n")

generate_data_nonlinear <- function(n, seed = NULL) {
  # Modified DGP: adds nonlinear terms to outcome equation
  # IPW should still work (doesn't model outcome)
  # OLS should break (outcome model is wrong)
  if (!is.null(seed)) set.seed(seed)

  U  <- rnorm(n); e1 <- rnorm(n); e2 <- rnorm(n)
  e3 <- rnorm(n); e4 <- rnorm(n)

  D <- as.integer(0.25 * U + e4 > 0)
  V <- 0.25 * U + e3
  M <- 0.5 * D + 0.5 * V + 0.25 * U + e2

  # Nonlinear outcome: add D*U and D*V interactions
  Y <- 0.5*D + M + 0.5*D*M + V + 0.25*U +
       0.3*D*U +          # extra nonlinearity (misspecifies OLS)
       0.2*D*V + e1       # extra nonlinearity

  data.frame(Y=Y, D=D, M=M, V=V, U=U)
}

run_misspec <- function(method = "probit", n_sims = 500) {
  results <- vector("list", n_sims)
  set.seed(42)
  for (i in seq_len(n_sims)) {
    dat <- generate_data_nonlinear(N_OBS)
    results[[i]] <- tryCatch({
      Y <- dat$Y; D <- dat$D; M <- dat$M
      V <- dat$V; U <- dat$U

      p_X  <- estimate_pscore(D, cbind(V, U),    method = method)
      p_MX <- estimate_pscore(D, cbind(M, V, U), method = method)

      ipw_d <- ipw_direct_1(Y, D, p_X, p_MX)
      ipw_i <- ipw_indirect_1(Y, D, p_X, p_MX)
      ols_  <- ols_effects(Y, D, M, V, U, include_interaction = TRUE)

      list(
        theta_1_ipw   = ipw_d["theta_1"],
        theta_1_ols   = ols_["theta_1"],
        delta_1_ipw   = ipw_i["delta_1"],
        delta_1_ols   = ols_["delta_1"]
      )
    }, error = function(e) NULL)
  }
  results <- Filter(Negate(is.null), results)

  # True value under this nonlinear DGP is approximately 0.75 for theta(1)
  # (the direct interaction terms add to the bias of OLS but not IPW)
  extract <- function(name) sapply(results, function(r) r[[name]])

  data.frame(
    method     = method,
    theta1_ipw_bias = mean(extract("theta_1_ipw"), na.rm=TRUE) - 0.75,
    theta1_ols_bias = mean(extract("theta_1_ols"), na.rm=TRUE) - 0.75,
    delta1_ipw_bias = mean(extract("delta_1_ipw"), na.rm=TRUE) - 0.75,
    delta1_ols_bias = mean(extract("delta_1_ols"), na.rm=TRUE) - 0.75
  )
}

misspec_results <- do.call(rbind, lapply(METHODS, run_misspec))


# =============================================================================
# SAVE ALL EXTENSION RESULTS
# =============================================================================

save(ml_results_f1, ml_results_f2, size_results, misspec_results,
     file = "output/ml_extension_results.RData")

cat("\n✓ ML extension complete. Results saved to output/ml_extension_results.RData\n")
