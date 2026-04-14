
# =============================================================================
# FILE 05: EMPIRICAL APPLICATION (Job Corps Dataset)
# Analyzes real-world data to compare Probit vs LASSO
# =============================================================================

library(causalweight)
library(boot)


source("simulation/02_estimators.R")

data(JC)
cat("JC dataset loaded:", nrow(JC), "rows,", ncol(JC), "columns\n")

# ── 1. Variables & Data Prep ───────────────────────────────────────────────
D <- JC$assignment                       # Treatment: Job Corps assignment
M <- as.integer(JC$pworky2 > 0)          # Mediator: employed in year 2
Y <- as.integer(JC$health30 == 1)        # Outcome: very good health at 30 months
female <- JC$female

# Pre-treatment covariates X (Framework 1 & 2)
X_vars <- c("age","white","black","hispanic","educ","haschild",
            "everwkd","hhsize","health","smoke","alcohol","welfarechild")

# Post-treatment covariates W (Required for Framework 2)
W_vars <- c("everwkdy1","pworky1","health12","trainy1")

# Remove rows with NA in key variables
keep <- complete.cases(D, M, Y, JC[, X_vars], JC[, W_vars])
cat("Complete cases:", sum(keep), "out of", nrow(JC), "\n")

D <- D[keep]; M <- M[keep]; Y <- Y[keep]
X <- as.matrix(JC[keep, X_vars])
W <- as.matrix(JC[keep, W_vars])
female_keep <- female[keep]

# Split by gender
f_idx <- which(female_keep == 1)
m_idx <- which(female_keep == 0)
cat("Females:", length(f_idx), "| Males:", length(m_idx), "\n\n")

# ── 2. Core Estimation Functions ───────────────────────────────────────────
estimate_effects <- function(Y, D, M, X, W, method="probit") {
  p_X   <- estimate_pscore(D, X,              method=method)
  p_WX  <- estimate_pscore(D, cbind(W, X),    method=method)
  p_MWX <- estimate_pscore(D, cbind(M, W, X), method=method)
  
  ate      <- mean(Y[D==1]) - mean(Y[D==0])
  dir      <- ipw_direct_2(Y, D, p_X, p_MWX)
  part_ind <- ipw_partial_indirect(Y, D, p_X, p_WX, p_MWX)
  tot_ind  <- ipw_total_indirect(Y, D, M, W, X, p_X)
  
  c(ATE      = ate,
    theta_1  = unname(dir["theta_1"]),
    theta_0  = unname(dir["theta_0"]),
    delta_t1 = unname(tot_ind["delta_t1"]),
    delta_t0 = unname(tot_ind["delta_t0"]),
    delta_p1 = unname(part_ind["delta_p1"]),
    delta_p0 = unname(part_ind["delta_p0"]))
}

# Bootstrapping for standard errors
# NOTE: Set B=1999 for final run. If testing, change to B=199 to save time.
bootstrap_se <- function(Y, D, M, X, W, method="probit", B=1999) {
  n <- length(Y)
  boot_mat <- matrix(NA, nrow=B, ncol=7)
  for (b in 1:B) {
    idx <- sample(n, n, replace=TRUE)
    boot_mat[b,] <- tryCatch(
      estimate_effects(Y[idx], D[idx], M[idx], X[idx,], W[idx,], method),
      error=function(e) rep(NA, 7)
    )
  }
  apply(boot_mat, 2, sd, na.rm=TRUE)
}

# ── 3. Baseline Replication (Probit) ───────────────────────────────────────
cat("=================================================================\n")
cat("BASELINE REPLICATION (Tables 6: Probit)\n")
cat("=================================================================\n")

# Females
cat("\n--- FEMALES ---\n")
est_f <- estimate_effects(Y[f_idx], D[f_idx], M[f_idx], X[f_idx,], W[f_idx,])
cat("Bootstrapping SEs (B=1999)... This will take a few minutes.\n")
se_f  <- bootstrap_se(Y[f_idx], D[f_idx], M[f_idx], X[f_idx,], W[f_idx,])
pval_f <- 2*pnorm(-abs(est_f/se_f))

cat("Female Results:\n")
print(round(data.frame(Estimate=est_f, SE=se_f, pvalue=pval_f), 4))

# Males
cat("\n--- MALES ---\n")
est_m <- estimate_effects(Y[m_idx], D[m_idx], M[m_idx], X[m_idx,], W[m_idx,])
cat("Bootstrapping SEs (B=1999)... This will take a few minutes.\n")
se_m  <- bootstrap_se(Y[m_idx], D[m_idx], M[m_idx], X[m_idx,], W[m_idx,])
pval_m <- 2*pnorm(-abs(est_m/se_m))

cat("Male Results:\n")
print(round(data.frame(Estimate=est_m, SE=se_m, pvalue=pval_m), 4))

# ── 4. ML Extension (LASSO vs Probit) ──────────────────────────────────────
cat("\n=================================================================\n")
cat("ML METHOD COMPARISON (LASSO vs PROBIT)\n")
cat("=================================================================\n")

methods <- c("probit", "lasso")

ml_f <- do.call(rbind, lapply(methods, function(m) {
  cat("Running Female |", m, "\n")
  est <- tryCatch(
    estimate_effects(Y[f_idx], D[f_idx], M[f_idx], X[f_idx,], W[f_idx,], m),
    error=function(e) rep(NA,7)
  )
  data.frame(method=m, t(est))
}))

ml_m <- do.call(rbind, lapply(methods, function(m) {
  cat("Running Male   |", m, "\n")
  est <- tryCatch(
    estimate_effects(Y[m_idx], D[m_idx], M[m_idx], X[m_idx,], W[m_idx,], m),
    error=function(e) rep(NA,7)
  )
  data.frame(method=m, t(est))
}))

cat("\n--- Female ML Comparison ---\n")
print(ml_f)
cat("\n--- Male ML Comparison ---\n")
print(ml_m)

# ── 5. Save Final Results ──────────────────────────────────────────────────
dir.create("output", showWarnings=FALSE)
write.csv(data.frame(Estimate=est_f, SE=se_f, pvalue=pval_f), "output/table6_females.csv")
write.csv(data.frame(Estimate=est_m, SE=se_m, pvalue=pval_m), "output/table6_males.csv")
write.csv(ml_f, "output/ml_comparison_females.csv", row.names=FALSE)
write.csv(ml_m, "output/ml_comparison_males.csv", row.names=FALSE)

cat("\n✓ All empirical results successfully saved to output/ folder.\n")