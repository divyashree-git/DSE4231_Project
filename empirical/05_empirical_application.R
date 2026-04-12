
library(causalweight)
library(boot)

source("simulation/02_estimators.R")

data(JC)
cat("JC dataset loaded:", nrow(JC), "rows,", ncol(JC), "columns
")

# ── Key variables ──────────────────────────────────────────────────────────
D <- JC$assignment                        # Treatment: Job Corps assignment
M <- as.integer(JC$pworky2 > 0)          # Mediator: employed in year 2
Y <- as.integer(JC$health30 == 1)        # Outcome: very good health at 30 months
female <- JC$female

# ── Pre-treatment covariates X ─────────────────────────────────────────────
X_vars <- c("age","white","black","hispanic","educ","haschild",
            "everwkd","hhsize","health","smoke","alcohol","welfarechild")

# ── Post-treatment covariates W (measured after treatment, before mediator) ─
W_vars <- c("everwkdy1","pworky1","health12","trainy1")

cat("X vars:", paste(X_vars, collapse=", "), "
")
cat("W vars:", paste(W_vars, collapse=", "), "
")

# ── Remove rows with NA in key variables ───────────────────────────────────
keep <- complete.cases(D, M, Y, JC[, X_vars], JC[, W_vars])
cat("Complete cases:", sum(keep), "out of", nrow(JC), "
")

D <- D[keep]; M <- M[keep]; Y <- Y[keep]
X <- as.matrix(JC[keep, X_vars])
W <- as.matrix(JC[keep, W_vars])
female_keep <- female[keep]

# ── Split by gender ────────────────────────────────────────────────────────
f_idx <- which(female_keep == 1)
m_idx <- which(female_keep == 0)
cat("Females:", length(f_idx), "| Males:", length(m_idx), "
")

# ── Estimation function ────────────────────────────────────────────────────
estimate_effects <- function(Y, D, M, X, W, method="probit") {
  p_X   <- estimate_pscore(D, X,              method=method)
  p_WX  <- estimate_pscore(D, cbind(W, X),    method=method)
  p_MWX <- estimate_pscore(D, cbind(M, W, X), method=method)

  ate      <- mean(Y[D==1]) - mean(Y[D==0])
  dir      <- ipw_direct_2(Y, D, p_X, p_MWX)
  part_ind <- ipw_partial_indirect(Y, D, p_X, p_WX, p_MWX)
  tot_ind  <- ipw_total_indirect(Y, D, M, W, X, p_X)

  p_MX_f1 <- estimate_pscore(D, cbind(M, X), method=method)
  ind_f1  <- ipw_indirect_1(Y, D, p_X, p_MX_f1)

  c(ATE      = ate,
    theta_1  = dir["theta_1"],
    theta_0  = dir["theta_0"],
    delta_t1 = tot_ind["delta_t1"],
    delta_t0 = tot_ind["delta_t0"],
    delta_p1 = part_ind["delta_p1"],
    delta_p0 = part_ind["delta_p0"],
    delta_f1 = ind_f1["delta_1"])
}

# ── Bootstrap SE function ──────────────────────────────────────────────────
bootstrap_se <- function(Y, D, M, X, W, method="probit", B=199) {
  n <- length(Y)
  boot_mat <- matrix(NA, nrow=B, ncol=8)
  for (b in 1:B) {
    idx <- sample(n, n, replace=TRUE)
    boot_mat[b,] <- tryCatch(
      estimate_effects(Y[idx], D[idx], M[idx], X[idx,], W[idx,], method),
      error=function(e) rep(NA, 8)
    )
  }
  apply(boot_mat, 2, sd, na.rm=TRUE)
}

# ── Run for females ────────────────────────────────────────────────────────
cat("
===== FEMALES =====
")
est_f <- estimate_effects(Y[f_idx], D[f_idx], M[f_idx], X[f_idx,], W[f_idx,])
cat("Point estimates:
"); print(round(est_f, 4))
cat("Bootstrapping SEs (B=199)...
")
se_f  <- bootstrap_se(Y[f_idx], D[f_idx], M[f_idx], X[f_idx,], W[f_idx,])
pval_f <- 2*pnorm(-abs(est_f/se_f))

cat("
Female Results:
")
print(round(data.frame(Estimate=est_f, SE=se_f, pvalue=pval_f), 4))

# ── Run for males ──────────────────────────────────────────────────────────
cat("
===== MALES =====
")
est_m <- estimate_effects(Y[m_idx], D[m_idx], M[m_idx], X[m_idx,], W[m_idx,])
cat("Point estimates:
"); print(round(est_m, 4))
cat("Bootstrapping SEs (B=199)...
")
se_m  <- bootstrap_se(Y[m_idx], D[m_idx], M[m_idx], X[m_idx,], W[m_idx,])
pval_m <- 2*pnorm(-abs(est_m/se_m))

cat("
Male Results:
")
print(round(data.frame(Estimate=est_m, SE=se_m, pvalue=pval_m), 4))

# ── ML comparison ─────────────────────────────────────────────────────────
cat("
===== ML METHOD COMPARISON =====
")
methods <- c("probit","lasso","rf","gbm")

ml_f <- do.call(rbind, lapply(methods, function(m) {
  cat("Female |", m, "
")
  est <- tryCatch(
    estimate_effects(Y[f_idx], D[f_idx], M[f_idx], X[f_idx,], W[f_idx,], m),
    error=function(e) rep(NA,8)
  )
  data.frame(method=m, t(est))
}))

ml_m <- do.call(rbind, lapply(methods, function(m) {
  cat("Male |", m, "
")
  est <- tryCatch(
    estimate_effects(Y[m_idx], D[m_idx], M[m_idx], X[m_idx,], W[m_idx,], m),
    error=function(e) rep(NA,8)
  )
  data.frame(method=m, t(est))
}))

cat("
Female ML Comparison:
"); print(ml_f)
cat("
Male ML Comparison:
"); print(ml_m)

# ── Save results ───────────────────────────────────────────────────────────
dir.create("output", showWarnings=FALSE)
write.csv(data.frame(Estimate=est_f, SE=se_f, pvalue=pval_f), "output/table6_females.csv")
write.csv(data.frame(Estimate=est_m, SE=se_m, pvalue=pval_m), "output/table6_males.csv")
write.csv(ml_f, "output/ml_comparison_females.csv")
write.csv(ml_m, "output/ml_comparison_males.csv")

cat("
✓ Empirical results saved to output/ folder
")

