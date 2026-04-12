
library(ggplot2)
library(dplyr)
library(tidyr)

# =============================================================================
# TABLE FORMATTING FUNCTION
# =============================================================================

format_sim_table <- function(results_b0, results_b05, effect = "theta_1") {
  tbl_b0  <- results_b0[[effect]]
  tbl_b05 <- results_b05[[effect]]
  cat(sprintf("
=== %s ===
", effect))
  cat(sprintf("%-14s | %6s %6s %6s | %6s %6s %6s
",
              "Estimator", "bias", "var", "MSE", "bias", "var", "MSE"))
  cat(sprintf("%-14s | %18s | %18s
", "", "  beta=0", "  beta=0.5"))
  cat(strrep("-", 64), "
")
  for (est in rownames(tbl_b0)) {
    cat(sprintf("%-14s | %6.3f %6.3f %6.3f | %6.3f %6.3f %6.3f
",
      est,
      tbl_b0[est,"bias"],  tbl_b0[est,"var"],  tbl_b0[est,"MSE"],
      tbl_b05[est,"bias"], tbl_b05[est,"var"], tbl_b05[est,"MSE"]
    ))
  }
}

# =============================================================================
# PLOT 1: Empirical ML comparison (from script 05)
# =============================================================================

plot_empirical_comparison <- function(ml_f, ml_m, save_path = NULL) {
  colnames(ml_f) <- c("method","ATE","theta_1","theta_0",
                       "delta_t1","delta_t0","delta_p1","delta_p0","delta_f1")
  colnames(ml_m) <- c("method","ATE","theta_1","theta_0",
                       "delta_t1","delta_t0","delta_p1","delta_p0","delta_f1")
  ml_f$gender <- "Female"
  ml_m$gender <- "Male"
  df <- rbind(ml_f, ml_m)
  df_long <- pivot_longer(df,
    cols = c("ATE","theta_1","theta_0","delta_t1","delta_p1"),
    names_to = "effect", values_to = "estimate"
  )
  df_long$method <- factor(df_long$method,
    levels = c("probit","lasso","rf","gbm"),
    labels = c("Probit","LASSO","Random Forest","GBM"))
  df_long$effect <- factor(df_long$effect,
    levels = c("ATE","theta_1","theta_0","delta_t1","delta_p1"),
    labels = c("ATE","Direct theta(1)","Direct theta(0)",
               "Total Indirect","Partial Indirect"))
  p <- ggplot(df_long, aes(x=effect, y=estimate, colour=method, shape=method)) +
    geom_hline(yintercept=0, linetype="dashed", colour="grey50") +
    geom_point(position=position_dodge(width=0.5), size=3) +
    facet_wrap(~ gender) +
    scale_colour_manual(values=c("#2C7BB6","#D7191C","#1A9641","#F46D43")) +
    labs(title="Empirical Estimates by Propensity Score Method",
         subtitle="Job Corps: direct and indirect health effects by gender",
         x="", y="Estimated Effect",
         colour="PS Method", shape="PS Method") +
    theme_bw(base_size=12) +
    theme(axis.text.x=element_text(angle=30, hjust=1),
          panel.grid.minor=element_blank())
  if (!is.null(save_path)) {
    ggsave(save_path, p, width=11, height=5, dpi=300)
    cat("Plot saved:", save_path, "
")
  }
  return(p)
}

# =============================================================================
# PLOT 2: Bias comparison across ML methods (from script 04)
# =============================================================================

plot_ml_bias <- function(ml_results_f1, save_path = NULL) {
  METHODS <- names(ml_results_f1)
  ml_bias <- data.frame(
    method = rep(METHODS, 2),
    effect = rep(c("Direct theta(1)","Indirect delta(1)"), each=length(METHODS)),
    bias   = c(
      sapply(METHODS, function(m) ml_results_f1[[m]]$theta_1["IPW","bias"]),
      sapply(METHODS, function(m) ml_results_f1[[m]]$delta_1["IPW","bias"])
    )
  )
  ml_bias$method <- factor(ml_bias$method,
    levels=c("probit","lasso","rf","gbm"),
    labels=c("Probit","LASSO","Random Forest","GBM"))
  p <- ggplot(ml_bias, aes(x=method, y=bias, fill=method)) +
    geom_col(width=0.6, alpha=0.85) +
    geom_hline(yintercept=0, linetype="dashed") +
    facet_wrap(~effect) +
    scale_fill_manual(values=c("#2C7BB6","#D7191C","#1A9641","#F46D43")) +
    labs(title="IPW Bias by Propensity Score Method",
         subtitle="Framework 1 (gamma=0, beta=0.5) — your ML extension",
         x="Propensity Score Method", y="Bias") +
    theme_bw(base_size=12) +
    theme(legend.position="none", panel.grid.minor=element_blank())
  if (!is.null(save_path)) {
    ggsave(save_path, p, width=10, height=5, dpi=300)
    cat("Plot saved:", save_path, "
")
  }
  return(p)
}

# =============================================================================
# PLOT 3: Wrong vs correct assumptions (key paper finding)
# =============================================================================

plot_wrong_vs_correct <- function(table3_b0, table3_b05, save_path = NULL) {
  assumption_df <- data.frame(
    Estimator = rownames(table3_b0$delta_t),
    bias_b0   = table3_b0$delta_t[,"bias"],
    bias_b05  = table3_b05$delta_t[,"bias"]
  )
  assumption_long <- pivot_longer(assumption_df,
    cols=c("bias_b0","bias_b05"),
    names_to="scenario", values_to="bias")
  assumption_long$scenario <- ifelse(
    assumption_long$scenario=="bias_b0","beta=0","beta=0.5")
  assumption_long$Estimator <- factor(assumption_long$Estimator,
    levels=rownames(table3_b0$delta_t))
  p <- ggplot(assumption_long, aes(x=Estimator, y=bias, fill=scenario)) +
    geom_col(position="dodge", width=0.6, alpha=0.85) +
    geom_hline(yintercept=0, linetype="dashed") +
    scale_fill_manual(values=c("#2C7BB6","#D7191C")) +
    labs(title="Bias of Total Indirect Effect: Correct vs Wrong Assumptions",
         subtitle="Framework 2 (gamma=0.2) — key finding of the paper",
         x="", y="Bias", fill="Scenario") +
    theme_bw(base_size=12) +
    theme(axis.text.x=element_text(angle=20, hjust=1),
          panel.grid.minor=element_blank())
  if (!is.null(save_path)) {
    ggsave(save_path, p, width=10, height=5, dpi=300)
    cat("Plot saved:", save_path, "
")
  }
  return(p)
}

# =============================================================================
# PLOT 4: MSE comparison across all estimators, both frameworks
# =============================================================================

plot_mse_comparison <- function(table2_b05, table3_b05, save_path = NULL) {
  df_f1 <- data.frame(
    Estimator = rownames(table2_b05$theta_1),
    MSE       = table2_b05$theta_1[,"MSE"],
    Framework = "Framework 1 (gamma=0)"
  )
  df_f2 <- data.frame(
    Estimator = rownames(table3_b05$theta),
    MSE       = table3_b05$theta[,"MSE"],
    Framework = "Framework 2 (gamma=0.2)"
  )
  df <- rbind(df_f1, df_f2)
  df$Estimator <- factor(df$Estimator, levels=unique(df$Estimator))
  p <- ggplot(df, aes(x=Estimator, y=MSE, fill=Framework)) +
    geom_col(position="dodge", width=0.6, alpha=0.85) +
    scale_fill_manual(values=c("#2C7BB6","#D7191C")) +
    labs(title="MSE of Direct Effect Estimators: Framework 1 vs Framework 2",
         subtitle="beta=0.5 scenario",
         x="", y="MSE", fill="") +
    theme_bw(base_size=12) +
    theme(panel.grid.minor=element_blank())
  if (!is.null(save_path)) {
    ggsave(save_path, p, width=10, height=5, dpi=300)
    cat("Plot saved:", save_path, "
")
  }
  return(p)
}

# =============================================================================
# GENERATE ALL OUTPUT
# =============================================================================

generate_all_output <- function() {
  dir.create("output/figures", showWarnings=FALSE, recursive=TRUE)

  # ── Tables 2 & 3 ─────────────────────────────────────────────────────────
  cat("
====== TABLE 2 REPLICATION (Framework 1, gamma=0) ======
")
  format_sim_table(table2_b0, table2_b05, "theta_1")
  format_sim_table(table2_b0, table2_b05, "theta_0")
  format_sim_table(table2_b0, table2_b05, "delta_1")
  format_sim_table(table2_b0, table2_b05, "delta_0")

  cat("
====== TABLE 3 REPLICATION (Framework 2, gamma=0.2) ======
")
  format_sim_table(table3_b0, table3_b05, "theta")
  format_sim_table(table3_b0, table3_b05, "delta_t")
  format_sim_table(table3_b0, table3_b05, "delta_p")

  # ── MSE comparison table (CSV) ────────────────────────────────────────────
  mse_table <- data.frame(
    Estimator    = rownames(table2_b05$theta_1),
    F1_theta_b0  = table2_b0$theta_1[,"MSE"],
    F1_theta_b05 = table2_b05$theta_1[,"MSE"],
    F2_theta_b0  = table3_b0$theta[1:4,"MSE"],
    F2_theta_b05 = table3_b05$theta[,"MSE"]
  )
  write.csv(mse_table, "output/mse_comparison.csv", row.names=FALSE)
  cat("
MSE comparison table saved
")

  # ── Plot: wrong vs correct assumptions ────────────────────────────────────
  plot_wrong_vs_correct(table3_b0, table3_b05,
    "output/figures/wrong_vs_correct.png")

  # ── Plot: MSE comparison ──────────────────────────────────────────────────
  plot_mse_comparison(table2_b05, table3_b05,
    "output/figures/mse_comparison.png")

  # ── Plot: ML bias comparison (needs 04 results) ───────────────────────────
  if (exists("ml_results_f1")) {
    plot_ml_bias(ml_results_f1,
      "output/figures/ml_bias_comparison.png")
  } else {
    cat("Skipping ML bias plot: run 04_ml_extension.R first
")
  }

  # ── Plot: empirical ML comparison (needs 05 results) ─────────────────────
  if (exists("ml_f") && exists("ml_m")) {
    plot_empirical_comparison(ml_f, ml_m,
      "output/figures/empirical_ml_comparison.png")
  } else {
    cat("Skipping empirical plot: run 05_empirical_application.R first
")
  }

  cat("
✓ All output generated in output/ and output/figures/
")
}

generate_all_output()

