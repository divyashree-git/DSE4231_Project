# =============================================================================
# FILE 1: Data Generating Process (DGP)
# Huber (2013), Section 4, Equations (17)-(20)
# =============================================================================
#
# The DGP generates fake data where we KNOW the true direct and indirect
# effects. This lets us test whether our estimators recover the right answer.
#
# The four equations are:
#   Y = 0.5*D + M + beta*D*M + V + 0.25*U + e1   (outcome)
#   M = 0.5*D + 0.5*V + 0.25*U + e2              (mediator)
#   V = gamma*D + 0.25*U + e3                     (post-treatment confounder)
#   D = I(0.25*U + e4 > 0)                        (treatment, prob = 0.5)
#
# Key parameters:
#   gamma = 0   => V is NOT affected by D => use Framework 1 (Assumptions 1&2)
#   gamma = 0.2 => V IS affected by D    => use Framework 2 (Assumptions 3-5)
#   beta  = 0   => no treatment-mediator interaction
#   beta  = 0.5 => treatment-mediator interaction exists
# =============================================================================

generate_data <- function(n, gamma, beta, seed = NULL) {
  # n     : sample size
  # gamma : controls whether V is affected by D (0 = no, 0.2 = yes)
  # beta  : controls D-M interaction in outcome equation
  # seed  : for reproducibility

  if (!is.null(seed)) set.seed(seed)

  # Draw all error terms independently from N(0,1)
  U  <- rnorm(n)   # observed confounder of D
  e1 <- rnorm(n)   # outcome error
  e2 <- rnorm(n)   # mediator error
  e3 <- rnorm(n)   # V error
  e4 <- rnorm(n)   # treatment error

  # Treatment: D = 1 if 0.25*U + e4 > 0
  # Because U and e4 are both N(0,1), Pr(D=1) = 0.5
  D <- as.integer(0.25 * U + e4 > 0)

  # Post-treatment variable V (confounder of M and Y)
  # If gamma=0: V does not depend on D => pre-treatment confounder
  # If gamma>0: V depends on D => post-treatment confounder (W in the paper)
  V <- gamma * D + 0.25 * U + e3

  # Mediator M (e.g. employment status)
  M <- 0.5 * D + 0.5 * V + 0.25 * U + e2

  # Outcome Y (e.g. health status)
  Y <- 0.5 * D + M + beta * D * M + V + 0.25 * U + e1

  data.frame(Y = Y, D = D, M = M, V = V, U = U)
}


# =============================================================================
# TRUE PARAMETER VALUES (derived analytically in Appendix A.3 of the paper)
# =============================================================================
# These are the "right answers" our estimators should recover.
#
# Scenario 1: gamma=0, beta=0
#   theta(1) = theta(0) = 0.5
#   delta(1) = delta(0) = 0.5
#
# Scenario 2: gamma=0, beta=0.5
#   theta(1) = 0.75, theta(0) = 0.5
#   delta(1) = 0.75, delta(0) = 0.5
#
# Scenario 3: gamma=0.2, beta=0
#   theta(1) = theta(0) = 0.5
#   delta_t(1) = delta_t(0) = 0.6
#   delta_p(1) = delta_p(0) = 0.5
#
# Scenario 4: gamma=0.2, beta=0.5
#   theta(1) = 0.8, theta(0) = 0.5
#   delta_t(1) = 0.9, delta_t(0) = 0.6
#   delta_p(1) = 0.75, delta_p(0) = 0.5

true_values <- list(
  # gamma=0, beta=0
  g0_b0 = list(
    theta_1 = 0.5,  theta_0 = 0.5,
    delta_1 = 0.5,  delta_0 = 0.5
  ),
  # gamma=0, beta=0.5
  g0_b05 = list(
    theta_1 = 0.75, theta_0 = 0.5,
    delta_1 = 0.75, delta_0 = 0.5
  ),
  # gamma=0.2, beta=0
  g02_b0 = list(
    theta_1  = 0.5,  theta_0  = 0.5,
    delta_t1 = 0.6,  delta_t0 = 0.6,
    delta_p1 = 0.5,  delta_p0 = 0.5
  ),
  # gamma=0.2, beta=0.5
  g02_b05 = list(
    theta_1  = 0.8,  theta_0  = 0.5,
    delta_t1 = 0.9,  delta_t0 = 0.6,
    delta_p1 = 0.75, delta_p0 = 0.5
  )
)
