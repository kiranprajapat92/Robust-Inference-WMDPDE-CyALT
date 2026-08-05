rm(list = ls())
library(optimx)
library(MASS)

source("C:/Users/Kiran/Downloads/WMDPDE_CyALT_lognormal/Main_codes/Testing/cyalt_lognormal_WMDPDE1.R")

##### reconstruct theta_true and data_obs exactly as in real_data_analysis1.R

failure_times <- c(26763, 31959, 32887, 33069, 34019, 34924, 36754,
                   37054, 37385, 38045, 41033, 41755, 42333,
                   42818, 44638, 44867, 48364, 49767)

fit   <- fitdistr(failure_times, "log-normal")
sigma <- round(fit$estimate["sdlog"], 3)
p_h   <- 0.90
p_u   <- 0.001
t_c   <- 50000
tau   <- 0.5

v_0F_phys <- 0.10; v_0C_phys <- 0.25; v_hC_phys <- 3.00
s_im  <- function(v) (log(v) - log(v_0F_phys)) / (log(v_hC_phys) - log(v_0F_phys))

s_0F  <- 0.00;  s_0C <- 0.27
s_hF  <- 0.50;  s_hC <- 1.00
s_1C  <- 0.70;  s_2C <- 1.00;  s_F <- 0.30

B_fun       <- function(a1, sC, sF) tau * exp(-a1 * sC) + (1 - tau) * exp(-a1 * sF)
mu_h        <- log(t_c) - sigma * qnorm(p_h)
mu_u        <- log(t_c) - sigma * qnorm(p_u)
eq_a1       <- function(a1) (mu_h - mu_u) - (log(B_fun(a1, s_0C, s_0F)) - log(B_fun(a1, s_hC, s_hF)))
alpha1_true <- uniroot(eq_a1, c(-30, -1e-4))$root
alpha0_true <- mu_h + log(B_fun(alpha1_true, s_hC, s_hF))
theta_true  <- c(alpha0_true, alpha1_true, sigma)

ITs        <- c(25000, 35000, 45000, 50000, 60000, 65000)
Kvec       <- c(140, 60)
stress_mat <- matrix(c(s_F, s_1C, s_F, s_2C), nrow = 2, byrow = TRUE)

data_obs <- sim_cyalt_data_continuous(theta_true, Kvec, stress_mat, tau, ITs, seed = 125)

theta0 <- theta_true   # null hypothesis: does the CyALT data support the pilot-based values?

cat(sprintf("\nNull hypothesis theta0 (from pilot experiment): alpha0=%.4f, alpha1=%.4f, sigma=%.4f\n",
            theta0[1], theta0[2], theta0[3]))


# =============================================================================
# WALD-TYPE TEST STATISTICS  (same definitions as the simulation study)
# =============================================================================

wald_H1 <- function(theta_hat, theta0, Kvec, stress_mat, tau, ITs, beta) {
  Sig  <- Sigma_hat(theta0, Kvec, stress_mat, tau, ITs, beta)
  Sinv <- tryCatch(solve(Sig), error = function(e) matrix(NA, 3, 3))
  diff <- theta_hat - theta0
  drop(t(diff) %*% Sinv %*% diff)
}

wald_H2a <- function(theta_hat, theta0, Kvec, stress_mat, tau, ITs, beta) {
  Sig <- Sigma_hat(theta0, Kvec, stress_mat, tau, ITs, beta)
  (theta_hat[1] - theta0[1])^2 / Sig[1,1]
}

wald_H2b <- function(theta_hat, theta0, Kvec, stress_mat, tau, ITs, beta) {
  Sig <- Sigma_hat(theta0, Kvec, stress_mat, tau, ITs, beta)
  (theta_hat[2] - theta0[2])^2 / Sig[2,2]
}

wald_H3 <- function(theta_hat, theta0, Kvec, stress_mat, tau, ITs, beta) {
  Sig   <- Sigma_hat(theta0, Kvec, stress_mat, tau, ITs, beta)
  Sig22 <- Sig[1:2, 1:2]
  Sinv2 <- tryCatch(solve(Sig22), error = function(e) matrix(NA, 2, 2))
  diff  <- theta_hat[1:2] - theta0[1:2]
  drop(t(diff) %*% Sinv2 %*% diff)
}


# =============================================================================
# RAO-TYPE TEST STATISTICS  (same definitions as the simulation study)
# =============================================================================

precompute_rao_null <- function(theta0, Kvec, stress_mat, tau, ITs, beta) {
  R <- nrow(stress_mat)
  p_list <- vector("list", R)
  W_list <- vector("list", R)
  D_list <- vector("list", R)
  for (i in seq_len(R)) {
    p_list[[i]] <- p_i_theta(theta0[1], theta0[2], theta0[3],
                             sC = stress_mat[i,2], sF = stress_mat[i,1], tau, ITs)
    W_list[[i]] <- W_i_theta_matrix(theta0[1], theta0[2], theta0[3],
                                    sC = stress_mat[i,2], sF = stress_mat[i,1], tau, ITs)
    D_list[[i]] <- diag(p_list[[i]]^(beta - 1))
  }
  Kmat <- K_beta_mat(theta0, Kvec, stress_mat, tau, ITs, beta)
  Kinv <- tryCatch(solve(Kmat), error = function(e) matrix(NA, 3, 3))
  list(p = p_list, W = W_list, D = D_list, Kmat = Kmat, Kinv = Kinv)
}

compute_U_beta <- function(counts_list, Kvec, precomp) {
  R <- length(counts_list)
  K <- sum(Kvec)
  U <- numeric(3)
  for (i in seq_len(R)) {
    phat_i <- counts_list[[i]] / Kvec[i]
    U <- U + (Kvec[i]/K) * drop(t(precomp$W[[i]]) %*% precomp$D[[i]] %*%
                                  (precomp$p[[i]] - phat_i))
  }
  U
}

rao_H1  <- function(U, precomp, K) K * drop(t(U) %*% precomp$Kinv %*% U)
rao_H2a <- function(U, precomp, K) K * U[1]^2 / precomp$Kmat[1,1]
rao_H2b <- function(U, precomp, K) K * U[2]^2 / precomp$Kmat[2,2]
rao_H3  <- function(U, precomp, K) {
  Ksub  <- precomp$Kmat[1:2, 1:2]
  Kinv2 <- solve(Ksub)
  U12   <- U[1:2]
  K * drop(t(U12) %*% Kinv2 %*% U12)
}


# =============================================================================
# APPLY BOTH TESTS TO THE REAL DATA, FOR EACH BETA
# =============================================================================

beta_vec    <- c(0, 0.2, 0.4, 0.6, 0.8, 1.0)
beta_labels <- c("MLE", "0.2", "0.4", "0.6", "0.8", "1.0")

cat("\nFitting WMDPDE for all beta (needed for the Wald-type tests)...\n")
estimates_all <- fit_all_betas(data_obs, Kvec, stress_mat, tau, ITs,
                               beta_vec = beta_vec, init = c(10, -1, 0.2))

K <- sum(Kvec)

results_wald <- data.frame()
results_rao  <- data.frame()

for (bi in seq_along(beta_vec)) {
  
  b   <- beta_vec[bi]
  nm  <- if (b == 0) "MLE" else paste0("MDPDE_", b)
  est <- estimates_all[[nm]]
  
  # ---- Wald-type ----
  # ---- Wald-type ---- (corrected — no external * K)
  W1  <- wald_H1(est,  theta0, Kvec, stress_mat, tau, ITs, b)
  W2a <- wald_H2a(est, theta0, Kvec, stress_mat, tau, ITs, b)
  W2b <- wald_H2b(est, theta0, Kvec, stress_mat, tau, ITs, b)
  W3  <- wald_H3(est,  theta0, Kvec, stress_mat, tau, ITs, b)
  
  results_wald <- rbind(results_wald, data.frame(
    beta = beta_labels[bi],
    H1  = W1,  p_H1  = 1 - pchisq(W1,  df = 3),
    H2a = W2a, p_H2a = 1 - pchisq(W2a, df = 1),
    H2b = W2b, p_H2b = 1 - pchisq(W2b, df = 1),
    H3  = W3,  p_H3  = 1 - pchisq(W3,  df = 2)
  ))
  
  # ---- Rao-type ----
  precomp <- precompute_rao_null(theta0, Kvec, stress_mat, tau, ITs, b)
  U       <- compute_U_beta(data_obs, Kvec, precomp)
  
  R1  <- rao_H1(U,  precomp, K)
  R2a <- rao_H2a(U, precomp, K)
  R2b <- rao_H2b(U, precomp, K)
  R3  <- rao_H3(U,  precomp, K)
  
  results_rao <- rbind(results_rao, data.frame(
    beta = beta_labels[bi],
    H1  = R1,  p_H1  = 1 - pchisq(R1,  df = 3),
    H2a = R2a, p_H2a = 1 - pchisq(R2a, df = 1),
    H2b = R2b, p_H2b = 1 - pchisq(R2b, df = 1),
    H3  = R3,  p_H3  = 1 - pchisq(R3,  df = 2)
  ))
}


# =============================================================================
# PRINT TABLES
# =============================================================================

cat("\n=== Wald-type test statistics and p-values ===\n")
print(results_wald, digits = 4, row.names = FALSE)

cat("\n=== Rao-type test statistics and p-values ===\n")
print(results_rao, digits = 4, row.names = FALSE)

save(results_wald, results_rao, theta0, theta_true, beta_vec, beta_labels,
     file = "C:/Users/Kiran/Downloads/WMDPDE_CyALT_lognormal/real_data_tests.RData")
cat("\nSaved: real_data_tests.RData\n")