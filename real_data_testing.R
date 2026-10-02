rm(list = ls())
library(optimx)
library(MASS)

# source("/mnt/nfs/home/nkp117/largefiles/KP_Leandro_Maria_1/Testing/arxiv/cyalt_lognormal_WMDPDE1.R")
source("C:/Users/Kiran/WMDPDE_CyALT_lognormal/Testing/cyalt_lognormal_WMDPDE1.R")

H_beta_objective_restricted_alpha1 <- function(par2, counts_list, Kvec, stress_mat,
                                               tau, ITs, beta, alpha1_fixed) {
    alpha0 <- par2[1]; sigma <- par2[2]
    if (sigma <= 0) return(1e10)
    if (alpha1_fixed > 0) return(1e10)   # relaxed from >= 0 to allow alpha1 = 0
    R <- length(counts_list)
    K <- sum(Kvec)
    obj <- 0
    for (i in seq_len(R)) {
        pvec   <- p_i_theta(alpha0, alpha1_fixed, sigma,
                            sC = stress_mat[i, 2], sF = stress_mat[i, 1], tau, ITs)
        phat_i <- counts_list[[i]] / Kvec[i]
        if (beta == 0) {
            obj <- obj - (Kvec[i] / K) * sum(phat_i * log(pvec))
        } else {
            obj <- obj + (Kvec[i] / K) * (sum(pvec^(1 + beta)) -
                                              (1 + 1 / beta) * sum(phat_i * pvec^beta))
        }
    }
    return(obj)
}

fit_mdpde_restricted_alpha1 <- function(counts_list, Kvec, stress_mat, tau, ITs,
                                        beta, alpha1_fixed, init = NULL) {
    if (is.null(init)) init <- c(10, 0.2)
    result <- tryCatch(
        optimx(par = init, fn = H_beta_objective_restricted_alpha1,
               counts_list = counts_list, Kvec = Kvec, stress_mat = stress_mat,
               tau = tau, ITs = ITs, beta = beta, alpha1_fixed = alpha1_fixed,
               method = "Nelder-Mead", control = list(maxit = 5000, reltol = 1e-10)),
        error = function(e) NULL)
    if (is.null(result) || result$convcode[1] != 0) return(rep(NA, 3))
    par2_out <- as.numeric(result[1, 1:2])
    if (any(is.na(par2_out))) return(rep(NA, 3))
    return(c(par2_out[1], alpha1_fixed, par2_out[2]))
}

fit_all_betas_restricted_alpha1 <- function(counts_list, Kvec, stress_mat, tau, ITs,
                                            alpha1_fixed,
                                            beta_vec = c(0, 0.2, 0.4, 0.6, 0.8, 1),
                                            init = NULL) {
    estimates    <- list()
    current_init <- init
    for (b in beta_vec) {
        label <- if (b == 0) "MLE" else paste0("MDPDE_", b)
        est   <- fit_mdpde_restricted_alpha1(counts_list, Kvec, stress_mat, tau, ITs,
                                             beta = b, alpha1_fixed = alpha1_fixed,
                                             init = current_init)
        estimates[[label]] <- est
        if (!any(is.na(est)) && !any(abs(est) > 1e6)) current_init <- c(est[1], est[3])
    }
    return(estimates)
}


############### reconstruct theta_true and data_obs ##############

      
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

theta0      <- theta_true   # H(A): does the CyALT data support the pilot-based values?
sigma0      <- theta0[3]
alpha1_null <- 0             # H(B): stress effect exactly zero

cat(sprintf("\nNull hypothesis theta0 (from pilot experiment): alpha0=%.4f, alpha1=%.4f, sigma=%.4f\n",
            theta0[1], theta0[2], theta0[3]))


###############  H(A): SIMPLE NULL  theta = theta0 ############### 

      
wald_HA <- function(theta_hat, theta0, Kvec, stress_mat, tau, ITs, beta) {
    Sig  <- Sigma_hat(theta0, Kvec, stress_mat, tau, ITs, beta)
    Sinv <- tryCatch(solve(Sig), error = function(e) matrix(NA, 3, 3))
    diff <- theta_hat - theta0
    drop(t(diff) %*% Sinv %*% diff)
}

rao_HA <- function(U, Kmat, K) {
    Kinv <- tryCatch(solve(Kmat), error = function(e) matrix(NA, 3, 3))
    if (any(is.na(Kinv))) return(NA)
    K * drop(t(U) %*% Kinv %*% U)
}


                     
###############  H(B) and H(C): COMPOSITE, SCALAR CONSTRAINTS ############### 
# H(B): alpha1 = 0    -> j = 2   (no stress effect)
# H(C): sigma  = sigma0  -> j = 3   (shape-parameter consistency)


                     
wald_composite_scalar <- function(theta_hat_unrestricted, j, null_val,
                                  Kvec, stress_mat, tau, ITs, beta) {
    Sig <- Sigma_hat(theta_hat_unrestricted, Kvec, stress_mat, tau, ITs, beta)
    if (any(is.na(Sig))) return(NA)
    sjj <- Sig[j, j]
    if (is.na(sjj) || sjj <= 0) return(NA)
    (theta_hat_unrestricted[j] - null_val)^2 / sjj
}

rao_composite_scalar <- function(theta_tilde_restricted, j,
                                 counts_list, Kvec, stress_mat, tau, ITs, beta) {
    K    <- sum(Kvec)
    Jmat <- J_beta_mat(theta_tilde_restricted, Kvec, stress_mat, tau, ITs, beta)
    Jinv <- tryCatch(solve(Jmat), error = function(e) matrix(NA, 3, 3))
    if (any(is.na(Jinv))) return(NA)
    
    Hvec <- rep(0, 3); Hvec[j] <- 1
    Q    <- Jinv %*% Hvec / drop(t(Hvec) %*% Jinv %*% Hvec)   # 3x1, r=1 case
    
    U    <- U_beta_vec(theta_tilde_restricted, counts_list, Kvec, stress_mat, tau, ITs, beta)
    Kmat <- K_beta_mat(theta_tilde_restricted, Kvec, stress_mat, tau, ITs, beta)
    
    num <- drop(t(Q) %*% U)^2
    den <- drop(t(Q) %*% Kmat %*% Q)
    if (den <= 0) return(NA)
    K * num / den
}


                     
############### FIT: unrestricted + restricted under H(B) and H(C), for all beta  ############### 
                     

beta_vec    <- c(0, 0.2, 0.4, 0.6, 0.8, 1.0)
beta_labels <- c("MLE", "0.2", "0.4", "0.6", "0.8", "1.0")

cat("\nFitting unrestricted WMDPDE for all beta...\n")
estimates_unrestricted <- fit_all_betas(data_obs, Kvec, stress_mat, tau, ITs,
                                        beta_vec = beta_vec, init = c(10, -1, 0.2))

cat("Fitting restricted WMDPDE under H(B): alpha1 = 0 (local override in effect) ...\n")
estimates_restricted_B <- fit_all_betas_restricted_alpha1(
    data_obs, Kvec, stress_mat, tau, ITs,
    alpha1_fixed = alpha1_null, beta_vec = beta_vec, init = c(10, 0.2))

cat("Fitting restricted WMDPDE under H(C): sigma = sigma0 ...\n")
estimates_restricted_C <- fit_all_betas_restricted_sigma(
    data_obs, Kvec, stress_mat, tau, ITs,
    sigma_fixed = sigma0, beta_vec = beta_vec, init = c(10, -1))


                     
###############  DIAGNOSTIC CHECK: confirm restricted fit under H(B) actually moved ############### 
                     

cat("\n=== Diagnostic: restricted fit under H(B): alpha1 = 0 ===\n")
for (bi in seq_along(beta_vec)) {
    b  <- beta_vec[bi]
    nm <- if (b == 0) "MLE" else paste0("MDPDE_", b)
    tt <- estimates_restricted_B[[nm]]
    cat(sprintf("beta=%.1f  theta_tilde=(%.4f, %.4f, %.4f)\n", b, tt[1], tt[2], tt[3]))
}


                     
################  APPLY ALL THREE TESTS, FOR EACH BETA ############### 

                     
K <- sum(Kvec)
results <- data.frame()

for (bi in seq_along(beta_vec)) {
    
    b  <- beta_vec[bi]
    nm <- if (b == 0) "MLE" else paste0("MDPDE_", b)
    
    theta_hat_unr <- estimates_unrestricted[[nm]]
    theta_tilde_B <- estimates_restricted_B[[nm]]
    theta_tilde_C <- estimates_restricted_C[[nm]]
    
    # ---- H(A): theta = theta0 (simple) ----
    WA     <- wald_HA(theta_hat_unr, theta0, Kvec, stress_mat, tau, ITs, b)
    U_A    <- U_beta_vec(theta0, data_obs, Kvec, stress_mat, tau, ITs, b)
    Kmat_A <- K_beta_mat(theta0, Kvec, stress_mat, tau, ITs, b)
    RA     <- rao_HA(U_A, Kmat_A, K)
    
    # ---- H(B): alpha1 = 0 (composite) ----
    WB <- wald_composite_scalar(theta_hat_unr, j = 2, null_val = alpha1_null,
                                Kvec, stress_mat, tau, ITs, b)
    RB <- if (any(is.na(theta_tilde_B))) NA else
        rao_composite_scalar(theta_tilde_B, j = 2,
                             data_obs, Kvec, stress_mat, tau, ITs, b)
    
    # ---- H(C): sigma = sigma0 (composite) ----
    WC <- wald_composite_scalar(theta_hat_unr, j = 3, null_val = sigma0,
                                Kvec, stress_mat, tau, ITs, b)
    RC <- rao_composite_scalar(theta_tilde_C, j = 3,
                               data_obs, Kvec, stress_mat, tau, ITs, b)
    
    results <- rbind(results, data.frame(
        beta = beta_labels[bi],
        WA = WA, p_WA = 1 - pchisq(WA, df = 3),
        RA = RA, p_RA = 1 - pchisq(RA, df = 3),
        WB = WB, p_WB = 1 - pchisq(WB, df = 1),
        RB = RB, p_RB = 1 - pchisq(RB, df = 1),
        WC = WC, p_WC = 1 - pchisq(WC, df = 1),
        RC = RC, p_RC = 1 - pchisq(RC, df = 1)
    ))
}



cat("\n=== Test statistics and p-values for H(A), H(B), H(C) ===\n")
print(results, digits = 3, row.names = FALSE)

save(results, theta0, theta_true, sigma0, alpha1_null, beta_vec, beta_labels,
     file = "C:/Users/Kiran/WMDPDE_CyALT_lognormal/Testing/real_data_tests_v3_alpha1zero.RData")
cat("\nSaved: real_data_tests_v3_alpha1zero.RData\n")


fmt_cell <- function(stat_w, stat_r, digits = 2) {
    sprintf("%.*f (%.*f)", digits, stat_w, digits, stat_r)
}
fmt_p <- function(p_w, p_r, digits = 2) {
    sprintf("%.*f (%.*f)", digits, p_w, digits, p_r)
}

results_fmt <- data.frame(
    beta         = results$beta,
    HA_Statistic = mapply(fmt_cell, results$WA, results$RA),
    HA_p         = mapply(fmt_p,    results$p_WA, results$p_RA),
    HB_Statistic = mapply(fmt_cell, results$WB, results$RB),
    HB_p         = mapply(fmt_p,    results$p_WB, results$p_RB, digits = 3),
    HC_Statistic = mapply(fmt_cell, results$WC, results$RC),
    HC_p         = mapply(fmt_p,    results$p_WC, results$p_RC)
)

print(results_fmt, row.names = FALSE)
