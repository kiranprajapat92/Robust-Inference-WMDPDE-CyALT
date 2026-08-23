library(optimx)
library(MASS)

# =============================================================================
# MODEL FUNCTIONS
# =============================================================================

eta_fun <- function(alpha0, alpha1, s) {
  exp(-alpha0 - alpha1 * s)
}

mu_fun <- function(alpha0, alpha1, sC, sF, tau) {
  etaC <- eta_fun(alpha0, alpha1, sC)
  etaF <- eta_fun(alpha0, alpha1, sF)
  -log(tau * etaC + (1 - tau) * etaF)
}

sbar_fun <- function(alpha0, alpha1, sC, sF, tau) {
  etaC <- eta_fun(alpha0, alpha1, sC)
  etaF <- eta_fun(alpha0, alpha1, sF)
  denom <- tau * etaC + (1 - tau) * etaF
  (tau * sC * etaC + (1 - tau) * sF * etaF) / denom
}

aij_fun <- function(IT, mu, sigma) {
  if (IT <= 0) return(-Inf)
  (log(IT) - mu) / sigma
}

interval_prob <- function(mu, sigma, IT_lo, IT_hi) {
  pnorm(aij_fun(IT_hi, mu, sigma)) - pnorm(aij_fun(IT_lo, mu, sigma))
}

surv_prob <- function(mu, sigma, t_c) {
  1 - pnorm(aij_fun(t_c, mu, sigma))
}

p_i_theta <- function(alpha0, alpha1, sigma, sC, sF, tau, ITs) {
  mu  <- mu_fun(alpha0, alpha1, sC, sF, tau)
  L   <- length(ITs)
  t_c <- ITs[L]
  pvec      <- numeric(L + 1)
  IT_bounds <- c(0, ITs)
  for (j in seq_len(L))
    pvec[j] <- interval_prob(mu, sigma, IT_bounds[j], IT_bounds[j + 1])
  pvec[L + 1] <- surv_prob(mu, sigma, t_c)
  pvec[pvec < 1e-8] <- 1e-8
  pvec[pvec > 1 - 1e-8] <- 1 - 1e-8
  pvec / sum(pvec)
}

prob_list <- function(theta, stress_mat, tau, ITs) {
  alpha0 <- theta[1]; alpha1 <- theta[2]; sigma <- theta[3]
  R <- nrow(stress_mat)
  lapply(seq_len(R), function(i)
    p_i_theta(alpha0, alpha1, sigma,
              sC = stress_mat[i, 2], sF = stress_mat[i, 1], tau, ITs))
}

W_i_theta_matrix <- function(alpha0, alpha1, sigma, sC, sF, tau, ITs) {
  mu   <- mu_fun(alpha0, alpha1, sC, sF, tau)
  sbar <- sbar_fun(alpha0, alpha1, sC, sF, tau)
  L    <- length(ITs)
  t_c  <- ITs[L]
  IT_bounds <- c(0, ITs)
  W <- matrix(0, nrow = L + 1, ncol = 3)
  for (j in seq_len(L)) {
    a_hi <- aij_fun(IT_bounds[j + 1], mu, sigma)
    a_lo <- aij_fun(IT_bounds[j],     mu, sigma)
    phi_hi    <- dnorm(a_hi)
    phi_lo    <- if (is.infinite(a_lo)) 0 else dnorm(a_lo)
    adphi_lo  <- if (is.infinite(a_lo)) 0 else a_lo * phi_lo
    dphi_diff  <- phi_hi - phi_lo
    adphi_diff <- a_hi * phi_hi - adphi_lo
    W[j, 1] <- -(1 / sigma) * dphi_diff
    W[j, 2] <- -(sbar / sigma) * dphi_diff
    W[j, 3] <- -(1 / sigma) * adphi_diff
  }
  a_L   <- aij_fun(t_c, mu, sigma)
  phi_L <- dnorm(a_L)
  W[L + 1, 1] <- phi_L / sigma
  W[L + 1, 2] <- sbar * phi_L / sigma
  W[L + 1, 3] <- a_L * phi_L / sigma
  return(W)
}

# =============================================================================
# OBJECTIVE AND ESTIMATION
# =============================================================================

H_beta_objective <- function(theta, counts_list, Kvec, stress_mat, tau, ITs, beta) {
  alpha0 <- theta[1]; alpha1 <- theta[2]; sigma <- theta[3]
  if (sigma <= 0) return(1e10)
  if (alpha1 >= 0) return(1e10)
  R <- length(counts_list)
  K <- sum(Kvec)
  obj <- 0
  for (i in seq_len(R)) {
    pvec   <- p_i_theta(alpha0, alpha1, sigma,
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

fit_mdpde <- function(counts_list, Kvec, stress_mat, tau, ITs, beta, init = NULL) {
  if (is.null(init)) init <- c(10, -1, 0.2)
  result <- tryCatch(
    optimx(par = init, fn = H_beta_objective,
           counts_list = counts_list, Kvec = Kvec, stress_mat = stress_mat,
           tau = tau, ITs = ITs, beta = beta,
           method = "Nelder-Mead", control = list(maxit = 5000, reltol = 1e-10)),
    error = function(e) NULL)
  if (is.null(result) || result$convcode[1] == 9999) return(rep(NA, 3))
  par_out <- as.numeric(result[1, 1:3])
  if (any(is.na(par_out))) return(rep(NA, 3))
  return(par_out)
}

fit_all_betas <- function(counts_list, Kvec, stress_mat, tau, ITs,
                          beta_vec = c(0, 0.2, 0.4, 0.6, 0.8, 1),
                          init = NULL) {
  estimates    <- list()
  current_init <- init
  for (b in beta_vec) {
    label <- if (b == 0) "MLE" else paste0("MDPDE_", b)
    est   <- fit_mdpde(counts_list, Kvec, stress_mat, tau, ITs,
                       beta = b, init = current_init)
    estimates[[label]] <- est
    if (!any(is.na(est)) && !any(abs(est) > 1e6)) current_init <- est
  }
  return(estimates)
}

# ---- H2a composite: H0: alpha0 = alpha0,0, alpha1 and sigma unknown ----

H_beta_objective_restricted_alpha0 <- function(par2, counts_list, Kvec, stress_mat,
                                                tau, ITs, beta, alpha0_fixed) {
  theta_full <- c(alpha0_fixed, par2[1], par2[2])
  H_beta_objective(theta_full, counts_list, Kvec, stress_mat, tau, ITs, beta)
}

fit_mdpde_restricted_alpha0 <- function(counts_list, Kvec, stress_mat, tau, ITs,
                                        beta, alpha0_fixed, init = NULL) {
  if (is.null(init)) init <- c(-1, 0.2)
  result <- tryCatch(
    optimx(par = init, fn = H_beta_objective_restricted_alpha0,
           counts_list = counts_list, Kvec = Kvec, stress_mat = stress_mat,
           tau = tau, ITs = ITs, beta = beta, alpha0_fixed = alpha0_fixed,
           method = "Nelder-Mead", control = list(maxit = 5000, reltol = 1e-10)),
    error = function(e) NULL)
  if (is.null(result) || result$convcode[1] == 9999) return(rep(NA, 3))
  par2_out <- as.numeric(result[1, 1:2])
  if (any(is.na(par2_out))) return(rep(NA, 3))
  return(c(alpha0_fixed, par2_out[1], par2_out[2]))
}

fit_all_betas_restricted_alpha0 <- function(counts_list, Kvec, stress_mat, tau, ITs,
                                            alpha0_fixed,
                                            beta_vec = c(0, 0.2, 0.4, 0.6, 0.8, 1),
                                            init = NULL) {
  estimates    <- list()
  current_init <- init
  for (b in beta_vec) {
    label <- if (b == 0) "MLE" else paste0("MDPDE_", b)
    est   <- fit_mdpde_restricted_alpha0(counts_list, Kvec, stress_mat, tau, ITs,
                                         beta = b, alpha0_fixed = alpha0_fixed,
                                         init = current_init)
    estimates[[label]] <- est
    if (!any(is.na(est)) && !any(abs(est) > 1e6)) current_init <- c(est[2], est[3])
  }
  return(estimates)
}


# ---- H2b composite: H0: alpha1 = alpha1,0, alpha0 and sigma unknown ----

H_beta_objective_restricted_alpha1 <- function(par2, counts_list, Kvec, stress_mat,
                                                tau, ITs, beta, alpha1_fixed) {
  theta_full <- c(par2[1], alpha1_fixed, par2[2])
  H_beta_objective(theta_full, counts_list, Kvec, stress_mat, tau, ITs, beta)
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
  if (is.null(result) || result$convcode[1] == 9999) return(rep(NA, 3))
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

# ---- H(C) composite: H0: sigma = sigma0, alpha0 and alpha1 unknown ----

H_beta_objective_restricted_sigma <- function(par2, counts_list, Kvec, stress_mat,
                                               tau, ITs, beta, sigma_fixed) {
  theta_full <- c(par2[1], par2[2], sigma_fixed)
  H_beta_objective(theta_full, counts_list, Kvec, stress_mat, tau, ITs, beta)
}

fit_mdpde_restricted_sigma <- function(counts_list, Kvec, stress_mat, tau, ITs,
                                        beta, sigma_fixed, init = NULL) {
  if (is.null(init)) init <- c(10, -1)
  result <- tryCatch(
    optimx(par = init, fn = H_beta_objective_restricted_sigma,
           counts_list = counts_list, Kvec = Kvec, stress_mat = stress_mat,
           tau = tau, ITs = ITs, beta = beta, sigma_fixed = sigma_fixed,
           method = "Nelder-Mead", control = list(maxit = 5000, reltol = 1e-10)),
    error = function(e) NULL)
  if (is.null(result) || result$convcode[1] == 9999) return(rep(NA, 3))
  par2_out <- as.numeric(result[1, 1:2])
  if (any(is.na(par2_out))) return(rep(NA, 3))
  return(c(par2_out[1], par2_out[2], sigma_fixed))
}

fit_all_betas_restricted_sigma <- function(counts_list, Kvec, stress_mat, tau, ITs,
                                            sigma_fixed,
                                            beta_vec = c(0, 0.2, 0.4, 0.6, 0.8, 1),
                                            init = NULL) {
  estimates    <- list()
  current_init <- init
  for (b in beta_vec) {
    label <- if (b == 0) "MLE" else paste0("MDPDE_", b)
    est   <- fit_mdpde_restricted_sigma(counts_list, Kvec, stress_mat, tau, ITs,
                                        beta = b, sigma_fixed = sigma_fixed,
                                        init = current_init)
    estimates[[label]] <- est
    if (!any(is.na(est)) && !any(abs(est) > 1e6)) current_init <- c(est[1], est[2])
  }
  return(estimates)
}

# =============================================================================
# ASYMPTOTIC COVARIANCE
# =============================================================================

J_beta_mat <- function(theta, Kvec, stress_mat, tau, ITs, beta) {
  alpha0 <- theta[1]; alpha1 <- theta[2]; sigma <- theta[3]
  K <- sum(Kvec)
  J <- matrix(0, 3, 3)
  for (i in seq_len(nrow(stress_mat))) {
    pvec <- p_i_theta(alpha0, alpha1, sigma,
                      sC = stress_mat[i, 2], sF = stress_mat[i, 1], tau, ITs)
    Wi   <- W_i_theta_matrix(alpha0, alpha1, sigma,
                             sC = stress_mat[i, 2], sF = stress_mat[i, 1], tau, ITs)
    J    <- J + (Kvec[i] / K) * t(Wi) %*% diag(pvec^(beta - 1)) %*% Wi
  }
  return(J)
}

K_beta_mat <- function(theta, Kvec, stress_mat, tau, ITs, beta) {
  alpha0 <- theta[1]; alpha1 <- theta[2]; sigma <- theta[3]
  K  <- sum(Kvec)
  Km <- matrix(0, 3, 3)
  for (i in seq_len(nrow(stress_mat))) {
    pvec <- p_i_theta(alpha0, alpha1, sigma,
                      sC = stress_mat[i, 2], sF = stress_mat[i, 1], tau, ITs)
    Wi   <- W_i_theta_matrix(alpha0, alpha1, sigma,
                             sC = stress_mat[i, 2], sF = stress_mat[i, 1], tau, ITs)
    pb   <- pvec^beta
    Km   <- Km + (Kvec[i] / K) * t(Wi) %*% (diag(pvec^(2 * beta - 1)) - pb %o% pb) %*% Wi
  }
  return(Km)
}

Sigma_hat <- function(theta_hat, Kvec, stress_mat, tau, ITs, beta) {
  K    <- sum(Kvec)
  J    <- J_beta_mat(theta_hat, Kvec, stress_mat, tau, ITs, beta)
  Kmat <- K_beta_mat(theta_hat, Kvec, stress_mat, tau, ITs, beta)
  Jinv <- tryCatch(solve(J), error = function(e) matrix(NA, 3, 3))
  if (any(is.na(Jinv))) return(matrix(NA, 3, 3))
  (1 / K) * Jinv %*% Kmat %*% Jinv
}

U_beta_vec <- function(theta, counts_list, Kvec, stress_mat, tau, ITs, beta) {
  R <- nrow(stress_mat)
  K <- sum(Kvec)
  U <- numeric(3)

  for (i in seq_len(R)) {
    p_i    <- p_i_theta(theta[1], theta[2], theta[3],
                        sC = stress_mat[i,2], sF = stress_mat[i,1], tau, ITs)
    W_i    <- W_i_theta_matrix(theta[1], theta[2], theta[3],
                               sC = stress_mat[i,2], sF = stress_mat[i,1], tau, ITs)
    D_i    <- diag(p_i^(beta - 1))
    phat_i <- counts_list[[i]] / Kvec[i]

    U <- U + (Kvec[i] / K) * drop(t(W_i) %*% D_i %*% (p_i - phat_i))
  }
  U
}

theta_ci <- function(theta_hat, Kvec, stress_mat, tau, ITs, beta, level = 0.95) {
  Sig <- Sigma_hat(theta_hat, Kvec, stress_mat, tau, ITs, beta)
  if (any(is.na(Sig))) return(matrix(NA, 3, 2))
  se <- sqrt(diag(Sig))
  z  <- qnorm(1 - (1 - level) / 2)
  cbind(theta_hat - z * se, theta_hat + z * se)
}

# =============================================================================
# LIFETIME CHARACTERISTICS AT USE CONDITION
# =============================================================================

mu_use <- function(theta, s0C, s0F, tau) {
  mu_fun(theta[1], theta[2], s0C, s0F, tau)
}

quantile_use <- function(theta, s0C, s0F, tau, q = 0.5) {
  exp(mu_use(theta, s0C, s0F, tau) + theta[3] * qnorm(q))
}

quantile_use_grad <- function(theta, s0C, s0F, tau, q = 0.5) {
  tq    <- quantile_use(theta, s0C, s0F, tau, q)
  sbar0 <- sbar_fun(theta[1], theta[2], s0C, s0F, tau)
  c(tq, tq * sbar0, tq * qnorm(q))
}

quantile_use_se <- function(theta_hat, Kvec, stress_mat, tau, ITs, beta,
                            s0C, s0F, q = 0.5) {
  Sig  <- Sigma_hat(theta_hat, Kvec, stress_mat, tau, ITs, beta)
  if (any(is.na(Sig))) return(NA)
  grad <- quantile_use_grad(theta_hat, s0C, s0F, tau, q)
  sqrt(drop(t(grad) %*% Sig %*% grad))
}

quantile_use_ci <- function(theta_hat, Kvec, stress_mat, tau, ITs, beta,
                            s0C, s0F, q = 0.5, level = 0.95) {
  tq <- quantile_use(theta_hat, s0C, s0F, tau, q)
  se <- quantile_use_se(theta_hat, Kvec, stress_mat, tau, ITs, beta, s0C, s0F, q)
  if (is.na(se)) return(c(NA, NA))
  z  <- qnorm(1 - (1 - level) / 2)
  c(tq * exp(-z * se / tq), tq * exp(z * se / tq))
}

mttf_use <- function(theta, s0C, s0F, tau) {
  exp(mu_fun(theta[1], theta[2], s0C, s0F, tau) + theta[3]^2 / 2)
}

mttf_use_grad <- function(theta, s0C, s0F, tau) {
  m     <- mttf_use(theta, s0C, s0F, tau)
  sbar0 <- sbar_fun(theta[1], theta[2], s0C, s0F, tau)
  c(m, m * sbar0, m * theta[3])
}

mttf_use_se <- function(theta_hat, Kvec, stress_mat, tau, ITs, beta, s0C, s0F) {
  Sig  <- Sigma_hat(theta_hat, Kvec, stress_mat, tau, ITs, beta)
  if (any(is.na(Sig))) return(NA)
  grad <- mttf_use_grad(theta_hat, s0C, s0F, tau)
  sqrt(drop(t(grad) %*% Sig %*% grad))
}

mttf_use_ci <- function(theta_hat, Kvec, stress_mat, tau, ITs, beta,
                        s0C, s0F, level = 0.95) {
  m  <- mttf_use(theta_hat, s0C, s0F, tau)
  se <- mttf_use_se(theta_hat, Kvec, stress_mat, tau, ITs, beta, s0C, s0F)
  if (is.na(se)) return(c(NA, NA))
  z  <- qnorm(1 - (1 - level) / 2)
  c(m * exp(-z * se / m), m * exp(z * se / m))
}

reliability_use <- function(theta, s0C, s0F, tau, t0) {
  mu0 <- mu_use(theta, s0C, s0F, tau)
  1 - pnorm((log(t0) - mu0) / theta[3])
}

reliability_use_grad <- function(theta, s0C, s0F, tau, t0) {
  mu0   <- mu_use(theta, s0C, s0F, tau)
  sbar0 <- sbar_fun(theta[1], theta[2], s0C, s0F, tau)
  a0    <- (log(t0) - mu0) / theta[3]
  phi0  <- dnorm(a0)
  c(phi0 / theta[3], phi0 * sbar0 / theta[3], phi0 * a0 / theta[3])
}

reliability_use_se <- function(theta_hat, Kvec, stress_mat, tau, ITs, beta,
                               s0C, s0F, t0) {
  Sig  <- Sigma_hat(theta_hat, Kvec, stress_mat, tau, ITs, beta)
  if (any(is.na(Sig))) return(NA)
  grad <- reliability_use_grad(theta_hat, s0C, s0F, tau, t0)
  sqrt(drop(t(grad) %*% Sig %*% grad))
}

reliability_use_ci <- function(theta_hat, Kvec, stress_mat, tau, ITs, beta,
                               s0C, s0F, t0, level = 0.95) {
  r  <- reliability_use(theta_hat, s0C, s0F, tau, t0)
  se <- reliability_use_se(theta_hat, Kvec, stress_mat, tau, ITs, beta, s0C, s0F, t0)
  if (is.na(se)) return(c(NA, NA))
  z  <- qnorm(1 - (1 - level) / 2)
  S  <- exp(z * se / (r * (1 - r)))
  c(r / (r + (1 - r) * S), r / (r + (1 - r) / S))
}

# =============================================================================
# DATA SIMULATION
# =============================================================================

sim_cyalt_data <- function(theta, Kvec, stress_mat, tau, ITs, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  alpha0 <- theta[1]; alpha1 <- theta[2]; sigma <- theta[3]
  R <- nrow(stress_mat)
  counts_list <- vector("list", R)
  for (i in seq_len(R)) {
    pvec <- p_i_theta(alpha0, alpha1, sigma,
                      sC = stress_mat[i, 2], sF = stress_mat[i, 1], tau, ITs)
    counts_list[[i]] <- drop(rmultinom(1, size = Kvec[i], prob = pvec))
  }
  return(counts_list)
}

sim_cyalt_data_continuous <- function(theta, Kvec, stress_mat, tau, ITs,
                                      seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  alpha0 <- theta[1]; alpha1 <- theta[2]; sigma <- theta[3]
  R      <- nrow(stress_mat)
  L      <- length(ITs)
  breaks <- c(0, ITs, Inf)
  counts_list <- vector("list", R)
  for (i in seq_len(R)) {
    mu_i             <- mu_fun(alpha0, alpha1, stress_mat[i, 2], stress_mat[i, 1], tau)
    T_i              <- exp(mu_i + sigma * rnorm(Kvec[i]))
    bins             <- cut(T_i, breaks = breaks, right = TRUE, labels = FALSE)
    counts_list[[i]] <- as.integer(tabulate(bins, nbins = L + 1))
  }
  return(counts_list)
}

contaminate_data <- function(theta, Kvec, stress_mat, tau, ITs,
                             eps, cont_grp, cont_cell, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  alpha0  <- theta[1]; alpha1 <- theta[2]; sigma <- theta[3]
  R       <- nrow(stress_mat)
  L_plus1 <- length(ITs) + 1
  counts_list <- vector("list", R)
  for (i in seq_len(R)) {
    pvec <- p_i_theta(alpha0, alpha1, sigma,
                      sC = stress_mat[i, 2], sF = stress_mat[i, 1], tau, ITs)
    K_i  <- Kvec[i]
    if (i == cont_grp && eps > 0) {
      n_clean               <- floor(K_i * (1 - eps))
      n_cont                <- K_i - n_clean
      clean_counts          <- drop(rmultinom(1, size = n_clean, prob = pvec))
      cont_counts           <- integer(L_plus1)
      cont_counts[cont_cell] <- n_cont
      counts_list[[i]]      <- clean_counts + cont_counts
    } else {
      counts_list[[i]] <- drop(rmultinom(1, size = K_i, prob = pvec))
    }
  }
  return(counts_list)
}

# MULTI-CELL CONTAMINATION FUNCTION
# Replaced the contaminate_data() for the multi-cell case.
gen_data_multicell <- function(theta, Kvec, stress_mat, tau_val, ITs,
                               eps, cont_grp, cont_cells, seed) {
  set.seed(seed)
  R <- nrow(stress_mat)
  L <- length(ITs)
  
  counts_list <- vector("list", R)
  
  for (i in seq_len(R)) {
    sC   <- stress_mat[i, 2]
    sF   <- stress_mat[i, 1]
    mu_i <- theta[1] - log(tau_val*exp(-theta[2]*sC) + (1-tau_val)*exp(-theta[2]*sF))
    
    # interval probabilities
    pvec <- numeric(L+1)
    for (j in seq_len(L)) {
      lo <- if (j==1) -Inf else (log(ITs[j-1]) - mu_i) / theta[3]
      hi <- (log(ITs[j]) - mu_i) / theta[3]
      pvec[j] <- pnorm(hi) - pnorm(lo)
    }
    pvec[L+1] <- 1 - pnorm((log(ITs[L]) - mu_i) / theta[3])
    pvec <- pmax(pvec, 1e-10)
    pvec <- pvec / sum(pvec)
    
    K_i <- Kvec[i]
    
    if (i == cont_grp && eps > 0) {
      n_clean <- floor(K_i * (1 - eps))
      n_cont  <- K_i - n_clean
      
      # clean portion: drawn from true model
      clean_counts <- drop(rmultinom(1, n_clean, pvec))
      
      # contaminated portion: drawn uniformly across cont_cells
      q_cont              <- numeric(L+1)
      q_cont[cont_cells]  <- 1 / length(cont_cells)
      cont_counts         <- drop(rmultinom(1, n_cont, q_cont))
      
      counts_list[[i]] <- clean_counts + cont_counts
      
    } else {
      counts_list[[i]] <- drop(rmultinom(1, K_i, pvec))
    }
  }
  counts_list
}


is_valid_dataset <- function(counts_list) {
  all(sapply(counts_list, function(n_i) sum(n_i[-length(n_i)]) >= 1))
}

# =============================================================================
# BCa BOOTSTRAP
# =============================================================================

bootstrap_estimates_cyalt <- function(theta_hat, Kvec, stress_mat, tau, ITs,
                                      beta, char_fun, B = 500, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  bestimates <- numeric(B)
  for (b in seq_len(B)) {
    valid <- FALSE
    while (!valid) {
      dat <- sim_cyalt_data_continuous(theta_hat, Kvec, stress_mat, tau, ITs)
      if (is_valid_dataset(dat)) valid <- TRUE
    }
    est          <- fit_mdpde(dat, Kvec, stress_mat, tau, ITs, beta = beta, init = theta_hat)
    bestimates[b] <- char_fun(est)
  }
  return(bestimates)
}

jackknife_cyalt <- function(theta_hat, data_obs, Kvec, stress_mat, tau, ITs,
                            beta, char_fun) {
  R       <- length(data_obs)
  L_plus1 <- length(ITs) + 1
  jack_vals   <- c()
  cell_counts <- c()
  for (i in seq_len(R)) {
    for (j in seq_len(L_plus1 - 1)) {
      n_ij <- data_obs[[i]][j]
      if (n_ij == 0) next
      dat_jack               <- data_obs
      dat_jack[[i]][j]       <- n_ij - 1
      dat_jack[[i]][L_plus1] <- dat_jack[[i]][L_plus1] + 1
      est_jack     <- fit_mdpde(dat_jack, Kvec, stress_mat, tau, ITs,
                                beta = beta, init = theta_hat)
      jack_vals   <- c(jack_vals,   char_fun(est_jack))
      cell_counts <- c(cell_counts, n_ij)
    }
  }
  if (length(jack_vals) < 2) return(0)
  w        <- cell_counts
  jack_bar <- sum(w * jack_vals) / sum(w)
  num      <- sum(w * (jack_vals - jack_bar)^3)
  denom    <- sum(w * (jack_vals - jack_bar)^2)
  if (denom < 1e-12) return(0)
  (1 / 6) * num * denom^(-3/2)
}

bca_ci <- function(bestimates, char_est, gamma_hat, alpha = 0.05) {
  B    <- length(bestimates)
  prop <- mean(bestimates <= char_est)
  if (prop == 0) prop <- 0.0001
  if (prop == 1) prop <- 0.9999
  z0      <- qnorm(prop)
  z_alpha <- qnorm(1 - alpha / 2)
  q_lo <- pnorm(z0 + (z0 - z_alpha) / (1 - gamma_hat * (z0 - z_alpha)))
  q_hi <- pnorm(z0 + (z0 + z_alpha) / (1 - gamma_hat * (z0 + z_alpha)))
  bs   <- sort(bestimates)
  c(bs[max(floor(q_lo * B), 1)], bs[min(floor(q_hi * B), B)])
}

ci_all_three <- function(theta_hat, data_obs, Kvec, stress_mat, tau, ITs,
                         beta, char_fun, se_fun,
                         type = c("log", "logit", "none"),
                         B = 500, alpha = 0.05, seed = NULL) {
  type     <- match.arg(type)
  z        <- qnorm(1 - alpha / 2)
  char_est <- char_fun(theta_hat)
  se       <- se_fun(theta_hat)
  ci_direct <- c(char_est - z * se, char_est + z * se)
  if (type == "log") {
    ci_trans <- c(char_est * exp(-z * se / char_est),
                  char_est * exp( z * se / char_est))
  } else if (type == "logit") {
    S        <- exp(z * se / (char_est * (1 - char_est)))
    ci_trans <- c(char_est / (char_est + (1 - char_est) * S),
                  char_est / (char_est + (1 - char_est) / S))
  } else {
    ci_trans <- ci_direct        # unbounded parameter — no transformation
  }
  bestimates <- bootstrap_estimates_cyalt(theta_hat, Kvec, stress_mat, tau, ITs,
                                          beta, char_fun, B = B, seed = seed)
  gamma_hat  <- jackknife_cyalt(theta_hat, data_obs, Kvec, stress_mat, tau, ITs,
                                beta, char_fun)
  ci_bca     <- bca_ci(bestimates, char_est, gamma_hat, alpha = alpha)
  list(estimate       = char_est,
       ci_direct      = ci_direct,
       ci_transformed = ci_trans,
       ci_bca         = ci_bca,
       n_boot         = B)
}

# =============================================================================
# SUMMARY TABLE
# =============================================================================

summary_table <- function(estimates_list, Kvec, stress_mat, tau, ITs,
                          beta_vec = c(0, 0.2, 0.4, 0.6, 0.8, 1),
                          level = 0.95) {
  beta_labels <- ifelse(beta_vec == 0, "MLE", paste0("MDPDE_", beta_vec))
  z    <- qnorm(1 - (1 - level) / 2)
  rows <- list()
  for (bi in seq_along(beta_vec)) {
    lbl   <- beta_labels[bi]
    theta <- estimates_list[[lbl]]
    if (is.null(theta) || any(is.na(theta))) next
    Sig <- Sigma_hat(theta, Kvec, stress_mat, tau, ITs, beta_vec[bi])
    se  <- if (any(is.na(Sig))) rep(NA, 3) else sqrt(diag(Sig))
    rows[[lbl]] <- data.frame(
      beta   = beta_vec[bi],
      a0_est = theta[1], a0_lo = theta[1] - z * se[1], a0_hi = theta[1] + z * se[1],
      a1_est = theta[2], a1_lo = theta[2] - z * se[2], a1_hi = theta[2] + z * se[2],
      sg_est = theta[3], sg_lo = theta[3] - z * se[3], sg_hi = theta[3] + z * se[3]
    )
  }
  do.call(rbind, rows)
}

# =============================================================================
# SIMULATION STUDY
# =============================================================================

run_simulation <- function(theta_true, Kvec, stress_mat, tau, ITs,
                           eps_vec   = c(0, 0.05, 0.10, 0.20, 0.30, 0.40),
                           beta_vec  = c(0, 0.2, 0.4, 0.6, 0.8, 1),
                           cont_grp  = 1, cont_cell = 3,
                           n_iter    = 2000, init = NULL, verbose = TRUE) {
  if (is.null(init)) init <- theta_true
  beta_labels <- ifelse(beta_vec == 0, "MLE", paste0("MDPDE_", beta_vec))
  n_eps   <- length(eps_vec)
  n_beta  <- length(beta_vec)
  mse_arr <- array(0, dim = c(n_eps, n_beta, 3),
                   dimnames = list(as.character(eps_vec), beta_labels,
                                   c("alpha0", "alpha1", "sigma")))
  count_arr <- matrix(0, n_eps, n_beta)
  for (ei in seq_along(eps_vec)) {
    eps <- eps_vec[ei]
    if (verbose) cat(sprintf("  eps = %.2f ...\n", eps))
    for (rep in seq_len(n_iter)) {
      dat <- contaminate_data(theta_true, Kvec, stress_mat, tau, ITs,
                              eps = eps, cont_grp = cont_grp, cont_cell = cont_cell,
                              seed = rep * 17 + ei * 1000)
      if (!is_valid_dataset(dat)) next
      ests <- fit_all_betas(dat, Kvec, stress_mat, tau, ITs,
                            beta_vec = beta_vec, init = init)
      for (bi in seq_along(beta_vec)) {
        est <- ests[[beta_labels[bi]]]
        if (any(is.na(est))) next
        mse_arr[ei, bi, ] <- mse_arr[ei, bi, ] + (est - theta_true)^2
        count_arr[ei, bi] <- count_arr[ei, bi] + 1
      }
    }
    for (bi in seq_along(beta_vec)) {
      cnt <- count_arr[ei, bi]
      if (cnt > 0) mse_arr[ei, bi, ] <- mse_arr[ei, bi, ] / cnt
    }
  }
  return(list(mse = mse_arr, counts = count_arr))
}

print_mse_table <- function(sim_result, param_idx = 1, param_name = "alpha0") {
  mse <- sim_result$mse[, , param_idx]
  cat(sprintf("\nMSE table for %s:\n", param_name))
  cat(sprintf("%-8s", "eps"))
  cat(paste(sprintf("%-12s", colnames(mse)), collapse = ""), "\n")
  for (i in seq_len(nrow(mse))) {
    cat(sprintf("%-8s", rownames(mse)[i]))
    cat(paste(sprintf("%-12.6f", mse[i, ]), collapse = ""), "\n")
  }
}