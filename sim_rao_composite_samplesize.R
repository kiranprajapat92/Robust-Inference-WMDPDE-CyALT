rm(list = ls())

library(foreach)
library(doParallel)
library(optimx)
library(MASS)

source("/mnt/nfs/home/nkp117/largefiles/KP_Leandro_Maria_1/Testing/cyalt_lognormal_WMDPDE1.R")


# =============================================================================
# SETUP
# =============================================================================

theta0        <- c(5.0, -2.0, 0.5)
theta_alt_H1  <- c(5.2, -2.2, 0.48)
theta_alt_H2a <- c(5.2, -2.0, 0.5)
theta_alt_H2b <- c(5.0, -2.2, 0.5)
theta_alt_H3  <- c(5.2, -2.2, 0.5)

tau        <- 0.40
s_F        <- 0.40; s_1C <- 0.65; s_2C <- 1.00
stress_mat <- matrix(c(s_F, s_1C, s_F, s_2C), nrow = 2, byrow = TRUE)
ITs        <- c(15, 25, 35, 50, 65, 80)

beta_vec    <- c(0, 0.2, 0.4, 0.6, 0.8, 1.0)
beta_labels <- c("MLE", "0.2", "0.4", "0.6", "0.8", "1.0")
nb          <- length(beta_vec)

alpha_level <- 0.05
cv_H1 <- qchisq(1 - alpha_level, df = 3)
cv_H2 <- qchisq(1 - alpha_level, df = 1)
cv_H3 <- qchisq(1 - alpha_level, df = 2)

n_iter <- 1000

K_list      <- list(c(60,40), c(90,60), c(120,80),
                    c(150,100), c(180,120), c(210,140), c(240,160))
K_total_vec <- sapply(K_list, sum)
nK          <- length(K_list)

base_path <- "/mnt/nfs/home/nkp117/largefiles/KP_Leandro_Maria_1/Testing/"

H_alpha0 <- matrix(c(1, 0, 0), ncol = 1)
H_alpha1 <- matrix(c(0, 1, 0), ncol = 1)


# =============================================================================
# H1, H3 — SIMPLE, unchanged
# =============================================================================

rao_H1 <- function(U, Kmat, K) {
  Kinv <- tryCatch(solve(Kmat), error = function(e) matrix(NA, 3, 3))
  if (any(is.na(Kinv))) return(NA)
  K * drop(t(U) %*% Kinv %*% U)
}

rao_H3 <- function(U, Kmat, K) {
  Ksub  <- Kmat[1:2, 1:2]
  Kinv2 <- tryCatch(solve(Ksub), error = function(e) matrix(NA, 2, 2))
  if (any(is.na(Kinv2))) return(NA)
  U12 <- U[1:2]
  K * drop(t(U12) %*% Kinv2 %*% U12)
}


# =============================================================================
# H2a, H2b — COMPOSITE, everything at theta_tilde
# =============================================================================

rao_composite <- function(theta_tilde, Kvec, stress_mat, tau, ITs,
                          beta, counts_list, K, H_constraint) {
  Jmat <- J_beta_mat(theta_tilde, Kvec, stress_mat, tau, ITs, beta)
  Jinv <- tryCatch(solve(Jmat), error = function(e) matrix(NA, 3, 3))
  if (any(is.na(Jinv))) return(NA)
  c_vec <- Jinv %*% H_constraint
  denom <- drop(t(H_constraint) %*% c_vec)
  if (is.na(denom) || denom == 0) return(NA)
  Q <- c_vec / denom
  U_tilde    <- U_beta_vec(theta_tilde, counts_list, Kvec, stress_mat, tau, ITs, beta)
  Kmat_tilde <- K_beta_mat(theta_tilde, Kvec, stress_mat, tau, ITs, beta)
  num <- drop(t(Q) %*% U_tilde)
  var <- drop(t(Q) %*% Kmat_tilde %*% Q)
  if (is.na(var) || var <= 0) return(NA)
  K * num^2 / var
}


# =============================================================================
# SINGLE REPLICATION — no contamination, varying sample size
# =============================================================================

one_rep_rao_K <- function(rep_id, Kvec) {

  dat0    <- gen_data_multicell(theta0,        Kvec, stress_mat, tau, ITs,
                                0, 1, c(1,2,3), NULL)
  dat_H1  <- gen_data_multicell(theta_alt_H1,  Kvec, stress_mat, tau, ITs,
                                0, 1, c(1,2,3), NULL)
  dat_H2a <- gen_data_multicell(theta_alt_H2a, Kvec, stress_mat, tau, ITs,
                                0, 1, c(1,2,3), NULL)
  dat_H2b <- gen_data_multicell(theta_alt_H2b, Kvec, stress_mat, tau, ITs,
                                0, 1, c(1,2,3), NULL)
  dat_H3  <- gen_data_multicell(theta_alt_H3,  Kvec, stress_mat, tau, ITs,
                                0, 1, c(1,2,3), NULL)

  if (!is_valid_dataset(dat0)) return(NULL)

  K   <- sum(Kvec)
  out <- numeric(nb * 8)

  tilde0_H2a_all <- fit_all_betas_restricted_alpha0(dat0, Kvec, stress_mat, tau, ITs,
                                                     alpha0_fixed = theta0[1],
                                                     beta_vec = beta_vec,
                                                     init = c(theta0[2], theta0[3]))
  tildeA_H2a_all <- fit_all_betas_restricted_alpha0(dat_H2a, Kvec, stress_mat, tau, ITs,
                                                     alpha0_fixed = theta0[1],
                                                     beta_vec = beta_vec,
                                                     init = c(theta0[2], theta0[3]))

  tilde0_H2b_all <- fit_all_betas_restricted_alpha1(dat0, Kvec, stress_mat, tau, ITs,
                                                     alpha1_fixed = theta0[2],
                                                     beta_vec = beta_vec,
                                                     init = c(theta0[1], theta0[3]))
  tildeA_H2b_all <- fit_all_betas_restricted_alpha1(dat_H2b, Kvec, stress_mat, tau, ITs,
                                                     alpha1_fixed = theta0[2],
                                                     beta_vec = beta_vec,
                                                     init = c(theta0[1], theta0[3]))

  for (bi in seq_along(beta_vec)) {

    b   <- beta_vec[bi]
    key <- if (b == 0) "MLE" else paste0("MDPDE_", b)

    Kmat0 <- K_beta_mat(theta0, Kvec, stress_mat, tau, ITs, b)
    U0    <- U_beta_vec(theta0, dat0,   Kvec, stress_mat, tau, ITs, b)
    UH1   <- U_beta_vec(theta0, dat_H1, Kvec, stress_mat, tau, ITs, b)
    UH3   <- U_beta_vec(theta0, dat_H3, Kvec, stress_mat, tau, ITs, b)

    lev_H1 <- as.numeric(rao_H1(U0,  Kmat0, K) > cv_H1)
    lev_H3 <- as.numeric(rao_H3(U0,  Kmat0, K) > cv_H3)
    pow_H1 <- as.numeric(rao_H1(UH1, Kmat0, K) > cv_H1)
    pow_H3 <- as.numeric(rao_H3(UH3, Kmat0, K) > cv_H3)

    tilde0_H2a <- tilde0_H2a_all[[key]]; tildeA_H2a <- tildeA_H2a_all[[key]]
    lev_H2a <- if (any(is.na(tilde0_H2a))) NA else
      as.numeric(rao_composite(tilde0_H2a, Kvec, stress_mat, tau, ITs, b, dat0,    K, H_alpha0) > cv_H2)
    pow_H2a <- if (any(is.na(tildeA_H2a))) NA else
      as.numeric(rao_composite(tildeA_H2a, Kvec, stress_mat, tau, ITs, b, dat_H2a, K, H_alpha0) > cv_H2)

    tilde0_H2b <- tilde0_H2b_all[[key]]; tildeA_H2b <- tildeA_H2b_all[[key]]
    lev_H2b <- if (any(is.na(tilde0_H2b))) NA else
      as.numeric(rao_composite(tilde0_H2b, Kvec, stress_mat, tau, ITs, b, dat0,    K, H_alpha1) > cv_H2)
    pow_H2b <- if (any(is.na(tildeA_H2b))) NA else
      as.numeric(rao_composite(tildeA_H2b, Kvec, stress_mat, tau, ITs, b, dat_H2b, K, H_alpha1) > cv_H2)

    idx <- (bi - 1) * 8
    out[idx+1] <- lev_H1;  out[idx+2] <- lev_H2a
    out[idx+3] <- lev_H2b; out[idx+4] <- lev_H3
    out[idx+5] <- pow_H1;  out[idx+6] <- pow_H2a
    out[idx+7] <- pow_H2b; out[idx+8] <- pow_H3
  }
  out
}


# =============================================================================
# PARALLEL SETUP
# =============================================================================

n_cores <- 48
cl      <- makeCluster(n_cores)
registerDoParallel(cl)
cat(sprintf("Using %d cores\n", n_cores))

clusterExport(cl, varlist = c(
  "theta0","theta_alt_H1","theta_alt_H2a","theta_alt_H2b","theta_alt_H3",
  "tau","stress_mat","ITs","beta_vec","nb",
  "cv_H1","cv_H2","cv_H3",
  "H_alpha0","H_alpha1",
  "rao_H1","rao_H3","rao_composite","one_rep_rao_K",
  "gen_data_multicell","is_valid_dataset",
  "fit_all_betas_restricted_alpha0","fit_mdpde_restricted_alpha0",
  "H_beta_objective_restricted_alpha0",
  "fit_all_betas_restricted_alpha1","fit_mdpde_restricted_alpha1",
  "H_beta_objective_restricted_alpha1",
  "H_beta_objective","U_beta_vec",
  "p_i_theta","W_i_theta_matrix","mu_fun","eta_fun",
  "sbar_fun","aij_fun","interval_prob","surv_prob",
  "J_beta_mat","K_beta_mat"
))
clusterEvalQ(cl, { library(optimx); library(MASS) })


# =============================================================================
# RUN SIMULATION
# =============================================================================

level_K <- array(NA, dim = c(nK, nb, 4),
                 dimnames = list(as.character(K_total_vec),
                                 beta_labels, c("H1","H2a","H2b","H3")))
power_K <- array(NA, dim = c(nK, nb, 4),
                 dimnames = list(as.character(K_total_vec),
                                 beta_labels, c("H1","H2a","H2b","H3")))

t_start <- proc.time()
cat("Running Rao-type simulation: varying sample size, no contamination\n\n")

for (ki in seq_along(K_list)) {

  Kvec <- K_list[[ki]]
  cat(sprintf("  K = %d ...\n", sum(Kvec)))

  reps <- foreach(
    rep            = seq_len(n_iter),
    .combine       = "rbind",
    .packages      = c("optimx","MASS"),
    .errorhandling = "remove"
  ) %dopar% {
    suppressWarnings(one_rep_rao_K(rep, Kvec))
  }

  if (!is.null(reps) && nrow(reps) > 0) {
    reps <- reps[complete.cases(reps), , drop = FALSE]
    for (bi in seq_along(beta_vec)) {
      idx <- (bi - 1) * 8
      for (h in 1:4) {
        level_K[ki, bi, h] <- mean(reps[, idx + h],     na.rm = TRUE)
        power_K[ki, bi, h] <- mean(reps[, idx + 4 + h], na.rm = TRUE)
      }
    }
  }

  save(level_K, power_K, K_total_vec, beta_vec, beta_labels,
       theta0, theta_alt_H1, theta_alt_H2a, theta_alt_H2b, theta_alt_H3,
       alpha_level,
       file = paste0(base_path, "sim_rao_mixed_samplesize.RData"))
  cat(sprintf("    saved after K = %d\n", sum(Kvec)))
}

cat(sprintf("\nDone in %.1f minutes\n",
            (proc.time() - t_start)["elapsed"] / 60))

save(level_K, power_K, K_total_vec, beta_vec, beta_labels,
     theta0, theta_alt_H1, theta_alt_H2a, theta_alt_H2b, theta_alt_H3,
     alpha_level,
     file = paste0(base_path, "results_rao_mixed_samplesize.RData"))
cat("Final results saved.\n")

stopCluster(cl)


# =============================================================================
# PLOTTING — set 1: mixed, all four hypotheses (H1, H2a composite, H2b composite, H3)
# =============================================================================

cols <- c("#000000","#0072B2","#009E73","#D55E00","#7B2FBE","#CC79A7")
ltys <- rep(1, 6); lwds <- rep(1.3, 6); pchs <- c(16, 1, 2, 5, 6, 0)

leg_expr <- c(expression(MLE ~ (beta==0)), expression(beta==0.2), expression(beta==0.4),
             expression(beta==0.6), expression(beta==0.8), expression(beta=="1.0"))

hyp_labels <- list(
  expression(H[0]:~alpha[0]==alpha[list(0,0)]~","~alpha[1]==alpha[list(1,0)]~","~sigma==sigma[0]),
  expression(H[0]:~alpha[0]==alpha[list(0,0)]~(alpha[1]*","~sigma~"unknown")),
  expression(H[0]:~alpha[1]==alpha[list(1,0)]~(alpha[0]*","~sigma~"unknown")),
  expression(H[0]:~alpha[0]==alpha[list(0,0)]~","~alpha[1]==alpha[list(1,0)]~(sigma~"known"))
)

one_panel <- function(xv, ymat, xlab, main_lab, ylab_txt,
                      ref_line = FALSE, ref_val = 0.05, leg_pos = "topright") {
  ymin <- max(0, min(ymat, na.rm = TRUE) - 0.05)
  ymax <- min(1, max(ymat, na.rm = TRUE) + 0.05)
  plot(NA, xlim = range(xv), ylim = c(ymin, ymax), xlab = xlab, ylab = "", main = main_lab,
       cex.main = 2.3, cex.lab = 2.2, cex.axis = 2.0, font.main = 2, font.lab = 1, las = 1, bty = "o")
  title(ylab = ylab_txt, line = 5.2, cex.lab = 2.2, font.lab = 1)
  grid(nx = NA, ny = 5, col = "grey93", lty = 1, lwd = 0.8)
  if (ref_line) abline(h = ref_val, lty = 2, col = "grey40", lwd = 1.2)
  box()
  for (bi in seq_along(beta_vec)) {
    lines(xv, ymat[, bi], col = cols[bi], lty = ltys[bi], lwd = lwds[bi], type = "b")
    points(xv, ymat[, bi], col = cols[bi], pch = pchs[bi], cex = 1.3, bg = cols[bi])
  }
  legend(leg_pos, legend = leg_expr, col = cols, lty = ltys, lwd = lwds,
         pch = pchs, pt.bg = cols, pt.cex = 1.2, bty = "n", cex = 1.7, y.intersp = 1.0)
}

plot_level_K <- function() {
  par(mfrow = c(1,4), mar = c(7, 8.5, 5.5, 2), font.axis = 1)
  for (h in 1:4)
    one_panel(K_total_vec, level_K[,,h], xlab = "Sample size K",
              main_lab = hyp_labels[[h]], ylab_txt = "Empirical level",
              ref_line = TRUE, ref_val = alpha_level, leg_pos = "topright")
}

plot_power_K <- function() {
  par(mfrow = c(1,4), mar = c(7, 8.5, 5.5, 2), font.axis = 1)
  for (h in 1:4)
    one_panel(K_total_vec, power_K[,,h], xlab = "Sample size K",
              main_lab = hyp_labels[[h]], ylab_txt = "Empirical power",
              ref_line = FALSE, leg_pos = "bottomright")
}

pdf(paste0(base_path,"fig_level_rao_mixed_samplesize.pdf"), width = 18, height = 5)
plot_level_K(); dev.off()

pdf(paste0(base_path,"fig_power_rao_mixed_samplesize.pdf"), width = 18, height = 5)
plot_power_K(); dev.off()


# =============================================================================
# PLOTTING — set 2: composite-only, H2a and H2b, all in one row
# =============================================================================

hyp_labels_composite <- list(
  expression(H[0]:~alpha[0]==alpha[list(0,0)]~(alpha[1]*","~sigma~"unknown")),
  expression(H[0]:~alpha[1]==alpha[list(1,0)]~(alpha[0]*","~sigma~"unknown"))
)

plot_composite_combined_K <- function() {
  par(mfrow = c(1,4), mar = c(7, 8.5, 5.5, 2), font.axis = 1)
  one_panel(K_total_vec, level_K[,,"H2a"], xlab = "Sample size K",
            main_lab = hyp_labels_composite[[1]], ylab_txt = "Empirical level",
            ref_line = TRUE, ref_val = alpha_level, leg_pos = "top")
  one_panel(K_total_vec, level_K[,,"H2b"], xlab = "Sample size K",
            main_lab = hyp_labels_composite[[2]], ylab_txt = "Empirical level",
            ref_line = TRUE, ref_val = alpha_level, leg_pos = "top")
  one_panel(K_total_vec, power_K[,,"H2a"], xlab = "Sample size K",
            main_lab = hyp_labels_composite[[1]], ylab_txt = "Empirical power",
            ref_line = FALSE, leg_pos = "topleft")
  one_panel(K_total_vec, power_K[,,"H2b"], xlab = "Sample size K",
            main_lab = hyp_labels_composite[[2]], ylab_txt = "Empirical power",
            ref_line = FALSE, leg_pos = "topleft")
}

pdf(paste0(base_path,"fig_rao_composite_samplesize.pdf"), width = 18, height = 5)
plot_composite_combined_K(); dev.off()

cat("All figures saved.\n")