rm(list = ls())

library(foreach)
library(doParallel)
library(optimx)
library(MASS)

source("/mnt/nfs/home/nkp117/largefiles/KP_Leandro_Maria_1/Testing/cyalt_lognormal_WMDPDE1.R")
# source("C:/Users/Kiran/Downloads/WMDPDE_CyALT_lognormal/Main_codes/Testing/cyalt_lognormal_WMDPDE1.R")


n_cores <- 48


theta0       <- c(5.0, -2.0, 0.5)
theta_alt_H3 <- c(5.2, -2.2, 0.5)   # sigma unchanged -- H3 assumes sigma known/fixed

tau        <- 0.40
s_F        <- 0.40;  s_1C <- 0.65;  s_2C <- 1.00
s_0F       <- 0.00;  s_0C <- 0.20
stress_mat <- matrix(c(s_F, s_1C, s_F, s_2C), nrow = 2, byrow = TRUE)
ITs        <- c(15, 25, 35, 50, 65, 80)

beta_vec    <- c(0, 0.2, 0.4, 0.6, 0.8, 1.0)
beta_labels <- c("MLE", "0.2", "0.4", "0.6", "0.8", "1.0")
nb          <- length(beta_vec)

alpha_level <- 0.05
cv_H3 <- qchisq(1 - alpha_level, df = 2)

n_iter     <- 1000
Kvec       <- c(120, 80)
eps_vec    <- seq(0, 0.12, by = 0.02)   # 0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12
ne         <- length(eps_vec)
cont_grp   <- 1
cont_cells <- c(1, 2, 3)

base_path <- "/mnt/nfs/home/nkp117/largefiles/KP_Leandro_Maria_1/Testing/"
# base_path <- "C:/Users/Kiran/Downloads/WMDPDE_CyALT_lognormal/Main_codes/Testing/"

############# SIX CANDIDATE H3 STATISTICS #############
# theta_hat below is always the FULL 3-vector (alpha0, alpha1, sigma) --
# for the "restricted" numerator (stats 2, 4) it is c(alpha0_hat, alpha1_hat,
# sigma0), the output of fit_mdpde_restricted_sigma, which already keeps the
# 3-vector shape (sigma fixed at sigma0 in the third slot).
#
# Panel 5 has no formula -- there is no individual-entry (J_beta, K_beta)
# version of the Remark for H0^(3); that case was left blank by Leandro.

# 1. Section 7.2 -- joint estimation, submatrix of Sigma THEN inverted
stat1_sec72_joint <- function(theta_hat, theta0, Kvec, stress_mat, tau, ITs, beta) {
  Sig <- Sigma_hat(theta0, Kvec, stress_mat, tau, ITs, beta)
  if (any(is.na(Sig))) return(NA)
  Sig12  <- Sig[1:2, 1:2]
  Sinv12 <- tryCatch(solve(Sig12), error = function(e) matrix(NA, 2, 2))
  if (any(is.na(Sinv12))) return(NA)
  diff <- theta_hat[1:2] - theta0[1:2]
  drop(t(diff) %*% Sinv12 %*% diff)
}

# 2. Section 7.2 formula -- same denominator, restricted (sigma fixed) 2D numerator
stat2_sec72_plugin <- function(theta_hat_restricted, theta0, Kvec, stress_mat, tau, ITs, beta) {
  stat1_sec72_joint(theta_hat_restricted, theta0, Kvec, stress_mat, tau, ITs, beta)
}

# 3. Remark quadratic form (Omega submatrix, NOT re-inverted) -- joint estimation
stat3_omega_joint <- function(theta_hat, theta0, Kvec, stress_mat, tau, ITs, beta) {
  Sig <- Sigma_hat(theta0, Kvec, stress_mat, tau, ITs, beta)
  if (any(is.na(Sig))) return(NA)
  Omega <- tryCatch(solve(Sig), error = function(e) matrix(NA, 3, 3))
  if (any(is.na(Omega))) return(NA)
  Omega12 <- Omega[1:2, 1:2]
  diff    <- theta_hat[1:2] - theta0[1:2]
  drop(t(diff) %*% Omega12 %*% diff)
}

# 4. -- no formula, panel left blank (see note above)

# 5. Omega submatrix (same denominator as stat 3) -- restricted (sigma fixed) 2D plug-in numerator
stat5_omega_plugin <- function(theta_hat_restricted, theta0, Kvec, stress_mat, tau, ITs, beta) {
  stat3_omega_joint(theta_hat_restricted, theta0, Kvec, stress_mat, tau, ITs, beta)
}

# 6. Rao-type test for H3 -- no estimation, score evaluated at theta0 using the given data
stat6_rao <- function(dat, theta0, Kvec, stress_mat, tau, ITs, beta) {
  Kmat <- K_beta_mat(theta0, Kvec, stress_mat, tau, ITs, beta)
  K12  <- Kmat[1:2, 1:2]
  Kinv12 <- tryCatch(solve(K12), error = function(e) matrix(NA, 2, 2))
  if (any(is.na(Kinv12))) return(NA)
  U <- U_beta_vec(theta0, dat, Kvec, stress_mat, tau, ITs, beta)
  if (any(is.na(U))) return(NA)
  K <- sum(Kvec)
  U12 <- U[1:2]
  K * drop(t(U12) %*% Kinv12 %*% U12)
}


############ SINGLE REPLICATION ##########

one_rep_eps <- function(rep_id, eps) {

  dat0   <- gen_data_multicell(theta0,       Kvec, stress_mat, tau, ITs,
                               eps, cont_grp, cont_cells, NULL)
  dat_H3 <- gen_data_multicell(theta_alt_H3, Kvec, stress_mat, tau, ITs,
                               eps, cont_grp, cont_cells, NULL)

  if (!is_valid_dataset(dat0)) return(NULL)

  # 12 values per beta: level x 6 statistics, power x 6 statistics
  # (slots 5 and 11 are always NA -- panel 5 has no formula)
  out <- numeric(nb * 12)

  for (bi in seq_along(beta_vec)) {

    b   <- beta_vec[bi]
    key <- if (b == 0) "MLE" else paste0("MDPDE_", b)

    # joint 3D fits (needed for stats 1, 3)
    est0_joint  <- fit_all_betas(dat0,   Kvec, stress_mat, tau, ITs,
                                 beta_vec = b, init = theta0)[[key]]
    estH3_joint <- fit_all_betas(dat_H3, Kvec, stress_mat, tau, ITs,
                                 beta_vec = b, init = theta_alt_H3)[[key]]

    # restricted (sigma fixed) 2D fits, using the real source-file function
    # (needed for stats 2, 4)
    est0_restr  <- fit_mdpde_restricted_sigma(dat0,   Kvec, stress_mat, tau, ITs,
                                              beta = b, sigma_fixed = theta0[3],
                                              init = theta0[1:2])
    estH3_restr <- fit_mdpde_restricted_sigma(dat_H3, Kvec, stress_mat, tau, ITs,
                                              beta = b, sigma_fixed = theta0[3],
                                              init = theta_alt_H3[1:2])

    # --- level (null data) ---
    lev1 <- if (any(is.na(est0_joint))) NA else
      as.numeric(stat1_sec72_joint(est0_joint, theta0,Kvec,stress_mat,tau,ITs,b) > cv_H3)
    lev2 <- if (any(is.na(est0_restr))) NA else
      as.numeric(stat2_sec72_plugin(est0_restr, theta0,Kvec,stress_mat,tau,ITs,b) > cv_H3)
    lev3 <- if (any(is.na(est0_joint))) NA else
      as.numeric(stat3_omega_joint(est0_joint, theta0,Kvec,stress_mat,tau,ITs,b) > cv_H3)
    lev4 <- NA   # panel 4 left blank -- no formula
    lev5 <- if (any(is.na(est0_restr))) NA else
      as.numeric(stat5_omega_plugin(est0_restr, theta0,Kvec,stress_mat,tau,ITs,b) > cv_H3)
    lev6 <- {
      s <- stat6_rao(dat0, theta0, Kvec, stress_mat, tau, ITs, b)
      if (is.na(s)) NA else as.numeric(s > cv_H3)
    }

    # --- power (alternative data) ---
    pow1 <- if (any(is.na(estH3_joint))) NA else
      as.numeric(stat1_sec72_joint(estH3_joint, theta0,Kvec,stress_mat,tau,ITs,b) > cv_H3)
    pow2 <- if (any(is.na(estH3_restr))) NA else
      as.numeric(stat2_sec72_plugin(estH3_restr, theta0,Kvec,stress_mat,tau,ITs,b) > cv_H3)
    pow3 <- if (any(is.na(estH3_joint))) NA else
      as.numeric(stat3_omega_joint(estH3_joint, theta0,Kvec,stress_mat,tau,ITs,b) > cv_H3)
    pow4 <- NA   # panel 4 left blank -- no formula
    pow5 <- if (any(is.na(estH3_restr))) NA else
      as.numeric(stat5_omega_plugin(estH3_restr, theta0,Kvec,stress_mat,tau,ITs,b) > cv_H3)
    pow6 <- {
      s <- stat6_rao(dat_H3, theta0, Kvec, stress_mat, tau, ITs, b)
      if (is.na(s)) NA else as.numeric(s > cv_H3)
    }

    idx <- (bi - 1) * 12
    out[idx+1]  <- lev1; out[idx+2]  <- lev2; out[idx+3]  <- lev3
    out[idx+4]  <- lev4; out[idx+5]  <- lev5; out[idx+6]  <- lev6
    out[idx+7]  <- pow1; out[idx+8]  <- pow2; out[idx+9]  <- pow3
    out[idx+10] <- pow4; out[idx+11] <- pow5; out[idx+12] <- pow6
  }
  out
}


############# PARALLEL SETUP ######################

cl <- makeCluster(n_cores)
registerDoParallel(cl)
cat(sprintf("Using %d cores\n", n_cores))

clusterExport(cl, varlist = c(
  "theta0","theta_alt_H3","tau","stress_mat","ITs","beta_vec","nb","Kvec",
  "cv_H3","cont_grp","cont_cells",
  "one_rep_eps","gen_data_multicell","is_valid_dataset",
  "fit_all_betas","fit_mdpde","H_beta_objective",
  "fit_mdpde_restricted_sigma","H_beta_objective_restricted_sigma",
  "p_i_theta","W_i_theta_matrix","mu_fun","eta_fun",
  "sbar_fun","aij_fun","interval_prob","surv_prob",
  "J_beta_mat","K_beta_mat","Sigma_hat","U_beta_vec",
  "stat1_sec72_joint","stat2_sec72_plugin","stat3_omega_joint",
  "stat5_omega_plugin","stat6_rao"
))
clusterEvalQ(cl, { library(optimx); library(MASS) })


# RUN SIMULATION

stat_names <- c("Sec7.2_joint","Sec7.2_plugin","Omega_joint",
                "Blank","Omega_plugin","Rao")

level_eps <- array(NA, dim = c(ne, nb, 6),
                   dimnames = list(as.character(eps_vec), beta_labels, stat_names))
power_eps <- array(NA, dim = c(ne, nb, 6),
                   dimnames = list(as.character(eps_vec), beta_labels, stat_names))

t_start <- proc.time()
cat("Running simulation: all 6 H3 formulas (panel 5 blank), K = 200\n\n")

for (ei in seq_along(eps_vec)) {

  eps <- eps_vec[ei]
  cat(sprintf("  eps = %.2f ...\n", eps))

  reps <- foreach(
    rep            = seq_len(n_iter),
    .combine       = "rbind",
    .packages      = c("optimx","MASS"),
    .errorhandling = "remove"
  ) %dopar% {
    suppressWarnings(one_rep_eps(rep, eps))
  }

  if (!is.null(reps) && nrow(reps) > 0) {
    # complete.cases check now excludes columns 4 and 10 (always NA, panel 4)
    keep_cols <- setdiff(seq_len(nb * 12), unlist(lapply(0:(nb-1), function(bi) (bi*12) + c(4, 10))))
    reps <- reps[complete.cases(reps[, keep_cols]), , drop = FALSE]
    for (bi in seq_along(beta_vec)) {
      idx <- (bi - 1) * 12
      for (s in c(1,2,3,5,6)) {
        level_eps[ei, bi, s] <- mean(reps[, idx + s],     na.rm = TRUE)
        power_eps[ei, bi, s] <- mean(reps[, idx + 6 + s], na.rm = TRUE)
      }
      # s = 4 stays NA in both arrays
    }
  }

  save(level_eps, power_eps, eps_vec, beta_vec, beta_labels, Kvec,
       theta0, theta_alt_H3, alpha_level, stat_names,
       file = paste0(base_path, "sim_H3_all_formulas.RData"))
  cat(sprintf("    saved after eps = %.2f\n", eps))
}

save(level_eps, power_eps, eps_vec, beta_vec, beta_labels, Kvec,
     theta0, theta_alt_H3, alpha_level, stat_names,
     file = paste0(base_path, "results_H3_all_formulas.RData"))
cat("Final results saved: results_H3_all_formulas.RData\n")

cat(sprintf("\nDone in %.1f minutes\n",
            (proc.time() - t_start)["elapsed"] / 60))
stopCluster(cl)


#################### Plots: all 6 formulas, level and power ####################

cols <- c("#000000","#0072B2","#009E73","#D55E00","#7B2FBE","#CC79A7")
pchs <- c(16, 1, 2, 5, 6, 0)

leg_expr <- c(
  expression(MLE ~ (beta==0)), expression(beta==0.2), expression(beta==0.4),
  expression(beta==0.6), expression(beta==0.8), expression(beta=="1.0")
)

panel_titles <- c(
  "1: Sec 7.2 (joint, submatrix then inv.)",
  "2: Sec 7.2 formula + 2D plug-in num.",
  "3: Omega submatrix (joint, no re-inv.)",
  "4: left blank",
  "5: Omega submatrix (2D plug-in num.)",
  "6: Rao-type"
)

one_panel <- function(xv, ymat, xlab, main_lab, ylab_txt,
                      ref_line = FALSE, ref_val = 0.05, leg_pos = "topright") {
  ymin <- max(0, min(ymat, na.rm = TRUE) - 0.02)
  ymax <- min(1, max(ymat, na.rm = TRUE) + 0.02)
  plot(NA, xlim = range(xv), ylim = c(ymin, ymax), xlab = xlab, ylab = "",
       main = main_lab, cex.main = 1.35, cex.lab = 1.3, cex.axis = 1.1,
       font.main = 2, las = 1, bty = "o")
  title(ylab = ylab_txt, line = 3.2, cex.lab = 1.3)
  grid(nx = NA, ny = 5, col = "grey93", lty = 1, lwd = 0.8)
  if (ref_line) abline(h = ref_val, lty = 2, col = "grey40", lwd = 1.2)
  box()
  for (bi in seq_along(beta_vec)) {
    lines(xv, ymat[, bi], col = cols[bi], lwd = 1.2, type = "b")
    points(xv, ymat[, bi], col = cols[bi], pch = pchs[bi], cex = 0.9)
  }
  legend(leg_pos, legend = leg_expr, col = cols, lwd = 1.2, pch = pchs,
         pt.cex = 0.9, bty = "n", cex = 0.85)
}

plot_all_level <- function() {
  par(mfrow = c(2, 3), mar = c(4.5, 5, 3.5, 1.5))
  for (s in 1:6) {
    if (s == 4) {
      plot.new()
      title(panel_titles[s], cex.main = 1.35, font.main = 2)
      text(0.5, 0.5, "not included", cex = 1.2, col = "grey50")
    } else {
      one_panel(eps_vec, level_eps[,,s], expression(epsilon), panel_titles[s],
                "Empirical level", TRUE, alpha_level, "topleft")
    }
  }
}

plot_all_power <- function() {
  par(mfrow = c(2, 3), mar = c(4.5, 5, 3.5, 1.5))
  for (s in 1:6) {
    if (s == 4) {
      plot.new()
      title(panel_titles[s], cex.main = 1.35, font.main = 2)
      text(0.5, 0.5, "not included", cex = 1.2, col = "grey50")
    } else {
      one_panel(eps_vec, power_eps[,,s], expression(epsilon), panel_titles[s],
                "Empirical power", FALSE, leg_pos = "topleft")
    }
  }
}

if (interactive() && capabilities("X11")) {
  dev.new(width = 15, height = 9); plot_all_level()
  dev.new(width = 15, height = 9); plot_all_power()
}

pdf(paste0(base_path,"fig_H3_all_formulas_level.pdf"), width = 15, height = 9)
plot_all_level(); dev.off()

pdf(paste0(base_path,"fig_H3_all_formulas_power.pdf"), width = 15, height = 9)
plot_all_power(); dev.off()

cat("Figures saved: fig_H3_all_formulas_level.pdf, fig_H3_all_formulas_power.pdf\n")