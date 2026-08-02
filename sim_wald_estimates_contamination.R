rm(list = ls())

library(foreach)
library(doParallel)
library(optimx)
library(MASS)

source("/mnt/nfs/home/nkp117/largefiles/KP_Leandro_Maria_1/Testing/cyalt_lognormal_WMDPDE1.R")
# source("C:/Users/Kiran/Downloads/WMDPDE_CyALT_lognormal/cyalt_lognormal_WMDPDE1.R")


n_cores <- 48


theta0        <- c(5.0, -2.0, 0.5)
theta_alt_H1  <- c(5.2, -2.2, 0.48)
theta_alt_H2a <- c(5.2,  -2.0,  0.5)
theta_alt_H2b <- c(5.0,  -2.2,  0.5)
theta_alt_H3  <- c(5.2, -2.2, 0.5)

tau        <- 0.40
s_F        <- 0.40; s_1C <- 0.65; s_2C <- 1.00
stress_mat <- matrix(c(s_F, s_1C, s_F, s_2C), nrow = 2, byrow = TRUE)
ITs        <- c(15, 25, 35, 50, 65, 80)

beta_vec    <- c(0, 0.2, 0.4, 0.6, 0.8, 1.0)
beta_labels <- c("MLE", "0.2", "0.4", "0.6", "0.8", "1.0")
nb          <- length(beta_vec)

Kvec       <- c(120, 80)
eps_vec <- c(0, 0.025, 0.05, 0.075, 0.10, 0.125, 0.15)
ne      <- length(eps_vec)
cont_grp   <- 1
cont_cells <- c(1, 2, 3)

n_iter <- 1000
base_path <- "/mnt/nfs/home/nkp117/largefiles/KP_Leandro_Maria_1/Testing/"


# =============================================================================
# ONE REPLICATION — record estimates only, no testing
# =============================================================================

one_rep_estimates <- function(rep_id, eps) {

  dat_H1  <- gen_data_multicell(theta_alt_H1,  Kvec, stress_mat, tau, ITs,
                                eps, cont_grp, cont_cells, NULL)
  dat_H2a <- gen_data_multicell(theta_alt_H2a, Kvec, stress_mat, tau, ITs,
                                eps, cont_grp, cont_cells, NULL)
  dat_H2b <- gen_data_multicell(theta_alt_H2b, Kvec, stress_mat, tau, ITs,
                                eps, cont_grp, cont_cells, NULL)
  dat_H3  <- gen_data_multicell(theta_alt_H3,  Kvec, stress_mat, tau, ITs,
                                eps, cont_grp, cont_cells, NULL)

  if (!is_valid_dataset(dat_H1)) return(NULL)

  # 4 hypotheses x 3 parameters x nb betas
  out <- numeric(nb * 12)

  for (bi in seq_along(beta_vec)) {

    b   <- beta_vec[bi]
    key <- if (b == 0) "MLE" else paste0("MDPDE_", b)

    eH1  <- fit_all_betas(dat_H1,  Kvec, stress_mat, tau, ITs,
                          beta_vec = b, init = theta_alt_H1)[[key]]
    eH2a <- fit_all_betas(dat_H2a, Kvec, stress_mat, tau, ITs,
                          beta_vec = b, init = theta_alt_H2a)[[key]]
    eH2b <- fit_all_betas(dat_H2b, Kvec, stress_mat, tau, ITs,
                          beta_vec = b, init = theta_alt_H2b)[[key]]
    eH3  <- fit_all_betas(dat_H3,  Kvec, stress_mat, tau, ITs,
                          beta_vec = b, init = theta_alt_H3)[[key]]

    idx <- (bi - 1) * 12
    out[idx+1]  <- if (any(is.na(eH1)))  NA else eH1[1]
    out[idx+2]  <- if (any(is.na(eH1)))  NA else eH1[2]
    out[idx+3]  <- if (any(is.na(eH1)))  NA else eH1[3]
    out[idx+4]  <- if (any(is.na(eH2a))) NA else eH2a[1]
    out[idx+5]  <- if (any(is.na(eH2a))) NA else eH2a[2]
    out[idx+6]  <- if (any(is.na(eH2a))) NA else eH2a[3]
    out[idx+7]  <- if (any(is.na(eH2b))) NA else eH2b[1]
    out[idx+8]  <- if (any(is.na(eH2b))) NA else eH2b[2]
    out[idx+9]  <- if (any(is.na(eH2b))) NA else eH2b[3]
    out[idx+10] <- if (any(is.na(eH3)))  NA else eH3[1]
    out[idx+11] <- if (any(is.na(eH3)))  NA else eH3[2]
    out[idx+12] <- if (any(is.na(eH3)))  NA else eH3[3]
  }
  out
}


# =============================================================================
# PARALLEL SETUP
# =============================================================================


cl      <- makeCluster(n_cores)
registerDoParallel(cl)
cat(sprintf("Using %d cores\n", n_cores))

clusterExport(cl, varlist = c(
  "theta_alt_H1","theta_alt_H2a","theta_alt_H2b","theta_alt_H3",
  "tau","stress_mat","ITs","beta_vec","nb","Kvec",
  "cont_grp","cont_cells",
  "one_rep_estimates","gen_data_multicell","is_valid_dataset",
  "fit_all_betas","fit_mdpde","H_beta_objective",
  "p_i_theta","W_i_theta_matrix","mu_fun","eta_fun",
  "sbar_fun","aij_fun","interval_prob","surv_prob"
))
clusterEvalQ(cl, { library(optimx); library(MASS) })


# =============================================================================
# RUN — mean estimate per hypothesis, per parameter, per beta, per eps
# =============================================================================

# dims: eps x beta x hypothesis x parameter
est_mean <- array(NA, dim = c(ne, nb, 4, 3),
                  dimnames = list(as.character(eps_vec), beta_labels,
                                  c("H1","H2a","H2b","H3"),
                                  c("alpha0","alpha1","sigma")))

for (ei in seq_along(eps_vec)) {

  eps <- eps_vec[ei]
  cat(sprintf("eps = %.2f ...\n", eps))

  reps <- foreach(
    rep            = seq_len(n_iter),
    .combine       = "rbind",
    .packages      = c("optimx","MASS"),
    .errorhandling = "remove"
  ) %dopar% {
    suppressWarnings(one_rep_estimates(rep, eps))
  }

  if (!is.null(reps) && nrow(reps) > 0) {
    for (bi in seq_along(beta_vec)) {
      idx <- (bi - 1) * 12
      est_mean[ei, bi, "H1",  ] <- colMeans(reps[, idx+1:3,  drop=FALSE], na.rm = TRUE)
      est_mean[ei, bi, "H2a", ] <- colMeans(reps[, idx+4:6,  drop=FALSE], na.rm = TRUE)
      est_mean[ei, bi, "H2b", ] <- colMeans(reps[, idx+7:9,  drop=FALSE], na.rm = TRUE)
      est_mean[ei, bi, "H3",  ] <- colMeans(reps[, idx+10:12,drop=FALSE], na.rm = TRUE)
    }
  }
}

stopCluster(cl)

save(est_mean, eps_vec, beta_vec, beta_labels,
     theta0, theta_alt_H1, theta_alt_H2a, theta_alt_H2b, theta_alt_H3,
     file = paste0(base_path, "estimate_bias_check.RData"))
cat("Saved estimate_bias_check.RData\n")


# =============================================================================
# PLOT — for each hypothesis, show the relevant parameter's mean estimate
# vs contamination, relative to theta0 and theta_alt
# =============================================================================

cols <- c("#000000","#0072B2","#009E73","#D55E00","#7B2FBE","#CC79A7")
pchs <- c(16, 1, 2, 5, 6, 0)

one_bias_panel <- function(eps_vec, ymat, null_val, alt_val, ylab_txt, main_lab) {
    plot(NA, xlim = range(eps_vec),
         ylim = range(c(ymat, null_val, alt_val), na.rm = TRUE),
         xlab = expression(epsilon), ylab = "",
         main = main_lab, cex.main = 2.4, cex.lab = 2.0, cex.axis = 1.8,
         font.main = 2, font.lab = 1, las = 1)
    title(ylab = ylab_txt, line = 4.5, cex.lab = 2.0, font.lab = 1)
    abline(h = null_val, lty = 2, col = "grey30", lwd = 1.3)
    abline(h = alt_val,  lty = 3, col = "red",    lwd = 1.3)
    for (bi in seq_along(beta_vec)) {
        lines(eps_vec, ymat[, bi], type = "b",
              col = cols[bi], pch = pchs[bi], cex = 1.2, lwd = 1.3)
    }
}


# =============================================================================
# FULL PLOT — one panel per parameter tested in each hypothesis
# =============================================================================

plot_all_bias <- function() {
    
    par(mfrow = c(2, 4), mar = c(6, 7.0, 5, 2), font.axis = 1)
    
    # H1: alpha0, alpha1, sigma
    one_bias_panel(eps_vec, est_mean[, , "H1", "alpha0"],
                   theta0[1], theta_alt_H1[1],
                   expression(hat(alpha)[0]), "H1")
    one_bias_panel(eps_vec, est_mean[, , "H1", "alpha1"],
                   theta0[2], theta_alt_H1[2],
                   expression(hat(alpha)[1]), "H1")
    one_bias_panel(eps_vec, est_mean[, , "H1", "sigma"],
                   theta0[3], theta_alt_H1[3],
                   expression(hat(sigma)), "H1")
    
    # H2a: alpha0 only
    one_bias_panel(eps_vec, est_mean[, , "H2a", "alpha0"],
                   theta0[1], theta_alt_H2a[1],
                   expression(hat(alpha)[0]), "H2a")
    
    # H2b: alpha1 only
    one_bias_panel(eps_vec, est_mean[, , "H2b", "alpha1"],
                   theta0[2], theta_alt_H2b[2],
                   expression(hat(alpha)[1]), "H2b")
    
    # H3: alpha0, alpha1
    one_bias_panel(eps_vec, est_mean[, , "H3", "alpha0"],
                   theta0[1], theta_alt_H3[1],
                   expression(hat(alpha)[0]), "H3")
    one_bias_panel(eps_vec, est_mean[, , "H3", "alpha1"],
                   theta0[2], theta_alt_H3[2],
                   expression(hat(alpha)[1]), "H3")
    
    # legend in the last empty panel
    plot.new()
    legend("center",
           legend = c("null value", "alt value", beta_labels),
           col    = c("grey30", "red", cols),
           lty    = c(2, 3, rep(1, 6)),
           pch    = c(NA, NA, pchs),
           lwd    = 1.3,
           cex    = 1.6,
           bty    = "n")
}

pdf(paste0(base_path, "fig_estimate_bias_full.pdf"), width = 18, height = 9)
plot_all_bias()
dev.off()

cat("Saved fig_estimate_bias_full.pdf\n")



