rm(list = ls())

library(foreach)
library(doParallel)
library(optimx)
library(MASS)

source("/mnt/nfs/home/nkp117/largefiles/KP_Leandro_Maria_1/Testing/arxiv/cyalt_lognormal_WMDPDE1.R")
#source("C:/Users/Kiran/Downloads/WMDPDE_CyALT_lognormal/cyalt_lognormal_WMDPDE1.R")

n_cores <- 48

theta0        <- c(5.0, -2.0, 0.5)

theta_alt_H1  <- c(5.2, -2.2, 0.48)
theta_alt_H2a <- c(5.2,  -2.0,  0.5)
theta_alt_H2b <- c(5.0,  -2.2,  0.5)
theta_alt_H3  <- c(5.2, -2.2, 0.5)

tau        <- 0.40
s_F        <- 0.40;  s_1C <- 0.65;  s_2C <- 1.00
s_0F       <- 0.00;  s_0C <- 0.20
stress_mat <- matrix(c(s_F, s_1C, s_F, s_2C), nrow = 2, byrow = TRUE)
ITs        <- c(15, 25, 35, 50, 65, 80)

beta_vec    <- c(0, 0.2, 0.4, 0.6, 0.8, 1.0)
beta_labels <- c("MLE", "0.2", "0.4", "0.6", "0.8", "1.0")
nb          <- length(beta_vec)

n_iter     <- 1000
Kvec       <- c(120, 80)

# contamination grid, capped at 0.10, matching the other final scripts
eps_vec <- seq(0, 0.10, length.out = 7)
ne      <- length(eps_vec)
cont_grp   <- 1
cont_cells <- c(1, 2, 3)

base_path <- "/mnt/nfs/home/nkp117/largefiles/KP_Leandro_Maria_1/Testing/arxiv/"
#base_path <- "C:/Users/Kiran/Downloads/WMDPDE_CyALT_lognormal/"


################# 1D / 2D PLUG-IN FITTING FOR H2a, H2b, H3 ###############

# Parameters stated as "known" are held fixed at their true value and
# never estimated -- a genuine reduced-dimension fit, matching the final
# agreed formula used in sim_wald_contamination.R / sim_wald_samplesize.R.

fit_alpha0_1D <- function(dat, Kvec, stress_mat, tau, ITs, beta,
                          alpha1_known, sigma_known, init_alpha0) {
  obj_1d <- function(a0) {
    theta_try <- c(a0, alpha1_known, sigma_known)
    H_beta_objective(theta_try, dat, Kvec, stress_mat, tau, ITs, beta)
  }
  fit <- tryCatch(optimize(obj_1d, interval = init_alpha0 + c(-3, 3)),
                  error = function(e) NULL)
  if (is.null(fit)) return(NA)
  fit$minimum
}

fit_alpha1_1D <- function(dat, Kvec, stress_mat, tau, ITs, beta,
                          alpha0_known, sigma_known, init_alpha1) {
  obj_1d <- function(a1) {
    theta_try <- c(alpha0_known, a1, sigma_known)
    H_beta_objective(theta_try, dat, Kvec, stress_mat, tau, ITs, beta)
  }
  fit <- tryCatch(optimize(obj_1d, interval = init_alpha1 + c(-3, 3)),
                  error = function(e) NULL)
  if (is.null(fit)) return(NA)
  fit$minimum
}

fit_alpha0_alpha1_2D <- function(dat, Kvec, stress_mat, tau, ITs, beta,
                                 sigma_known, init_alpha0_alpha1) {
  obj_2d <- function(par) {
    theta_try <- c(par[1], par[2], sigma_known)
    H_beta_objective(theta_try, dat, Kvec, stress_mat, tau, ITs, beta)
  }
  fit <- tryCatch(optim(par = init_alpha0_alpha1, fn = obj_2d,
                        method = "Nelder-Mead"),
                  error = function(e) NULL)
  if (is.null(fit) || fit$convergence != 0) return(c(NA, NA))
  fit$par
}

##################### SINGLE REPLICATION ######################

one_rep_est <- function(rep_id, eps) {

  dat_H1  <- gen_data_multicell(theta_alt_H1,  Kvec, stress_mat, tau, ITs,
                                eps, cont_grp, cont_cells, NULL)
  dat_H2a <- gen_data_multicell(theta_alt_H2a, Kvec, stress_mat, tau, ITs,
                                eps, cont_grp, cont_cells, NULL)
  dat_H2b <- gen_data_multicell(theta_alt_H2b, Kvec, stress_mat, tau, ITs,
                                eps, cont_grp, cont_cells, NULL)
  dat_H3  <- gen_data_multicell(theta_alt_H3,  Kvec, stress_mat, tau, ITs,
                                eps, cont_grp, cont_cells, NULL)

  if (!is_valid_dataset(dat_H1)) return(NULL)

  # per beta: H1 (alpha0,alpha1,sigma) + H2a (alpha0) + H2b (alpha1)
  #         + H3 (alpha0,alpha1)  =  3 + 1 + 1 + 2 = 7 values
  out <- numeric(nb * 7)

  for (bi in seq_along(beta_vec)) {

    b   <- beta_vec[bi]
    key <- if (b == 0) "MLE" else paste0("MDPDE_", b)

    # H1 -- full joint 3-parameter fit
    estH1 <- fit_all_betas(dat_H1, Kvec, stress_mat, tau, ITs,
                           beta_vec = b, init = theta_alt_H1)[[key]]

    # H2a -- 1D plug-in fit
    a0_H2a <- fit_alpha0_1D(dat_H2a, Kvec, stress_mat, tau, ITs, b,
                            alpha1_known = theta0[2], sigma_known = theta0[3],
                            init_alpha0  = theta_alt_H2a[1])

    # H2b -- 1D plug-in fit
    a1_H2b <- fit_alpha1_1D(dat_H2b, Kvec, stress_mat, tau, ITs, b,
                            alpha0_known = theta0[1], sigma_known = theta0[3],
                            init_alpha1  = theta_alt_H2b[2])

    # H3 -- 2D plug-in fit
    a01_H3 <- fit_alpha0_alpha1_2D(dat_H3, Kvec, stress_mat, tau, ITs, b,
                                   sigma_known        = theta0[3],
                                   init_alpha0_alpha1 = theta_alt_H3[1:2])

    idx <- (bi - 1) * 7
    out[idx+1] <- if (any(is.na(estH1))) NA else estH1[1]   # H1 alpha0
    out[idx+2] <- if (any(is.na(estH1))) NA else estH1[2]   # H1 alpha1
    out[idx+3] <- if (any(is.na(estH1))) NA else estH1[3]   # H1 sigma
    out[idx+4] <- a0_H2a                                     # H2a alpha0
    out[idx+5] <- a1_H2b                                     # H2b alpha1
    out[idx+6] <- a01_H3[1]                                  # H3 alpha0
    out[idx+7] <- a01_H3[2]                                  # H3 alpha1
  }
  out
}

################### PARALLEL SETUP ######################

cl <- makeCluster(n_cores)
registerDoParallel(cl)
cat(sprintf("Using %d cores\n", n_cores))

clusterExport(cl, varlist = c(
  "theta0","theta_alt_H1","theta_alt_H2a","theta_alt_H2b","theta_alt_H3",
  "tau","stress_mat","ITs","beta_vec","nb","Kvec",
  "cont_grp","cont_cells",
  "one_rep_est","gen_data_multicell","is_valid_dataset",
  "fit_all_betas","fit_mdpde","H_beta_objective",
  "p_i_theta","W_i_theta_matrix","mu_fun","eta_fun",
  "sbar_fun","aij_fun","interval_prob","surv_prob",
  "J_beta_mat","K_beta_mat","Sigma_hat",
  "fit_alpha0_1D","fit_alpha1_1D","fit_alpha0_alpha1_2D"
))
clusterEvalQ(cl, { library(optimx); library(MASS) })

################# RUN : mean estimate per (eps, beta, hypothesis, parameter) ###################

param_names <- c("alpha0","alpha1","sigma")

est_mean <- array(NA, dim = c(ne, nb, 4, 3),
                  dimnames = list(as.character(round(eps_vec, 4)),
                                  beta_labels,
                                  c("H1","H2a","H2b","H3"),
                                  param_names))

t_start <- proc.time()
cat("Running estimate-bias simulation: varying contamination, K = 200\n\n")

for (ei in seq_along(eps_vec)) {

  eps <- eps_vec[ei]
  cat(sprintf("  eps = %.4f ...\n", eps))

  reps <- foreach(
    rep            = seq_len(n_iter),
    .combine       = "rbind",
    .packages      = c("optimx","MASS"),
    .errorhandling = "remove"
  ) %dopar% {
    suppressWarnings(one_rep_est(rep, eps))
  }

  if (!is.null(reps) && nrow(reps) > 0) {
    for (bi in seq_along(beta_vec)) {
      idx <- (bi - 1) * 7
      est_mean[ei, bi, "H1",  "alpha0"] <- mean(reps[, idx+1], na.rm = TRUE)
      est_mean[ei, bi, "H1",  "alpha1"] <- mean(reps[, idx+2], na.rm = TRUE)
      est_mean[ei, bi, "H1",  "sigma"]  <- mean(reps[, idx+3], na.rm = TRUE)
      est_mean[ei, bi, "H2a", "alpha0"] <- mean(reps[, idx+4], na.rm = TRUE)
      est_mean[ei, bi, "H2b", "alpha1"] <- mean(reps[, idx+5], na.rm = TRUE)
      est_mean[ei, bi, "H3",  "alpha0"] <- mean(reps[, idx+6], na.rm = TRUE)
      est_mean[ei, bi, "H3",  "alpha1"] <- mean(reps[, idx+7], na.rm = TRUE)
    }
  }

  save(est_mean, eps_vec, beta_vec, beta_labels, Kvec,
       theta0, theta_alt_H1, theta_alt_H2a, theta_alt_H2b, theta_alt_H3,
       file = paste0(base_path, "sim_estimate_bias_contamination.RData"))
  cat(sprintf("    saved after eps = %.4f\n", eps))
}

save(est_mean, eps_vec, beta_vec, beta_labels, Kvec,
     theta0, theta_alt_H1, theta_alt_H2a, theta_alt_H2b, theta_alt_H3,
     file = paste0(base_path, "results_estimate_bias_contamination.RData"))
cat("Final results saved: results_estimate_bias_contamination.RData\n")

cat(sprintf("\nDone in %.1f minutes\n",
            (proc.time() - t_start)["elapsed"] / 60))
stopCluster(cl)

############## PLOT : for each hypothesis ###################

cols <- c("#000000","#0072B2","#009E73","#D55E00","#7B2FBE","#CC79A7")
pchs <- c(16, 1, 2, 5, 6, 0)

one_bias_panel <- function(eps_vec, ymat, null_val, alt_val, ylab_txt, main_lab) {
  plot(NA, xlim = range(eps_vec),
       ylim = range(c(ymat, null_val, alt_val), na.rm = TRUE),
       xlab = expression(epsilon), ylab = "",
       main = main_lab, cex.main = 2.4, cex.lab = 2.2, cex.axis = 2.0,
       font.main = 2, font.lab = 1, las = 1)
  title(ylab = ylab_txt, line = 5.2, cex.lab = 2.2, font.lab = 1)
  abline(h = null_val, lty = 2, col = "grey30", lwd = 1.3)
  abline(h = alt_val,  lty = 3, col = "red",    lwd = 1.3)
  for (bi in seq_along(beta_vec)) {
    lines(eps_vec, ymat[, bi], type = "b",
          col = cols[bi], pch = pchs[bi], cex = 1.3, lwd = 1.3)
  }
}

########### FULL PLOT : one panel per parameter tested in each hypothesis ############

                  
plot_all_bias <- function() {

  par(mfrow = c(2, 4), mar = c(7, 8.5, 5.5, 2), font.axis = 1)

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
         cex    = 1.7,
         bty    = "n")
}

pdf(paste0(base_path, "fig_estimate_bias_full.pdf"), width = 18, height = 9)
plot_all_bias()
dev.off()
cat("Saved fig_estimate_bias_full.pdf\n")
