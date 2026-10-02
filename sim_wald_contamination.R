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

alpha_level <- 0.05
cv_H1 <- qchisq(1 - alpha_level, df = 3)
cv_H2 <- qchisq(1 - alpha_level, df = 1)
cv_H3 <- qchisq(1 - alpha_level, df = 2)

n_iter     <- 1000
Kvec       <- c(120, 80)
eps_vec <- seq(0, 0.10, length.out = 7)
ne      <- length(eps_vec)
cont_grp   <- 1
cont_cells <- c(1, 2, 3)

base_path <- "/mnt/nfs/home/nkp117/largefiles/KP_Leandro_Maria_1/Testing/arxiv/"
#base_path <- "C:/Users/Kiran/Downloads/WMDPDE_CyALT_lognormal/"

############# WALD-TYPE TEST STATISTICS #############

# H1 unchanged — full 3-parameter joint estimation

wald_H1 <- function(theta_hat, theta0, Kvec, stress_mat, tau, ITs, beta) {
  Sig  <- Sigma_hat(theta0, Kvec, stress_mat, tau, ITs, beta)
  if (any(is.na(Sig))) return(NA)
  Sinv <- tryCatch(solve(Sig), error = function(e) matrix(NA, 3, 3))
  if (any(is.na(Sinv))) return(NA)
  diff <- theta_hat - theta0
  drop(t(diff) %*% Sinv %*% diff)
}

# H2a, H2b, H3 — FINAL AGREED FORMULA (Maria/Leandro, confirmed)

wald_H2a <- function(alpha0_hat, theta0, Kvec, stress_mat, tau, ITs, beta) {
  K  <- sum(Kvec)
  J1 <- J_beta_mat(theta0, Kvec, stress_mat, tau, ITs, beta)[1, 1]
  K1 <- K_beta_mat(theta0, Kvec, stress_mat, tau, ITs, beta)[1, 1]
  if (is.na(K1) || K1 <= 0) return(NA)
  Sigma2a <- K1 / (K * J1^2)
  (alpha0_hat - theta0[1])^2 / Sigma2a
}

wald_H2b <- function(alpha1_hat, theta0, Kvec, stress_mat, tau, ITs, beta) {
  K  <- sum(Kvec)
  J2 <- J_beta_mat(theta0, Kvec, stress_mat, tau, ITs, beta)[2, 2]
  K2 <- K_beta_mat(theta0, Kvec, stress_mat, tau, ITs, beta)[2, 2]
  if (is.na(K2) || K2 <= 0) return(NA)
  Sigma2b <- K2 / (K * J2^2)
  (alpha1_hat - theta0[2])^2 / Sigma2b
}

wald_H3 <- function(alpha01_hat, theta0, Kvec, stress_mat, tau, ITs, beta) {
  K    <- sum(Kvec)
  Jmat <- J_beta_mat(theta0, Kvec, stress_mat, tau, ITs, beta)
  Kmat <- K_beta_mat(theta0, Kvec, stress_mat, tau, ITs, beta)
  J12  <- Jmat[1:2, 1:2]
  K12  <- Kmat[1:2, 1:2]
  Jinv12 <- tryCatch(solve(J12), error = function(e) matrix(NA, 2, 2))
  if (any(is.na(Jinv12))) return(NA)
  Sig3  <- (1 / K) * Jinv12 %*% K12 %*% Jinv12
  Sinv  <- tryCatch(solve(Sig3), error = function(e) matrix(NA, 2, 2))
  if (any(is.na(Sinv))) return(NA)
  diff  <- alpha01_hat - theta0[1:2]
  drop(t(diff) %*% Sinv %*% diff)
}

############# 1D / 2D PLUG-IN FITTING FOR H2a, H2b, H3 #############

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

################## SINGLE REPLICATION — varying contamination, fixed sample size ###########3

one_rep_eps <- function(rep_id, eps) {

  dat0    <- gen_data_multicell(theta0,        Kvec, stress_mat, tau, ITs,
                                eps, cont_grp, cont_cells, NULL)
  dat_H1  <- gen_data_multicell(theta_alt_H1,  Kvec, stress_mat, tau, ITs,
                                eps, cont_grp, cont_cells, NULL)
  dat_H2a <- gen_data_multicell(theta_alt_H2a, Kvec, stress_mat, tau, ITs,
                                eps, cont_grp, cont_cells, NULL)
  dat_H2b <- gen_data_multicell(theta_alt_H2b, Kvec, stress_mat, tau, ITs,
                                eps, cont_grp, cont_cells, NULL)
  dat_H3  <- gen_data_multicell(theta_alt_H3,  Kvec, stress_mat, tau, ITs,
                                eps, cont_grp, cont_cells, NULL)

  if (!is_valid_dataset(dat0)) return(NULL)

  out <- numeric(nb * 8)

  for (bi in seq_along(beta_vec)) {

    b   <- beta_vec[bi]
    key <- if (b == 0) "MLE" else paste0("MDPDE_", b)

    # full 3-parameter joint fit — H1 only
    est0   <- fit_all_betas(dat0,   Kvec, stress_mat, tau, ITs,
                            beta_vec = b, init = theta0)[[key]]
    estH1  <- fit_all_betas(dat_H1, Kvec, stress_mat, tau, ITs,
                            beta_vec = b, init = theta_alt_H1)[[key]]

    # 1D plug-in fits — H2a, H2b
    alpha0_est0   <- fit_alpha0_1D(dat0,    Kvec, stress_mat, tau, ITs, b,
                                   alpha1_known = theta0[2], sigma_known = theta0[3],
                                   init_alpha0  = theta0[1])
    alpha0_estH2a <- fit_alpha0_1D(dat_H2a, Kvec, stress_mat, tau, ITs, b,
                                   alpha1_known = theta0[2], sigma_known = theta0[3],
                                   init_alpha0  = theta_alt_H2a[1])

    alpha1_est0   <- fit_alpha1_1D(dat0,    Kvec, stress_mat, tau, ITs, b,
                                   alpha0_known = theta0[1], sigma_known = theta0[3],
                                   init_alpha1  = theta0[2])
    alpha1_estH2b <- fit_alpha1_1D(dat_H2b, Kvec, stress_mat, tau, ITs, b,
                                   alpha0_known = theta0[1], sigma_known = theta0[3],
                                   init_alpha1  = theta_alt_H2b[2])

    # 2D plug-in fit — H3
    alpha01_est0  <- fit_alpha0_alpha1_2D(dat0,   Kvec, stress_mat, tau, ITs, b,
                                          sigma_known        = theta0[3],
                                          init_alpha0_alpha1 = theta0[1:2])
    alpha01_estH3 <- fit_alpha0_alpha1_2D(dat_H3, Kvec, stress_mat, tau, ITs, b,
                                          sigma_known        = theta0[3],
                                          init_alpha0_alpha1 = theta_alt_H3[1:2])

    # empirical level
    lev_H1  <- if (any(is.na(est0)))        NA else
      as.numeric(wald_H1(est0,  theta0,Kvec,stress_mat,tau,ITs,b) > cv_H1)
    lev_H2a <- if (is.na(alpha0_est0))       NA else
      as.numeric(wald_H2a(alpha0_est0,  theta0,Kvec,stress_mat,tau,ITs,b) > cv_H2)
    lev_H2b <- if (is.na(alpha1_est0))       NA else
      as.numeric(wald_H2b(alpha1_est0,  theta0,Kvec,stress_mat,tau,ITs,b) > cv_H2)
    lev_H3  <- if (any(is.na(alpha01_est0))) NA else
      as.numeric(wald_H3(alpha01_est0,  theta0,Kvec,stress_mat,tau,ITs,b) > cv_H3)

    # empirical power
    pow_H1  <- if (any(is.na(estH1)))         NA else
      as.numeric(wald_H1(estH1,  theta0,Kvec,stress_mat,tau,ITs,b) > cv_H1)
    pow_H2a <- if (is.na(alpha0_estH2a))       NA else
      as.numeric(wald_H2a(alpha0_estH2a,theta0,Kvec,stress_mat,tau,ITs,b) > cv_H2)
    pow_H2b <- if (is.na(alpha1_estH2b))       NA else
      as.numeric(wald_H2b(alpha1_estH2b,theta0,Kvec,stress_mat,tau,ITs,b) > cv_H2)
    pow_H3  <- if (any(is.na(alpha01_estH3)))  NA else
      as.numeric(wald_H3(alpha01_estH3, theta0,Kvec,stress_mat,tau,ITs,b) > cv_H3)

    idx <- (bi - 1) * 8
    out[idx+1] <- lev_H1;  out[idx+2] <- lev_H2a
    out[idx+3] <- lev_H2b; out[idx+4] <- lev_H3
    out[idx+5] <- pow_H1;  out[idx+6] <- pow_H2a
    out[idx+7] <- pow_H2b; out[idx+8] <- pow_H3
  }
  out
}

############# PARALLEL SETUP ######################

cl <- makeCluster(n_cores)
registerDoParallel(cl)
cat(sprintf("Using %d cores\n", n_cores))

clusterExport(cl, varlist = c(
  "theta0","theta_alt_H1","theta_alt_H2a","theta_alt_H2b","theta_alt_H3",
  "tau","stress_mat","ITs","beta_vec","nb","Kvec",
  "cv_H1","cv_H2","cv_H3","cont_grp","cont_cells",
  "one_rep_eps","gen_data_multicell","is_valid_dataset",
  "fit_all_betas","fit_mdpde","H_beta_objective",
  "p_i_theta","W_i_theta_matrix","mu_fun","eta_fun",
  "sbar_fun","aij_fun","interval_prob","surv_prob",
  "J_beta_mat","K_beta_mat",
  "Sigma_hat","wald_H1","wald_H2a","wald_H2b","wald_H3",
  "fit_alpha0_1D","fit_alpha1_1D","fit_alpha0_alpha1_2D"
))
clusterEvalQ(cl, { library(optimx); library(MASS) })

# RUN SIMULATION

level_eps <- array(NA, dim = c(ne, nb, 4),
                   dimnames = list(as.character(round(eps_vec, 4)),
                                   beta_labels, c("H1","H2a","H2b","H3")))
power_eps <- array(NA, dim = c(ne, nb, 4),
                   dimnames = list(as.character(round(eps_vec, 4)),
                                   beta_labels, c("H1","H2a","H2b","H3")))

t_start <- proc.time()
cat("Running simulation: varying contamination, K = 200\n\n")

for (ei in seq_along(eps_vec)) {

  eps <- eps_vec[ei]
  cat(sprintf("  eps = %.4f ...\n", eps))

  reps <- foreach(
    rep            = seq_len(n_iter),
    .combine       = "rbind",
    .packages      = c("optimx","MASS"),
    .errorhandling = "remove"
  ) %dopar% {
    suppressWarnings(one_rep_eps(rep, eps))
  }

  if (!is.null(reps) && nrow(reps) > 0) {
    reps <- reps[complete.cases(reps), , drop = FALSE]
    for (bi in seq_along(beta_vec)) {
      idx <- (bi - 1) * 8
      for (h in 1:4) {
        level_eps[ei, bi, h] <- mean(reps[, idx + h],     na.rm = TRUE)
        power_eps[ei, bi, h] <- mean(reps[, idx + 4 + h], na.rm = TRUE)
      }
    }
  }

  save(level_eps, power_eps, eps_vec, beta_vec, beta_labels, Kvec,
       theta0, theta_alt_H1, theta_alt_H2a, theta_alt_H2b, theta_alt_H3,
       alpha_level,
       file = paste0(base_path, "sim_wald_contamination.RData"))
  cat(sprintf("    saved after eps = %.4f\n", eps))
}

save(level_eps, power_eps, eps_vec, beta_vec, beta_labels, Kvec,
     theta0, theta_alt_H1, theta_alt_H2a, theta_alt_H2b, theta_alt_H3,
     alpha_level,
     file = paste0(base_path, "results_wald_contamination.RData"))
cat("Final results saved: results_wald_contamination.RData\n")

cat(sprintf("\nDone in %.1f minutes\n",
            (proc.time() - t_start)["elapsed"] / 60))
stopCluster(cl)

#################### Plots of Empirical power and type-I error/levels ####################

cols <- c("#000000","#0072B2","#009E73","#D55E00","#7B2FBE","#CC79A7")
ltys <- rep(1, 6)
lwds <- rep(1.3, 6)
pchs <- c(16, 1, 2, 5, 6, 0)

leg_expr <- c(
  expression(MLE ~ (beta==0)),
  expression(beta==0.2),
  expression(beta==0.4),
  expression(beta==0.6),
  expression(beta==0.8),
  expression(beta=="1.0")
)

hyp_labels <- list(
  expression(H[0]:~alpha[0]==alpha[list(0,0)]~","~alpha[1]==alpha[list(1,0)]~","~sigma==sigma[0]),
  expression(H[0]:~alpha[0]==alpha[list(0,0)]~(alpha[1]*","~sigma~"known")),
  expression(H[0]:~alpha[1]==alpha[list(1,0)]~(alpha[0]*","~sigma~"known")),
  expression(H[0]:~alpha[0]==alpha[list(0,0)]~","~alpha[1]==alpha[list(1,0)]~(sigma~"known"))
)

# x-axis tick labels rounded to 2 decimals for display only;
# the underlying eps_vec used for simulation stays at full precision
eps_axis_labels <- sprintf("%.2f", eps_vec)

one_panel <- function(xv, ymat, xlab, main_lab, ylab_txt,
                      ref_line = FALSE, ref_val = 0.05,
                      leg_pos  = "topright", xtick_labels = NULL) {
  ymin <- max(0, min(ymat, na.rm = TRUE) - 0.05)
  ymax <- min(1, max(ymat, na.rm = TRUE) + 0.05)
  plot(NA, xlim = range(xv), ylim = c(ymin, ymax),
       xlab = xlab, ylab = "",
       main = main_lab,
       cex.main = 2.3, cex.lab = 2.2, cex.axis = 2.0,
       font.main = 2, font.lab = 1, las = 1, bty = "o", xaxt = "n")
  if (!is.null(xtick_labels)) {
    axis(1, at = xv, labels = xtick_labels, cex.axis = 1.6)
  } else {
    axis(1, cex.axis = 1.8)
  }
  title(ylab = ylab_txt, line = 5.2, cex.lab = 2.2, font.lab = 1)
  grid(nx = NA, ny = 5, col = "grey93", lty = 1, lwd = 0.8)
  if (ref_line)
    abline(h = ref_val, lty = 2, col = "grey40", lwd = 1.2)
  box()
  for (bi in seq_along(beta_vec)) {
    lines(xv, ymat[, bi],
          col = cols[bi], lty = ltys[bi], lwd = lwds[bi], type = "b")
    points(xv, ymat[, bi],
           col = cols[bi], pch = pchs[bi], cex = 1.3, bg = cols[bi])
  }
  legend(leg_pos, legend = leg_expr,
         col = cols, lty = ltys, lwd = lwds,
         pch = pchs, pt.bg = cols, pt.cex = 1.2,
         bty = "n", cex = 1.7, y.intersp = 1.0)
}

plot_level_eps <- function() {
  par(mfrow = c(1,4), mar = c(7, 8.5, 5.5, 2), font.axis = 1)
  for (h in 1:4)
    one_panel(eps_vec, level_eps[,,h],
              xlab     = expression(epsilon),
              main_lab = hyp_labels[[h]],
              ylab_txt = "Empirical level",
              ref_line = TRUE, ref_val = alpha_level,
              leg_pos  = "topleft", xtick_labels = eps_axis_labels)
}

plot_power_eps <- function() {
  par(mfrow = c(1,4), mar = c(7, 8.5, 5.5, 2), font.axis = 1)
  one_panel(eps_vec, power_eps[,,1],
            xlab = expression(epsilon), main_lab = hyp_labels[[1]],
            ylab_txt = "Empirical power", ref_line = FALSE,
            leg_pos = "topright", xtick_labels = eps_axis_labels)
  one_panel(eps_vec, power_eps[,,2],
            xlab = expression(epsilon), main_lab = hyp_labels[[2]],
            ylab_txt = "Empirical power", ref_line = FALSE,
            leg_pos = "bottomleft", xtick_labels = eps_axis_labels)
  one_panel(eps_vec, power_eps[,,3],
            xlab = expression(epsilon), main_lab = hyp_labels[[3]],
            ylab_txt = "Empirical power", ref_line = FALSE,
            leg_pos = "bottomright", xtick_labels = eps_axis_labels)
  one_panel(eps_vec, power_eps[,,4],
            xlab = expression(epsilon), main_lab = hyp_labels[[4]],
            ylab_txt = "Empirical power", ref_line = FALSE,
            leg_pos = "topright", xtick_labels = eps_axis_labels)
}

if (interactive() && capabilities("X11")) {
  dev.new(width = 18, height = 5); plot_level_eps()
  dev.new(width = 18, height = 5); plot_power_eps()
}

pdf(paste0(base_path,"fig_level_contamination.pdf"), width = 18, height = 5)
plot_level_eps(); dev.off()

pdf(paste0(base_path,"fig_power_contamination.pdf"), width = 18, height = 5)
plot_power_eps(); dev.off()

cat("Figures saved.\n")
