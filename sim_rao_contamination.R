rm(list = ls())

library(foreach)
library(doParallel)
library(optimx)
library(MASS)

source("/mnt/nfs/home/nkp117/largefiles/KP_Leandro_Maria_1/Testing/cyalt_lognormal_WMDPDE1.R")
#source("C:/Users/Kiran/Downloads/WMDPDE_CyALT_lognormal/cyalt_lognormal_WMDPDE1.R")


n_cores <- 48


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

n_iter     <- 1000
Kvec       <- c(120, 80)
eps_vec <- seq(0, 0.10, length.out = 7)
ne         <- length(eps_vec)
cont_grp   <- 1
cont_cells <- c(1, 2, 3)


base_path <- "/mnt/nfs/home/nkp117/largefiles/KP_Leandro_Maria_1/Testing/"
#base_path <- "C:/Users/Kiran/Downloads/WMDPDE_CyALT_lognormal/"

# =============================================================================
# RAO-TYPE TEST STATISTICS
# =============================================================================
# All evaluated at theta0 — simple null hypotheses. U_beta_vec and
# K_beta_mat are sourced from the source file. No numerical
# optimisation is required for the Rao-type test.

rao_H1 <- function(U, Kmat, K) {
  Kinv <- tryCatch(solve(Kmat), error = function(e) matrix(NA, 3, 3))
  if (any(is.na(Kinv))) return(NA)
  K * drop(t(U) %*% Kinv %*% U)
}

rao_H2a <- function(U, Kmat, K) {
  k11 <- Kmat[1, 1]
  if (is.na(k11) || k11 <= 0) return(NA)
  K * U[1]^2 / k11
}

rao_H2b <- function(U, Kmat, K) {
  k22 <- Kmat[2, 2]
  if (is.na(k22) || k22 <= 0) return(NA)
  K * U[2]^2 / k22
}

rao_H3 <- function(U, Kmat, K) {
  Ksub  <- Kmat[1:2, 1:2]
  Kinv2 <- tryCatch(solve(Ksub), error = function(e) matrix(NA, 2, 2))
  if (any(is.na(Kinv2))) return(NA)
  U12 <- U[1:2]
  K * drop(t(U12) %*% Kinv2 %*% U12)
}


# =============================================================================
# SINGLE REPLICATION — varying contamination, fixed sample size
# =============================================================================

one_rep_rao_eps <- function(rep_id, eps) {

  # generate data under null — used for empirical level
  dat0    <- gen_data_multicell(theta0,        Kvec, stress_mat, tau, ITs,
                                eps, cont_grp, cont_cells, NULL)

  # generate data under each alternative — used for empirical power
  dat_H1  <- gen_data_multicell(theta_alt_H1,  Kvec, stress_mat, tau, ITs,
                                eps, cont_grp, cont_cells, NULL)
  dat_H2a <- gen_data_multicell(theta_alt_H2a, Kvec, stress_mat, tau, ITs,
                                eps, cont_grp, cont_cells, NULL)
  dat_H2b <- gen_data_multicell(theta_alt_H2b, Kvec, stress_mat, tau, ITs,
                                eps, cont_grp, cont_cells, NULL)
  dat_H3  <- gen_data_multicell(theta_alt_H3,  Kvec, stress_mat, tau, ITs,
                                eps, cont_grp, cont_cells, NULL)

  if (!is_valid_dataset(dat0)) return(NULL)

  K   <- sum(Kvec)
  out <- numeric(nb * 8)

  for (bi in seq_along(beta_vec)) {

    b    <- beta_vec[bi]
    Kmat <- K_beta_mat(theta0, Kvec, stress_mat, tau, ITs, b)

    U0   <- U_beta_vec(theta0, dat0,    Kvec, stress_mat, tau, ITs, b)
    UH1  <- U_beta_vec(theta0, dat_H1,  Kvec, stress_mat, tau, ITs, b)
    UH2a <- U_beta_vec(theta0, dat_H2a, Kvec, stress_mat, tau, ITs, b)
    UH2b <- U_beta_vec(theta0, dat_H2b, Kvec, stress_mat, tau, ITs, b)
    UH3  <- U_beta_vec(theta0, dat_H3,  Kvec, stress_mat, tau, ITs, b)

    lev_H1  <- as.numeric(rao_H1(U0,  Kmat, K) > cv_H1)
    lev_H2a <- as.numeric(rao_H2a(U0, Kmat, K) > cv_H2)
    lev_H2b <- as.numeric(rao_H2b(U0, Kmat, K) > cv_H2)
    lev_H3  <- as.numeric(rao_H3(U0,  Kmat, K) > cv_H3)

    pow_H1  <- as.numeric(rao_H1(UH1,  Kmat, K) > cv_H1)
    pow_H2a <- as.numeric(rao_H2a(UH2a,Kmat, K) > cv_H2)
    pow_H2b <- as.numeric(rao_H2b(UH2b,Kmat, K) > cv_H2)
    pow_H3  <- as.numeric(rao_H3(UH3,  Kmat, K) > cv_H3)

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
  "rao_H1","rao_H2a","rao_H2b","rao_H3",
  "one_rep_rao_eps","gen_data_multicell","is_valid_dataset",
  "U_beta_vec","K_beta_mat",
  "p_i_theta","W_i_theta_matrix","mu_fun","eta_fun",
  "sbar_fun","aij_fun","interval_prob","surv_prob"
))
clusterEvalQ(cl, { library(optimx); library(MASS) })


# =============================================================================
# RUN SIMULATION
# =============================================================================

level_rao_eps <- array(NA, dim = c(ne, nb, 4),
                       dimnames = list(as.character(eps_vec),
                                       beta_labels, c("H1","H2a","H2b","H3")))
power_rao_eps <- array(NA, dim = c(ne, nb, 4),
                       dimnames = list(as.character(eps_vec),
                                       beta_labels, c("H1","H2a","H2b","H3")))

t_start <- proc.time()
cat("Running Rao-type simulation: varying contamination, K = 200\n\n")

for (ei in seq_along(eps_vec)) {

  eps <- eps_vec[ei]
  cat(sprintf("  eps = %.3f ...\n", eps))

  reps <- foreach(
    rep            = seq_len(n_iter),
    .combine       = "rbind",
    .packages      = c("optimx","MASS"),
    .errorhandling = "remove"
  ) %dopar% {
    suppressWarnings(one_rep_rao_eps(rep, eps))
  }

  if (!is.null(reps) && nrow(reps) > 0) {
    reps <- reps[complete.cases(reps), , drop = FALSE]
    for (bi in seq_along(beta_vec)) {
      idx <- (bi - 1) * 8
      for (h in 1:4) {
        level_rao_eps[ei, bi, h] <- mean(reps[, idx + h],     na.rm = TRUE)
        power_rao_eps[ei, bi, h] <- mean(reps[, idx + 4 + h], na.rm = TRUE)
      }
    }
  }

  # save after each contamination level in case of crash
  save(level_rao_eps, power_rao_eps, eps_vec, beta_vec, beta_labels, Kvec,
       theta0, theta_alt_H1, theta_alt_H2a, theta_alt_H2b, theta_alt_H3,
       alpha_level,
       file = paste0(base_path, "sim_rao_contamination.RData"))
  cat(sprintf("    saved after eps = %.3f\n", eps))
}

cat(sprintf("\nDone in %.1f minutes\n",
            (proc.time() - t_start)["elapsed"] / 60))

# final clean save for plotting on laptop
save(level_rao_eps, power_rao_eps, eps_vec, beta_vec, beta_labels, Kvec,
     theta0, theta_alt_H1, theta_alt_H2a, theta_alt_H2b, theta_alt_H3,
     alpha_level,
     file = paste0(base_path, "results_rao_contamination.RData"))
cat("Final results saved.\n")

stopCluster(cl)


# =============================================================================
# PLOTTING
# =============================================================================

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

one_panel <- function(xv, ymat, xlab, main_lab, ylab_txt,
                      ref_line = FALSE, ref_val = 0.05,
                      leg_pos  = "topright") {
  ymin <- max(0, min(ymat, na.rm = TRUE) - 0.05)
  ymax <- min(1, max(ymat, na.rm = TRUE) + 0.05)
  plot(NA, xlim = range(xv), ylim = c(ymin, ymax),
       xlab = xlab, ylab = "",
       main = main_lab,
       cex.main = 2.3, cex.lab = 2.2, cex.axis = 2.0,
       font.main = 2, font.lab = 1, las = 1, bty = "o")
  title(ylab = ylab_txt, line = 5.2, cex.lab = 2.2, font.lab = 1)
  grid(nx = NA, ny = 5, col = "grey93", lty = 1, lwd = 0.8)
  if (ref_line)
    abline(h = ref_val, lty = 2, col = "grey40", lwd = 1.2)
  box()
  for (bi in seq_along(beta_vec)) {
    lines(xv, ymat[, bi],
          col  = cols[bi],
          lty  = ltys[bi],
          lwd  = lwds[bi],
          type = "b")
    points(xv, ymat[, bi],
           col = cols[bi],
           pch = pchs[bi],
           cex = 1.3,
           bg  = cols[bi])
  }
  legend(leg_pos, legend = leg_expr,
         col = cols, lty = ltys, lwd = lwds,
         pch = pchs, pt.bg = cols, pt.cex = 1.2,
         bty = "n", cex = 1.7, y.intersp = 1.0)
}

plot_level_rao_eps <- function() {
  par(mfrow = c(1,4), mar = c(7, 8.5, 5.5, 2), font.axis = 1)
  for (h in 1:4)
    one_panel(eps_vec, level_rao_eps[,,h],
              xlab     = expression(epsilon),
              main_lab = hyp_labels[[h]],
              ylab_txt = "Empirical level",
              ref_line = TRUE, ref_val = alpha_level,
              leg_pos  = "topleft")
}

plot_power_rao_eps <- function() {
  par(mfrow = c(1,4), mar = c(7, 8.5, 5.5, 2), font.axis = 1)

  one_panel(eps_vec, power_rao_eps[,,1],
            xlab = expression(epsilon), main_lab = hyp_labels[[1]],
            ylab_txt = "Empirical power", ref_line = FALSE,
            leg_pos = "topleft")

  one_panel(eps_vec, power_rao_eps[,,2],
            xlab = expression(epsilon), main_lab = hyp_labels[[2]],
            ylab_txt = "Empirical power", ref_line = FALSE,
            leg_pos = "bottomleft")

  one_panel(eps_vec, power_rao_eps[,,3],
            xlab = expression(epsilon), main_lab = hyp_labels[[3]],
            ylab_txt = "Empirical power", ref_line = FALSE,
            leg_pos = "bottomright")

  one_panel(eps_vec, power_rao_eps[,,4],
            xlab = expression(epsilon), main_lab = hyp_labels[[4]],
            ylab_txt = "Empirical power", ref_line = FALSE,
            leg_pos = "topright")
}

pdf(paste0(base_path,"fig_level_rao_contamination.pdf"), width = 18, height = 5)
plot_level_rao_eps(); dev.off()

pdf(paste0(base_path,"fig_power_rao_contamination.pdf"), width = 18, height = 5)
plot_power_rao_eps(); dev.off()

cat("Figures saved.\n")
