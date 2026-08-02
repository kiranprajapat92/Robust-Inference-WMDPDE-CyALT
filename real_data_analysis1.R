rm(list = ls())
library(optimx)
library(MASS)

source("C:/Users/Kiran/Downloads/WMDPDE_CyALT_lognormal/cyalt_lognormal_WMDPDE1.R")

##### preliminary data

failure_times <- c(26763, 31959, 32887, 33069, 34019, 34924, 36754,
                   37054, 37385, 38045, 41033, 41755, 42333,
                   42818, 44638, 44867, 48364, 49767)

fitdistr(failure_times, "log-normal")

fit   <- fitdistr(failure_times, "log-normal")
sigma <- round(fit$estimate["sdlog"], 3)   # 0.156
p_h   <- 0.90
p_u   <- 0.001
t_c   <- 50000
tau   <- 0.5


v_0F_phys <- 0.10; v_0C_phys <- 0.25; v_hC_phys <- 3.00
s_im  <- function(v) (log(v) - log(v_0F_phys)) / (log(v_hC_phys) - log(v_0F_phys))

s_0F  <- 0.00;  s_0C <- 0.27
s_hF  <- 0.50;  s_hC <- 1.00
s_1C  <- 0.70;  s_2C <- 1.00;  s_F <- 0.30

cat(sprintf("Standardised use-condition stresses: s_0F=%.4f, s_0C=%.4f\n", s_im(v_0F_phys), s_im(v_0C_phys)))

B_fun       <- function(a1, sC, sF) tau * exp(-a1 * sC) + (1 - tau) * exp(-a1 * sF)
mu_h        <- log(t_c) - sigma * qnorm(p_h)
mu_u        <- log(t_c) - sigma * qnorm(p_u)
eq_a1       <- function(a1) (mu_h - mu_u) - (log(B_fun(a1, s_0C, s_0F)) - log(B_fun(a1, s_hC, s_hF)))
alpha1_true <- uniroot(eq_a1, c(-30, -1e-4))$root
alpha0_true <- mu_h + log(B_fun(alpha1_true, s_hC, s_hF))
theta_true  <- c(alpha0_true, alpha1_true, sigma)

cat(sprintf("theta_true: alpha0=%.4f, alpha1=%.4f, sigma=%.4f\n",
            theta_true[1], theta_true[2], theta_true[3]))
cat(sprintf("Verification: F_h(t_c)=%.4f (target 0.90) | F_u(t_c)=%.6f (target 0.001)\n",
            pnorm((log(t_c) - mu_h) / sigma), pnorm((log(t_c) - mu_u) / sigma)))

# =============================================================================
# CyALT EXPERIMENT DESIGN AND DATA SIMULATION
# =============================================================================

ITs        <- c(25000, 35000, 45000, 50000, 60000, 65000)
Kvec       <- c(140, 60)
stress_mat <- matrix(c(s_F, s_1C, s_F, s_2C), nrow = 2, byrow = TRUE)
t0         <- 75000

data_obs <- sim_cyalt_data_continuous(theta_true, Kvec, stress_mat, tau, ITs, seed = 125)

L <- length(ITs)
intervals <- c(sprintf("(0, %d]", ITs[1]),
               sapply(seq_len(L - 1), function(j) sprintf("(%d, %d]", ITs[j], ITs[j + 1])),
               "Survivors")
cat("\nObserved interval-censored counts:\n")
cat(sprintf("  %-24s  %8s  %8s\n", "Interval", "Group 1", "Group 2"))
cat(sprintf("  %s\n", strrep("-", 44)))
for (j in seq_len(L + 1))
  cat(sprintf("  %-24s  %8d  %8d\n", intervals[j], data_obs[[1]][j], data_obs[[2]][j]))
cat(sprintf("  %-24s  %8d  %8d\n", "Total", Kvec[1], Kvec[2]))

# =============================================================================
# FIT WMDPDE
# =============================================================================

beta_vec    <- c(0, 0.2, 0.4, 0.6, 0.8, 1.0)
beta_labels <- c("MLE", "0.2", "0.4", "0.6", "0.8", "1.0")
nb          <- length(beta_vec)

cat("\nFitting WMDPDE for all beta...\n")
estimates_all <- fit_all_betas(data_obs, Kvec, stress_mat, tau, ITs,
                               beta_vec = beta_vec, init = c(10, -1, 0.2))
cat("Done.\n")

# =============================================================================
# TABLE 2: PARAMETER ESTIMATES AND 95% DIRECT CIs
# =============================================================================

z95 <- qnorm(0.975)
cat("\nTable 2 — Parameter estimates and 95% CIs:\n")
cat(sprintf("  %-5s  %8s  %18s  %8s  %18s  %7s  %14s\n",
            "beta", "a0", "CI(a0)", "a1", "CI(a1)", "sigma", "CI(sigma)"))
cat(sprintf("  %s\n", strrep("-", 82)))
for (bi in seq_along(beta_vec)) {
  nm  <- if (beta_vec[bi] == 0) "MLE" else paste0("MDPDE_", beta_vec[bi])
  est <- estimates_all[[nm]]
  Sig <- Sigma_hat(est, Kvec, stress_mat, tau, ITs, beta_vec[bi])
  se  <- sqrt(diag(Sig))
  cat(sprintf("  %-5s  %8.4f  (%7.4f,%7.4f)  %8.4f  (%7.4f,%7.4f)  %7.4f  (%6.4f,%6.4f)\n",
              beta_labels[bi],
              est[1], est[1] - z95 * se[1], est[1] + z95 * se[1],
              est[2], est[2] - z95 * se[2], est[2] + z95 * se[2],
              est[3], est[3] - z95 * se[3], est[3] + z95 * se[3]))
}
cat(sprintf("  True  %8.4f  %20s  %8.4f  %20s  %7.4f\n",
            theta_true[1], "", theta_true[2], "", theta_true[3]))

# =============================================================================
# TABLE 3: LIFETIME CHARACTERISTICS — ESTIMATE + DIRECT + TRANSFORMED + BCa
# =============================================================================

B_boot <- 500

cat(sprintf("\nComputing BCa bootstrap CIs (B=%d)...\n", B_boot))
results <- list()

for (bi in seq_along(beta_vec)) {
  nm   <- if (beta_vec[bi] == 0) "MLE" else paste0("MDPDE_", beta_vec[bi])
  est  <- estimates_all[[nm]]
  beta <- beta_vec[bi]
  
  cat(sprintf("  beta = %s\n", beta_labels[bi]))
  
  med_char  <- function(th) quantile_use(th, s_0C, s_0F, tau, q = 0.5)
  mttf_char <- function(th) mttf_use(th, s_0C, s_0F, tau)
  rel_char  <- function(th) reliability_use(th, s_0C, s_0F, tau, t0 = t0)
  
  med_se    <- function(th) quantile_use_se(th, Kvec, stress_mat, tau, ITs, beta, s_0C, s_0F, 0.5)
  mttf_se   <- function(th) mttf_use_se(th, Kvec, stress_mat, tau, ITs, beta, s_0C, s_0F)
  rel_se    <- function(th) reliability_use_se(th, Kvec, stress_mat, tau, ITs, beta, s_0C, s_0F, t0)
  
  a0_char <- function(th) th[1]
  a1_char <- function(th) th[2]
  sg_char <- function(th) th[3]
  
  a0_se   <- function(th) sqrt(Sigma_hat(th, Kvec, stress_mat, tau, ITs, beta)[1, 1])
  a1_se   <- function(th) sqrt(Sigma_hat(th, Kvec, stress_mat, tau, ITs, beta)[2, 2])
  sg_se   <- function(th) sqrt(Sigma_hat(th, Kvec, stress_mat, tau, ITs, beta)[3, 3])
  
  results[[beta_labels[bi]]] <- list(
    alpha0 = ci_all_three(est, data_obs, Kvec, stress_mat, tau, ITs,
                          beta, a0_char, a0_se, type = "none", B = B_boot, seed = bi * 10),
    alpha1 = ci_all_three(est, data_obs, Kvec, stress_mat, tau, ITs,
                          beta, a1_char, a1_se, type = "none", B = B_boot, seed = bi * 10 + 1),
    sigma  = ci_all_three(est, data_obs, Kvec, stress_mat, tau, ITs,
                          beta, sg_char, sg_se,  type = "log",  B = B_boot, seed = bi * 10 + 2),
    median = ci_all_three(est, data_obs, Kvec, stress_mat, tau, ITs,
                          beta, med_char, med_se,   type = "log",   B = B_boot, seed = bi * 100),
    mttf   = ci_all_three(est, data_obs, Kvec, stress_mat, tau, ITs,
                          beta, mttf_char, mttf_se, type = "log",   B = B_boot, seed = bi * 100 + 1),
    rel    = ci_all_three(est, data_obs, Kvec, stress_mat, tau, ITs,
                          beta, rel_char, rel_se,   type = "logit", B = B_boot, seed = bi * 100 + 2)
  )
}

med_true  <- quantile_use(theta_true, s_0C, s_0F, tau, 0.5)
mttf_true <- mttf_use(theta_true, s_0C, s_0F, tau)
rel_true  <- reliability_use(theta_true, s_0C, s_0F, tau, t0 = t0)

print_char_table <- function(title, key, true_val, fmt = "%10.0f") {
  cat(sprintf("\n%s  [true = %s]:\n", title, sprintf(fmt, true_val)))
  cat(sprintf("  %-5s  %10s  %24s  %24s  %24s\n",
              "beta", "Estimate", "Direct 95% CI", "Transformed 95% CI", "BCa 95% CI"))
  cat(sprintf("  %s\n", strrep("-", 92)))
  for (lbl in beta_labels) {
    r <- results[[lbl]][[key]]
    cat(sprintf(paste0("  %-5s  ", fmt, "  (", fmt, ", ", fmt, ")  (",
                       fmt, ", ", fmt, ")  (", fmt, ", ", fmt, ")\n"),
                lbl,
                r$estimate,
                r$ci_direct[1],      r$ci_direct[2],
                r$ci_transformed[1], r$ci_transformed[2],
                r$ci_bca[1],         r$ci_bca[2]))
  }
}

cat("\n\nTable 3 — Lifetime characteristics at use condition:\n")
print_char_table("alpha0", "alpha0", alpha0_true, fmt = "%10.4f")
print_char_table("alpha1", "alpha1", alpha1_true, fmt = "%10.4f")
print_char_table("sigma",  "sigma",  sigma,        fmt = "%10.4f")
print_char_table("Median lifetime (cycles)", "median", med_true, "%10.0f")
print_char_table("MTTF (cycles)",            "mttf",   mttf_true, "%10.0f")
print_char_table(sprintf("Reliability at t0=%d", t0), "rel", rel_true, "%10.4f")

# =============================================================================
# FIGURE: ESTIMATES AND CIs (TRANSFORMED + BCa) FOR ALL 6 QUANTITIES
# =============================================================================

# ---- extract arrays ----------------------------------------------------------
nm_est  <- function(b) if (b == 0) "MLE" else paste0("MDPDE_", b)  # for estimates_all
nm_res  <- function(b) beta_labels[match(b, beta_vec)]
get_est <- function(b, key) results[[nm_res(b)]][[key]]$estimate
get_tlo <- function(b, key) results[[nm_res(b)]][[key]]$ci_transformed[1]
get_thi <- function(b, key) results[[nm_res(b)]][[key]]$ci_transformed[2]
get_blo <- function(b, key) results[[nm_res(b)]][[key]]$ci_bca[1]
get_bhi <- function(b, key) results[[nm_res(b)]][[key]]$ci_bca[2]
get_dlo <- function(b, key) results[[nm_res(b)]][[key]]$ci_direct[1]
get_dhi <- function(b, key) results[[nm_res(b)]][[key]]$ci_direct[2]

a0_est <- sapply(beta_vec, function(b) { nm <- nm_est(b); estimates_all[[nm]][1] })
a1_est <- sapply(beta_vec, function(b) { nm <- nm_est(b); estimates_all[[nm]][2] })
sg_est <- sapply(beta_vec, function(b) { nm <- nm_est(b); estimates_all[[nm]][3] })

a0_lo <- sapply(beta_vec, function(b) {
  nm <- nm_est(b); est <- estimates_all[[nm]]
  Sig <- Sigma_hat(est, Kvec, stress_mat, tau, ITs, b)
  est[1] - z95 * sqrt(diag(Sig))[1]
})
a0_hi <- sapply(beta_vec, function(b) {
  nm <- nm_est(b); est <- estimates_all[[nm]]
  Sig <- Sigma_hat(est, Kvec, stress_mat, tau, ITs, b)
  est[1] + z95 * sqrt(diag(Sig))[1]
})
a1_lo <- sapply(beta_vec, function(b) {
  nm <- nm_est(b); est <- estimates_all[[nm]]
  Sig <- Sigma_hat(est, Kvec, stress_mat, tau, ITs, b)
  est[2] - z95 * sqrt(diag(Sig))[2]
})
a1_hi <- sapply(beta_vec, function(b) {
  nm <- nm_est(b); est <- estimates_all[[nm]]
  Sig <- Sigma_hat(est, Kvec, stress_mat, tau, ITs, b)
  est[2] + z95 * sqrt(diag(Sig))[2]
})
sg_lo <- sapply(beta_vec, function(b) {
  nm <- nm_est(b); est <- estimates_all[[nm]]
  Sig <- Sigma_hat(est, Kvec, stress_mat, tau, ITs, b)
  est[3] - z95 * sqrt(diag(Sig))[3]
})
sg_hi <- sapply(beta_vec, function(b) {
  nm <- nm_est(b); est <- estimates_all[[nm]]
  Sig <- Sigma_hat(est, Kvec, stress_mat, tau, ITs, b)
  est[3] + z95 * sqrt(diag(Sig))[3]
})

a0_blo <- sapply(beta_vec, function(b) get_blo(b, "alpha0"))
a0_bhi <- sapply(beta_vec, function(b) get_bhi(b, "alpha0"))
a1_blo <- sapply(beta_vec, function(b) get_blo(b, "alpha1"))
a1_bhi <- sapply(beta_vec, function(b) get_bhi(b, "alpha1"))
sg_blo <- sapply(beta_vec, function(b) get_blo(b, "sigma"))
sg_bhi <- sapply(beta_vec, function(b) get_bhi(b, "sigma"))

med_est  <- sapply(beta_vec, function(b) get_est(b, "median"))
med_tlo  <- sapply(beta_vec, function(b) get_tlo(b, "median"))
med_thi  <- sapply(beta_vec, function(b) get_thi(b, "median"))
med_blo  <- sapply(beta_vec, function(b) get_blo(b, "median"))
med_bhi  <- sapply(beta_vec, function(b) get_bhi(b, "median"))

mttf_est <- sapply(beta_vec, function(b) get_est(b, "mttf"))
mttf_tlo <- sapply(beta_vec, function(b) get_tlo(b, "mttf"))
mttf_thi <- sapply(beta_vec, function(b) get_thi(b, "mttf"))
mttf_blo <- sapply(beta_vec, function(b) get_blo(b, "mttf"))
mttf_bhi <- sapply(beta_vec, function(b) get_bhi(b, "mttf"))

rel_est  <- sapply(beta_vec, function(b) get_est(b, "rel"))
rel_tlo  <- sapply(beta_vec, function(b) get_tlo(b, "rel"))
rel_thi  <- sapply(beta_vec, function(b) get_thi(b, "rel"))
rel_blo  <- sapply(beta_vec, function(b) get_blo(b, "rel"))
rel_bhi  <- sapply(beta_vec, function(b) get_bhi(b, "rel"))

med_dlo  <- sapply(beta_vec, function(b) get_dlo(b, "median"))
med_dhi  <- sapply(beta_vec, function(b) get_dhi(b, "median"))
mttf_dlo <- sapply(beta_vec, function(b) get_dlo(b, "mttf"))
mttf_dhi <- sapply(beta_vec, function(b) get_dhi(b, "mttf"))
rel_dlo  <- sapply(beta_vec, function(b) get_dlo(b, "rel"))
rel_dhi  <- sapply(beta_vec, function(b) get_dhi(b, "rel"))

# # ---- panel function ----------------------------------------------------------
# # Shows transformed CI (solid bars) and BCa CI (dashed bars) together
# one_panel <- function(x, est, tlo, thi, blo, bhi, true_val,
#                       main_lab, ylab = "Estimate", ylim = NULL) {
#   all_vals <- c(tlo, thi, blo, bhi)
#   if (is.null(ylim)) {
#     rng  <- range(all_vals, na.rm = TRUE)
#     pad  <- diff(rng) * 0.18
#     ylim <- c(rng[1] - pad, rng[2] + pad)
#   }
#   plot(x, est,
#        type = "b", pch = 19, lwd = 3.0, col = "black",
#        xaxt = "n", xlab = expression(beta), ylab = ylab,
#        main = main_lab, ylim = ylim,
#        cex = 1.80,  
#        cex.main = 2.0, cex.lab = 1.9, cex.axis = 1.7,
#        font.main = 2, font.lab = 2)                 # bold title and axis labels
#   axis(1, at = x, labels = beta_labels, cex.axis = 1.7, font = 2)  # bold tick labels
#   
#   cap <- 0.2    # wider end-caps (was 0.08)
#   off <- 0.11    # larger offset so caps do not overlap (was 0.08)
#   
#   for (i in x) {
#     lines(c(i - cap, i + cap),           c(tlo[i], tlo[i]), lwd = 2.4)
#     lines(c(i - cap, i + cap),           c(thi[i], thi[i]), lwd = 2.4)
#     lines(c(i, i),                       c(tlo[i], thi[i]), lwd = 2.4)
#     lines(c(i+off-cap, i+off+cap),       c(blo[i], blo[i]), lwd = 2.4, lty = 2, col = "#555555")
#     lines(c(i+off-cap, i+off+cap),       c(bhi[i], bhi[i]), lwd = 2.4, lty = 2, col = "#555555")
#     lines(c(i+off, i+off),               c(blo[i], bhi[i]), lwd = 2.4, lty = 2, col = "#555555")
#   }
#   abline(h = true_val, lty = 2, lwd = 2.4, col = "red")
# }
# 
# make_combined <- function() {
#   par(mfrow = c(2, 3), mar = c(4.5, 5.5, 3.5, 1.5), oma = c(0, 0, 3.5, 0),
#       font.axis = 2, font.lab = 2)
#   x <- seq_along(beta_vec)
#   
#   one_panel(x, a0_est, a0_lo, a0_hi, a0_blo, a0_bhi, alpha0_true,
#             main_lab = expression(alpha[0]), ylab = "Estimate")
#   one_panel(x, a1_est, a1_lo, a1_hi, a1_blo, a1_bhi, alpha1_true,
#             main_lab = expression(alpha[1]), ylab = "Estimate")
#   one_panel(x, sg_est, sg_lo, sg_hi, sg_blo, sg_bhi, sigma,
#             main_lab = expression(sigma),   ylab = "Estimate")
#   one_panel(x, med_est,  med_tlo,  med_thi,  med_blo,  med_bhi,  med_true,
#             main_lab = "Median (cycles)",   ylab = "Estimate")
#   one_panel(x, mttf_est, mttf_tlo, mttf_thi, mttf_blo, mttf_bhi, mttf_true,
#             main_lab = "MTTF (cycles)",     ylab = "Estimate")
#   one_panel(x, rel_est,  rel_tlo,  rel_thi,  rel_blo,  rel_bhi,  rel_true,
#             main_lab = bquote(R[0](t[0])),  ylab = "Estimate",
#             ylim = c(max(0, min(rel_blo, rel_tlo, na.rm = TRUE) - 0.04),
#                      min(1, max(rel_bhi, rel_thi, na.rm = TRUE) + 0.04)))
#   
#   par(fig = c(0,1,0,1), oma = c(0,0,0,0), mar = c(0,0,0,0), new = TRUE)
#   legend(x = 0.5, y = 1.03,
#          legend = c("Estimate", "Transformed CI", "BCa CI", "True value"),
#          lty    = c(1, 1, 2, 2),
#          lwd    = c(3.0, 2.4, 2.4, 2.4),
#          pch    = c(19, NA, NA, NA),
#          col    = c("black", "black", "#555555", "red"),
#          bty = "n", cex = 2, horiz = TRUE,
#          xjust = -0.25, yjust = -0.25, xpd = NA)
# }
# 
# dev.new(width = 14, height = 9)
# make_combined()
# 
# pdf("C:/Users/Kiran/Downloads/WMDPDE_CyALT_lognormal/fig_combined.pdf", width = 14, height = 9)
# make_combined()
# dev.off()
# cat("\nSaved: fig_combined.pdf\n")


one_panel <- function(x, est, tlo, thi, blo, bhi, true_val,
                      main_lab, ylab = "Estimate", ylim = NULL,
                      dlo = NULL, dhi = NULL) {
  
  all_vals <- c(tlo, thi, blo, bhi)
  if (!is.null(dlo)) all_vals <- c(all_vals, dlo, dhi)
  if (is.null(ylim)) {
    rng  <- range(all_vals, na.rm = TRUE)
    pad  <- diff(rng) * 0.18
    ylim <- c(rng[1] - pad, rng[2] + pad)
  }
  
  plot(x, est,
       type  = "b", pch = 19, lwd = 2.5, col = "black",
       xaxt  = "n",
       xlab  = expression(beta),
       ylab  = ylab,
       main  = main_lab,
       xlim  = c(min(x) - 0.5, max(x) + 0.5),   # ← side padding
       ylim  = ylim,
       cex   = 1.60,
       cex.main = 2.0, cex.lab = 1.9, cex.axis = 1.7,
       font.main = 2, font.lab = 2)
  
  axis(1, at = x, labels = beta_labels, cex.axis = 1.7, font = 2)
  
  cap <- 0.10   # ← reduced from 0.18
  off <- 0.09   # ← reduced from 0.13
  
  for (i in x) {
    # Direct CI — blue solid, offset left
    if (!is.null(dlo)) {
      lines(c(i-off-cap, i-off+cap), c(dlo[i], dlo[i]), lwd = 2.0, col = "#0072B2")
      lines(c(i-off-cap, i-off+cap), c(dhi[i], dhi[i]), lwd = 2.0, col = "#0072B2")
      lines(c(i-off, i-off),         c(dlo[i], dhi[i]), lwd = 2.0, col = "#0072B2")
    }
    # Transformed CI — black solid, centred
    lines(c(i-cap, i+cap), c(tlo[i], tlo[i]), lwd = 2.0)
    lines(c(i-cap, i+cap), c(thi[i], thi[i]), lwd = 2.0)
    lines(c(i, i),         c(tlo[i], thi[i]), lwd = 2.0)
    # BCa CI — grey dashed, offset right
    # BCa CI — grey, SOLID end caps, dashed vertical bar only
    lines(c(i+off-cap, i+off+cap), c(blo[i], blo[i]), lwd = 2.0, lty = 1, col = "#555555")
    lines(c(i+off-cap, i+off+cap), c(bhi[i], bhi[i]), lwd = 2.0, lty = 1, col = "#555555")
    lines(c(i+off, i+off),         c(blo[i], bhi[i]), lwd = 2.0, lty = 2, col = "#555555")
  }
  abline(h = true_val, lty = 2, lwd = 2.0, col = "red")
}

make_combined <- function() {
  par(mfrow = c(2, 3), mar = c(4.5, 5.5, 3.5, 1.5), oma = c(0, 0, 3.5, 0),
      font.axis = 2, font.lab = 2)
  x <- seq_along(beta_vec)
  
  # parameters — no direct CI shown (direct = transformed for alpha0, alpha1;
  # for sigma type="log" so they differ slightly — pass if desired)
  one_panel(x, a0_est, a0_lo, a0_hi, a0_blo, a0_bhi, alpha0_true,
            main_lab = expression(alpha[0]), ylab = "Estimate")
  one_panel(x, a1_est, a1_lo, a1_hi, a1_blo, a1_bhi, alpha1_true,
            main_lab = expression(alpha[1]), ylab = "Estimate")
  one_panel(x, sg_est, sg_lo, sg_hi, sg_blo, sg_bhi, sigma,
            main_lab = expression(sigma),   ylab = "Estimate")
  
  # characteristics — direct CI added (blue bars, offset left)
  one_panel(x, med_est,  med_tlo,  med_thi,  med_blo,  med_bhi,  med_true,
            main_lab = "Median (cycles)", ylab = "Estimate",
            dlo = med_dlo,  dhi = med_dhi)
  one_panel(x, mttf_est, mttf_tlo, mttf_thi, mttf_blo, mttf_bhi, mttf_true,
            main_lab = "MTTF (cycles)",   ylab = "Estimate",
            dlo = mttf_dlo, dhi = mttf_dhi)
  one_panel(x, rel_est,  rel_tlo,  rel_thi,  rel_blo,  rel_bhi,  rel_true,
            main_lab = bquote(R[0](t[0])), ylab = "Estimate",
            ylim = c(max(0, min(c(rel_blo, rel_tlo, rel_dlo), na.rm = TRUE) - 0.04),
                     min(1, max(c(rel_bhi, rel_thi, rel_dhi), na.rm = TRUE) + 0.04)),
            dlo = rel_dlo,  dhi = rel_dhi)
  
  # global legend — now includes direct CI
  par(fig = c(0,1,0,1), oma = c(0,0,0,0), mar = c(0,0,0,0), new = TRUE)
  legend(x = 0.5, y = 1.03,
         legend = c("Estimate", "Direct CI", "Transformed CI      ", "BCa CI", "True value"),
         lty    = c(1,   1,         1,              2,         2),
         lwd    = c(2.5, 2.0,      2.0,            2.0,       2.0),
         pch    = c(19,  NA,       NA,             NA,        NA),
         col    = c("black", "#0072B2", "black", "#555555", "red"),
         bty = "n", cex = 1.8, horiz = TRUE, lwd = 2.0,
         xjust = 0.025, yjust = -1.5, xpd = NA)
}

dev.new(width = 14, height = 9)
make_combined()

pdf("C:/Users/Kiran/Downloads/WMDPDE_CyALT_lognormal/fig_combined.pdf",
    width = 14, height = 9)
make_combined()
dev.off()
cat("\nSaved: fig_combined.pdf\n")