rm(list = ls())

base_path <- "/mnt/nfs/home/nkp117/largefiles/KP_Leandro_Maria_1/"

# =============================================================================
# CI SIMULATION — PARAMETERS + LIFETIME CHARACTERISTICS
# Quantities : alpha0, alpha1, sigma, median, MTTF, R0(t0)
# CI types   : Direct + Transformed + BCa  (all three, all six quantities)
# Saves      : everything in RData for later use
# Reports    : paper-format table (direct + transformed + BCa for characteristics)
#
# RUN ONE eps PER HPC JOB:
#   Rscript sim_CIs.R 0.00
#   Rscript sim_CIs.R 0.05
#   Rscript sim_CIs.R 0.10
#   Rscript sim_CIs.R 0.20
#   Rscript sim_CIs.R 0.25
#
# CHANGES FROM PREVIOUS VERSION:
#   1. All 6 beta values: 0, 0.2, 0.4, 0.6, 0.8, 1.0
#   2. Direct CI added alongside Transformed and BCa in LaTeX table
#   3. B_boot default = 500
# =============================================================================

library(foreach)
library(doParallel)
library(optimx)
library(MASS)

source(paste0(base_path, "cyalt_lognormal_WMDPDE1.R"))

# =============================================================================
# READ eps FROM COMMAND LINE
# =============================================================================

args    <- commandArgs(trailingOnly = TRUE)
eps_val <- if (length(args) >= 1) as.numeric(args[1]) else 0.00
n_iter  <- if (length(args) >= 2) as.integer(args[2]) else 1000
B_boot  <- if (length(args) >= 3) as.integer(args[3]) else 500
cat(sprintf("eps=%.3f  n_iter=%d  B_boot=%d\n\n", eps_val, n_iter, B_boot))

# =============================================================================
# DESIGN — must match MSE simulation exactly
# =============================================================================

alpha0_true <- 5.0;  alpha1_true <- -2.0;  sigma_true <- 0.5
theta_true  <- c(alpha0_true, alpha1_true, sigma_true)

tau        <- 0.40
Kvec       <- c(120, 80)
s_0F       <- 0.00;  s_0C <- 0.20
s_F        <- 0.40;  s_1C <- 0.65;  s_2C <- 1.00
stress_mat <- matrix(c(s_F, s_1C, s_F, s_2C), nrow = 2, byrow = TRUE)
ITs        <- c(15, 25, 35, 50, 65, 80)
t0         <- 100

med_true  <- quantile_use(theta_true, s_0C, s_0F, tau, q = 0.5)
mttf_true <- mttf_use(theta_true, s_0C, s_0F, tau)
rel_true  <- reliability_use(theta_true, s_0C, s_0F, tau, t0 = t0)

cat(sprintf("True theta: (%.2f, %.2f, %.3f)\n", theta_true[1], theta_true[2], theta_true[3]))
cat(sprintf("True chars: median=%.4f, MTTF=%.4f, R0=%.4f\n\n",
            med_true, mttf_true, rel_true))

# =============================================================================
# CHANGE 1: All 6 beta values
# =============================================================================
beta_vec    <- c(0, 0.2, 0.4, 0.6, 0.8, 1.0)
beta_labels <- c("MLE", "0.2", "0.4", "0.6", "0.8", "1.0")

cont_grp    <- 1
cont_cells  <- c(1, 2, 3)

alpha       <- 0.05
z_val       <- qnorm(1 - alpha / 2)
nb          <- length(beta_vec)

# =============================================================================
# QUANTITY DEFINITIONS
# 6 quantities total:
#   params [1:3] : alpha0, alpha1, sigma  → direct CI (proposed)
#   chars  [4:6] : median, MTTF, R0      → transformed CI (proposed)
# All quantities also get BCa CI
# =============================================================================

qty_names  <- c("alpha0","alpha1","sigma","median","mttf","reliability")
true_vals  <- c(theta_true, med_true, mttf_true, rel_true)
ci_types   <- c("none","none","none","log","log","logit")

# =============================================================================
# gen_data_multicell
# =============================================================================

gen_data_multicell <- function(theta, Kvec, stress_mat, tau_val, ITs,
                               eps, cont_grp, cont_cells, seed) {
  set.seed(seed)
  R <- nrow(stress_mat);  L <- length(ITs)
  counts_list <- vector("list", R)
  for (i in seq_len(R)) {
    sC   <- stress_mat[i, 2];  sF <- stress_mat[i, 1]
    mu_i <- theta[1] - log(tau_val * exp(-theta[2]*sC) +
                            (1-tau_val) * exp(-theta[2]*sF))
    pvec <- numeric(L + 1)
    for (j in seq_len(L)) {
      lo <- if (j == 1) -Inf else (log(ITs[j-1]) - mu_i) / theta[3]
      pvec[j] <- pnorm((log(ITs[j]) - mu_i) / theta[3]) - pnorm(lo)
    }
    pvec[L+1] <- 1 - pnorm((log(ITs[L]) - mu_i) / theta[3])
    pvec      <- pmax(pvec, 1e-10);  pvec <- pvec / sum(pvec)
    K_i       <- Kvec[i]
    if (i == cont_grp && eps > 0) {
      n_clean          <- floor(K_i * (1 - eps));  n_cont <- K_i - n_clean
      clean_c          <- drop(rmultinom(1, n_clean, pvec))
      q_c              <- numeric(L+1);  q_c[cont_cells] <- 1/length(cont_cells)
      counts_list[[i]] <- clean_c + drop(rmultinom(1, n_cont, q_c))
    } else {
      counts_list[[i]] <- drop(rmultinom(1, K_i, pvec))
    }
  }
  counts_list
}

# =============================================================================
# PARALLEL SETUP
# =============================================================================

n_cores <- as.integer(Sys.getenv("SLURM_NTASKS",
                     unset = as.character(max(1, detectCores() - 1))))
cl      <- makeCluster(n_cores)
registerDoParallel(cl)
cat(sprintf("Using %d cores\n\n", n_cores))

clusterExport(cl, "base_path")
clusterEvalQ(cl, {
  library(optimx);  library(MASS)
  source(paste0(base_path, "cyalt_lognormal_WMDPDE1.R"))
})

clusterExport(cl, varlist = c(
  "theta_true","Kvec","stress_mat","tau","ITs","t0",
  "s_0C","s_0F","beta_vec","beta_labels","nb",
  "cont_grp","cont_cells","eps_val",
  "true_vals","ci_types","qty_names",
  "med_true","mttf_true","rel_true",
  "alpha","z_val","B_boot",
  "gen_data_multicell"
))

# =============================================================================
# HELPER: extract all 6 quantities from one theta vector
# =============================================================================

get_all_quantities <- function(th, s_0C, s_0F, tau, t0) {
  c(th[1],
    th[2],
    th[3],
    quantile_use(th, s_0C, s_0F, tau, 0.5),
    mttf_use(th, s_0C, s_0F, tau),
    reliability_use(th, s_0C, s_0F, tau, t0))
}

clusterExport(cl, "get_all_quantities")

# =============================================================================
# PARALLEL SIMULATION
# =============================================================================

cat("Running CI simulation...\n")
t_start <- proc.time()

rep_results <- foreach(
  rep            = seq_len(n_iter),
  .combine       = "rbind",
  .packages      = c("optimx","MASS"),
  .errorhandling = "remove"
) %dopar% {

  # ----------------------------------------------------------------
  # 1. Generate contaminated dataset
  # ----------------------------------------------------------------
  dat <- gen_data_multicell(
    theta = theta_true, Kvec = Kvec, stress_mat = stress_mat,
    tau_val = tau, ITs = ITs,
    eps = eps_val, cont_grp = cont_grp, cont_cells = cont_cells,
    seed = rep * 31 + round(eps_val * 10000)
  )
  if (!is_valid_dataset(dat)) return(NULL)

  # ----------------------------------------------------------------
  # 2. Fit all beta values
  # ----------------------------------------------------------------
  ests <- suppressWarnings(
    fit_all_betas(dat, Kvec, stress_mat, tau, ITs,
                  beta_vec = beta_vec, init = theta_true)
  )

  row_vals <- numeric(nb * 36)

  for (bi in seq_along(beta_vec)) {

    b   <- beta_vec[bi]
    key <- if (b == 0) "MLE" else paste0("MDPDE_", b)
    est <- ests[[key]]

    if (is.null(est) || any(is.na(est)) ||
        est[2] >= 0 || est[3] <= 0 || any(abs(est) > 1e6)) {
      row_vals[((bi-1)*36 + 1):(bi*36)] <- NA
      next
    }

    # Asymptotic covariance
    Sig <- tryCatch(Sigma_hat(est, Kvec, stress_mat, tau, ITs, b),
                    error = function(e) matrix(NA, 3, 3))
    if (any(is.na(Sig))) {
      row_vals[((bi-1)*36 + 1):(bi*36)] <- NA
      next
    }
    se_par <- sqrt(diag(Sig))

    se_med  <- tryCatch(quantile_use_se(est, Kvec, stress_mat, tau, ITs,
                                         b, s_0C, s_0F, 0.5),
                        error = function(e) NA)
    se_mttf <- tryCatch(mttf_use_se(est, Kvec, stress_mat, tau, ITs,
                                     b, s_0C, s_0F),
                        error = function(e) NA)
    se_rel  <- tryCatch(reliability_use_se(est, Kvec, stress_mat, tau, ITs,
                                            b, s_0C, s_0F, t0),
                        error = function(e) NA)

    qty_hat <- get_all_quantities(est, s_0C, s_0F, tau, t0)
    qty_se  <- c(se_par[1], se_par[2], se_par[3], se_med, se_mttf, se_rel)

    # Bootstrap — one pass per beta, all 6 quantities
    boot_mat  <- matrix(NA, B_boot, 6)
    b_success <- 0;  attempt <- 0

    while (b_success < B_boot && attempt < 5 * B_boot) {
      attempt <- attempt + 1
      dat_b   <- sim_cyalt_data_continuous(est, Kvec, stress_mat, tau, ITs)
      if (!is_valid_dataset(dat_b)) next
      est_b <- tryCatch(
        fit_mdpde(dat_b, Kvec, stress_mat, tau, ITs, beta = b, init = est),
        error = function(e) rep(NA, 3))
      if (any(is.na(est_b)) || est_b[2] >= 0 || est_b[3] <= 0) next
      b_success <- b_success + 1
      boot_mat[b_success, ] <- get_all_quantities(est_b, s_0C, s_0F, tau, t0)
    }
    boot_mat <- boot_mat[seq_len(b_success), , drop = FALSE]

    # Jackknife — one pass per beta, all 6 quantities
    R_grp <- length(dat);  Lp1 <- length(ITs) + 1
    jack_list <- list();   jack_w <- integer(0)

    for (i in seq_len(R_grp)) {
      for (j in seq_len(Lp1 - 1)) {
        n_ij <- dat[[i]][j]
        if (n_ij == 0) next
        dat_j          <- dat
        dat_j[[i]][j]  <- n_ij - 1
        dat_j[[i]][Lp1] <- dat_j[[i]][Lp1] + 1
        est_j <- tryCatch(
          fit_mdpde(dat_j, Kvec, stress_mat, tau, ITs, beta = b, init = est),
          error = function(e) rep(NA, 3))
        psi_j <- if (any(is.na(est_j)) || est_j[2] >= 0 || est_j[3] <= 0) {
          qty_hat
        } else {
          get_all_quantities(est_j, s_0C, s_0F, tau, t0)
        }
        jack_list <- c(jack_list, list(psi_j))
        jack_w    <- c(jack_w, n_ij)
      }
    }
    jack_mat <- do.call(rbind, jack_list)

    # Compute all 3 CI types for all 6 quantities
    qty_block <- numeric(36)

    for (q in 1:6) {
      qhat  <- qty_hat[q]
      qtrue <- true_vals[q]
      se_q  <- qty_se[q]
      type_q <- ci_types[q]

      if (is.na(se_q) || se_q <= 0 || is.na(qhat)) {
        qty_block[((q-1)*6 + 1):(q*6)] <- NA
        next
      }

      # Direct CI
      ci_dir <- c(qhat - z_val * se_q, qhat + z_val * se_q)
      cp_dir <- as.integer(ci_dir[1] <= qtrue && qtrue <= ci_dir[2])
      aw_dir <- ci_dir[2] - ci_dir[1]

      # Transformed CI
      if (type_q == "log") {
        ci_tr <- c(qhat * exp(-z_val * se_q / qhat),
                   qhat * exp( z_val * se_q / qhat))
      } else if (type_q == "logit") {
        S     <- exp(z_val * se_q / (qhat * (1 - qhat)))
        ci_tr <- c(qhat / (qhat + (1 - qhat) * S),
                   qhat / (qhat + (1 - qhat) / S))
      } else {
        ci_tr <- ci_dir
      }
      cp_tr <- as.integer(ci_tr[1] <= qtrue && qtrue <= ci_tr[2])
      aw_tr <- ci_tr[2] - ci_tr[1]

      # BCa CI
      if (b_success < 10) {
        cp_bca <- NA;  aw_bca <- NA
      } else {
        boot_q <- boot_mat[, q]

        if (!is.null(jack_mat) && nrow(jack_mat) >= 2) {
          jq   <- jack_mat[, q]
          jbar <- sum(jack_w * jq) / sum(jack_w)
          num_j <- sum(jack_w * (jq - jbar)^3)
          den_j <- sum(jack_w * (jq - jbar)^2)
          gam_q <- if (den_j < 1e-12) 0 else (1/6) * num_j * den_j^(-3/2)
        } else {
          gam_q <- 0
        }

        prop <- mean(boot_q <= qhat)
        if (prop <= 0) prop <- 0.0001
        if (prop >= 1) prop <- 0.9999
        z0      <- qnorm(prop)
        z_alpha <- z_val
        q_lo    <- pnorm(z0 + (z0 - z_alpha) / (1 - gam_q * (z0 - z_alpha)))
        q_hi    <- pnorm(z0 + (z0 + z_alpha) / (1 - gam_q * (z0 + z_alpha)))
        bs      <- sort(boot_q)
        B_eff   <- length(bs)
        ci_bca  <- c(bs[max(floor(q_lo * B_eff), 1)],
                     bs[min(floor(q_hi * B_eff), B_eff)])
        cp_bca  <- as.integer(ci_bca[1] <= qtrue && qtrue <= ci_bca[2])
        aw_bca  <- ci_bca[2] - ci_bca[1]
      }

      idx <- (q - 1) * 6
      qty_block[idx + 1] <- cp_dir
      qty_block[idx + 2] <- aw_dir
      qty_block[idx + 3] <- cp_tr
      qty_block[idx + 4] <- aw_tr
      qty_block[idx + 5] <- cp_bca
      qty_block[idx + 6] <- aw_bca
    }

    row_vals[((bi-1)*36 + 1):(bi*36)] <- qty_block
  }

  row_vals
}

t_elapsed <- proc.time() - t_start
cat(sprintf("\nTotal time: %.1f minutes\n", t_elapsed["elapsed"] / 60))
stopCluster(cl)

# =============================================================================
# AGGREGATE
# res: [nb, 6 quantities, 6 metrics]
# =============================================================================

metric_names <- c("CP_dir","AW_dir","CP_trans","AW_trans","CP_bca","AW_bca")

res <- array(
  NA,
  dim      = c(nb, 6, 6),
  dimnames = list(beta_labels, qty_names, metric_names)
)
n_valid <- integer(nb)

if (!is.null(rep_results) && nrow(rep_results) > 0) {
  mat <- rep_results[complete.cases(rep_results), , drop = FALSE]

  for (bi in seq_len(nb)) {
    cols_bi <- ((bi-1)*36 + 1):(bi*36)
    block   <- mat[, cols_bi, drop = FALSE]
    ok      <- rowSums(is.na(block)) == 0
    if (!any(ok)) next
    n_valid[bi] <- sum(ok)
    blk <- block[ok, , drop = FALSE]

    for (q in 1:6) {
      idx        <- (q - 1) * 6
      res[bi, q, ] <- colMeans(blk[, (idx+1):(idx+6), drop=FALSE], na.rm=TRUE)
    }
  }
}

cat(sprintf("\nValid replications per beta:\n"))
cat(paste(beta_labels, n_valid, sep = "=", collapse = "  "), "\n")

# =============================================================================
# SAVE
# =============================================================================

save_path <- sprintf(
  "%sci_full_eps%s.RData",
  base_path,
  gsub("\\.", "", sprintf("%.3f", eps_val))
)
save(
  res, n_valid,
  eps_val, beta_vec, beta_labels,
  qty_names, metric_names, ci_types, true_vals,
  theta_true, med_true, mttf_true, rel_true, t0,
  Kvec, stress_mat, tau, ITs, s_0C, s_0F,
  n_iter, B_boot, alpha,
  file = save_path
)
cat(sprintf("\nSaved: %s\n", save_path))

# =============================================================================
# PRINT CONSOLE SUMMARY
# =============================================================================

par_idx  <- 1:3
char_idx <- 4:6

cat("\n\n=== PARAMETERS (proposed: direct CI) ===\n")
for (q in par_idx) {
  cat(sprintf("\n--- %s (true = %.5f) ---\n", qty_names[q], true_vals[q]))
  cat(sprintf("%-6s  %7s %8s  %7s %8s\n",
              "beta", "CP_dir", "AW_dir", "CP_bca", "AW_bca"))
  for (bi in seq_len(nb)) {
    r <- res[bi, q, ]
    cat(sprintf("%-6s  %7.4f %8.4f  %7.4f %8.4f\n",
                beta_labels[bi],
                r["CP_dir"],  r["AW_dir"],
                r["CP_bca"],  r["AW_bca"]))
  }
}

cat("\n\n=== CHARACTERISTICS (proposed: transformed CI) ===\n")
for (q in char_idx) {
  cat(sprintf("\n--- %s (true = %.5f) ---\n", qty_names[q], true_vals[q]))
  cat(sprintf("%-6s  %9s %8s  %9s %8s  %7s %8s\n",
              "beta", "CP_dir", "AW_dir", "CP_trans", "AW_trans",
              "CP_bca", "AW_bca"))
  for (bi in seq_len(nb)) {
    r <- res[bi, q, ]
    cat(sprintf("%-6s  %9.4f %8.4f  %9.4f %8.4f  %7.4f %8.4f\n",
                beta_labels[bi],
                r["CP_dir"],   r["AW_dir"],
                r["CP_trans"], r["AW_trans"],
                r["CP_bca"],   r["AW_bca"]))
  }
}

# =============================================================================
# CHANGE 2: LaTeX table — Direct + Transformed + BCa for characteristics
# =============================================================================

make_ci_latex <- function(res, eps_val, beta_vec, beta_labels,
                          n_iter, B_boot,
                          char_idx = 4:6,
                          char_tex = c("$t_{0.5,0}$",
                                       "$\\mathrm{MTTF}_0$",
                                       "$R_0(t_0)$")) {
  nb      <- length(beta_vec)
  eps_lbl <- gsub("\\.", "", sprintf("%.3f", eps_val))

  rows <- character(0)
  for (bi in seq_len(nb)) {
    blab <- if (beta_vec[bi] == 0) "$0$" else sprintf("$%.1f$", beta_vec[bi])

    # Direct CI columns
    dir_str <- paste(sapply(char_idx, function(q)
      sprintf("%.4f & %.3f",
              res[bi, q, "CP_dir"],
              res[bi, q, "AW_dir"])),
      collapse = " & ")

    # Transformed CI columns
    tr_str <- paste(sapply(char_idx, function(q)
      sprintf("%.4f & %.3f",
              res[bi, q, "CP_trans"],
              res[bi, q, "AW_trans"])),
      collapse = " & ")

    # BCa CI columns
    bc_str <- paste(sapply(char_idx, function(q)
      sprintf("%.4f & %.3f",
              res[bi, q, "CP_bca"],
              res[bi, q, "AW_bca"])),
      collapse = " & ")

    rows <- c(rows, sprintf("%s & %s & %s & %s \\\\",
                            blab, dir_str, tr_str, bc_str))
  }

  # Header: char names repeated under Direct, Transformed, BCa
  char_hdr <- paste(
    paste0("\\multicolumn{2}{c}{", char_tex, "}", collapse = " & ")
  )

  cat(paste(c(
    "",
    "% ============================================================",
    sprintf("%% CI TABLE: Direct + Transformed + BCa, eps = %.3f", eps_val),
    "% ============================================================",
    "\\begin{table}[htbp]",
    "\\centering\\small",
    "\\setlength{\\tabcolsep}{3pt}",
    sprintf(
      "\\caption{Coverage probability (CP) and average width (AW) of $95\\%%$ direct, transformed and BCa confidence intervals for lifetime characteristics at $\\varepsilon = %.3f$ ($n = %d$ replications, $B = %d$ bootstrap replicates).}",
      eps_val, n_iter, B_boot),
    sprintf("\\label{tab:ci_chars_eps%s}", eps_lbl),
    "\\begin{tabular}{l rr rr rr rr rr rr rr rr rr}\\toprule",
    "& \\multicolumn{6}{c}{Direct} & \\multicolumn{6}{c}{Transformed} & \\multicolumn{6}{c}{BCa} \\\\",
    "\\cmidrule(lr){2-7}\\cmidrule(lr){8-13}\\cmidrule(lr){14-19}",
    paste0("& ", char_hdr, " & ", char_hdr, " & ", char_hdr, " \\\\"),
    paste0("$\\beta$ & ", paste(rep("CP & AW", 9), collapse = " & "), " \\\\"),
    "\\midrule",
    paste(rows, collapse = "\n"),
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}",
    ""
  ), collapse = "\n"))
}

make_ci_latex(res, eps_val, beta_vec, beta_labels, n_iter, B_boot)

# In Rstudio after loading RData file.

# make_ci_latex <- function(res, eps_val, beta_vec, beta_labels,
#                           n_iter, B_boot,
#                           char_idx = 4:6,
#                           char_tex = c("$t_{0.5,0}$",
#                                        "$\\mathrm{MTTF}_0$",
#                                        "$R_0(t_0)$")) {
#   nb      <- length(beta_vec)
#   eps_lbl <- gsub("\\.", "", sprintf("%.3f", eps_val))
#   rows <- character(0)
#   for (bi in seq_len(nb)) {
#     blab <- if (beta_vec[bi] == 0) "$0$" else sprintf("$%.1f$", beta_vec[bi])
#     dir_str <- paste(sapply(char_idx, function(q)
#       sprintf("%.2f & %.2f", res[bi,q,"CP_dir"], res[bi,q,"AW_dir"])),
#       collapse = " & ")
#     tr_str <- paste(sapply(char_idx, function(q)
#       sprintf("%.2f & %.2f", res[bi,q,"CP_trans"], res[bi,q,"AW_trans"])),
#       collapse = " & ")
#     bc_str <- paste(sapply(char_idx, function(q)
#       sprintf("%.2f & %.2f", res[bi,q,"CP_bca"], res[bi,q,"AW_bca"])),
#       collapse = " & ")
#     rows <- c(rows, sprintf("%s & %s & %s & %s \\\\",
#                             blab, dir_str, tr_str, bc_str))
#   }
#   char_hdr <- paste(
#     paste0("\\multicolumn{2}{c}{", char_tex, "}", collapse = " & "))
#   cat(paste(c(
#     "\\begin{table}[htbp]",
#     "\\centering\\small",
#     "\\setlength{\\tabcolsep}{3pt}",
#     sprintf("\\caption{CP and AW of $95\\%%$ direct, transformed and BCa CIs at $\\varepsilon = %.3f$ ($n = %d$, $B = %d$).}",
#             eps_val, n_iter, B_boot),
#     sprintf("\\label{tab:ci_chars_eps%s}", eps_lbl),
#     "\\begin{tabular}{l rr rr rr rr rr rr rr rr rr}\\toprule",
#     "& \\multicolumn{6}{c}{Direct} & \\multicolumn{6}{c}{Transformed} & \\multicolumn{6}{c}{BCa} \\\\",
#     "\\cmidrule(lr){2-7}\\cmidrule(lr){8-13}\\cmidrule(lr){14-19}",
#     paste0("& ", char_hdr, " & ", char_hdr, " & ", char_hdr, " \\\\"),
#     paste0("$\\beta$ & ", paste(rep("CP & AW", 9), collapse = " & "), " \\\\"),
#     "\\midrule",
#     paste(rows, collapse = "\n"),
#     "\\bottomrule",
#     "\\end{tabular}",
#     "\\end{table}"
#   ), collapse = "\n"))
# }
# 
# # Call it
# make_ci_latex(res, eps_val, beta_vec, beta_labels, n_iter, B_boot)