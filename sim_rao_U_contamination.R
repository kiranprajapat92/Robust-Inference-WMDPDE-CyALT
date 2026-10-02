rm(list = ls())

library(foreach)
library(doParallel)
library(optimx)
library(MASS)

source("/mnt/nfs/home/nkp117/largefiles/KP_Leandro_Maria_1/Testing/arxiv/cyalt_lognormal_WMDPDE1.R")
#source("C:/Users/Kiran/Downloads/WMDPDE_CyALT_lognormal/Testing/cyalt_lognormal_WMDPDE1.R")

n_cores <- 48


########################  SETUP  ########################  


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

Kvec       <- c(120, 80)
eps_vec <- seq(0, 0.10, length.out = 7)
ne         <- length(eps_vec)
cont_grp   <- 1
cont_cells <- c(1, 2, 3)

n_iter <- 1000
base_path <- "/mnt/nfs/home/nkp117/largefiles/KP_Leandro_Maria_1/Testing/arxiv/"


########################   ONE REPLICATION — record U_beta(theta0) components under each alternative ########################  

one_rep_U_check <- function(rep_id, eps) {
    
    dat_H1  <- gen_data_multicell(theta_alt_H1,  Kvec, stress_mat, tau, ITs,
                                  eps, cont_grp, cont_cells, NULL)
    dat_H2a <- gen_data_multicell(theta_alt_H2a, Kvec, stress_mat, tau, ITs,
                                  eps, cont_grp, cont_cells, NULL)
    dat_H2b <- gen_data_multicell(theta_alt_H2b, Kvec, stress_mat, tau, ITs,
                                  eps, cont_grp, cont_cells, NULL)
    dat_H3  <- gen_data_multicell(theta_alt_H3,  Kvec, stress_mat, tau, ITs,
                                  eps, cont_grp, cont_cells, NULL)
    
    if (!is_valid_dataset(dat_H1)) return(NULL)
    
    # 4 hypotheses x 3 components x nb betas
    out <- numeric(nb * 12)
    
    for (bi in seq_along(beta_vec)) {
        
        b <- beta_vec[bi]
        
        UH1  <- U_beta_vec(theta0, dat_H1,  Kvec, stress_mat, tau, ITs, b)
        UH2a <- U_beta_vec(theta0, dat_H2a, Kvec, stress_mat, tau, ITs, b)
        UH2b <- U_beta_vec(theta0, dat_H2b, Kvec, stress_mat, tau, ITs, b)
        UH3  <- U_beta_vec(theta0, dat_H3,  Kvec, stress_mat, tau, ITs, b)
        
        idx <- (bi - 1) * 12
        out[idx+1]  <- UH1[1]
        out[idx+2]  <- UH1[2]
        out[idx+3]  <- UH1[3]
        out[idx+4]  <- UH2a[1]
        out[idx+5]  <- UH2a[2]
        out[idx+6]  <- UH2a[3]
        out[idx+7]  <- UH2b[1]
        out[idx+8]  <- UH2b[2]
        out[idx+9]  <- UH2b[3]
        out[idx+10] <- UH3[1]
        out[idx+11] <- UH3[2]
        out[idx+12] <- UH3[3]
    }
    out
}


########################  PARALLEL SETUP ########################  


cl      <- makeCluster(n_cores)
registerDoParallel(cl)
cat(sprintf("Using %d cores\n", n_cores))

clusterExport(cl, varlist = c(
    "theta0","theta_alt_H1","theta_alt_H2a","theta_alt_H2b","theta_alt_H3",
    "tau","stress_mat","ITs","beta_vec","nb","Kvec",
    "cont_grp","cont_cells",
    "one_rep_U_check","gen_data_multicell","is_valid_dataset",
    "U_beta_vec",
    "p_i_theta","W_i_theta_matrix","mu_fun","eta_fun",
    "sbar_fun","aij_fun","interval_prob","surv_prob"
))
clusterEvalQ(cl, { library(optimx); library(MASS) })


########################  RUN  ########################  

# dims: eps x beta x hypothesis x component

U_mean <- array(NA, dim = c(ne, nb, 4, 3),
                dimnames = list(as.character(eps_vec), beta_labels,
                                c("H1","H2a","H2b","H3"),
                                c("U1","U2","U3")))

for (ei in seq_along(eps_vec)) {
    
    eps <- eps_vec[ei]
    cat(sprintf("eps = %.3f ...\n", eps))
    
    reps <- foreach(
        rep            = seq_len(n_iter),
        .combine       = "rbind",
        .packages      = c("optimx","MASS"),
        .errorhandling = "remove"
    ) %dopar% {
        suppressWarnings(one_rep_U_check(rep, eps))
    }
    
    if (!is.null(reps) && nrow(reps) > 0) {
        for (bi in seq_along(beta_vec)) {
            idx <- (bi - 1) * 12
            U_mean[ei, bi, "H1",  ] <- colMeans(reps[, idx+1:3,  drop=FALSE], na.rm = TRUE)
            U_mean[ei, bi, "H2a", ] <- colMeans(reps[, idx+4:6,  drop=FALSE], na.rm = TRUE)
            U_mean[ei, bi, "H2b", ] <- colMeans(reps[, idx+7:9,  drop=FALSE], na.rm = TRUE)
            U_mean[ei, bi, "H3",  ] <- colMeans(reps[, idx+10:12,drop=FALSE], na.rm = TRUE)
        }
    }
}

stopCluster(cl)

save(U_mean, eps_vec, beta_vec, beta_labels,
     theta0, theta_alt_H1, theta_alt_H2a, theta_alt_H2b, theta_alt_H3,
     file = paste0(base_path, "U_bias_check.RData"))
cat("Saved U_bias_check.RData\n")


########################  PLOT  ########################  

cols <- c("#000000","#0072B2","#009E73","#D55E00","#7B2FBE","#CC79A7")
pchs <- c(16, 1, 2, 5, 6, 0)

one_U_panel <- function(eps_vec, ymat, ylab_txt, main_lab) {
    plot(NA, xlim = range(eps_vec),
         ylim = range(c(ymat, 0), na.rm = TRUE),
         xlab = expression(epsilon), ylab = "",
         main = main_lab, cex.main = 2.3, cex.lab = 2.2, cex.axis = 2.0,
         font.main = 2, font.lab = 1, las = 1)
    title(ylab = ylab_txt, line = 5.2, cex.lab = 2.2, font.lab = 1)
    abline(h = 0, lty = 2, col = "grey30", lwd = 1.3)
    for (bi in seq_along(beta_vec)) {
        lines(eps_vec, ymat[, bi], type = "b",
              col = cols[bi], pch = pchs[bi], cex = 1.3, lwd = 1.3)
    }
}

plot_all_U_bias <- function() {
    
    par(mfrow = c(2, 4), mar = c(7, 8.5, 5.5, 2), font.axis = 1)
    
    # H1: U1, U2, U3
    one_U_panel(eps_vec, U_mean[, , "H1", "U1"],
                expression(U[1]^beta*(theta[0])), "H1")
    one_U_panel(eps_vec, U_mean[, , "H1", "U2"],
                expression(U[2]^beta*(theta[0])), "H1")
    one_U_panel(eps_vec, U_mean[, , "H1", "U3"],
                expression(U[3]^beta*(theta[0])), "H1")
    
    # H2a: U1 only
    one_U_panel(eps_vec, U_mean[, , "H2a", "U1"],
                expression(U[1]^beta*(theta[0])), "H2a")
    
    # H2b: U2 only
    one_U_panel(eps_vec, U_mean[, , "H2b", "U2"],
                expression(U[2]^beta*(theta[0])), "H2b")
    
    # H3: U1, U2
    one_U_panel(eps_vec, U_mean[, , "H3", "U1"],
                expression(U[1]^beta*(theta[0])), "H3")
    one_U_panel(eps_vec, U_mean[, , "H3", "U2"],
                expression(U[2]^beta*(theta[0])), "H3")
    
    # legend
    plot.new()
    legend("center",
           legend = c("null value (0)", beta_labels),
           col    = c("grey30", cols),
           lty    = c(2, rep(1, 6)),
           pch    = c(NA, pchs),
           lwd    = 1.3,
           cex    = 1.7,
           bty    = "n")
}

pdf(paste0(base_path, "fig_U_bias_full.pdf"), width = 18, height = 9)
plot_all_U_bias()
dev.off()

cat("Saved fig_U_bias_full.pdf\n")
