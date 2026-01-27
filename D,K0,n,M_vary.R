# ==============================================================================
# SCRIPT: Sensitivity Analysis of Online Misdiscovery Rate Control (OMDRC)
# Description: This script evaluates the robustness of OMDRC-DD under varying 
#              hyperparameters (n, K0, D, M) compared to the Oracle baseline.
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. Environment Setup
# ------------------------------------------------------------------------------
library(Matrix)
library(REBayes)
library(foreach)
library(doParallel)
library(ggplot2)
library(dplyr)
library(patchwork)
library(RColorBrewer) 
library(ggpubr)

# The user must set the working directory to the location of 'OMDRC.R'.
# A placeholder path is provided below.
# Example: setwd("/path/to/your/project/folder")
setwd("/path/to/your/project/folder") # <-- USER: PLEASE MODIFY THIS PATH
source('OMDRC.R')

# ------------------------------------------------------------------------------
# 1. Core Simulation Engine (Final Version)
# ------------------------------------------------------------------------------
# This version calculates final MDR and FDR values directly within the function.
# MDR = E[num_mdr] / E[den_mdr]  (Ratio of Expectations)
# FDR = E[num_fdr / den_fdr]   (Expectation of Ratio)
run_simulation_ratio <- function(m, ini, n, pi, alpha, reps, D, M, methods_to_run = c("OR", "DD")) {
  
  # For MDR: Store Numerators and Denominators for averaging
  num_mdr_or <- matrix(0, reps, length(m)); den_mdr_or <- matrix(0, reps, length(m))
  num_mdr_dd <- matrix(0, reps, length(m)); den_mdr_dd <- matrix(0, reps, length(m))
  
  # For FDR: Store the final ratio of each replication for averaging
  fdr_or <- matrix(0, reps, length(m))
  fdr_dd <- matrix(0, reps, length(m))
  
  # Ground truth parameters
  mu1 <- 2; mu2 <- -3; sd1 <- 0.7
  
  # Parallel Backend Setup
  cl <- makeCluster(min(parallel::detectCores() - 1, 10))
  registerDoParallel(cl)
  clusterExport(cl, c("m", "ini", "n", "pi", "reps", "alpha", "D", "M", "mu1", "mu2", "sd1", "methods_to_run"), envir = environment())
  
  result <- foreach(r = 1:reps, .packages = c("kedd", "REBayes")) %dopar% {
    source('OMDRC.R')
    set.seed(r)
    
    # Data Generation
    N <- max(m) + ini
    theta <- rbinom(N, size = 1, prob = pi)
    z0 <- rnorm(N, 0, 1)
    mu_alt <- ifelse(rbinom(N, 1, 0.5) == 1, mu1, mu2)
    z_alt <- rnorm(N, mu_alt, sd1)
    z <- ifelse(theta == 0, z0, z_alt)
    f0 <- dnorm(z, 0, 1)
    f1 <- 0.5 * dnorm(z, mu1, sd1) + 0.5 * dnorm(z, mu2, sd1)
    z.lmdr <- (pi * f1) / ((1 - pi) * f0 + pi * f1)
    z_stream <- z[(ini + 1):N]
    z.lmdr_stream <- z.lmdr[(ini + 1):N]
    theta_stream <- theta[(ini + 1):N]
    z1 <- rnorm(n, ifelse(rbinom(n, 1, 0.5) == 1, mu1, mu2), sd1)
    
    res <- list()
    # Execute Oracle OMDRC
    if ("OR" %in% methods_to_run) {
      de <- OMDRC_OR(z.lmdr_stream, alpha)$de
      res$num_mdr_or <- sapply(m, function(k) sum(theta_stream[1:k] * (1 - de[1:k])))
      res$den_mdr_or <- sapply(m, function(k) sum(theta_stream[1:k]))
      res$fdr_or     <- sapply(m, function(k) sum((1-theta_stream[1:k]) * de[1:k]) / max(sum(de[1:k]), 1))
    }
    # Execute Data-Driven OMDRC
    if ("DD" %in% methods_to_run) {
      de <- OMDRC_DD(z_stream, z[1:ini], z1, alpha, D = D, M = M)$de
      res$num_mdr_dd <- sapply(m, function(k) sum(theta_stream[1:k] * (1 - de[1:k])))
      res$den_mdr_dd <- sapply(m, function(k) sum(theta_stream[1:k]))
      res$fdr_dd     <- sapply(m, function(k) sum((1-theta_stream[1:k]) * de[1:k]) / max(sum(de[1:k]), 1))
    }
    res
  }
  stopCluster(cl)
  
  # Aggregate results
  for (r in 1:reps) {
    if ("OR" %in% methods_to_run && !is.null(result[[r]]$num_mdr_or)) {
      num_mdr_or[r,] <- result[[r]]$num_mdr_or
      den_mdr_or[r,] <- result[[r]]$den_mdr_or
      fdr_or[r,]     <- result[[r]]$fdr_or
    }
    if ("DD" %in% methods_to_run && !is.null(result[[r]]$num_mdr_dd)) {
      num_mdr_dd[r,] <- result[[r]]$num_mdr_dd
      den_mdr_dd[r,] <- result[[r]]$den_mdr_dd
      fdr_dd[r,]     <- result[[r]]$fdr_dd
    }
  }
  
  # Calculate final metrics here and return them
  final_res <- list()
  if ("OR" %in% methods_to_run) {
    final_res$mdr_or <- colMeans(num_mdr_or) / pmax(colMeans(den_mdr_or), 1)
    final_res$fdr_or <- colMeans(fdr_or)
  }
  if ("DD" %in% methods_to_run) {
    final_res$mdr_dd <- colMeans(num_mdr_dd) / pmax(colMeans(den_mdr_dd), 1)
    final_res$fdr_dd <- colMeans(fdr_dd)
  }
  
  return(final_res)
}

# ------------------------------------------------------------------------------
# 2. Experimental Configuration and Execution
# ------------------------------------------------------------------------------
# --- Fixed Parameters ---
m_seq = seq(100, 1000, 50); pi_v = 0.1; alpha_v = 0.1; reps_v = 500;

# --- Parameter Lists for Sensitivity Analysis ---
n_list   <- c(5, 50, 100, 200, 500, 1000, 2000)   
ini_list <- c(5, 50, 100, 200, 500, 1000, 2000) 
D_list   <- c(10, 50, 100, 200, 500, 1000, 2000)
M_list   <- c(2, 5, 10, 20, 50, 100, 1000)

# --- Default values when a parameter is NOT being varied ---
default_n   <- 500
default_ini <- 1000
default_D   <- 1000
default_M   <- 10

# --- Baseline: OMDRC-OR (Calculated Once) ---
base_res <- run_simulation_ratio(m = m_seq, ini = default_ini, n = default_n, pi = pi_v, alpha = alpha_v, reps = reps_v, D = default_D, M = default_M, methods_to_run = "OR")
mdr_or_vec <- base_res$mdr_or
fdr_or_vec <- base_res$fdr_or

# --- Experiments ---
res_n <- lapply(n_list, function(nv) {
  run_simulation_ratio(m = m_seq, ini = default_ini, n = nv, pi = pi_v, alpha = alpha_v, reps = reps_v, D = default_D, M = default_M, methods_to_run = "DD")
})
res_ini <- lapply(ini_list, function(iv) {
  run_simulation_ratio(m = m_seq, ini = iv, n = default_n, pi = pi_v, alpha = alpha_v, reps = reps_v, D = default_D, M = default_M, methods_to_run = "DD")
})
res_D <- lapply(D_list, function(dv) {
  run_simulation_ratio(m = m_seq, ini = default_ini, n = default_n, pi = pi_v, alpha = alpha_v, reps = reps_v, D = dv, M = default_M, methods_to_run = "DD")
})
res_M <- lapply(M_list, function(mv) {
  run_simulation_ratio(m = m_seq, ini = default_ini, n = default_n, pi = pi_v, alpha = alpha_v, reps = reps_v, D = default_D, M = mv, methods_to_run = "DD")
})

# ------------------------------------------------------------------------------
# 3. Data Processing and Visualization (Simplified Version)
# ------------------------------------------------------------------------------
get_plot_df <- function(res_list, param_name, param_values, or_mdr, or_fdr) {
  df_or <- data.frame(t=m_seq, MDR=or_mdr, FDR=or_fdr, Method="OMDRC.OR")
  
  df_dd <- do.call(rbind, lapply(1:length(param_values), function(i) {
    res <- res_list[[i]]
    data.frame(t = m_seq, 
               MDR = res$mdr_dd,
               FDR = res$fdr_dd,
               Method = paste0("OMDRC.DD (", param_name, "=", param_values[i], ")"))
  }))
  
  df <- rbind(df_or, df_dd)
  df$Method <- factor(df$Method, levels = c("OMDRC.OR", unique(df_dd$Method)))
  return(df)
}

# --- Visualization Aesthetics ---
custom_theme <- theme_bw() + theme(legend.position='bottom', legend.title=element_blank(), panel.border=element_rect(linewidth=0.8))
my_shapes <- c(16, 17, 15, 3, 4, 8, 2, 7)

# --- Plot A: Sensitivity to n (Greens) ---
df_n <- get_plot_df(res_n, "n", n_list, mdr_or_vec, fdr_or_vec)
colors_n <- c("#F8766D", brewer.pal(9, "Greens")[3:9])
p_mdr_n <- ggplot(df_n, aes(t, MDR, color=Method, shape=Method)) + geom_line() + geom_point() + 
  geom_hline(yintercept=alpha_v, linetype="dashed") + scale_color_manual(values=colors_n) + 
  scale_shape_manual(values=my_shapes) + labs(subtitle="(a.1)", y="MDR") + custom_theme
p_fdr_n <- ggplot(df_n, aes(t, FDR, color=Method, shape=Method)) + geom_line() + geom_point() + 
  scale_color_manual(values=colors_n) + scale_shape_manual(values=my_shapes) + labs(subtitle="(a.2)", y="FDR") + custom_theme

# --- Plot B: Sensitivity to K0 (Blues) ---
df_ini <- get_plot_df(res_ini, "K0", ini_list, mdr_or_vec, fdr_or_vec)
colors_ini <- c("#F8766D", brewer.pal(9, "Blues")[3:9])
p_mdr_ini <- ggplot(df_ini, aes(t, MDR, color=Method, shape=Method)) + geom_line() + geom_point() + 
  geom_hline(yintercept=alpha_v, linetype="dashed") + scale_color_manual(values=colors_ini) + 
  scale_shape_manual(values=my_shapes) + labs(subtitle="(b.1)", y="MDR") + custom_theme
p_fdr_ini <- ggplot(df_ini, aes(t, FDR, color=Method, shape=Method)) + geom_line() + geom_point() + 
  scale_color_manual(values=colors_ini) + scale_shape_manual(values=my_shapes) + labs(subtitle="(b.2)", y="FDR") + custom_theme

# --- Plot C: Sensitivity to D (Purples) ---
df_D <- get_plot_df(res_D, "D", D_list, mdr_or_vec, fdr_or_vec)
colors_D <- c("#F8766D", brewer.pal(9, "Purples")[3:9])
p_mdr_D <- ggplot(df_D, aes(t, MDR, color=Method, shape=Method)) + geom_line() + geom_point() + 
  geom_hline(yintercept=alpha_v, linetype="dashed") + scale_color_manual(values=colors_D) + 
  scale_shape_manual(values=my_shapes) + labs(subtitle="(c.1)", y="MDR") + custom_theme
p_fdr_D <- ggplot(df_D, aes(t, FDR, color=Method, shape=Method)) + geom_line() + geom_point() + 
  scale_color_manual(values=colors_D) + scale_shape_manual(values=my_shapes) + labs(subtitle="(c.2)", y="FDR") + custom_theme

# --- Plot D: Sensitivity to M (Oranges/Reds) ---
df_M <- get_plot_df(res_M, "M", M_list, mdr_or_vec, fdr_or_vec)
colors_M <- c("#F8766D", brewer.pal(8, "Oranges")[2:8])
p_mdr_M <- ggplot(df_M, aes(t, MDR, color=Method, shape=Method)) + geom_line() + geom_point() + 
  geom_hline(yintercept=alpha_v, linetype="dashed") + scale_color_manual(values=colors_M) + 
  scale_shape_manual(values=my_shapes[1:(length(M_list)+1)]) + labs(subtitle="(d.1)", y="MDR") + custom_theme
p_fdr_M <- ggplot(df_M, aes(t, FDR, color=Method, shape=Method)) + geom_line() + geom_point() + 
  scale_color_manual(values=colors_M) + scale_shape_manual(values=my_shapes[1:(length(M_list)+1)]) + labs(subtitle="(d.2)", y="FDR") + custom_theme

# ------------------------------------------------------------------------------
# 4. Final Integration and Layout
# ------------------------------------------------------------------------------
fig_n   <- (p_mdr_n | p_fdr_n) + plot_layout(guides="collect") & theme(legend.position="bottom")
fig_ini <- (p_mdr_ini | p_fdr_ini) + plot_layout(guides="collect") & theme(legend.position="bottom")
fig_D   <- (p_mdr_D | p_fdr_D) + plot_layout(guides="collect") & theme(legend.position="bottom")
fig_M   <- (p_mdr_M | p_fdr_M) + plot_layout(guides="collect") & theme(legend.position="bottom")

# Vertical stack of all sensitivity experiments
final_plot <- ggarrange(fig_n, fig_ini, fig_D, fig_M, ncol=1, nrow=4)

print(final_plot)
# ggsave("Sensitivity_Analysis_Final.pdf", final_plot, width=12, height=20)
