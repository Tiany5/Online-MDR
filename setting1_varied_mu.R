# ==============================================================================
# SCRIPT: Sensitivity Analysis for Signal Strength (mu) - Setting 1
# Description: This script evaluates the robustness of the Online MDR Control 
#              framework under varying signal strengths (mu). Performance 
#              metrics are computed using the Ratio of Expectations formula.
# Submission: Anonymous for Review (ICML)
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Initialization and Environment Setup
# ------------------------------------------------------------------------------
library(Matrix)
library(REBayes)
library(foreach)
library(doParallel)
library(ggplot2)
library(dplyr)
library(kedd)
library(onlineFDR)
library(patchwork)

# NOTE: Please set your working directory to the folder containing 'OMDRC.R'
# setwd("path/to/reproducible/code")
source('OMDRC.R')

# ------------------------------------------------------------------------------
# 2. Simulation Logic for Varying Signal Strength (mu)
# ------------------------------------------------------------------------------
run_sim_mu <- function(mu_val, m, n, ini, alpha, reps, D) {
  
  # Initialize storage matrices for metrics across replications
  # MDR = E[Numerator] / E[Denominator]
  n_mdr_or      <- matrix(0, reps, length(m)); d_mdr_or      <- matrix(0, reps, length(m))
  n_mdr_dd      <- matrix(0, reps, length(m)); d_mdr_dd      <- matrix(0, reps, length(m))
  n_mdr_off     <- matrix(0, reps, length(m)); d_mdr_off     <- matrix(0, reps, length(m))
  n_mdr_saffron <- matrix(0, reps, length(m)); d_mdr_saffron <- matrix(0, reps, length(m))
  
  # FDR = E[V/R]
  n_fdr_or      <- matrix(0, reps, length(m)); d_fdr_or      <- matrix(0, reps, length(m))
  n_fdr_dd      <- matrix(0, reps, length(m)); d_fdr_dd      <- matrix(0, reps, length(m))
  n_fdr_off     <- matrix(0, reps, length(m)); d_fdr_off     <- matrix(0, reps, length(m))
  n_fdr_saffron <- matrix(0, reps, length(m)); d_fdr_saffron <- matrix(0, reps, length(m))
  
  # Distribution parameters
  mu1 <- mu_val; mu2 <- -3; sd1 <- 0.7
  M_ref <- 200000; set.seed(202401)
  z_ref_global <- c(rnorm(M_ref/2, mu1, sd1), rnorm(M_ref/2, mu2, sd1))
  
  # Parallel Backend Configuration
  cl <- makeCluster(min(parallel::detectCores() - 1, 10))
  registerDoParallel(cl)
  on.exit({ stopCluster(cl); registerDoSEQ() }) 
  clusterExport(cl, c("m", "n", "ini", "reps", "alpha", "D", "z_ref_global", "mu1", "mu2", "sd1"), envir = environment())
  
  # Parallel Monte Carlo Simulations
  results_list <- foreach(r = 1:reps, .packages = c("kedd", "onlineFDR")) %dopar% {
    source('OMDRC.R')
    set.seed(r)
    
    # --- Data Generation ---
    N <- max(m) + ini; theta <- rbinom(N, 1, 0.1) 
    z <- ifelse(theta == 0, rnorm(N, 0, 1), 
                rnorm(N, ifelse(rbinom(N, 1, 0.5) == 1, mu1, mu2), sd1))
    
    # --- Likelihood Ratio (Local MDR) Calculation ---
    f0 <- dnorm(z, 0, 1); f1 <- 0.5 * dnorm(z, mu1, sd1) + 0.5 * dnorm(z, mu2, sd1)
    z_lmdr_all <- (0.1 * f1) / (0.9 * f0 + 0.1 * f1)
    
    # --- Stream Partitioning ---
    z_stream      <- z[(ini + 1):N]
    z_lmdr        <- z_lmdr_all[(ini + 1):N]
    theta_stream  <- theta[(ini + 1):N]
    
    # P-value calculation for baseline (SAFFRON)
    p_values <- pmin(pmax(ecdf(abs(z_ref_global))(abs(z_stream)), 0), 1)
    z1 <- rnorm(n, ifelse(rbinom(n, 1, 0.5) == 1, mu1, mu2), sd1)
    
    # --- Algorithm Execution ---
    dec_or      <- OMDRC_OR(z_lmdr, alpha)$de
    dec_dd      <- OMDRC_DD(z_stream, z[1:ini], z1, alpha, D)$de
    dec_off     <- OMDRC_OFF(z_lmdr, alpha)$de
    dec_saffron <- 1 - onlineFDR::SAFFRON(p_values, alpha = alpha)$R
    
    # --- Performance Metrics Logic (Numerator/Denominator extraction) ---
    get_n_mdr <- function(th, de) sapply(m, function(k) sum(th[1:k] * (1 - de[1:k])))
    get_d_mdr <- function(th)     sapply(m, function(k) max(sum(th[1:k]), 1))
    get_n_fdr <- function(th, de) sapply(m, function(k) sum((1 - th[1:k]) * de[1:k]))
    get_d_fdr <- function(de)     sapply(m, function(k) max(sum(de[1:k]), 1))
    
    list(n_or = get_n_mdr(theta_stream, dec_or),   d_or = get_d_mdr(theta_stream),
         nf_or = get_n_fdr(theta_stream, dec_or),  df_or = get_d_fdr(dec_or),
         n_dd = get_n_mdr(theta_stream, dec_dd),   d_dd = get_d_mdr(theta_stream),
         nf_dd = get_n_fdr(theta_stream, dec_dd),  df_dd = get_d_fdr(dec_dd),
         n_off = get_n_mdr(theta_stream, dec_off), d_off = get_d_mdr(theta_stream),
         nf_off = get_n_fdr(theta_stream, dec_off), df_off = get_d_fdr(dec_off),
         n_saffron = get_n_mdr(theta_stream, dec_saffron), d_saffron = get_d_mdr(theta_stream),
         nf_saffron = get_n_fdr(theta_stream, dec_saffron), df_saffron = get_d_fdr(dec_saffron))
  }
  
  # Aggregating results from parallel workers
  for (r in 1:reps) {
    n_mdr_or[r,]  <- results_list[[r]]$n_or;  d_mdr_or[r,]  <- results_list[[r]]$d_or
    n_fdr_or[r,]  <- results_list[[r]]$nf_or; d_fdr_or[r,]  <- results_list[[r]]$df_or
    n_mdr_dd[r,]  <- results_list[[r]]$n_dd;  d_mdr_dd[r,]  <- results_list[[r]]$d_dd
    n_fdr_dd[r,]  <- results_list[[r]]$nf_dd; d_fdr_dd[r,]  <- results_list[[r]]$df_dd
    n_mdr_off[r,] <- results_list[[r]]$n_off; d_mdr_off[r,] <- results_list[[r]]$d_off
    n_fdr_off[r,] <- results_list[[r]]$nf_off;d_fdr_off[r,] <- results_list[[r]]$df_off
    n_mdr_saffron[r,] <- results_list[[r]]$n_saffron; d_mdr_saffron[r,] <- results_list[[r]]$d_saffron
    n_fdr_saffron[r,] <- results_list[[r]]$nf_saffron;d_fdr_saffron[r,] <- results_list[[r]]$df_saffron
  }
  
  # Metric consolidation
  return(list(
    mdr = rbind(colMeans(n_mdr_or) / colMeans(d_mdr_or),
                colMeans(n_mdr_dd) / colMeans(d_mdr_dd),
                colMeans(n_mdr_off) / colMeans(d_mdr_off),
                colMeans(n_mdr_saffron) / colMeans(d_mdr_saffron)),
    fdr = rbind(colMeans(n_fdr_or / d_fdr_or),
                colMeans(n_fdr_dd / d_fdr_dd),
                colMeans(n_fdr_off / d_fdr_off),
                colMeans(n_fdr_saffron / d_fdr_saffron))
  ))
}

# ------------------------------------------------------------------------------
# 3. Visualization Configuration
# ------------------------------------------------------------------------------
my_colors <- c("OMDRC.OR"="#F8766D", "OMDRC.DD"="#7CAE00", "FT"="#00BFC4", "Adj-SAFFRON"="#C77CFF")
my_shapes <- c("OMDRC.OR"=16, "OMDRC.DD"=17, "FT"=15, "Adj-SAFFRON"=3)

custom_theme <- theme_bw() + 
  theme(plot.subtitle    = element_text(size = 15, hjust = 0.5, margin = margin(b = 5)),
        legend.title     = element_blank(), 
        legend.text      = element_text(size = 15),
        axis.text        = element_text(size = 15, colour = "black"), 
        axis.title       = element_text(size = 15),
        panel.grid.major = element_line(colour = "grey92", linewidth = 0.4),
        panel.border     = element_rect(colour = "black", fill = NA, linewidth = 0.8))

# ------------------------------------------------------------------------------
# 4. Visualization Helper (Standardized Panels)
# ------------------------------------------------------------------------------
plot_panel <- function(res, mu_val, type = "MDR") {
  method_labels <- c('OMDRC.OR', 'OMDRC.DD', 'FT', 'Adj-SAFFRON')
  mat <- if(type == "MDR") res$mdr else res$fdr
  
  # Map matrix rows to their respective labels for tidy visualization
  df <- data.frame(
    t     = rep(params$m, times = 4),
    value = as.vector(t(mat)), 
    type  = factor(rep(method_labels, each = length(params$m)), levels = method_labels)
  )
  
  p_main <- ggplot(df, aes(x = t, y = value, color = type, shape = type)) +
    geom_line(linewidth = 0.7) + geom_point(size = 1.5) +
    scale_color_manual(values = my_colors) + scale_shape_manual(values = my_shapes) +
    labs(subtitle = bquote(mu == .(mu_val)), x = "Time (t)", y = type) + custom_theme
  
  if (type == "MDR") {
    p_main <- p_main + geom_hline(yintercept = params$alpha, linetype = 'dashed', color = 'black') +
      scale_y_continuous(breaks = seq(0, 0.5, by = 0.1))
    
    # --- Inset detailed view ---
    df_zoom <- df %>% filter(type %in% c("OMDRC.OR", "FT"))
    limits  <- df_zoom %>% summarise(ymin = min(value), ymax = max(value))
    y_pad   <- (limits$ymax - limits$ymin) * 0.15
    y_breaks <- seq(round(limits$ymin, 3), round(limits$ymax, 3), length.out = 2)
    
    p_inset <- ggplot(df %>% filter(type != "Adj-SAFFRON"), aes(x = t, y = value, color = type, shape = type)) +
      geom_line(linewidth = 0.7, show.legend = FALSE) + geom_point(size = 0.8, show.legend = FALSE) +
      geom_hline(yintercept = params$alpha, linetype = 'dashed', linewidth = 0.6, color = 'black') +
      coord_cartesian(ylim = c(limits$ymin - y_pad, limits$ymax + y_pad)) + 
      scale_y_continuous(breaks = y_breaks) + scale_x_continuous(breaks = c(250, 750)) + 
      scale_color_manual(values = my_colors) + scale_shape_manual(values = my_shapes) +
      theme_void() + theme(panel.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
                           axis.text = element_text(size = 9, color = "black"))
    
    p_final <- p_main + inset_element(p_inset, 0.35, 0.3, 0.98, 0.75)
  } else {
    p_final <- p_main + scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
  }
  return(p_final)
}

# ------------------------------------------------------------------------------
# 5. Execution and Plot Integration
# ------------------------------------------------------------------------------
mu_list <- c(1.25, 1.5, 1.75, 2, 2.25, 2.5)
params <- list(ini = 500, m = seq(100, 1000, 50), n = 100, alpha = 0.1, reps = 1000, D = 1000)

cat("Running Sensitivity Analysis for mu sequence...\n")
all_results_mu <- lapply(mu_list, function(p) run_sim_mu(p, params$m, params$n, params$ini, params$alpha, params$reps, params$D))

# Generate Panels for MDR and FDR
mdr_plots <- lapply(1:6, function(i) plot_panel(all_results_mu[[i]], mu_list[i], "MDR"))
fdr_plots <- lapply(1:6, function(i) plot_panel(all_results_mu[[i]], mu_list[i], "FDR"))

# Consolidation of results into final figure
final_layout <- (wrap_plots(mdr_plots, ncol=2) | wrap_plots(fdr_plots, ncol=2)) + 
  plot_layout(guides = 'collect') & theme(legend.position = 'bottom')

print(final_layout)

# Save the workspace for reproduction
save.image(file = "Sensitivity_Analysis_Mu_Results.RData")
