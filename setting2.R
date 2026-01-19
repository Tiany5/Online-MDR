#' Simulation and Visualization for Online MDR Control (Setting 2: Skewed Signals)
#' 
#' This script reproduces Figure (b.1) and (b.2) using the Ratio of Expectations 
#' formula for MDR and FDR.
#' 
#' Author: [Anonymous]
#' Date: January 2025

# --- 1. Load Required Libraries ---
library(Matrix)
library(REBayes)
library(foreach)
library(doParallel)
library(ggplot2)
library(dplyr)
library(kedd)
library(onlineFDR)
library(patchwork)

# --- 2. Core Simulation Function (Ratio of Expectations Logic) ---

run_simulation_setting2 <- function(m, n, ini, alpha, reps, D) {
  # Matrices to store Numerators and Denominators for each replication
  n_mdr_or <- matrix(0, reps, length(m)); d_mdr_or <- matrix(0, reps, length(m))
  n_mdr_dd <- matrix(0, reps, length(m)); d_mdr_dd <- matrix(0, reps, length(m))
  n_mdr_off <- matrix(0, reps, length(m)); d_mdr_off <- matrix(0, reps, length(m))
  n_mdr_saffron <- matrix(0, reps, length(m)); d_mdr_saffron <- matrix(0, reps, length(m))
  
  n_fdr_or <- matrix(0, reps, length(m)); d_fdr_or <- matrix(0, reps, length(m))
  n_fdr_dd <- matrix(0, reps, length(m)); d_fdr_dd <- matrix(0, reps, length(m))
  n_fdr_off <- matrix(0, reps, length(m)); d_fdr_off <- matrix(0, reps, length(m))
  n_fdr_saffron <- matrix(0, reps, length(m)); d_fdr_saffron <- matrix(0, reps, length(m))
  
  ## --- Optimization: Pre-generate Gamma reference for p-values ---
  shape_alt <- 3; scale_alt <- 1
  M_ref <- 200000
  set.seed(202402)
  z_ref_global <- rgamma(M_ref, shape = shape_alt, scale = scale_alt)
  
  cores <- min(parallel::detectCores() - 1, 10)
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  on.exit({ stopCluster(cl); registerDoSEQ() }) 
  
  clusterExport(cl, c("m", "n", "ini", "reps", "alpha", "D", "z_ref_global", "shape_alt", "scale_alt"), envir = environment())
  
  results_list <- foreach(r = 1:reps, .packages = c("kedd", "onlineFDR")) %dopar% {
    source('OMDRC.R')
    set.seed(r)
    
    # Data Generation (Skewed: Exp vs Gamma)
    N <- max(m) + ini
    pi_t <- rep(0.1, N)
    theta <- rbinom(N, 1, pi_t)
    
    z0 <- rexp(N, rate = 1) # Null ~ Exp(1)
    z_alt <- rgamma(N, shape = shape_alt, scale = scale_alt) # Alternative ~ Gamma(3, 1)
    z <- ifelse(theta == 0, z0, z_alt)
    
    # Oracle values
    f0 <- dexp(z, rate = 1)
    f1 <- dgamma(z, shape = shape_alt, scale = scale_alt)
    z_lmdr_all <- (pi_t * f1) / ((1 - pi_t) * f0 + pi_t * f1)
    
    z_ini <- z[1:ini]
    z_stream <- z[(ini + 1):N]
    z_lmdr <- z_lmdr_all[(ini + 1):N]
    theta_stream <- theta[(ini + 1):N]
    
    # P-values for SAFFRON (ECDF based on Gamma reference)
    p_values <- pmin(pmax(ecdf(z_ref_global)(z_stream), 0.00001), 0.99999)
    
    # Labeled set z1 for DD method
    z1 <- rgamma(n, shape = shape_alt, scale = scale_alt)
    
    # Calculate Numerators and Denominators
    get_n_mdr <- function(th, de) sapply(m, function(k) sum(th[1:k] * (1 - de[1:k])))
    get_d_mdr <- function(th)     sapply(m, function(k) max(sum(th[1:k]), 1))
    get_n_fdr <- function(th, de) sapply(m, function(k) sum((1 - th[1:k]) * de[1:k]))
    get_d_fdr <- function(de)     sapply(m, function(k) max(sum(de[1:k]), 1))
    
    list(
      n_or = get_n_mdr(theta_stream, dec_or <- OMDRC_OR(z_lmdr, alpha)$de),   d_or = get_d_mdr(theta_stream),
      nf_or = get_n_fdr(theta_stream, dec_or),                              df_or = get_d_fdr(dec_or),
      
      n_dd = get_n_mdr(theta_stream, dec_dd <- OMDRC_DD(z_stream, z_ini, z1, alpha, D)$de), d_dd = get_d_mdr(theta_stream),
      nf_dd = get_n_fdr(theta_stream, dec_dd),                              df_dd = get_d_fdr(dec_dd),
      
      n_off = get_n_mdr(theta_stream, dec_off <- OMDRC_OFF(z_lmdr, alpha)$de), d_off = get_d_mdr(theta_stream),
      nf_off = get_n_fdr(theta_stream, dec_off),                              df_off = get_d_fdr(dec_off),
      
      n_saffron = get_n_mdr(theta_stream, dec_saffron <- (1 - onlineFDR::SAFFRON(p_values, alpha = alpha)$R)), d_saffron = get_d_mdr(theta_stream),
      nf_saffron = get_n_fdr(theta_stream, dec_saffron),                             df_saffron = get_d_fdr(dec_saffron)
    )
  }
  
  for (r in 1:reps) {
    n_mdr_or[r,] <- results_list[[r]]$n_or; d_mdr_or[r,] <- results_list[[r]]$d_or
    n_fdr_or[r,] <- results_list[[r]]$nf_or; d_fdr_or[r,] <- results_list[[r]]$df_or
    n_mdr_dd[r,] <- results_list[[r]]$n_dd; d_mdr_dd[r,] <- results_list[[r]]$d_dd
    n_fdr_dd[r,] <- results_list[[r]]$nf_dd; d_fdr_dd[r,] <- results_list[[r]]$df_dd
    n_mdr_off[r,] <- results_list[[r]]$n_off; d_mdr_off[r,] <- results_list[[r]]$d_off
    n_fdr_off[r,] <- results_list[[r]]$nf_off; d_fdr_off[r,] <- results_list[[r]]$df_off
    n_mdr_saffron[r,] <- results_list[[r]]$n_saffron; d_mdr_saffron[r,] <- results_list[[r]]$d_saffron
    n_fdr_saffron[r,] <- results_list[[r]]$nf_saffron; d_fdr_saffron[r,] <- results_list[[r]]$df_saffron
  }
  
  return(list(
    mdr.or = colMeans(n_mdr_or) / colMeans(d_mdr_or),
    fdr.or = colMeans(n_fdr_or / d_fdr_or),
    mdr.dd = colMeans(n_mdr_dd) / colMeans(d_mdr_dd),
    fdr.dd = colMeans(n_fdr_dd / d_fdr_dd),
    mdr.off = colMeans(n_mdr_off) / colMeans(d_mdr_off),
    fdr.off = colMeans(n_fdr_off/ d_fdr_off),
    mdr.saffron = colMeans(n_mdr_saffron) / colMeans(d_mdr_saffron),
    fdr.saffron = colMeans(n_fdr_saffron / d_fdr_saffron)
  ))
}

# --- 3. Robust Data Preparation Helper ---
prepare_tidy_data <- function(res_list, metric_prefix, time_points, type_labels) {
  cols <- paste0(metric_prefix, c(".or", ".dd", ".off", ".saffron"))
  data_list <- lapply(res_list[cols], function(x) {
    if (is.matrix(x)) return(colMeans(x, na.rm = TRUE))
    return(x)
  })
  data.frame(
    value = unlist(data_list),
    type = factor(rep(type_labels, each = length(time_points)), levels = type_labels),
    t = rep(time_points, times = length(type_labels))
  )
}

# --- 4. Plotting Configuration ---
my_colors <- c("OMDRC.OR"="#F8766D", "OMDRC.DD"="#7CAE00", "FT"="#00BFC4", "Adj-SAFFRON"="#C77CFF")
my_shapes <- c("OMDRC.OR"=16, "OMDRC.DD"=17, "FT"=15, "Adj-SAFFRON"=3)
custom_theme <- theme_bw() + 
  theme(
    plot.subtitle = element_text(size = 15, hjust = 0.5, margin = margin(b = 5)),
    legend.title = element_blank(),
    legend.text = element_text(size = 15),
    axis.text = element_text(size = 15, colour = "black"),
    axis.title = element_text(size = 15),
    panel.grid.major = element_line(colour = "grey92", linewidth = 0.4),
    panel.border = element_rect(colour = "black", fill=NA, linewidth=0.8)
  )

# --- 5. Run and Visualize ---
params <- list(ini = 500, m = seq(100, 1000, 50), n = 100, alpha = 0.1, reps = 1000, D = 1000)
out_exp2 <- run_simulation_setting2(params$m, params$n, params$ini, params$alpha, params$reps, params$D)

method_labels <- c('OMDRC.OR', 'OMDRC.DD', 'FT', 'Adj-SAFFRON')
df_mdr2 <- prepare_tidy_data(out_exp2, "mdr", params$m, method_labels)
df_fdr2 <- prepare_tidy_data(out_exp2, "fdr", params$m, method_labels)
df_mdr2_zoom <- df_mdr2 %>% filter(type != "Adj-SAFFRON")

# (b.1) Main MDR Plot
g2.1_main <- ggplot(df_mdr2, aes(x = t, y = value, color = type, shape = type)) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.5) +
  geom_hline(yintercept = 0.1, linetype = 'dashed', linewidth = 0.8, col = 'black') + 
  scale_color_manual(values = my_colors) + scale_shape_manual(values = my_shapes) +
  labs(subtitle = "(b.1)", x = "Time (t)", y = "MDR") + custom_theme +
  scale_y_continuous(breaks = seq(0, 0.5, by = 0.1))

# Inset Calculation
inset_limits2 <- df_mdr2 %>% filter(type %in% c("OMDRC.OR", "FT")) %>%
  summarise(ymin = min(value), ymax = max(value))
y_range2 <- (inset_limits2$ymax - inset_limits2$ymin) * 0.1
y_breaks2 <- seq(round(inset_limits2$ymin, 3), round(inset_limits2$ymax, 3), length.out = 2)

# Inset Plot
g2.1_inset <- ggplot(df_mdr2_zoom, aes(x = t, y = value, color = type, shape = type)) +
  geom_line(linewidth = 0.7, show.legend = FALSE) + geom_point(size = 1, show.legend = FALSE) +
  geom_hline(yintercept = 0.1, linetype = 'dashed', linewidth = 0.7, col = 'black') +
  coord_cartesian(ylim = c(inset_limits2$ymin - y_range2, inset_limits2$ymax + y_range2)) + 
  scale_y_continuous(breaks = y_breaks2) + scale_color_manual(values = my_colors) +
  scale_shape_manual(values = my_shapes) + scale_x_continuous(breaks = c(250, 750)) + 
  theme_void() + theme(panel.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
                       axis.text = element_text(size = 9, colour = "black"))

g2.1_final <- g2_main <- g2.1_main + inset_element(g2.1_inset, 0.4, 0.2, 0.98, 0.6)

# (b.2) FDR Plot
g2.2 <- g2_fdr <- ggplot(df_fdr2, aes(x = t, y = value, color = type, shape = type)) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.5) +
  labs(subtitle = "(b.2)", x = "Time (t)", y = "FDR") +
  scale_color_manual(values = my_colors) + scale_shape_manual(values = my_shapes) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) + custom_theme

# --- 6. Final Combined Plot & Save ---
if (exists("g1.1_final") && exists("g1.2")) {
  combined_final_plot <- (
    (g1.1_final | g1.2) /
      (g2.1_final | g2.2)
  ) + 
    plot_layout(guides = 'collect') & 
    theme(legend.position = 'bottom', legend.text = element_text(size=15))
  
  print(combined_final_plot)
} else {
  final_plot2 <- (g2.1_final + g2.2) + plot_layout(guides = 'collect') & theme(legend.position = 'bottom')
  print(final_plot2)
}

# Save Workspace
save.image(file = "Setting1,2_pi=0.1.RData")
