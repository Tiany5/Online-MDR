# ==============================================================================
# SCRIPT: Online Misdiscovery Rate (MDR) Control - Simulation & Visualization
# Description: This script reproduces Figures (a.1) and (a.2) for Setting 1 
#              (Fixed Signal Proportion). It evaluates performance using the 
#              "Ratio of Expectations" formula.
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Environment Setup
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

# NOTE: Set working directory to the path containing 'OMDRC.R'
# setwd("path/to/reproducible/code") 

# ------------------------------------------------------------------------------
# 2. Core Simulation Engine (Ratio of Expectations Logic)
# ------------------------------------------------------------------------------
run_simulation_setting1 <- function(m, n, ini, alpha, reps, D) {
  # Initialization: Matrices to store Numerators and Denominators across replications
  # Goal: Compute MDR = E[Uncovered Signals] / E[Total Active Signals]
  n_mdr_or <- matrix(0, reps, length(m)); d_mdr_or <- matrix(0, reps, length(m))
  n_mdr_dd <- matrix(0, reps, length(m)); d_mdr_dd <- matrix(0, reps, length(m))
  n_mdr_off <- matrix(0, reps, length(m)); d_mdr_off <- matrix(0, reps, length(m))
  n_mdr_saffron <- matrix(0, reps, length(m)); d_mdr_saffron <- matrix(0, reps, length(m))
  
  # Matrices for FDR (False Discovery Rate) tracking
  n_fdr_or <- matrix(0, reps, length(m)); d_fdr_or <- matrix(0, reps, length(m))
  n_fdr_dd <- matrix(0, reps, length(m)); d_fdr_dd <- matrix(0, reps, length(m))
  n_fdr_off <- matrix(0, reps, length(m)); d_fdr_off <- matrix(0, reps, length(m))
  n_fdr_saffron <- matrix(0, reps, length(m)); d_fdr_saffron <- matrix(0, reps, length(m))
  
  # Data-generating parameters
  mu1 <- 2; mu2 <- -3; sd1 <- 0.7; M_ref <- 200000; pi=0.08
  set.seed(202401)
  z_ref_global <- c(rnorm(M_ref/2, mu1, sd1), rnorm(M_ref/2, mu2, sd1))
  
  # Parallel computing setup
  cores <- min(parallel::detectCores() - 1, 10)
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  on.exit({ stopCluster(cl); registerDoSEQ() }) 
  
  clusterExport(cl, c("m", "n", "ini", "reps", "alpha", "D", "z_ref_global", "mu1", "mu2", "sd1"), envir = environment())
  
  results_list <- foreach(r = 1:reps, .packages = c("kedd", "onlineFDR")) %dopar% {
    source('OMDRC.R')
    set.seed(r)
    
    # --- Data Generation Process ---
    N <- max(m) + ini
    theta <- rbinom(N, 1, pi) # 10% signal proportion
    z <- ifelse(theta == 0, rnorm(N, 0, 1), 
                rnorm(N, ifelse(rbinom(N, 1, 0.5) == 1, mu1, mu2), sd1))
    
    # --- Likelihood Ratio Calculation (Oracle) ---
    f0 <- dnorm(z, 0, 1)
    f1 <- 0.5 * dnorm(z, mu1, sd1) + 0.5 * dnorm(z, mu2, sd1)
    z_lmdr_all <- (pi * f1) / ((1-pi) * f0 + pi * f1)
    
    # Stream partitioning
    z_stream     <- z[(ini + 1):N]
    z_lmdr       <- z_lmdr_all[(ini + 1):N]
    theta_stream <- theta[(ini + 1):N]
    
    # P-value derivation for baseline comparison (SAFFRON)
    p_values <- pmin(pmax(ecdf(abs(z_ref_global))(abs(z_stream)), 0), 1)
    z1 <- rnorm(n, ifelse(rbinom(n, 1, 0.5) == 1, mu1, mu2), sd1)
    
    # --- Execution of Algorithms ---
    dec_or      <- OMDRC_OR(z_lmdr, alpha)
    dec_dd      <- OMDRC_DD(z_stream, z[1:ini], z1, alpha, D)
    dec_off     <- OMDRC_OFF(z_lmdr, alpha)
    dec_saffron <- onlineFDR::SAFFRON(p_values, alpha = alpha)
    
    # --- Metrics Computation: Numerators and Denominators ---
    get_n_mdr <- function(th, de) sapply(m, function(k) sum(th[1:k] * (1 - de[1:k])))
    get_d_mdr <- function(th)     sapply(m, function(k) max(sum(th[1:k]), 1))
    get_n_fdr <- function(th, de) sapply(m, function(k) sum((1 - th[1:k]) * de[1:k]))
    get_d_fdr <- function(de)     sapply(m, function(k) max(sum(de[1:k]), 1))
    
    list(
      n_or  = get_n_mdr(theta_stream, dec_or$de),    d_or = get_d_mdr(theta_stream),
      nf_or = get_n_fdr(theta_stream, dec_or$de),   df_or = get_d_fdr(dec_or$de),
      n_dd  = get_n_mdr(theta_stream, dec_dd$de),    d_dd = get_d_mdr(theta_stream),
      nf_dd = get_n_fdr(theta_stream, dec_dd$de),   df_dd = get_d_fdr(dec_dd$de),
      n_off = get_n_mdr(theta_stream, dec_off$de),   d_off = get_d_mdr(theta_stream),
      nf_off= get_n_fdr(theta_stream, dec_off$de),  df_off = get_d_fdr(dec_off$de),
      n_saffron  = get_n_mdr(theta_stream, 1 - dec_saffron$R), d_saffron = get_d_mdr(theta_stream),
      nf_saffron = get_n_fdr(theta_stream, 1 - dec_saffron$R), df_saffron = get_d_fdr(1 - dec_saffron$R)
    )
  }
  
  # Accumulating results from all workers
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
  
  # Final performance calculation: Ratio of Expectations
  return(list(
    mdr.or      = colMeans(n_mdr_or) / colMeans(d_mdr_or),
    fdr.or      = colMeans(n_fdr_or / d_fdr_or),
    mdr.dd      = colMeans(n_mdr_dd) / colMeans(d_mdr_dd),
    fdr.dd      = colMeans(n_fdr_dd / d_fdr_dd),
    mdr.off     = colMeans(n_mdr_off) / colMeans(d_mdr_off),
    fdr.off     = colMeans(n_fdr_off / d_fdr_off),
    mdr.saffron = colMeans(n_mdr_saffron) / colMeans(d_mdr_saffron),
    fdr.saffron = colMeans(n_fdr_saffron / d_fdr_saffron)
  ))
}

# ------------------------------------------------------------------------------
# 3. Data Transformation Helper
# ------------------------------------------------------------------------------
prepare_tidy_data <- function(res_list, metric_prefix, time_points, type_labels) {
  cols <- paste0(metric_prefix, c(".or", ".dd", ".off", ".saffron"))
  
  data_list <- lapply(res_list[cols], function(x) {
    if (is.matrix(x)) return(colMeans(x, na.rm = TRUE))
    return(x) 
  })
  
  data.frame(
    value = unlist(data_list),
    type  = factor(rep(type_labels, each = length(time_points)), levels = type_labels),
    t     = rep(time_points, times = length(type_labels))
  )
}

# ------------------------------------------------------------------------------
# 4. Visualization Settings
# ------------------------------------------------------------------------------
my_colors <- c("OMDRC.OR"="#F8766D", "OMDRC.DD"="#7CAE00", "FT"="#00BFC4", "Adj-SAFFRON"="#C77CFF")
my_shapes <- c("OMDRC.OR"=16, "OMDRC.DD"=17, "FT"=15, "Adj-SAFFRON"=3)

custom_theme <- theme_bw() + 
  theme(
    plot.subtitle    = element_text(size = 15, hjust = 0.5, margin = margin(b = 5)),
    legend.title     = element_blank(),
    legend.text      = element_text(size = 15),
    axis.text        = element_text(size = 15, colour = "black"),
    axis.title       = element_text(size = 15),
    panel.grid.major = element_line(colour = "grey92", linewidth = 0.4),
    panel.border     = element_rect(colour = "black", fill = NA, linewidth = 0.8)
  )

# ------------------------------------------------------------------------------
# 5. Simulation Execution and Plot Generation
# ------------------------------------------------------------------------------
params <- list(ini = 500, m = seq(100, 1000, 50), n = 100, alpha = 0.1, reps = 1000, D = 1000)

cat("Starting Setting 1 Simulation...\n")
out_exp1 <- run_simulation_setting1(params$m, params$n, params$ini, params$alpha, params$reps, params$D)

method_labels <- c('OMDRC.OR', 'OMDRC.DD', 'FT', 'Adj-SAFFRON')
df_mdr <- prepare_tidy_data(out_exp1, "mdr", params$m, method_labels)
df_fdr <- prepare_tidy_data(out_exp1, "fdr", params$m, method_labels)
df_mdr_zoom <- df_mdr %>% filter(type != "Adj-SAFFRON")

# (a.1) Main MDR Plot
g1.1_main <- ggplot(df_mdr, aes(x = t, y = value, color = type, shape = type)) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.5) +
  geom_hline(yintercept = 0.1, linetype = 'dashed', linewidth = 0.8, col = 'black') + 
  scale_color_manual(values = my_colors) + scale_shape_manual(values = my_shapes) +
  labs(subtitle = "(a.1)", x = "Time (t)", y = "MDR") + custom_theme +
  coord_cartesian(ylim = c(0, NA)) +
  scale_y_continuous(breaks = seq(0, 0.5, by = 0.1))

# --- Inset Calculation ---
inset_limits <- df_mdr %>% filter(type %in% c("OMDRC.OR", "FT")) %>%
  summarise(ymin = min(value), ymax = max(value))
y_range  <- (inset_limits$ymax - inset_limits$ymin) * 0.1
y_breaks <- seq(round(inset_limits$ymin, 3), round(inset_limits$ymax, 3), length.out = 2)

# --- Inset Plot Creation ---
g1.1_inset_clean <- ggplot(df_mdr_zoom, aes(x = t, y = value, color = type, shape = type)) +
  geom_line(linewidth = 0.7, show.legend = FALSE) + geom_point(size = 1, show.legend = FALSE) +
  geom_hline(yintercept = 0.1, linetype = 'dashed', linewidth = 0.7, col = 'black') +
  coord_cartesian(ylim = c(inset_limits$ymin - y_range, inset_limits$ymax + y_range)) + 
  scale_y_continuous(breaks = y_breaks) + scale_color_manual(values = my_colors) +
  scale_shape_manual(values = my_shapes) + scale_x_continuous(breaks = c(250, 750)) + 
  theme_void() + theme(panel.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
                       axis.text = element_text(size = 9, colour = "black"))

g1.1_final <- g1.1_main + inset_element(g1.1_inset_clean, 0.4, 0.4, 0.98, 0.8)

# (a.2) FDR Comparison Plot
g1.2 <- ggplot(df_fdr, aes(x = t, y = value, color = type, shape = type)) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.5) +
  labs(subtitle = "(a.2)", x = "Time (t)", y = "FDR") +
  scale_color_manual(values = my_colors) + scale_shape_manual(values = my_shapes) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) + custom_theme

# --- Final Layout Consolidation ---
final_plot1 <- (g1.1_final + g1.2) + plot_layout(guides = 'collect') & theme(legend.position = 'bottom')
print(final_plot1)

# Export workspace for reproducibility
save.image(file = "Setting1_simulation_results.RData")
