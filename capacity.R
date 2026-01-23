# ==============================================================================
# Reproduction Script for Online Misdiscovery Rate Control (OMDRC)
# Submitted to ICML (Anonymous Submission)
# ==============================================================================
# This script performs numerical simulations to evaluate the performance of
# OMDRC-OR versus Fixed Threshold (FT) baselines under different signal patterns.
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Dependencies and Environment Setup
# ------------------------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(patchwork)
library(ggpubr)
library(tidyr)

# NOTE: Set the working directory to the folder containing 'OMDRC.R'
# setwd("path/to/your/code/directory") 
source("OMDRC.R")

# ------------------------------------------------------------------------------
# 2. Core Simulation Function: MDR Performance Evaluation
# ------------------------------------------------------------------------------
# Evaluates Misdiscovery Rate (MDR) as the Ratio of Expectations: E[N] / E[D]
run_or_vs_off_sim <- function(pi_t, alpha, reps, N_total, mu_params) {
  cat(paste0("--- Starting simulation with ", reps, " iterations ---\n"))
  
  # Initialize matrices to store Numerator (uncovered signals) and Denominator (total signals)
  num_or  <- matrix(0, nrow = reps, ncol = N_total)
  den_or  <- matrix(0, nrow = reps, ncol = N_total)
  num_off <- matrix(0, nrow = reps, ncol = N_total)
  den_off <- matrix(0, nrow = reps, ncol = N_total)
  
  for (r in 1:reps) {
    if (r %% 100 == 0) cat(paste("  Iteration", r, "/", reps, "\n"))
    set.seed(r)
    
    # --- Data Generation Process ---
    theta <- rbinom(N_total, 1, pi_t)
    z0 <- rnorm(N_total, 0, 1) # Null distribution
    # Alternative distribution (Mixture Model)
    z_alt <- ifelse(rbinom(N_total, 1, 0.5) == 1, 
                    rnorm(N_total, mu_params$mu1, mu_params$sd1), 
                    rnorm(N_total, mu_params$mu2, mu_params$sd1))
    z_stream <- ifelse(theta == 0, z0, z_alt)
    
    # --- Likelihood Ratio Calculation ---
    f0 <- dnorm(z_stream, 0, 1)
    f1 <- 0.5 * dnorm(z_stream, mu_params$mu1, mu_params$sd1) + 0.5 * dnorm(z_stream, mu_params$mu2, mu_params$sd1)
    lmdr <- (pi_t * f1) / ((1 - pi_t) * f0 + pi_t * f1)
    
    # --- Algorithmic Decision Making ---
    dec_or  <- OMDRC_OR(lmdr, alpha) 
    dec_off <- OMDRC_OFF(lmdr, alpha)
    
    # --- Performance Metrics Accumulation ---
    # Numerator: Undiscovered signals (theta=1 AND decision=0)
    num_or[r, ]  <- cumsum(theta * (1 - dec_or$de))
    num_off[r, ] <- cumsum(theta * (1 - dec_off$de))
    
    # Denominator: Total number of active signals (theta=1)
    total_signals <- cumsum(theta)
    den_or[r, ]  <- total_signals
    den_off[r, ] <- total_signals
  }
  
  # Final MDR Calculation: E[Numerator] / E[Denominator]
  # Using pmax(..., 1) to prevent division by zero during early stages
  e_num_or  <- colMeans(num_or)
  e_den_or  <- pmax(colMeans(den_or), 1)
  
  e_num_off <- colMeans(num_off)
  e_den_off <- pmax(colMeans(den_off), 1)
  
  list(
    mdr_or_mean  = e_num_or / e_den_or, 
    mdr_off_mean = e_num_off / e_den_off
  )
}

# ------------------------------------------------------------------------------
# 3. Global Parameters and Visualization Settings
# ------------------------------------------------------------------------------
SIM_PARAMS <- list(
  alpha_val = 0.1,
  N_total   = 1000,
  reps_val  = 1000,
  mu_params = list(mu1 = 2, mu2 = -3, sd1 = 0.7),
  seed      = 2024
)

# Color palettes for plots
colors_wealth <- c("Cap" = "#466300", "Cost" = "#F8766D", "FT" = "#00BFC4")
colors_mdr    <- c("OMDRC.OR" = "#466300", "FT" = "#00BFC4")
shapes_mdr    <- c("OMDRC.OR" = 16, "FT" = 15)

# Standardized ggplot theme for academic publication
CUSTOM_THEME <- theme_bw() + 
  theme(
    plot.subtitle    = element_text(size = 14, hjust = 0.5, margin = margin(b = 5)),
    legend.title     = element_blank(),
    legend.text      = element_text(size = 13),
    axis.text        = element_text(size = 12, colour = "black"),
    axis.title       = element_text(size = 13),
    panel.grid.major = element_line(colour = "grey92", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(colour = "black", fill = NA, linewidth = 0.8)
  )

# ------------------------------------------------------------------------------
# 4. Evolution of Testing Capacity (Left Panel Visualization)
# ------------------------------------------------------------------------------
create_wealth_plot <- function(pi_t, subtitle_text) {
  set.seed(SIM_PARAMS$seed)
  
  # Generate sample realization for visualization
  theta <- rbinom(SIM_PARAMS$N_total, 1, pi_t)
  z0 <- rnorm(SIM_PARAMS$N_total, 0, 1)
  z_alt <- ifelse(rbinom(SIM_PARAMS$N_total, 1, 0.5) == 1, 
                  rnorm(SIM_PARAMS$N_total, SIM_PARAMS$mu_params$mu1, SIM_PARAMS$mu_params$sd1), 
                  rnorm(SIM_PARAMS$N_total, SIM_PARAMS$mu_params$mu2, SIM_PARAMS$mu_params$sd1))
  z_test <- ifelse(theta == 0, z0, z_alt)
  
  f0 <- dnorm(z_test, 0, 1)
  f1 <- 0.5 * dnorm(z_test, SIM_PARAMS$mu_params$mu1, SIM_PARAMS$mu_params$sd1) + 
        0.5 * dnorm(z_test, SIM_PARAMS$mu_params$mu2, SIM_PARAMS$mu_params$sd1)
  L_oracle <- (pi_t * f1) / ((1 - pi_t) * f0 + pi_t * f1)
  
  res_or <- OMDRC_OR(L_oracle, SIM_PARAMS$alpha_val)
  res_off <- OMDRC_OFF(L_oracle, SIM_PARAMS$alpha_val)
  
  df_plot <- data.frame(
    t = 1:SIM_PARAMS$N_total, 
    Capacity = res_or$Chat, 
    Cost = (1 - SIM_PARAMS$alpha_val) * res_or$Lmdr
  )
  
  ggplot(df_plot, aes(x = t)) +
    geom_area(aes(y = Capacity), fill = "#466300", alpha = 0.1) +
    geom_line(aes(y = Cost, color = "Cost"), alpha = 0.5, linewidth = 0.4) +
    geom_line(aes(y = Capacity, color = "Cap"), linewidth = 0.8) +
    geom_hline(aes(yintercept = (1 - SIM_PARAMS$alpha_val) * res_off$lambda, color = "FT"), 
               linetype = "dashed", linewidth = 0.8) +
    scale_color_manual(values = colors_wealth, 
                       labels = c("Cap" = bquote(C[t]*": Oracle Capacity"), 
                                  "Cost" = bquote((1-alpha)*Lmdr[t]*": Local Cost"), 
                                  "FT" = bquote((1-alpha)*lambda*": FT Baseline"))) +
    labs(x = "Time (t)", y = "Value", subtitle = subtitle_text) +
    coord_cartesian(ylim = c(0, 1)) + CUSTOM_THEME
}

# Scenario (a): Constant Signal Proportion
p1_cap <- create_wealth_plot(rep(0.1, SIM_PARAMS$N_total), 
                             bquote("(a.1) Fixed signal proportion: " * pi[t] == 0.1))

# Scenario (b): Signal Clustering (Transient Burst)
pi_clust_cap <- rep(0.1, SIM_PARAMS$N_total); pi_clust_cap[400:600] <- 0.5
p2_cap <- create_wealth_plot(pi_clust_cap, 
                             bquote("(b.1) Signal clustering: " * pi[t] == 0.5 ~ "for" ~ t %in% group("[", list(400, 600), "]")))

final_output1 <- ggarrange(p1_cap, p2_cap, ncol = 1, nrow = 2, common.legend = TRUE, legend = "bottom")

# ------------------------------------------------------------------------------
# 5. Comparative Performance Analysis (Right Panel Visualization)
# ------------------------------------------------------------------------------
# Run Monte Carlo Simulations for Setting (a) & (b)
pi_fixed_mdr <- rep(0.05, SIM_PARAMS$N_total)
res_fixed <- run_or_vs_off_sim(pi_fixed_mdr, SIM_PARAMS$alpha_val, SIM_PARAMS$reps_val, SIM_PARAMS$N_total, SIM_PARAMS$mu_params)

pi_clust_mdr <- rep(0.05, SIM_PARAMS$N_total); pi_clust_mdr[400:600] <- 0.5
res_clust <- run_or_vs_off_sim(pi_clust_mdr, SIM_PARAMS$alpha_val, SIM_PARAMS$reps_val, SIM_PARAMS$N_total, SIM_PARAMS$mu_params)

# Helper function to prepare data for plotting
prep_mdr_df <- function(res) {
  data.frame(t = 1:SIM_PARAMS$N_total, OMDRC.OR = res$mdr_or_mean, FT = res$mdr_off_mean) %>%
    pivot_longer(cols = -t, names_to = "Method", values_to = "MDR")
}

df_fixed_mdr <- prep_mdr_df(res_fixed)
df_clust_mdr <- prep_mdr_df(res_clust)
point_indices <- seq(50, 1000, 50) # Downsampling points for clarity

p1_mdr <- ggplot(df_fixed_mdr, aes(x = t, y = MDR, color = Method, shape = Method)) +
  geom_line(linewidth = 0.6) + 
  geom_point(data = . %>% filter(t %in% point_indices), size = 1.8) +
  geom_hline(yintercept = SIM_PARAMS$alpha_val, linetype = "dashed") +
  scale_color_manual(values = colors_mdr) + scale_shape_manual(values = shapes_mdr) +
  labs(subtitle = "(a.2) Fixed signal proportion", x = "Time (t)", y = "MDR") +
  coord_cartesian(ylim = c(0, 0.12)) + CUSTOM_THEME

p2_mdr <- ggplot(df_clust_mdr, aes(x = t, y = MDR, color = Method, shape = Method)) +
  geom_line(linewidth = 0.6) + 
  geom_point(data = . %>% filter(t %in% point_indices), size = 1.8) +
  geom_hline(yintercept = SIM_PARAMS$alpha_val, linetype = "dashed") +
  scale_color_manual(values = colors_mdr) + scale_shape_manual(values = shapes_mdr) +
  labs(subtitle = "(b.2) Signal clustering", x = "Time (t)", y = "MDR") +
  coord_cartesian(ylim = c(0, 0.3)) + CUSTOM_THEME

final_output2 <- ggarrange(p1_mdr, p2_mdr, ncol = 1, nrow = 2, common.legend = TRUE, legend = "bottom")

# ------------------------------------------------------------------------------
# 6. Final Figure Integration
# ------------------------------------------------------------------------------
final_plot <- ggarrange(
  final_output1, 
  final_output2, 
  ncol = 2, 
  nrow = 1, 
  widths = c(2, 1) # Maintain 2:1 ratio between Left and Right panels
)

# Output results
print(final_plot)
# ggsave("simulation_results.pdf", final_plot, width = 12, height = 8)
