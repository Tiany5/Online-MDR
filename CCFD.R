# ==============================================================================
# SCRIPT: Real-World Experiment on Credit Card Fraud Detection
# Description: This script evaluates the Online Misdiscovery Rate Control (OMDRC)
#              framework using the Isolation Forest anomaly detector.
# Submission: Anonymous for Review
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Initialization and Dependency Management
# ------------------------------------------------------------------------------
# NOTE: Please set your working directory to the folder containing 'OMDRC.R'
# and the dataset 'creditcard.csv' before running.
# setwd("path/to/reproducible/directory")

required_packages <- c("data.table", "isotree", "ggplot2", "dplyr", "patchwork", "onlineFDR", "scales")

install_if_missing <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

cat("Initializing environment and loading packages...\n")
invisible(sapply(required_packages, install_if_missing))

# Load core algorithmic functions
if (file.exists("OMDRC.R")) {
  cat("Source: OMDRC.R found. Loading core functions...\n")
  source("OMDRC.R")
} else {
  stop("Critical Error: 'OMDRC.R' not found. Ensure it is in the working directory.")
}

# ------------------------------------------------------------------------------
# 2. Data Preparation and Feature Engineering (Isolation Forest)
# ------------------------------------------------------------------------------
cat("Loading Credit Card Fraud dataset...\n")
if (!file.exists("creditcard.csv")) {
  stop("Dataset 'creditcard.csv' not found.")
}
credit_data <- fread("creditcard.csv")

# Select feature columns (V1-V28, Time, and Amount)
features_to_use <- c(paste0("V", 1:28), "Time", "Amount")
model_data <- as.data.frame(credit_data[, ..features_to_use])

cat("Training Isolation Forest (n_trees=100, sample_size=256)...\n")
iso_model <- isolation.forest(data = model_data, ntrees = 100, sample_size = 256, seed = 42)

cat("Generating 1D anomaly score stream...\n")
anomaly_scores <- predict(iso_model, model_data)
credit_data[, anomaly_score := anomaly_scores]

# ------------------------------------------------------------------------------
# 3. Experimental Configuration
# ------------------------------------------------------------------------------
set.seed(123)
ALPHA     <- 0.1     # Target MDR level
K0        <- 2000    # Initial unlabeled samples for density estimation
N_LABELED <- 100     # Prior knowledge: number of labeled fraud samples
D_WINDOW  <- 1000    # Sliding window size for data-driven estimation

normal_tx <- credit_data[Class == 0]
fraud_tx  <- credit_data[Class == 1]

# Sample prior knowledge set (z1)
z1_samples <- sample_n(fraud_tx, N_LABELED)
z1 <- z1_samples$anomaly_score

# Construct the online test stream (sorted by time)
online_stream_df <- anti_join(credit_data, z1_samples, by = names(z1_samples)) %>%
  arrange(Time)

z_full_stream <- online_stream_df$anomaly_score
theta_full_stream <- online_stream_df$Class

# Partition into initialization set and online testing set
z_ini <- z_full_stream[1:K0]
z     <- z_full_stream[(K0 + 1):length(z_full_stream)]
theta <- theta_full_stream[(K0 + 1):length(theta_full_stream)]

cat(sprintf("\nConfiguration Summary:\n - Target MDR (alpha): %.2f\n - Initial Burn-in (K0): %d\n - Labeled Frauds (n): %d\n - Online Stream Size: %d\n",
            ALPHA, K0, N_LABELED, length(z)))

# ------------------------------------------------------------------------------
# 4. Algorithmic Execution
# ------------------------------------------------------------------------------
cat("\nExecuting OMDRC.DD (Data-Driven OMDRC)...\n")
omdrc_results <- OMDRC_DD(z = z, z_ini = z_ini, z1 = z1, alpha = ALPHA, D = D_WINDOW)
omdrc_decisions <- omdrc_results$de

cat("Executing Adj-SAFFRON (Baseline FDR Control)...\n")
# P-value calculation using empirical distribution of initial scores
p_values_saffron <- sapply(z, function(score) {
  p_val <- (sum(z_ini >= score) + 1) / (length(z_ini) + 1)
  return(p_val)
})

# SAFFRON decision logic (map rejection indices to binary vector)
saffron_fdr_results <- SAFFRON(p_values_saffron, alpha = ALPHA)
saffron_decisions   <- rep(0, length(z))
saffron_decisions[saffron_fdr_results$R] <- 1 

cat(sprintf("Results:\n - OMDRC.DD discoveries: %d\n - Adj-SAFFRON discoveries: %d\n",
            sum(omdrc_decisions), sum(saffron_decisions)))

# ------------------------------------------------------------------------------
# 5. Performance Evaluation
# ------------------------------------------------------------------------------
calculate_metrics <- function(decisions, true_labels) {
  len <- length(decisions)
  true_labels <- true_labels[1:len]
  
  cumulative_tp <- cumsum(true_labels * decisions)
  cumulative_fp <- cumsum((1 - true_labels) * decisions)
  cumulative_discoveries <- cumsum(decisions)
  cumulative_signals <- cumsum(true_labels)
  
  # Calculate Cumulative Rejection Rate
  cumulative_rejection_rate <- cumulative_discoveries / seq_along(decisions)
  
  # Numerical stability for denominator
  safe_discoveries <- ifelse(cumulative_discoveries == 0, 1, cumulative_discoveries)
  safe_signals     <- ifelse(cumulative_signals == 0, 1, cumulative_signals)
  
  fdr_empirical  <- cumulative_fp / safe_discoveries
  missed_signals <- cumulative_signals - cumulative_tp
  mdr_empirical  <- missed_signals / safe_signals
  
  return(data.frame(
    mdr = mdr_empirical, 
    fdr = fdr_empirical,
    tp  = cumulative_tp,
    rej_rate = cumulative_rejection_rate
  ))
}

omdrc_metrics   <- calculate_metrics(omdrc_decisions, theta)
saffron_metrics <- calculate_metrics(saffron_decisions, theta)

# ------------------------------------------------------------------------------
# 6. Visualization (Four-Panel Layout)
# ------------------------------------------------------------------------------
plot_df <- data.frame(
  time             = 1:length(z),
  mdr_omdrc        = omdrc_metrics$mdr,
  mdr_saffron      = saffron_metrics$mdr,
  tp_omdrc         = omdrc_metrics$tp,
  tp_saffron       = saffron_metrics$tp,
  fdr_omdrc        = omdrc_metrics$fdr,
  fdr_saffron      = saffron_metrics$fdr,
  rej_rate_omdrc   = omdrc_metrics$rej_rate,
  rej_rate_saffron = saffron_metrics$rej_rate
)

my_colors    <- c("OMDRC.DD" = "dodgerblue", "Adj-SAFFRON" = "darkorange")
my_linetypes <- c("OMDRC.DD" = "solid", "Adj-SAFFRON" = "dashed")
plot_linewidth <- 1.2

# (a) MDR Control Plot
plot_a <- ggplot(plot_df, aes(x = time)) +
  geom_line(aes(y = mdr_omdrc, color = "OMDRC.DD", linetype = "OMDRC.DD"), linewidth = plot_linewidth) +
  geom_line(aes(y = mdr_saffron, color = "Adj-SAFFRON", linetype = "Adj-SAFFRON"), linewidth = plot_linewidth) +
  geom_hline(yintercept = ALPHA, linetype = "dotted", color = "red", linewidth = 1) +
  annotate("text", x = length(z) * 0.8, y = ALPHA + 0.01, label = paste("Target MDR =", ALPHA), color = "red") +
  labs(title = "(a) MDR Control", x = NULL, y = "Empirical MDR") +
  theme_minimal(base_size = 14) +
  coord_cartesian(ylim = c(0, 0.2))

# (b) Cumulative Discoveries Plot
total_fraud_count_in_stream <- sum(theta)
plot_b <- ggplot(plot_df, aes(x = time)) +
  geom_line(aes(y = tp_omdrc, color = "OMDRC.DD", linetype = "OMDRC.DD"), linewidth = plot_linewidth) +
  geom_line(aes(y = tp_saffron, color = "Adj-SAFFRON", linetype = "Adj-SAFFRON"), linewidth = plot_linewidth) +
  geom_hline(yintercept = total_fraud_count_in_stream, linetype = "dotted", color = "darkgreen", linewidth = 1) +
  annotate("text", x = length(z) * 0.8, y = total_fraud_count_in_stream * 0.9, 
           label = paste("Total Frauds:", total_fraud_count_in_stream), color = "darkgreen") +
  labs(title = "(b) Cumulative True Discoveries", x = NULL, y = "Count") +
  theme_minimal(base_size = 14)

# (c) Empirical FDR Plot
plot_c <- ggplot(plot_df, aes(x = time)) +
  geom_line(aes(y = fdr_omdrc, color = "OMDRC.DD", linetype = "OM
