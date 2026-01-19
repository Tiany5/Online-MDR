# =========================================================================
# SCRIPT FOR REPRODUCING THE CREDIT CARD FRAUD DETECTION EXPERIMENT
# (VERSION WITH 4 PLOTS: 2024-05-23)
# This version adds a "Cumulative Rejection Rate" plot and confirms
# the FDR calculation is correct. The layout is updated to a 2x2 grid.
# =========================================================================
setwd('/Users/tiany/Desktop/zju/online MDR/code-semi-github')

# --- 1. SETUP: Install and Load Packages ---
required_packages <- c("data.table", "isotree", "ggplot2", "dplyr", "patchwork", "onlineFDR", "scales")

install_if_missing <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}
cat("Loading required packages...\n")
sapply(required_packages, install_if_missing)

# 确保 OMDRC.R 在您的工作目录中
if (file.exists("OMDRC.R")) {
  cat("Loading OMDRC functions from OMDRC.R...\n")
  source("OMDRC.R")
} else {
  stop("OMDRC.R not found in the working directory. Please place it there.")
}


# --- 2. DATA PREPARATION & FEATURE ENGINEERING (ISOLATION FOREST) ---
cat("Loading creditcard.csv dataset...\n")
credit_data <- fread("creditcard.csv")

features_to_use <- c(paste0("V", 1:28), "Time", "Amount")
model_data <- as.data.frame(credit_data[, ..features_to_use])

cat("Building Isolation Forest model...\n")
iso_model <- isolation.forest(data = model_data, ntrees = 100, sample_size = 256, seed = 42)

cat("Predicting anomaly scores...\n")
anomaly_scores <- predict(iso_model, model_data)
credit_data[, anomaly_score := anomaly_scores]
cat("Successfully generated 1D anomaly score stream.\n")


# --- 3. EXPERIMENT SETUP ---
set.seed(123)
ALPHA <- 0.1
K0 <- 2000
N_LABELED <- 100
D_WINDOW <- 1000

normal_tx <- credit_data[Class == 0]
fraud_tx <- credit_data[Class == 1]

z1_samples <- sample_n(fraud_tx, N_LABELED)
z1 <- z1_samples$anomaly_score

online_stream_df <- anti_join(credit_data, z1_samples, by = names(z1_samples)) %>%
  arrange(Time)

z_full_stream <- online_stream_df$anomaly_score
theta_full_stream <- online_stream_df$Class

z_ini <- z_full_stream[1:K0]
z <- z_full_stream[(K0 + 1):length(z_full_stream)]
theta <- theta_full_stream[(K0 + 1):length(theta_full_stream)]

cat(sprintf("Experiment Setup:\n - Target MDR (alpha): %.2f\n - Initial unlabeled (K0): %d\n - Labeled frauds (n): %d\n - Online stream size: %d\n",
            ALPHA, K0, N_LABELED, length(z)))


# --- 4. RUN ALGORITHMS ---
cat("\nRunning OMDRC.DD...\n")
omdrc_results <- OMDRC_DD(z = z, z_ini = z_ini, z1 = z1, alpha = ALPHA, D = D_WINDOW)
omdrc_decisions <- omdrc_results$de

cat("Running Adj-SAFFRON with standard FDR-controlling setup...\n")
p_values_saffron <- sapply(z, function(score) {
  p_val <- (sum(z_ini >= score) + 1) / (length(z_ini) + 1)
  return(p_val)
})

# --- 核心修改：Adj-SAFFRON的决策逻辑 ---
# SAFFRON的R输出的是拒绝的索引，我们需要一个完整的决策向量
saffron_fdr_results <- SAFFRON(p_values_saffron, alpha = ALPHA)
saffron_decisions <- rep(0, length(z))
saffron_decisions[saffron_fdr_results$R] <- 1 # 1 表示拒绝 (发现)

cat(sprintf("OMDRC.DD made %d discoveries. Adj-SAFFRON made %d discoveries.\n",
            sum(omdrc_decisions), sum(saffron_decisions)))


# =========================================================================
# --- 5. EVALUATION & VISUALIZATION (FOUR-PANEL PLOT VERSION) ---
# =========================================================================

# --- 核心修改 1: 更新 calculate_metrics 函数以包含拒绝率 ---
calculate_metrics <- function(decisions, true_labels) {
  len <- length(decisions)
  true_labels <- true_labels[1:len]
  
  cumulative_tp <- cumsum(true_labels * decisions)
  cumulative_fp <- cumsum((1 - true_labels) * decisions)
  cumulative_discoveries <- cumsum(decisions)
  cumulative_signals <- cumsum(true_labels)
  
  # 累计拒绝率 (Cumulative Rejection Rate)
  cumulative_rejection_rate <- cumulative_discoveries / seq_along(decisions)
  
  safe_discoveries <- ifelse(cumulative_discoveries == 0, 1, cumulative_discoveries)
  safe_signals <- ifelse(cumulative_signals == 0, 1, cumulative_signals)
  
  fdr_empirical <- cumulative_fp / safe_discoveries
  missed_signals <- cumulative_signals - cumulative_tp
  mdr_empirical <- missed_signals / safe_signals
  
  return(data.frame(
    mdr = mdr_empirical, 
    fdr = fdr_empirical,
    tp = cumulative_tp,
    rej_rate = cumulative_rejection_rate
  ))
}

omdrc_metrics <- calculate_metrics(omdrc_decisions, theta)
saffron_metrics <- calculate_metrics(saffron_decisions, theta)

# --- 核心修改 2: 在 plot_df 中加入拒绝率数据 ---
plot_df <- data.frame(
  time = 1:length(z),
  mdr_omdrc = omdrc_metrics$mdr,
  mdr_saffron = saffron_metrics$mdr,
  tp_omdrc = omdrc_metrics$tp,
  tp_saffron = saffron_metrics$tp,
  fdr_omdrc = omdrc_metrics$fdr,
  fdr_saffron = saffron_metrics$fdr,
  rej_rate_omdrc = omdrc_metrics$rej_rate,
  rej_rate_saffron = saffron_metrics$rej_rate
)

my_colors <- c("OMDRC.DD" = "dodgerblue", "Adj-SAFFRON" = "darkorange")
my_linetypes <- c("OMDRC.DD" = "solid", "Adj-SAFFRON" = "dashed")
plot_linewidth <- 1.2

# --- 创建四个独立的图 ---
plot_a <- ggplot(plot_df, aes(x = time)) +
  geom_line(aes(y = mdr_omdrc, color = "OMDRC.DD", linetype = "OMDRC.DD"), linewidth = plot_linewidth) +
  geom_line(aes(y = mdr_saffron, color = "Adj-SAFFRON", linetype = "Adj-SAFFRON"), linewidth = plot_linewidth) +
  geom_hline(yintercept = ALPHA, linetype = "dotted", color = "red", linewidth = 1) +
  annotate("text", x = length(z) * 0.8, y = ALPHA + 0.01, label = paste("Target MDR =", ALPHA), color = "red") +
  labs(title = "(a) MDR Control", x = NULL, y = "Empirical MDR") +
  theme_minimal(base_size = 14)+
  coord_cartesian(ylim = c(0, 0.2))

total_fraud_count_in_stream <- sum(theta)
plot_b <- ggplot(plot_df, aes(x = time)) +
  geom_line(aes(y = tp_omdrc, color = "OMDRC.DD", linetype = "OMDRC.DD"), linewidth = plot_linewidth) +
  geom_line(aes(y = tp_saffron, color = "Adj-SAFFRON", linetype = "Adj-SAFFRON"), linewidth = plot_linewidth) +
  geom_hline(yintercept = total_fraud_count_in_stream, linetype = "dotted", color = "darkgreen", linewidth = 1) +
  annotate("text", x = length(z) * 0.8, y = total_fraud_count_in_stream * 0.9, label = paste("Total Frauds:", total_fraud_count_in_stream), color = "darkgreen") +
  labs(title = "(b) Cumulative True Discoveries", x = NULL, y = "True Frauds Detected") +
  theme_minimal(base_size = 14)

plot_c <- ggplot(plot_df, aes(x = time)) +
  geom_line(aes(y = fdr_omdrc, color = "OMDRC.DD", linetype = "OMDRC.DD"), linewidth = plot_linewidth) +
  geom_line(aes(y = fdr_saffron, color = "Adj-SAFFRON", linetype = "Adj-SAFFRON"), linewidth = plot_linewidth) +
  labs(title = "(c) Empirical FDR", x = "Time (t)", y = "Empirical FDR") +
  theme_minimal(base_size = 14)

# --- 核心修改 3: 创建新的拒绝率图 (plot_d) ---
plot_d <- ggplot(plot_df, aes(x = time)) +
  geom_line(aes(y = rej_rate_omdrc, color = "OMDRC.DD", linetype = "OMDRC.DD"), linewidth = plot_linewidth) +
  geom_line(aes(y = rej_rate_saffron, color = "Adj-SAFFRON", linetype = "Adj-SAFFRON"), linewidth = plot_linewidth) +
  labs(title = "(d) Cumulative Rejection Rate", x = "Time (t)", y = "Rejection Rate") +
  scale_y_continuous(labels = scales::percent) + # 将y轴格式化为百分比
  theme_minimal(base_size = 14)

# --- 核心修改 4: 使用 2x2 布局组合所有图 ---
final_plot <- (plot_a + plot_b) / (plot_c + plot_d) +
  plot_layout(guides = 'collect') &
  scale_color_manual(name = "Method", values = my_colors) &
  scale_linetype_manual(name = "Method", values = my_linetypes) &
  theme(legend.position = 'bottom')

cat("\nDisplaying final four-panel comparison plot...\n")
print(final_plot)

