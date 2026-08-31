# ==============================================================================
# ABLATION: high-dimensional density-ratio estimation for OMDRC.DD in Setting 3.
#
# Same 20-dim data-generating process as setting3.R (F0 = N_d(0, Sigma),
# F1 = N_d(mu, Sigma), Sigma_jk = 0.5^|j-k|, mu = (1.5 x 5, 0 x 15), multi-regime
# staircase pi_t). We compare THREE ways of scoring the stream:
#   OMDRC.OR      - oracle: exact 20-dim linear log density ratio x'Sigma^{-1}mu.
#   DD-iForest    - reduce 20-dim -> 1-dim Isolation-Forest anomaly score, then
#                   fit a GAM density ratio on that 1-D score (ablation arm).
#   DD-GAM20      - NO dimension reduction: additive GAM classifier density ratio
#                   directly on the 20-dim features (y ~ s(V1)+...+s(V20)).
#
# Goal: show DD-GAM20 recovers near-oracle FDR (the true log-ratio is additive/
# linear here, so the additive GAM is well specified), whereas DD-iForest pays a
# large FDR penalty from the lossy 20->1 compression. MDR stays controlled for all.
# ==============================================================================

suppressPackageStartupMessages({
  library(foreach)
  library(parallel)
  library(doParallel)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

code_dir <- file.path(Sys.getenv("OMDRC_ROOT", unset = getwd()), "code-semi")
omdrc_file <- file.path(code_dir, "OMDRC.R")
if (!file.exists(omdrc_file)) stop("Cannot find OMDRC.R in: ", code_dir)
source(omdrc_file)

# ---- Setting-3 configuration (kept identical to setting3.R) ----
ini   <- 500
m     <- seq(from = 100, to = 1000, by = 50)
n     <- 1000            # labeled reference size per class
alpha <- 0.1
reps  <- 100
D     <- 150
d            <- 20
mu_val       <- 1.5
n_signal_dims<- 5
rho          <- 0.5
n_train      <- 2000
pi_levels <- c(0.35, 0.18, 0.04)
pi_breaks <- c(0.30, 0.60)
pi_width  <- 0.03
pi_bounds <- c(0.01, 0.99)
iso_ntrees <- 100L
iso_sample <- 256L

make_pi_staircase <- function(N, levels = c(0.35, 0.18, 0.04),
                              breaks = c(0.30, 0.60), width = 0.03) {
  u <- seq(0, 1, length.out = as.integer(N))
  pi <- rep(levels[1], length(u))
  for (i in seq_along(breaks)) {
    pi <- pi - (levels[i] - levels[i + 1]) * stats::plogis((u - breaks[i]) / width)
  }
  pmin(pmax(pi, min(levels)), max(levels))
}

# Correlated-Gaussian components and exact linear log-ratio (oracle).
idx <- seq_len(d)
Sigma <- rho^abs(outer(idx, idx, "-"))
mu <- c(rep(mu_val, n_signal_dims), rep(0, d - n_signal_dims))
Sigma_chol <- chol(Sigma)
w_lin <- solve(Sigma, mu)
cc_lin <- 0.5 * sum(mu * w_lin)

N <- max(m) + ini
pi_path <- make_pi_staircase(N, levels = pi_levels, breaks = pi_breaks,
                             width = pi_width)

results_file <- file.path(code_dir, "setting1,2,3", "compare_gam_dim_results.rds")
reuse_cache <- TRUE

if (reuse_cache && file.exists(results_file)) {
  result <- readRDS(results_file)
} else {

n_cores <- max(1L, min(parallel::detectCores(logical = TRUE), reps))
cl <- makeCluster(n_cores)
registerDoParallel(cl)
clusterExport(cl, c("m", "n", "ini", "alpha", "D", "N", "d", "mu", "Sigma_chol",
                    "w_lin", "cc_lin", "n_train", "pi_path", "pi_bounds",
                    "iso_ntrees", "iso_sample", "omdrc_file"),
              envir = environment())

result <- foreach(r = seq_len(reps), .packages = "isotree") %dopar% {
  source(omdrc_file)
  set.seed(r)

  rmv <- function(nrow, mean_vec) {
    if (nrow <= 0L) return(matrix(0, 0, d))
    Z <- matrix(rnorm(nrow * d), nrow, d)
    sweep(Z %*% Sigma_chol, 2L, mean_vec, "+")
  }

  # ---- Evaluation stream (d-dimensional, time-varying pi) ----
  theta <- rbinom(N, size = 1, prob = pi_path)
  X <- matrix(0, N, d)
  if (any(theta == 0)) X[theta == 0, ] <- rmv(sum(theta == 0), rep(0, d))
  if (any(theta == 1)) X[theta == 1, ] <- rmv(sum(theta == 1), mu)

  # Exact 20-dim oracle Lmdr via the linear log-ratio.
  R_all <- exp(as.numeric(X %*% w_lin) - cc_lin)
  lmdr_all <- (pi_path * R_all) / ((1 - pi_path) + pi_path * R_all)

  # Labeled references + pure-null iForest training set.
  X0_ref <- rmv(n, rep(0, d))
  X1_ref <- rmv(n, mu)
  X_train <- rmv(n_train, rep(0, d))

  # Split into initial calibration batch and stream.
  X_ini <- X[seq_len(ini), , drop = FALSE]
  X_str <- X[(ini + 1):N, , drop = FALSE]
  lmdr  <- lmdr_all[(ini + 1):N]
  theta_s <- theta[(ini + 1):N]

  # ---- (1) Oracle ----
  d_or <- OMDRC_OR(lmdr, alpha, w = 100L, x.Lmdr_ini = lmdr_ini)

  # ---- (2) DD-iForest: 20-dim -> 1-dim score -> GAM ratio ----
  iso <- isotree::isolation.forest(
    as.data.frame(X_train), ntrees = iso_ntrees,
    sample_size = min(iso_sample, n_train), nthreads = 1L, seed = r
  )
  score <- function(mat) as.numeric(predict(iso, as.data.frame(mat)))
  s_ini  <- score(X_ini)
  s_str  <- score(X_str)
  s0_ref <- score(X0_ref)
  s1_ref <- score(X1_ref)
  d_iso <- OMDRC_DD(z = s_str, z_ini = s_ini, z0 = s0_ref, z1 = s1_ref,
                    alpha = alpha, ratio_method = "gam",
                    pi_bounds = pi_bounds, em_tol = 1e-6, em_max_iter = 50L,
                    D_mode = "growing", D_beta = 0.6, D_min = 10L, w = 100L)

  # ---- (3) DD-GAM20: additive GAM ratio directly on the 20-dim features ----
  d_g20 <- OMDRC_DD(z = X_str, z_ini = X_ini, z0 = X0_ref, z1 = X1_ref,
                    alpha = alpha, ratio_method = "gam",
                    pi_bounds = pi_bounds, em_tol = 1e-6, em_max_iter = 50L,
                    D_mode = "growing", D_beta = 0.6, D_min = 10L, w = 100L)

  expected_signals <- cumsum(pi_path[(ini + 1):N])[m]
  metrics <- function(de) {
    cmiss <- cumsum(theta_s * (1 - de)); cfp <- cumsum((1 - theta_s) * de)
    cde <- cumsum(de)
    list(mdr = cmiss[m] / pmax(expected_signals, .Machine$double.eps),
         fdr = cfp[m] / pmax(cde[m], 1))
  }
  list(or = metrics(d_or$de), iso = metrics(d_iso$de), g20 = metrics(d_g20$de))
}

stopCluster(cl); registerDoSEQ()
saveRDS(result, results_file)
}

# ---- Aggregate terminal (t = 1000) metrics ----
idxT <- length(m)
agg <- function(key, metric) {
  vals <- sapply(result, function(x) x[[key]][[metric]][idxT])
  c(mean = mean(vals, na.rm = TRUE),
    se = sd(vals, na.rm = TRUE) / sqrt(sum(is.finite(vals))))
}
tab <- data.frame(
  Method = c("OMDRC.OR (20-dim exact)",
             "DD-iForest (1-dim score + GAM)",
             "DD-GAM20 (direct 20-dim additive GAM)"),
  MDR = c(agg("or", "mdr")["mean"], agg("iso", "mdr")["mean"], agg("g20", "mdr")["mean"]),
  FDR = c(agg("or", "fdr")["mean"], agg("iso", "fdr")["mean"], agg("g20", "fdr")["mean"]),
  MDR_se = c(agg("or", "mdr")["se"], agg("iso", "mdr")["se"], agg("g20", "mdr")["se"]),
  FDR_se = c(agg("or", "fdr")["se"], agg("iso", "fdr")["se"], agg("g20", "fdr")["se"])
)

cat(sprintf("\n=== Setting 3 density-ratio ablation (t=%d, alpha=%.2f, reps=%d) ===\n",
            m[idxT], alpha, reps))
cat(sprintf("%-40s %8s %8s\n", "Method", "MDR", "FDR"))
for (i in seq_len(nrow(tab))) {
  tag <- if (tab$MDR[i] > alpha + 1e-9) " [MDR>alpha]" else ""
  cat(sprintf("%-40s %6.3f(%.3f) %6.3f(%.3f)%s\n",
              tab$Method[i], tab$MDR[i], tab$MDR_se[i],
              tab$FDR[i], tab$FDR_se[i], tag))
}
cat("\nFDR gap DD-iForest - OMDRC.OR : ",
    sprintf("%.3f", tab$FDR[2] - tab$FDR[1]), "\n")
cat("FDR gap DD-GAM20   - OMDRC.OR : ",
    sprintf("%.3f", tab$FDR[3] - tab$FDR[1]), "\n")

saveRDS(tab, file.path(code_dir, "setting1,2,3", "compare_gam_dim_summary.rds"))

# ---- Timeline figure (MDR / FDR over the full time grid) ----
method_keys   <- c("or", "iso", "g20")
method_labels <- c("OMDRC.OR", "OMDRC.DD (iForest)", "OMDRC.DD (PC-DRE)")
my_colors <- c("OMDRC.OR" = "#984ea3", "OMDRC.DD (iForest)" = "#E69F00",
               "OMDRC.DD (PC-DRE)" = "#009E73")
my_shapes <- c("OMDRC.OR" = 16, "OMDRC.DD (iForest)" = 17, "OMDRC.DD (PC-DRE)" = 15)

build_df <- function(metric) {
  do.call(rbind, lapply(seq_along(method_keys), function(j) {
    mat <- t(sapply(result, function(x) x[[method_keys[j]]][[metric]]))  # reps x |m|
    est <- colMeans(mat, na.rm = TRUE)
    se  <- apply(mat, 2, sd, na.rm = TRUE) / sqrt(nrow(mat))
    data.frame(t = m, value = est,
               lower = pmax(0, est - 1.96 * se),
               upper = pmin(1, est + 1.96 * se),
               type = factor(method_labels[j], levels = method_labels))
  }))
}
df_mdr <- build_df("mdr"); df_fdr <- build_df("fdr")

thm <- theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 18),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 13),
        panel.grid.minor = element_blank())

g_mdr <- ggplot(df_mdr, aes(t, value, color = type, shape = type, group = type)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = type),
              alpha = 0.16, colour = NA, show.legend = FALSE) +
  geom_line(size = 0.8) + geom_point(size = 1.8) +
  geom_hline(yintercept = alpha, linetype = "dashed", colour = "black") +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  scale_shape_manual(values = my_shapes) +
  labs(x = "Time (t)", y = "MDR") + thm

g_fdr <- ggplot(df_fdr, aes(t, value, color = type, shape = type, group = type)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = type),
              alpha = 0.16, colour = NA, show.legend = FALSE) +
  geom_line(size = 0.8) + geom_point(size = 1.8) +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  scale_shape_manual(values = my_shapes) +
  labs(x = "Time (t)", y = "FDR") + thm

g <- (g_mdr + g_fdr) + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

fig_dir <- file.path(code_dir, "setting1,2,3")
ggsave(file.path(fig_dir, "compare_gam_dim_timeline.pdf"), g, width = 10, height = 5)
ggsave(file.path(fig_dir, "compare_gam_dim_timeline.png"), g, width = 10, height = 5,
       dpi = 300)
cat("Saved figure: compare_gam_dim_timeline.{pdf,png}\n")
