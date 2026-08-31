# ==============================================================================
# SCRIPT: Online MDR Control Simulation - Setting 3 (High-dimensional Signal)
#
# High-dimensional Signal Detection. Observations are d = 20 dimensional with
# correlated features:
#   F0 ~ N_d(0, Sigma),  Sigma_{jk} = 0.5^|j-k|
#   F1 ~ N_d(mu, Sigma), mu = (1.5, 1.5, 1.5, 1.5, 1.5, 0, ..., 0)
# so the signal is concentrated in only 5 of the 20 correlated coordinates.
#
# OMDRC.DD estimates the 20-dimensional density ratio DIRECTLY with an additive
# GAM probabilistic classifier (y ~ s(V1) + ... + s(V20), no interactions): no
# Isolation-Forest dimension reduction and no KDE. Because F0 and F1 share
# Sigma, the true log density-ratio is additive/linear, so the additive GAM is
# well specified and uses all 20 coordinates. The ORACLE procedures (OMDRC.OR,
# FT) instead use the exact 20-dimensional Gaussian Lmdr as the gold standard.
#
# Signal prevalence pi_t follows a multi-regime, NON-MONOTONE "valley" STAIRCASE
# ("regime-shift") path: moderate during the initial calibration batch (0.20),
# then DROPS to a low plateau (0.04), then RISES to a high plateau (0.26) above
# the calibration level -- a shape distinct from Setting 1's S-drop and Setting 2's
# hump. This two-sided departure breaks the fixed-rule baselines even
# with an accurate (GAM) 20-dim score:
#   * Static Lmdr keeps the calibration threshold -> misses the low-pi (0.04)
#     segment -> MDR loses control;
#   * Rolling Top-k keeps a fixed discovery FRACTION calibrated at pi=0.20 ->
#     under-covers the later high-pi (0.26) segment -> MDR loses control;
#   * FT (offline fixed threshold at the nominal budget alpha) provides the
#     offline benchmark;
# while OMDRC.DD (additive GAM on the raw 20-dim features) tracks pi_t and keeps
# MDR <= alpha.
#
# Six procedures: OMDRC.OR, OMDRC.DD, FT, Static Lmdr, Rolling Top-k, Adj-SAFFRON.
# ==============================================================================

suppressPackageStartupMessages({
  library(foreach)
  library(parallel)
  library(doParallel)
  library(ggplot2)
  library(dplyr)
  library(onlineFDR)
  library(patchwork)
})

# Directory that CONTAINS OMDRC.R.
code_dir <- file.path(Sys.getenv("OMDRC_ROOT", unset = getwd()), "code-semi")
setwd(code_dir)
omdrc_file <- file.path(code_dir, "OMDRC.R")
if (!file.exists(omdrc_file)) stop("Cannot find OMDRC.R in: ", code_dir)
source(omdrc_file)

# Multi-regime, NON-MONOTONE "valley" STAIRCASE prevalence path (a regime-shift
# scenario), distinct from Setting 1's S-drop and Setting 2's hump:
# pi sits at a moderate level during the initial calibration batch, DROPS to a
# low plateau, then RISES to a high plateau above the calibration level. This
# two-sided departure breaks BOTH fixed-rule baselines even with an accurate
# (GAM) score: the low segment makes SCT's stale threshold miss low-pi signals,
# and the later high segment makes RTK's stale small discovery fraction
# under-cover high-pi signals -> both lose MDR control, while OMDRC.DD tracks pi_t.
make_pi_staircase <- function(N, levels = c(0.20, 0.04, 0.26),
                              breaks = c(0.40, 0.62), width = 0.06) {
  u <- seq(0, 1, length.out = as.integer(N))
  pi <- rep(levels[1], length(u))
  for (i in seq_along(breaks)) {
    pi <- pi - (levels[i] - levels[i + 1]) * stats::plogis((u - breaks[i]) / width)
  }
  pmin(pmax(pi, min(levels)), max(levels))
}

# ------------------------------------------------------------------------------
# Monte Carlo engine
# ------------------------------------------------------------------------------
exp3 <- function(m, n, ini, alpha, reps, D,
                 d = 20, mu_val = 1.5, n_signal_dims = 5, rho = 0.5,
                 pi_levels = c(0.20, 0.04, 0.26),
                 pi_breaks = c(0.40, 0.62), pi_width = 0.06,
                 pi_bounds = c(0.01, 0.99),
                 topk_window = ini,
                 n_cores = max(1L, parallel::detectCores(logical = TRUE) - 2L)) {

  method_suffixes <- c("or", "dd", "off", "static", "topk", "lord")
  metric_names <- as.vector(outer(c("mdr", "fdr"), method_suffixes,
                                   paste, sep = "."))
  metric_store <- setNames(
    lapply(metric_names, function(.) matrix(0, reps, length(m))),
    metric_names
  )
  static_threshold <- numeric(reps)
  topk_k <- integer(reps)
  topk_rho <- numeric(reps)

  # Correlated Gaussian components.
  idx <- seq_len(d)
  Sigma <- rho^abs(outer(idx, idx, "-"))
  mu <- c(rep(mu_val, n_signal_dims), rep(0, d - n_signal_dims))
  Sigma_chol <- chol(Sigma)   # upper-triangular; rows of Z %*% chol are N(0,Sigma)

  # Since F0 and F1 share Sigma, the exact log density-ratio is LINEAR:
  #   log(f1/f0)(x) = x' Sigma^{-1} mu - 0.5 mu' Sigma^{-1} mu = x'w - cc.
  # This gives the oracle Lmdr without any multivariate-density package.
  w_lin <- solve(Sigma, mu)
  cc_lin <- 0.5 * sum(mu * w_lin)

  # Fixed time-varying prevalence path shared by every replicate.
  N <- max(m) + ini
  pi_path <- make_pi_staircase(N, levels = pi_levels, breaks = pi_breaks,
                               width = pi_width)

  n_cores <- max(1L, min(as.integer(n_cores), reps))
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  on.exit({
    try(stopCluster(cl), silent = TRUE)
    registerDoSEQ()
  }, add = TRUE)

  clusterExport(
    cl,
    c("m", "n", "ini", "alpha", "D", "N", "d", "mu", "Sigma_chol",
      "w_lin", "cc_lin", "pi_path", "pi_bounds", "topk_window",
      "omdrc_file"),
    envir = environment()
  )

  result <- foreach(
    r = seq_len(reps),
    .packages = "onlineFDR"
  ) %dopar% {
    source(omdrc_file)
    set.seed(r)

    # Fast correlated-Gaussian sampler via the shared Cholesky factor.
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

    # Exact 20-dim Gaussian Lmdr (oracle gold standard) via the linear log-ratio.
    R_all <- exp(as.numeric(X %*% w_lin) - cc_lin)   # f1/f0
    lmdr_all <- (pi_path * R_all) / ((1 - pi_path) + pi_path * R_all)

    # ---- Labeled reference samples (F0 / F1) ----
    X0_ref <- rmv(n, rep(0, d))          # labeled null reference (F0)
    X1_ref <- rmv(n, mu)                 # labeled alternative reference (F1)

    # Split the d-dimensional stream into the initial calibration batch and the
    # online stream (no Isolation-Forest score; OMDRC.DD sees the raw features).
    X_ini <- X[seq_len(ini), , drop = FALSE]
    X_str <- X[(ini + 1):N, , drop = FALSE]

    lmdr <- lmdr_all[(ini + 1):N]
    theta_s <- theta[(ini + 1):N]

    # ---- Six procedures ----
    d_or <- OMDRC_OR(lmdr, alpha, w = 100L,
                     x.Lmdr_ini = lmdr_all[seq_len(ini)])
    # OMDRC.DD: additive-GAM density ratio estimated DIRECTLY on the 20-dim
    # features (no Isolation-Forest dimension reduction, no KDE).
    d_dd <- OMDRC_DD(z = X_str, z_ini = X_ini, z0 = X0_ref, z1 = X1_ref,
                     alpha = alpha, ratio_method = "gam",
                     pi_bounds = pi_bounds, em_tol = 1e-6, em_max_iter = 50L,
                     D_mode = "growing", D_beta = 0.6, D_min = 10L, w = 100L)

    # Adj-SAFFRON p-values from the GAM density-ratio score (higher ratio =>
    # more signal-like => high p => not rejected => discovered), using the
    # alternative-reference ECDF. Mirrors Settings 1-2 but with the GAM score.
    ratio_alt_ref <- predict_ratio(d_dd$ratio_model, X1_ref)
    r_sorted <- sort(ratio_alt_ref)
    p_value <- findInterval(d_dd$LR, r_sorted) / length(r_sorted)
    p_value <- pmin(pmax(p_value, 1e-6), 1 - 1e-6)
    # FT: offline fixed threshold on the true Lmdr at the nominal budget alpha.
    d_off <- OMDRC_OFF(lmdr, alpha)
    d_static <- STATIC_LMDR_DD(d_dd$DR, d_dd$DR_ini, alpha)
    d_topk <- ROLLING_TOPK_DD(d_dd$DR, d_dd$DR_ini, alpha,
                              window_size = min(topk_window,
                                                length(d_dd$DR_ini)))
    d_lord <- onlineFDR::SAFFRON(p_value, alpha = alpha)

    decision_list <- list(
      or = d_or$de, dd = d_dd$de, off = d_off$de,
      static = d_static$de, topk = d_topk$de, lord = 1L - d_lord$R
    )

    expected_signals <- cumsum(pi_path[(ini + 1):N])[m]
    compute_metrics <- function(de) {
      cmiss <- cumsum(theta_s * (1 - de))
      cfp <- cumsum((1 - theta_s) * de)
      cde <- cumsum(de)
      list(mdr = cmiss[m] / pmax(expected_signals, .Machine$double.eps),
           fdr = cfp[m] / pmax(cde[m], 1))
    }
    metrics <- lapply(decision_list, compute_metrics)

    list(
      mdr.or = metrics$or$mdr, fdr.or = metrics$or$fdr,
      mdr.dd = metrics$dd$mdr, fdr.dd = metrics$dd$fdr,
      mdr.off = metrics$off$mdr, fdr.off = metrics$off$fdr,
      mdr.static = metrics$static$mdr, fdr.static = metrics$static$fdr,
      mdr.topk = metrics$topk$mdr, fdr.topk = metrics$topk$fdr,
      mdr.lord = metrics$lord$mdr, fdr.lord = metrics$lord$fdr,
      static.threshold = d_static$threshold,
      topk.k = d_topk$k, topk.rho = d_topk$rho
    )
  }

  stopCluster(cl)
  registerDoSEQ()

  for (r in seq_len(reps)) {
    for (nm in metric_names) metric_store[[nm]][r, ] <- result[[r]][[nm]]
    static_threshold[r] <- result[[r]]$static.threshold
    topk_k[r] <- result[[r]]$topk.k
    topk_rho[r] <- result[[r]]$topk.rho
  }

  c(metric_store,
    list(static.threshold = static_threshold,
         topk.k = topk_k, topk.rho = topk_rho))
}

# ------------------------------------------------------------------------------
# Run the experiment
# ------------------------------------------------------------------------------
ini <- 500
m <- seq(from = 100, to = 1000, by = 50)
n <- 500           # paper default: labeled reference size per class
alpha <- 0.1
reps <- 200
D <- 90            # legacy fixed-window fallback (UNUSED: DD now uses the
                   # paper growing window D_t = min{K0+t, max{10, (K0+t)^0.6}})
topk_window <- ini

# Cache of the intermediate numerical results, so plot/label tweaks do not
# require rerunning the experiment. Set reuse_cache <- TRUE to skip the Monte
# Carlo run and reload the saved results before re-plotting.
fig_dir <- file.path(code_dir, "setting1,2,3")
results_file <- file.path(fig_dir, "Setting3_results.rds")
reuse_cache <- FALSE

if (reuse_cache && file.exists(results_file)) {
  cat("Reusing cached numerical results from", results_file, "\n")
  out_exp3 <- readRDS(results_file)
} else {
  cat("Starting Setting 3 Simulation (High-dimensional, direct 20-dim GAM)...\n")
  out_exp3 <- exp3(m = m, n = n, ini = ini, alpha = alpha, reps = reps, D = D,
                   d = 20, mu_val = 1.5, n_signal_dims = 5, rho = 0.5,
                   pi_levels = c(0.20, 0.04, 0.26),
                   pi_breaks = c(0.40, 0.62), pi_width = 0.06,
                   pi_bounds = c(0.01, 0.99), topk_window = topk_window)
  saveRDS(out_exp3, results_file)
  cat("Saved numerical results to", results_file, "\n")
}

cat(sprintf("Static Lmdr threshold: mean %.3f (SD %.3f)\n",
            mean(out_exp3$static.threshold), sd(out_exp3$static.threshold)))
cat(sprintf("Rolling Top-k: mean k %.1f of %d (mean review fraction %.3f)\n",
            mean(out_exp3$topk.k), topk_window, mean(out_exp3$topk.rho)))

# Diagnostic: terminal (t = max m) MDR / FDR per method, flagging failures.
.diag_terminal <- function(out, m, alpha) {
  suf <- c("or", "dd", "off", "static", "topk", "lord")
  lab <- c("OMDRC.OR", "OMDRC.DD", "FT", "SCT",
           "RTK", "Adj-SAFFRON")
  idx <- length(m)
  fdr_dd <- mean(out[["fdr.dd"]][, idx], na.rm = TRUE)
  cat(sprintf("\n--- Terminal metrics at t=%d (alpha=%.2f) ---\n", m[idx], alpha))
  cat(sprintf("%-14s %8s %8s\n", "Method", "MDR", "FDR"))
  for (j in seq_along(suf)) {
    mdr <- mean(out[[paste0("mdr.", suf[j])]][, idx], na.rm = TRUE)
    fdr <- mean(out[[paste0("fdr.", suf[j])]][, idx], na.rm = TRUE)
    tag <- ""
    if (mdr > alpha + 1e-9) tag <- paste0(tag, " [MDR>alpha]")
    if (fdr > fdr_dd + 1e-9) tag <- paste0(tag, " [FDR>DD]")
    cat(sprintf("%-14s %8.3f %8.3f%s\n", lab[j], mdr, fdr, tag))
  }
}
.diag_terminal(out_exp3, m, alpha)

# ------------------------------------------------------------------------------
# Plotting (pointwise 95% Monte Carlo confidence bands)
# ------------------------------------------------------------------------------
my_colors <- c(
  "OMDRC.OR" = "#F8766D",
  "OMDRC.DD" = "#7CAE00",
  "FT" = "#00BFC4",
  "SCT" = "#619CFF",
  "RTK" = "#FF61C3",
  "Adj-SAFFRON" = "#C77CFF"
)
my_shapes <- c(
  "OMDRC.OR" = 16, "OMDRC.DD" = 17, "FT" = 15,
  "SCT" = 18, "RTK" = 8, "Adj-SAFFRON" = 3
)

custom_theme <- theme_bw() +
  theme(
    plot.subtitle = element_text(size = 15, hjust = 0.5, margin = margin(b = 5)),
    legend.title = element_blank(),
    legend.text = element_text(size = 13),
    axis.text = element_text(size = 15, colour = "black"),
    axis.title = element_text(size = 15),
    panel.grid.major = element_line(colour = "grey92", size = 0.4),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.8)
  )

method_suffixes <- c("or", "dd", "static", "topk", "lord")
method_labels <- c("OMDRC.OR", "OMDRC.DD", "SCT",
                   "RTK", "Adj-SAFFRON")

prepare_tidy_data <- function(res_list, metric_prefix, time_points,
                              suffixes, labels) {
  pieces <- lapply(seq_along(suffixes), function(j) {
    mat <- res_list[[paste0(metric_prefix, ".", suffixes[j])]]
    estimate <- colMeans(mat, na.rm = TRUE)
    se <- apply(mat, 2, stats::sd, na.rm = TRUE) / sqrt(nrow(mat))
    data.frame(
      value = estimate,
      lower = pmax(0, estimate - 1.96 * se),
      upper = pmin(1, estimate + 1.96 * se),
      type = labels[j],
      t = time_points
    )
  })
  df <- bind_rows(pieces)
  df$type <- factor(df$type, levels = labels)
  df
}

df_mdr3 <- prepare_tidy_data(out_exp3, "mdr", m, method_suffixes, method_labels)
df_fdr3 <- prepare_tidy_data(out_exp3, "fdr", m, method_suffixes, method_labels)

g3.1 <- ggplot(df_mdr3,
               aes(x = t, y = value, color = type, shape = type, group = type)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = type),
              alpha = 0.16, colour = NA, show.legend = FALSE) +
  geom_line(size = 0.75) +
  geom_point(size = 1.4) +
  geom_hline(yintercept = alpha, linetype = "dashed", size = 0.8,
             colour = "black") +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  scale_shape_manual(values = my_shapes) +
  labs(subtitle = "(c.1)", x = "Time (t)", y = "MDR") +
  custom_theme +
  coord_cartesian(ylim = c(0, NA)) +
  scale_y_continuous(breaks = seq(0, 0.5, by = 0.1))

g3.2 <- ggplot(df_fdr3,
               aes(x = t, y = value, color = type, shape = type, group = type)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = type),
              alpha = 0.16, colour = NA, show.legend = FALSE) +
  geom_line(size = 0.75) +
  geom_point(size = 1.4) +
  labs(subtitle = "(c.2)", x = "Time (t)", y = "FDR") +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  scale_shape_manual(values = my_shapes) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.10))) +
  custom_theme

final_plot3 <- (g3.1 + g3.2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
if (interactive()) print(final_plot3)
ggsave(file.path(fig_dir, "Setting3_fig.pdf"),
       plot = final_plot3, width = 10, height = 5)
ggsave(file.path(fig_dir, "Setting3_fig.png"),
       plot = final_plot3, width = 10, height = 5, dpi = 300)
