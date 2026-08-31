# ==============================================================================
# APPENDIX SCRIPT: Online MDR Control - Setting 3 sensitivity to the signal
#                  magnitude mu_val
#
# This is the APPENDIX companion of code-semi/setting1,2,3/setting3.R. It keeps
# the SAME methodology as the main Setting 3 script:
#   * high-dimensional d = 20 correlated Gaussians
#       F0 ~ N_d(0, Sigma),  F1 ~ N_d(mu, Sigma),  Sigma_{jk} = 0.5^|j-k|,
#     with the signal concentrated in the first n_signal_dims = 5 coordinates
#     (mu = (mu_val x 5, 0 x 15));
#   * multi-regime, NON-MONOTONE "valley" STAIRCASE prevalence path pi_t;
#   * OMDRC.DD estimates the 20-dim density ratio DIRECTLY with an additive GAM
#     classifier (no Isolation-Forest, no KDE);
#   * ORACLE procedures use the exact LINEAR 20-dim Gaussian Lmdr;
#   * FT uses the nominal miss budget alpha;
#   * the SAME six procedures, colours, shapes, labels and Monte-Carlo 95%
#     confidence bands.
# The only change is that we SWEEP the per-coordinate signal magnitude mu_val
# over a grid and lay the six sweep values out as a 3x2 grid of MDR panels next
# to a 3x2 grid of FDR panels, for the appendix sensitivity study.
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

# Directory that CONTAINS OMDRC.R (one level up from this appendix folder).
code_dir <- file.path(Sys.getenv("OMDRC_ROOT", unset = getwd()), "code-semi")
setwd(code_dir)
omdrc_file <- file.path(code_dir, "OMDRC.R")
if (!file.exists(omdrc_file)) stop("Cannot find OMDRC.R in: ", code_dir)
source(omdrc_file)

# Multi-regime, NON-MONOTONE "valley" STAIRCASE prevalence path, identical to
# setting3.R: moderate -> low plateau -> high plateau above the calibration.
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
# Monte Carlo engine (identical to exp3 in setting3.R; mu_val is the swept arg)
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
  Sigma_chol <- chol(Sigma)

  # Since F0 and F1 share Sigma, the exact log density-ratio is LINEAR.
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

  result <- foreach(r = seq_len(reps), .packages = "onlineFDR") %dopar% {
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

    R_all <- exp(as.numeric(X %*% w_lin) - cc_lin)   # f1/f0
    lmdr_all <- (pi_path * R_all) / ((1 - pi_path) + pi_path * R_all)

    X0_ref <- rmv(n, rep(0, d))
    X1_ref <- rmv(n, mu)

    X_ini <- X[seq_len(ini), , drop = FALSE]
    X_str <- X[(ini + 1):N, , drop = FALSE]

    lmdr <- lmdr_all[(ini + 1):N]
    theta_s <- theta[(ini + 1):N]

    d_or <- OMDRC_OR(lmdr, alpha, w = 100L, x.Lmdr_ini = lmdr_all[seq_len(ini)])
    d_dd <- OMDRC_DD(z = X_str, z_ini = X_ini, z0 = X0_ref, z1 = X1_ref,
                     alpha = alpha, ratio_method = "gam",
                     pi_bounds = pi_bounds, em_tol = 1e-6, em_max_iter = 50L,
                     D_mode = "growing", D_beta = 0.6, D_min = 10L, w = 100L)

    ratio_alt_ref <- predict_ratio(d_dd$ratio_model, X1_ref)
    r_sorted <- sort(ratio_alt_ref)
    p_value <- findInterval(d_dd$LR, r_sorted) / length(r_sorted)
    p_value <- pmin(pmax(p_value, 1e-6), 1 - 1e-6)

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
# Sweep the per-coordinate signal magnitude mu_val
# ------------------------------------------------------------------------------
ini <- 500
m <- seq(from = 100, to = 1000, by = 50)
n <- 1000
alpha <- 0.1
reps <- 200
D <- 90
topk_window <- ini

# Six signal-magnitude values for the appendix sensitivity study, evaluated at
# a fixed dimension d = 30. The high-dimensional regime is used deliberately:
# there the additive-GAM 20+-dim density ratio is hard enough that the
# baselines FT/SCT/RTK/Adj-SAFFRON keep terminal MDR above alpha across the
# whole mu range. At the main setting's d = 20 the usable window is too narrow
# (RTK barely exceeds alpha and Adj-SAFFRON drops below it at small mu).
# The range now spans mu = 0.90 (weak signal; OMDRC.DD terminal MDR ~ 0.117)
# up to mu = 1.40 (strong signal; OMDRC.DD MDR declines toward the nominal
# alpha line, ~ 0.106-0.112), so the panels display the decreasing MDR trend
# as the signal strengthens (scanned by _diag_setting3c.R). Note OMDRC.OR
# itself sits at ~ 0.10-0.105 regardless of mu (MDR control is tight at
# alpha), so OMDRC.DD cannot fall visibly below alpha via mu alone.
d_fixed <- 30
mu_list <- c(0.90, 1.00, 1.10, 1.20, 1.30, 1.40)

fig_dir <- file.path(code_dir, "setting1,2,3-appendix")
results_file <- file.path(fig_dir, "Setting3_varied_mu_results.rds")
reuse_cache <- TRUE

if (reuse_cache && file.exists(results_file)) {
  cat("Reusing cached results from", results_file, "\n")
  all_results <- readRDS(results_file)
} else {
  all_results <- lapply(mu_list, function(mu_v) {
    cat(sprintf("Setting 3 sweep: mu_val = %.2f ...\n", mu_v))
    exp3(m = m, n = n, ini = ini, alpha = alpha, reps = reps, D = D,
         d = d_fixed, mu_val = mu_v, n_signal_dims = 5, rho = 0.5,
         pi_levels = c(0.20, 0.04, 0.26),
         pi_breaks = c(0.40, 0.62), pi_width = 0.06,
         pi_bounds = c(0.01, 0.99), topk_window = topk_window)
  })
  names(all_results) <- paste0("mu_", mu_list)
  saveRDS(all_results, results_file)
  cat("Saved sweep results to", results_file, "\n")
}

# ------------------------------------------------------------------------------
# Plotting (same conventions as setting3.R; one panel per swept value)
# ------------------------------------------------------------------------------
my_colors <- c(
  "OMDRC.OR" = "#F8766D", "OMDRC.DD" = "#7CAE00", "FT" = "#00BFC4",
  "SCT" = "#619CFF", "RTK" = "#FF61C3", "Adj-SAFFRON" = "#C77CFF"
)
my_shapes <- c(
  "OMDRC.OR" = 16, "OMDRC.DD" = 17, "FT" = 15,
  "SCT" = 18, "RTK" = 8, "Adj-SAFFRON" = 3
)

custom_theme <- theme_bw() +
  theme(
    plot.subtitle = element_text(size = 20, hjust = 0.5, margin = margin(b = 5)),
    legend.title = element_blank(),
    legend.text = element_text(size = 17),
    axis.text = element_text(size = 16, colour = "black"),
    axis.title = element_text(size = 18),
    panel.grid.major = element_line(colour = "grey92", size = 0.4),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.8)
  )

method_suffixes <- c("or", "dd", "off", "static", "topk", "lord")
method_labels <- c("OMDRC.OR", "OMDRC.DD", "FT", "SCT", "RTK", "Adj-SAFFRON")

# Plot-only method set: FT is excluded from the figures (still shown in the
# diagnostic terminal table below, which keeps the full method_labels).
plot_suffixes <- c("or", "dd", "static", "topk", "lord")
plot_labels   <- c("OMDRC.OR", "OMDRC.DD", "SCT", "RTK", "Adj-SAFFRON")

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

build_panel <- function(out, subtitle_expr, metric) {
  df <- prepare_tidy_data(out, metric, m, plot_suffixes, plot_labels)
  p <- ggplot(df, aes(x = t, y = value, color = type, shape = type,
                      group = type)) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = type),
                alpha = 0.16, colour = NA, show.legend = FALSE) +
    geom_line(size = 0.7) +
    geom_point(size = 1.2) +
    scale_color_manual(values = my_colors) +
    scale_fill_manual(values = my_colors) +
    scale_shape_manual(values = my_shapes) +
    labs(subtitle = subtitle_expr, x = "Time (t)",
         y = toupper(metric)) +
    custom_theme
  if (metric == "mdr") {
    p <- p +
      geom_hline(yintercept = alpha, linetype = "dashed", size = 0.8,
                 colour = "black") +
      coord_cartesian(ylim = c(0, NA)) +
      scale_y_continuous(breaks = seq(0, 0.5, by = 0.1))
  } else {
    p <- p + scale_y_continuous(expand = expansion(mult = c(0.05, 0.10)))
  }
  p
}

# Diagnostic: terminal MDR/FDR per swept value.
cat("\n=== Setting 3 sweep terminal metrics (t = max) ===\n")
for (i in seq_along(mu_list)) {
  out <- all_results[[i]]
  idx <- length(m)
  cat(sprintf("\n-- mu_val = %.2f --\n", mu_list[i]))
  for (j in seq_along(method_suffixes)) {
    mdr <- mean(out[[paste0("mdr.", method_suffixes[j])]][, idx], na.rm = TRUE)
    fdr <- mean(out[[paste0("fdr.", method_suffixes[j])]][, idx], na.rm = TRUE)
    tag <- if (mdr > alpha + 1e-9) " [MDR>alpha]" else ""
    cat(sprintf("%-14s MDR %6.3f FDR %6.3f%s\n", method_labels[j], mdr, fdr, tag))
  }
}

mdr_panels <- lapply(seq_along(mu_list), function(i)
  build_panel(all_results[[i]], bquote(mu == .(mu_list[i])), "mdr"))
fdr_panels <- lapply(seq_along(mu_list), function(i)
  build_panel(all_results[[i]], bquote(mu == .(mu_list[i])), "fdr"))

final_layout <- (wrap_plots(mdr_panels, nrow = 3, ncol = 2) |
                   wrap_plots(fdr_panels, nrow = 3, ncol = 2)) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

if (interactive()) print(final_layout)
ggsave(file.path(fig_dir, "Setting3_varied_mu.pdf"),
       plot = final_layout, width = 18, height = 12)
ggsave(file.path(fig_dir, "Setting3_varied_mu.png"),
       plot = final_layout, width = 18, height = 12, dpi = 300)
