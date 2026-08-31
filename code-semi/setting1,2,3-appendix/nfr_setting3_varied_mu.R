# ==============================================================================
# APPENDIX SCRIPT (NFR VARIANT): Online MDR Control - Setting 3 sensitivity to
#                                the signal magnitude mu_val, WITH the NFR
#                                (anti-free-riding) barrier on BOTH OMDRC arms
#                                -- credited ledger, D_NFR = 100.
#
# This is the NFR companion of setting3_varied_mu.R in this appendix folder. The
# DGP is kept VERBATIM from that script (Setting 3 was NOT re-tuned; only
# Setting 1 received updated parameters):
#   * high-dimensional correlated Gaussians F0 ~ N_d(0, Sigma),
#     F1 ~ N_d(mu, Sigma), Sigma_{jk} = 0.5^|j-k|, signal in the first
#     n_signal_dims = 5 coordinates (mu = (mu_val x 5, 0 x (d-5))), evaluated at
#     a FIXED dimension d = 30;
#   * multi-regime, NON-MONOTONE "valley" STAIRCASE prevalence pi_t;
#   * OMDRC.DD estimates the d-dim density ratio DIRECTLY with an additive GAM
#     classifier, with pi_bounds = c(0.01, 0.99), D_mode = "growing",
#     D_beta = 0.6, D_min = 10L and labeled reference size n = 1000 -- the
#     appendix folder's own values, kept so the NFR panels stay directly
#     comparable with setting3_varied_mu.R;
#   * ORACLE procedures use the exact LINEAR d-dim Gaussian Lmdr;
#   * FT uses the nominal miss budget alpha;
#   * the SAME six procedures, colours, shapes, labels and Monte-Carlo 95%
#     confidence bands; we SWEEP the per-coordinate signal magnitude mu_val over
#     the same grid as setting3_varied_mu.R (3x2 MDR panels next to 3x2 FDR
#     panels). The high-dimensional d = 30 regime is used deliberately: there
#     the additive-GAM density ratio is hard enough that the baselines
#     FT/SCT/RTK/Adj-SAFFRON keep terminal MDR above alpha across the whole mu
#     range.
#
# NOTE on the two different "d"s in this script: the FIXED feature dimension is
# `d_fixed = 30`, whereas the NFR barrier window is the SEPARATE constant
# `D_NFR = 100` (the trailing-window length of the anti-free-riding gate). They
# are unrelated.
#
# NFR change (identical to nfr_setting3_fig.R): both OMDRC arms carry the
# CREDITED local barrier at w = 100 (Algorithms 1 and 2):
#   OMDRC.OR      -> OMDRC.OR+NFR   = OMDRC_OR_NFR(true q,      alpha, D_NFR)
#   OMDRC.DD      -> OMDRC.DD+NFR   = OMDRC_OR_NFR(estimated q, alpha, D_NFR)
# "credited" = a barrier-forced rejection credits +alpha*q to the live ledger
# (_nfr_core.R::OMDRC_OR_NFR). SCT / RTK / Adj-SAFFRON / FT are untouched and
# still read the SAME estimated Lmdr score path, so the comparison is unchanged.
#
# WARNING (documented, expected): under the CREDITED ledger the DD arm can have
# alpha-breaching PATHS even where its mean MDR stays below alpha -- the forced
# rejections remove exactly the high-q_hat points whose skip-set bias
# compensated the GAM's inflation of the total signal mass. The oracle arm
# OR+NFR remains controlled pathwise. Read the ledger table printed below.
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

# Shared NFR core (gamma path + the three barrier ledgers).
nfr_dir <- file.path(code_dir, "nofreeride-setting1,2,3")
nfr_core_file <- file.path(nfr_dir, "_nfr_core.R")
if (!file.exists(nfr_core_file)) stop("Cannot find _nfr_core.R in: ", nfr_dir)
source(nfr_core_file)

# Barrier window, identical to nfr_setting{1,2,3}_fig.R. NOT the feature
# dimension (that is the swept `d` below).
D_NFR <- 100L

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
# Monte Carlo engine (setting3_varied_mu.R engine; mu_val is the swept arg;
# both OMDRC arms carry the credited NFR barrier)
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
      "omdrc_file", "nfr_core_file", "D_NFR"),
    envir = environment()
  )

  result <- foreach(r = seq_len(reps), .packages = "onlineFDR") %dopar% {
    source(omdrc_file)
    source(nfr_core_file)
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
    lmdr_ini <- lmdr_all[seq_len(ini)]   # fills the gamma window during warm-up
    theta_s <- theta[(ini + 1):N]

    # ---- Six procedures (both OMDRC arms carry the credited NFR barrier) ----
    d_or <- OMDRC_OR_NFR(lmdr, alpha, d = D_NFR, x.Lmdr_ini = lmdr_ini)
    # OMDRC_DD is run to produce the ESTIMATED score path (DR / DR_ini). We
    # replay the same barrier explicitly below so that both OMDRC arms expose
    # identical diagnostics; the score path does not depend on the decisions.
    d_dd <- OMDRC_DD(z = X_str, z_ini = X_ini, z0 = X0_ref, z1 = X1_ref,
                     alpha = alpha, ratio_method = "gam",
                     pi_bounds = pi_bounds, em_tol = 1e-6, em_max_iter = 50L,
                     D_mode = "growing", D_beta = 0.6, D_min = 10L)
    d_dd_nfr <- OMDRC_OR_NFR(d_dd$DR, alpha, d = D_NFR,
                             x.Lmdr_ini = d_dd$DR_ini)

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
      or = d_or$de, dd = d_dd_nfr$de, off = d_off$de,
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
      topk.k = d_topk$k, topk.rho = d_topk$rho,
      # Ledger check: MDR_q = sum_skip q / sum q in the TRUE q. For the ORACLE
      # arm the capacity identity pins this at <= alpha on EVERY path; for the
      # DD arm the ledger runs on q_hat, so the identity is broken and this is
      # exactly where the estimation error shows up.
      mdrq.or = sum(lmdr * (1 - d_or$de)) / sum(lmdr),
      mdrq.dd = sum(lmdr * (1 - d_dd_nfr$de)) / sum(lmdr),
      blocked.or = sum(d_or$blocked), blocked.dd = sum(d_dd_nfr$blocked)
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
         topk.k = topk_k, topk.rho = topk_rho,
         mdrq.or = vapply(result, `[[`, numeric(1), "mdrq.or"),
         mdrq.dd = vapply(result, `[[`, numeric(1), "mdrq.dd"),
         blocked.or = vapply(result, function(x) as.numeric(x$blocked.or), numeric(1)),
         blocked.dd = vapply(result, function(x) as.numeric(x$blocked.dd), numeric(1)),
         d_nfr = D_NFR))
}

# ------------------------------------------------------------------------------
# Sweep the per-coordinate signal magnitude mu_val (dimension fixed at d = 30)
# ------------------------------------------------------------------------------
ini <- 500
m <- seq(from = 100, to = 1000, by = 50)
n <- 1000
alpha <- 0.1
reps <- 200
D <- 90
topk_window <- ini

# Six signal-magnitude values, SAME grid as setting3_varied_mu.R, evaluated at
# the FIXED dimension d = 30. The range spans mu = 0.90 (weak signal) up to
# mu = 1.40 (strong signal), so the panels display the decreasing MDR trend as
# the signal strengthens.
d_fixed <- 30
mu_list <- c(0.90, 1.00, 1.10, 1.20, 1.30, 1.40)

fig_dir <- file.path(code_dir, "setting1,2,3-appendix")
results_file <- file.path(fig_dir, "NFR_Setting3_varied_mu_results.rds")
reuse_cache <- identical(Sys.getenv("NFR_REUSE_CACHE"), "1")

cat(sprintf("Starting Setting 3 mu_val-sweep + NFR credited barrier (D_NFR=%d)...\n",
            D_NFR))
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
# Plotting (same conventions as setting3_varied_mu.R; one panel per swept value)
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
method_labels <- c("OMDRC.OR+NFR", "OMDRC.DD+NFR", "FT", "SCT",
                   "RTK", "Adj-SAFFRON")

# Plot-only method set: FT is excluded from the figures (still shown in the
# diagnostic terminal table below). Legend labels stay PLAIN OMDRC.OR /
# OMDRC.DD on purpose (same convention as nfr_setting{1,2,3}_fig.R); the
# console diagnostic keeps the +NFR labels so logs stay unambiguous.
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
    scale_color_manual(values = my_colors, limits = plot_labels, breaks = plot_labels) +
    scale_fill_manual(values = my_colors, limits = plot_labels, breaks = plot_labels) +
    scale_shape_manual(values = my_shapes, limits = plot_labels, breaks = plot_labels) +
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

# Diagnostic: terminal MDR/FDR per swept value + barrier ledger check.
cat("\n=== Setting 3 sweep terminal metrics (t = max, NFR arms) ===\n")
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
  cat(sprintf("  ledger (D_NFR=%d, credited): OR mean MDR_q %.4f max %.4f [%s]",
              out$d_nfr, mean(out$mdrq.or), max(out$mdrq.or),
              if (max(out$mdrq.or) <= alpha + 1e-9) "PASS" else "FAIL"))
  cat(sprintf(" | DD mean MDR_q %.4f max %.4f [%s]\n",
              mean(out$mdrq.dd), max(out$mdrq.dd),
              if (max(out$mdrq.dd) <= alpha + 1e-9) "PASS" else "FAIL"))
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
ggsave(file.path(fig_dir, "NFR_Setting3_varied_mu.pdf"),
       plot = final_layout, width = 18, height = 12)
ggsave(file.path(fig_dir, "NFR_Setting3_varied_mu.png"),
       plot = final_layout, width = 18, height = 12, dpi = 300)
cat("Figures written to", file.path(fig_dir, "NFR_Setting3_varied_mu.pdf/png"), "\n")
