# ==============================================================================
# SCRIPT: Online MDR Control Simulation - Setting 2 (Gaussian location shift)
#
# Setting 2 uses a Gaussian location-shift alternative (F0 = N(0,1),
# F1 = N(gs_mean, 1)) under a TIME-VARYING signal prevalence pi_t following the
# non-monotone Gaussian SURGE path (low -> peak -> low, 0.05 -> 0.30 -> 0.05),
# identical to the gaussian_shift branch of
#   code-semi/compare_uLSIF/sample_size_experiment.R
# that produced samplesize_kde_ulsif_2x2_timeline.png. The hump-shaped
# prevalence produces the valley-shaped MDR/FDR curves of that figure. This
# is used by Setting 2.
#
# Six procedures are compared (same set as Setting 1):
#   OMDRC.OR       - oracle online MDR control (knows true Lmdr)
#   OMDRC.DD       - data-driven online MDR control (GAM classifier density ratio)
#   FT             - offline oracle fixed threshold
#   Static Lmdr    - threshold calibrated once on the initial batch (SCT)
#   Rolling Top-k  - causal rolling top-k policy (RTK)
#   Adj-SAFFRON    - online-FDR baseline after swapping null/alternative
#
# Speed-ups: absolute-path sourcing, findInterval ECDF p-values, vectorised
# cumsum metrics, baselines reuse the DD scores, warm-start EM (em_max_iter=50),
# all-but-two logical cores.
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

# Directory that CONTAINS OMDRC.R (with make_smooth_prior_drop, OMDRC_DD taking
# z0 & z1, STATIC_LMDR_DD, ROLLING_TOPK_DD, OMDRC_OR, OMDRC_OFF).
code_dir <- file.path(Sys.getenv("OMDRC_ROOT", unset = getwd()), "code-semi")
setwd(code_dir)
omdrc_file <- file.path(code_dir, "OMDRC.R")
if (!file.exists(omdrc_file)) stop("Cannot find OMDRC.R in: ", code_dir)
source(omdrc_file)

# Gaussian-surge prevalence path (low -> peak -> low), identical to the
# gaussian_shift branch of compare_uLSIF/sample_size_experiment.R that produced
# samplesize_kde_ulsif_2x2_timeline.png. The hump-shaped prevalence yields the
# valley-shaped MDR/FDR curves seen in that figure. This is used by Setting 2.
make_pi_surge <- function(N, pi_low = 0.05, pi_high = 0.30,
                          center = 0.60, width = 0.15) {
  u <- seq(0, 1, length.out = as.integer(N))
  pmin(pmax(pi_low + (pi_high - pi_low) *
              exp(-((u - center)^2) / (2 * width^2)), pi_low), pi_high)
}

# ------------------------------------------------------------------------------
# Monte Carlo engine
# ------------------------------------------------------------------------------
exp2 <- function(m, n, ini, alpha, reps, D,
                 gs_mean = 3,
                 pi_low = 0.05, pi_high = 0.30,
                 surge_center = 0.60, surge_width = 0.15,
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

  # Reference sample for Adj-SAFFRON p-values (alternative F1, sorted once).
  M_ref <- 200000L
  set.seed(202402)
  z_ref_sorted <- sort(rnorm(M_ref, mean = gs_mean, sd = 1))
  n_ref <- length(z_ref_sorted)

  # Fixed time-varying prevalence path shared by every replicate (Gaussian
  # surge: low -> peak -> low).
  N <- max(m) + ini
  pi_path <- make_pi_surge(N, pi_low = pi_low, pi_high = pi_high,
                           center = surge_center, width = surge_width)

  n_cores <- max(1L, min(as.integer(n_cores), reps))
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  on.exit({
    try(stopCluster(cl), silent = TRUE)
    registerDoSEQ()
  }, add = TRUE)

  clusterExport(
    cl,
    c("m", "n", "ini", "alpha", "D", "N", "pi_path", "pi_bounds",
      "topk_window", "z_ref_sorted", "n_ref", "gs_mean", "omdrc_file"),
    envir = environment()
  )

  result <- foreach(r = seq_len(reps), .packages = "onlineFDR") %dopar% {
    source(omdrc_file)
    set.seed(r)

    # ---- Evaluation stream (Gaussian location shift, time-varying pi) ----
    theta <- rbinom(N, size = 1, prob = pi_path)
    z0_stream <- rnorm(N, mean = 0, sd = 1)
    z_alt <- rnorm(N, mean = gs_mean, sd = 1)
    z_all <- ifelse(theta == 0, z0_stream, z_alt)

    f0 <- dnorm(z_all, mean = 0, sd = 1)
    f1 <- dnorm(z_all, mean = gs_mean, sd = 1)
    lmdr_all <- (pi_path * f1) / ((1 - pi_path) * f0 + pi_path * f1)

    z_ini <- z_all[seq_len(ini)]
    z <- z_all[(ini + 1):N]
    lmdr <- lmdr_all[(ini + 1):N]
    theta_s <- theta[(ini + 1):N]

    # Adj-SAFFRON p-values from the alternative-distribution ECDF.
    p_value <- findInterval(z, z_ref_sorted) / n_ref
    p_value <- pmin(pmax(p_value, 1e-6), 1 - 1e-6)

    # ---- Labeled reference samples for the data-driven method ----
    z1_ref <- rnorm(n, mean = gs_mean, sd = 1)
    z0_ref <- rnorm(n, mean = 0, sd = 1)

    # ---- Six procedures ----
    d_or <- OMDRC_OR(lmdr, alpha, w = 100L, x.Lmdr_ini = lmdr_all[seq_len(ini)])
    d_dd <- OMDRC_DD(z = z, z_ini = z_ini, z0 = z0_ref, z1 = z1_ref,
                     alpha = alpha, ratio_method = "gam",
                     pi_bounds = pi_bounds, em_tol = 1e-6, em_max_iter = 50L,
                     D_mode = "growing", D_beta = 0.6, D_min = 10L, w = 100L)
    d_off <- OMDRC_OFF(lmdr, alpha)
    d_static <- STATIC_LMDR_DD(d_dd$DR, d_dd$DR_ini, alpha)
    d_topk <- ROLLING_TOPK_DD(d_dd$DR, d_dd$DR_ini, alpha,
                              window_size = min(topk_window,
                                                length(d_dd$DR_ini)))
    d_lord <- onlineFDR::SAFFRON(p_value, alpha = alpha)

    # Adj-SAFFRON rejects the swapped null, so the anomaly decision is 1 - R.
    decision_list <- list(
      or = d_or$de, dd = d_dd$de, off = d_off$de,
      static = d_static$de, topk = d_topk$de, lord = 1L - d_lord$R
    )

    # Vectorised cumulative MDR / FDR read out at the time grid m.
    expected_signals <- cumsum(pi_path[(ini + 1):N])[m]
    compute_metrics <- function(de) {
      cmiss <- cumsum(theta_s * (1 - de))
      cfp <- cumsum((1 - theta_s) * de)
      cde <- cumsum(de)
      list(
        mdr = cmiss[m] / pmax(expected_signals, .Machine$double.eps),
        fdr = cfp[m] / pmax(cde[m], 1)
      )
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
n <- 500           # paper default: n0 = n1 = 500
alpha <- 0.1
reps <- 200
D <- 80            # legacy fixed-window fallback (UNUSED: DD now uses the
                   # paper growing window D_t = min{K0+t, max{10, (K0+t)^0.6}})
topk_window <- ini

# Gaussian location shift F1 = N(gs_mean, 1) under a non-monotone Gaussian surge
# prevalence (low -> peak -> low, 0.05 -> 0.30 -> 0.05), centred at 0.60 with
# width 0.15 as specified in the manuscript. Under the surge, RTK
# calibrates a small discovery fraction on the low-pi initial batch and then
# under-discovers at the pi peak -> its MDR loses control (~0.33), while
# OMDRC.DD tracks pi_t and keeps MDR <= alpha.
# Cache of the intermediate numerical results, so plot/label tweaks do not
# require rerunning the experiment. Set reuse_cache <- TRUE to skip the Monte
# Carlo run and reload the saved results before re-plotting.
fig_dir <- file.path(code_dir, "setting1,2,3")
results_file <- file.path(fig_dir, "Setting2_results.rds")
reuse_cache <- FALSE

if (reuse_cache && file.exists(results_file)) {
  cat("Reusing cached numerical results from", results_file, "\n")
  out_exp2 <- readRDS(results_file)
} else {
  out_exp2 <- exp2(m = m, n = n, ini = ini, alpha = alpha, reps = reps, D = D,
                   gs_mean = 3,
                   pi_low = 0.05, pi_high = 0.30,
                   surge_center = 0.60, surge_width = 0.15,
                   pi_bounds = c(0.01, 0.99), topk_window = topk_window)
  saveRDS(out_exp2, results_file)
  cat("Saved numerical results to", results_file, "\n")
}

cat(sprintf("Static Lmdr threshold: mean %.3f (SD %.3f)\n",
            mean(out_exp2$static.threshold), sd(out_exp2$static.threshold)))
cat(sprintf("Rolling Top-k: mean k %.1f of %d (mean review fraction %.3f)\n",
            mean(out_exp2$topk.k), topk_window, mean(out_exp2$topk.rho)))

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
    if (suf[j] == "topk" && mdr <= alpha + 1e-9 && fdr > fdr_dd)
      tag <- paste0(tag, " [FDR>DD]")
    cat(sprintf("%-14s %8.3f %8.3f%s\n", lab[j], mdr, fdr, tag))
  }
}
.diag_terminal(out_exp2, m, alpha)

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
    plot.title = element_text(hjust = 0.5, size = 18),
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

df_mdr <- prepare_tidy_data(out_exp2, "mdr", m, method_suffixes, method_labels)
df_fdr <- prepare_tidy_data(out_exp2, "fdr", m, method_suffixes, method_labels)

g2.1_main <- ggplot(df_mdr,
                    aes(x = t, y = value, color = type, shape = type,
                        group = type)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = type),
              alpha = 0.16, colour = NA, show.legend = FALSE) +
  geom_line(size = 0.75) +
  geom_point(size = 1.4) +
  geom_hline(yintercept = alpha, linetype = "dashed", size = 0.8,
             colour = "black") +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  scale_shape_manual(values = my_shapes) +
  labs(subtitle = "(b.1)", x = "Time (t)", y = "MDR") +
  custom_theme +
  scale_y_continuous(breaks = seq(0, 0.5, by = 0.1))

g2.1_final <- g2.1_main

g2.2 <- ggplot(df_fdr,
               aes(x = t, y = value, color = type, shape = type,
                   group = type)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = type),
              alpha = 0.16, colour = NA, show.legend = FALSE) +
  geom_line(size = 0.75) +
  geom_point(size = 1.4) +
  labs(subtitle = "(b.2)", x = "Time (t)", y = "FDR") +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  scale_shape_manual(values = my_shapes) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.10))) +
  custom_theme

final_plot2 <- (g2.1_final + g2.2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

if (interactive()) print(final_plot2)
ggsave(file.path(fig_dir, "Setting2_fig.pdf"),
       plot = final_plot2, width = 10, height = 5)
ggsave(file.path(fig_dir, "Setting2_fig.png"),
       plot = final_plot2, width = 10, height = 5, dpi = 300)
