# ==============================================================================
# APPENDIX SCRIPT: Online MDR Control - Setting 2 sensitivity to the surge PEAK
#                  height pi_high
#
# Companion of code-semi/setting1,2,3/setting2.R (same methodology):
#   * Gaussian location shift F0 = N(0,1), F1 = N(gs_mean = 2.4, 1);
#   * non-monotone Gaussian SURGE prevalence pi_t (low -> peak -> low); here we
#     SWEEP the peak height pi_high;
#   * OMDRC.DD uses the GAM classifier density ratio (z0 & z1 references);
#   * FT uses the nominal miss budget alpha;
#   * the SAME six procedures, colours, shapes, labels and Monte-Carlo 95% bands.
#
# Range note (requirement: baselines must "fail"): across the swept pi_high the
# terminal MDR of FT, RTK and Adj-SAFFRON stays above alpha, while SCT
# over-discovers under the surge so its MDR stays <= alpha but its later-segment
# FDR exceeds OMDRC.DD's (verified by _diag_setting1.R). gs_mean is fixed at a
# small value (2.4) for the same reason (large mu makes the baselines succeed).
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

# Gaussian-surge prevalence path (low -> peak -> low), identical to setting1.R.
make_pi_surge <- function(N, pi_low = 0.05, pi_high = 0.30,
                          center = 0.60, width = 0.15) {
  u <- seq(0, 1, length.out = as.integer(N))
  pmin(pmax(pi_low + (pi_high - pi_low) *
              exp(-((u - center)^2) / (2 * width^2)), pi_low), pi_high)
}

# ------------------------------------------------------------------------------
# Monte Carlo engine (identical to exp1 in setting1.R; pi_high is the swept arg)
# ------------------------------------------------------------------------------
exp1 <- function(m, n, ini, alpha, reps, D,
                 gs_mean = 2.4,
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
  set.seed(202401)
  z_ref_sorted <- sort(rnorm(M_ref, mean = gs_mean, sd = 1))
  n_ref <- length(z_ref_sorted)

  # Fixed time-varying prevalence path shared by every replicate.
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
      "topk_window", "z_ref_sorted", "n_ref", "gs_mean",
      "omdrc_file"),
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

    p_value <- findInterval(z, z_ref_sorted) / n_ref
    p_value <- pmin(pmax(p_value, 1e-6), 1 - 1e-6)

    z1_ref <- rnorm(n, mean = gs_mean, sd = 1)
    z0_ref <- rnorm(n, mean = 0, sd = 1)

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

    decision_list <- list(
      or = d_or$de, dd = d_dd$de, off = d_off$de,
      static = d_static$de, topk = d_topk$de, lord = 1L - d_lord$R
    )

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
# Sweep the surge peak height pi_high (gs_mean fixed at 2.4)
# ------------------------------------------------------------------------------
ini <- 500
m <- seq(from = 100, to = 1000, by = 50)
n <- 1000
alpha <- 0.1
reps <- 200
D <- 80
topk_window <- ini
gs_mean_fixed <- 2.4

# Six surge-peak values for the appendix sensitivity study. Verified by
# _diag_setting1.R: across this range FT/RTK/Adj-SAFFRON keep terminal MDR
# above alpha and SCT's later-segment FDR exceeds OMDRC.DD's.
pi_high_list <- c(0.18, 0.22, 0.26, 0.30, 0.34, 0.38)

fig_dir <- file.path(code_dir, "setting1,2,3-appendix")
results_file <- file.path(fig_dir, "Setting2_varied_pi_results.rds")
reuse_cache <- TRUE

if (reuse_cache && file.exists(results_file)) {
  cat("Reusing cached results from", results_file, "\n")
  all_results <- readRDS(results_file)
} else {
  all_results <- lapply(pi_high_list, function(ph) {
    cat(sprintf("Setting 2 sweep: pi_high = %.2f ...\n", ph))
    exp1(m = m, n = n, ini = ini, alpha = alpha, reps = reps, D = D,
         gs_mean = gs_mean_fixed,
         pi_low = 0.05, pi_high = ph,
         surge_center = 0.60, surge_width = 0.15,
         pi_bounds = c(0.01, 0.99), topk_window = topk_window)
  })
  names(all_results) <- paste0("pi_", pi_high_list)
  saveRDS(all_results, results_file)
  cat("Saved sweep results to", results_file, "\n")
}

# ------------------------------------------------------------------------------
# Plotting (same conventions as setting1.R; one panel per swept value)
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
cat("\n=== Setting 2 pi_high sweep terminal metrics (t = max) ===\n")
for (i in seq_along(pi_high_list)) {
  out <- all_results[[i]]
  idx <- length(m)
  cat(sprintf("\n-- pi_high = %.2f --\n", pi_high_list[i]))
  for (j in seq_along(method_suffixes)) {
    mdr <- mean(out[[paste0("mdr.", method_suffixes[j])]][, idx], na.rm = TRUE)
    fdr <- mean(out[[paste0("fdr.", method_suffixes[j])]][, idx], na.rm = TRUE)
    tag <- if (mdr > alpha + 1e-9) " [MDR>alpha]" else ""
    cat(sprintf("%-14s MDR %6.3f FDR %6.3f%s\n", method_labels[j], mdr, fdr, tag))
  }
}

mdr_panels <- lapply(seq_along(pi_high_list), function(i)
  build_panel(all_results[[i]], bquote(pi[max] == .(pi_high_list[i])), "mdr"))
fdr_panels <- lapply(seq_along(pi_high_list), function(i)
  build_panel(all_results[[i]], bquote(pi[max] == .(pi_high_list[i])), "fdr"))

final_layout <- (wrap_plots(mdr_panels, nrow = 3, ncol = 2) |
                   wrap_plots(fdr_panels, nrow = 3, ncol = 2)) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

if (interactive()) print(final_layout)
ggsave(file.path(fig_dir, "Setting2_varied_pi.pdf"),
       plot = final_layout, width = 18, height = 12)
ggsave(file.path(fig_dir, "Setting2_varied_pi.png"),
       plot = final_layout, width = 18, height = 12, dpi = 300)
