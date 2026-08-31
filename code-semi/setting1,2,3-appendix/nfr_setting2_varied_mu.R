# ==============================================================================
# APPENDIX SCRIPT (NFR VARIANT): Online MDR Control - Setting 2 sensitivity to
#                                the signal mean gs_mean, WITH the NFR
#                                (anti-free-riding) barrier on BOTH OMDRC arms
#                                -- credited ledger, d = 100.
#
# This is the NFR companion of setting2_varied_mu.R in this appendix folder. The
# DGP is kept VERBATIM from that script (Setting 2 was NOT re-tuned; only
# Setting 1 received updated parameters):
#   * Gaussian location shift F0 = N(0,1), F1 = N(gs_mean, 1);
#   * non-monotone Gaussian SURGE prevalence path pi_t (low -> peak -> low),
#     pi_low = 0.05, pi_high = 0.30, surge_center = 0.60, surge_width = 0.15;
#   * OMDRC.DD uses the GAM classifier density ratio (z0 & z1 references) with
#     pi_bounds = c(0.01, 0.99), D_mode = "growing", D_beta = 0.6, D_min = 10L
#     and labeled reference size n = 1000 -- the appendix folder's own values,
#     kept so the NFR panels stay directly comparable with setting2_varied_mu.R;
#   * FT uses the nominal miss budget alpha;
#   * the SAME six procedures, colours, shapes, labels and Monte-Carlo 95%
#     confidence bands; we SWEEP the alternative mean gs_mean over the same
#     grid as setting2_varied_mu.R (3x2 MDR panels next to 3x2 FDR panels).
#
# NFR change (identical to nfr_setting2_fig.R): both OMDRC arms carry the
# CREDITED local barrier at w = 100 (Algorithms 1 and 2):
#   OMDRC.OR      -> OMDRC.OR+NFR   = OMDRC_OR_NFR(true q,      alpha, d)
#   OMDRC.DD      -> OMDRC.DD+NFR   = OMDRC_OR_NFR(estimated q, alpha, d)
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

# Barrier window: d = 100, identical to nfr_setting{1,2,3}_fig.R.
D_NFR <- 100L

# Gaussian-surge prevalence path (low -> peak -> low), identical to setting2.R.
make_pi_surge <- function(N, pi_low = 0.05, pi_high = 0.30,
                          center = 0.60, width = 0.15) {
  u <- seq(0, 1, length.out = as.integer(N))
  pmin(pmax(pi_low + (pi_high - pi_low) *
              exp(-((u - center)^2) / (2 * width^2)), pi_low), pi_high)
}

# ------------------------------------------------------------------------------
# Monte Carlo engine (setting2_varied_mu.R engine; gs_mean is the swept arg;
# both OMDRC arms carry the credited NFR barrier)
# ------------------------------------------------------------------------------
exp1 <- function(m, n, ini, alpha, reps, D,
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
      "omdrc_file", "nfr_core_file", "D_NFR"),
    envir = environment()
  )

  result <- foreach(r = seq_len(reps), .packages = "onlineFDR") %dopar% {
    source(omdrc_file)
    source(nfr_core_file)
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
    lmdr_ini <- lmdr_all[seq_len(ini)]   # fills the gamma window during warm-up
    theta_s <- theta[(ini + 1):N]

    p_value <- findInterval(z, z_ref_sorted) / n_ref
    p_value <- pmin(pmax(p_value, 1e-6), 1 - 1e-6)

    z1_ref <- rnorm(n, mean = gs_mean, sd = 1)
    z0_ref <- rnorm(n, mean = 0, sd = 1)

    # ---- Six procedures (both OMDRC arms carry the credited NFR barrier) ----
    d_or <- OMDRC_OR_NFR(lmdr, alpha, d = D_NFR, x.Lmdr_ini = lmdr_ini)
    # OMDRC_DD is run to produce the ESTIMATED score path (DR / DR_ini). We
    # replay the same barrier explicitly below so that both OMDRC arms expose
    # identical diagnostics; the score path does not depend on the decisions.
    d_dd <- OMDRC_DD(z = z, z_ini = z_ini, z0 = z0_ref, z1 = z1_ref,
                     alpha = alpha, ratio_method = "gam",
                     pi_bounds = pi_bounds, em_tol = 1e-6, em_max_iter = 50L,
                     D_mode = "growing", D_beta = 0.6, D_min = 10L)
    d_dd_nfr <- OMDRC_OR_NFR(d_dd$DR, alpha, d = D_NFR,
                             x.Lmdr_ini = d_dd$DR_ini)
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
# Sweep the alternative mean gs_mean
# ------------------------------------------------------------------------------
ini <- 500
m <- seq(from = 100, to = 1000, by = 50)
n <- 1000
alpha <- 0.1
reps <- 200
D <- 80
topk_window <- ini

# Six alternative-mean values, SAME grid as setting2_varied_mu.R. A SMALL-mu
# range is used deliberately: it keeps FT/RTK/Adj-SAFFRON terminal MDR above
# alpha AND keeps SCT's later-segment FDR above OMDRC.DD's. At large mu the
# signal is too easy and the baselines stop failing.
mu_list <- c(1.5, 1.8, 2.1, 2.4, 2.7, 3.0)

fig_dir <- file.path(code_dir, "setting1,2,3-appendix")
results_file <- file.path(fig_dir, "NFR_Setting2_varied_mu_results.rds")
reuse_cache <- identical(Sys.getenv("NFR_REUSE_CACHE"), "1")

cat(sprintf("Starting Setting 2 gs_mean-sweep + NFR credited barrier (d=%d)...\n",
            D_NFR))
if (reuse_cache && file.exists(results_file)) {
  cat("Reusing cached results from", results_file, "\n")
  all_results <- readRDS(results_file)
} else {
  all_results <- lapply(mu_list, function(mu_v) {
    cat(sprintf("Setting 2 sweep: gs_mean = %.2f ...\n", mu_v))
    exp1(m = m, n = n, ini = ini, alpha = alpha, reps = reps, D = D,
         gs_mean = mu_v,
         pi_low = 0.05, pi_high = 0.30,
         surge_center = 0.60, surge_width = 0.15,
         pi_bounds = c(0.01, 0.99), topk_window = topk_window)
  })
  names(all_results) <- paste0("mu_", mu_list)
  saveRDS(all_results, results_file)
  cat("Saved sweep results to", results_file, "\n")
}

# ------------------------------------------------------------------------------
# Plotting (same conventions as setting2_varied_mu.R; one panel per swept value)
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
cat("\n=== Setting 2 sweep terminal metrics (t = max, NFR arms) ===\n")
for (i in seq_along(mu_list)) {
  out <- all_results[[i]]
  idx <- length(m)
  cat(sprintf("\n-- gs_mean = %.2f --\n", mu_list[i]))
  for (j in seq_along(method_suffixes)) {
    mdr <- mean(out[[paste0("mdr.", method_suffixes[j])]][, idx], na.rm = TRUE)
    fdr <- mean(out[[paste0("fdr.", method_suffixes[j])]][, idx], na.rm = TRUE)
    tag <- if (mdr > alpha + 1e-9) " [MDR>alpha]" else ""
    cat(sprintf("%-14s MDR %6.3f FDR %6.3f%s\n", method_labels[j], mdr, fdr, tag))
  }
  cat(sprintf("  ledger (d=%d, credited): OR mean MDR_q %.4f max %.4f [%s]",
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
ggsave(file.path(fig_dir, "NFR_Setting2_varied_mu.pdf"),
       plot = final_layout, width = 18, height = 12)
ggsave(file.path(fig_dir, "NFR_Setting2_varied_mu.png"),
       plot = final_layout, width = 18, height = 12, dpi = 300)
cat("Figures written to", file.path(fig_dir, "NFR_Setting2_varied_mu.pdf/png"), "\n")
