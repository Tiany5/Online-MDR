# ==============================================================================
# SCRIPT: Main-text Setting (a), including the local barrier in Algorithms 1
#         and 2 (w = 100) for both OMDRC arms.
#
#   OMDRC.OR      -> OMDRC.OR+NFR   = OMDRC_OR_NFR(true q,      alpha, d)
#   OMDRC.DD      -> OMDRC.DD+NFR   = OMDRC_OR_NFR(estimated q, alpha, d)
# "credited" = a barrier-forced rejection credits +alpha*q to the live ledger
# (_nfr_core.R::OMDRC_OR_NFR).  SCT / RTK / Adj-SAFFRON / FT are untouched and
# still read the SAME estimated Lmdr score path, so the comparison is unchanged.
#
# DGP: F0 = Exp(1), F1 = Gamma(shape = 3, scale = 2), with the smooth
# prevalence drop 0.35 -> 0.03 specified in the main text. Results use
# n0 = n1 = 500, K0 = 500 and 200 Monte Carlo replications.
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

# Shared NFR core (gamma path + the three barrier ledgers).
nfr_dir <- file.path(code_dir, "nofreeride")
nfr_core_file <- file.path(nfr_dir, "_nfr_core.R")
if (!file.exists(nfr_core_file)) stop("Cannot find _nfr_core.R in: ", nfr_dir)
source(nfr_core_file)

# Barrier window w = 100 (called d in the shared replay helper).
D_NFR <- 100L

# Smooth DECREASING prevalence path (high -> low) used by Setting 1.
make_pi_drop <- function(N, pi_low = 0.03, pi_high = 0.35,
                         drop_center = 0.50, drop_width = 0.10) {
  make_smooth_prior_drop(N = N, pi_low = pi_low, pi_high = pi_high,
                         drop_center = drop_center, drop_width = drop_width)
}

# ------------------------------------------------------------------------------
# Monte Carlo engine
# ------------------------------------------------------------------------------
exp1 <- function(m, n, ini, alpha, reps, D,
                 shape_alt = 3, scale_alt = 2,
                 pi_low = 0.03, pi_high = 0.35,
                 drop_center = 0.50, drop_width = 0.10,
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

  # Reference sample for Adj-SAFFRON p-values (sorted once for findInterval).
  M_ref <- 200000L
  set.seed(202401)
  z_ref_sorted <- sort(rgamma(M_ref, shape = shape_alt, scale = scale_alt))
  n_ref <- length(z_ref_sorted)

  # Fixed time-varying prevalence path shared by every replicate.
  N <- max(m) + ini
  pi_path <- make_pi_drop(N, pi_low = pi_low, pi_high = pi_high,
                          drop_center = drop_center, drop_width = drop_width)

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
      "topk_window", "z_ref_sorted", "n_ref", "shape_alt",
      "scale_alt", "omdrc_file", "nfr_core_file", "D_NFR"),
    envir = environment()
  )

  result <- foreach(r = seq_len(reps), .packages = "onlineFDR") %dopar% {
    source(omdrc_file)
    source(nfr_core_file)
    set.seed(r)

    # ---- Evaluation stream (skewed alternative, time-varying pi) ----
    theta <- rbinom(N, size = 1, prob = pi_path)
    z0_stream <- rexp(N, rate = 1)                                # F0 ~ Exp(1)
    z_alt <- rgamma(N, shape = shape_alt, scale = scale_alt)      # F1 ~ Gamma(3,2)
    z_all <- ifelse(theta == 0, z0_stream, z_alt)

    f0 <- dexp(z_all, rate = 1)
    f1 <- dgamma(z_all, shape = shape_alt, scale = scale_alt)
    lmdr_all <- (pi_path * f1) / ((1 - pi_path) * f0 + pi_path * f1)

    z_ini <- z_all[seq_len(ini)]
    z <- z_all[(ini + 1):N]
    lmdr <- lmdr_all[(ini + 1):N]
    lmdr_ini <- lmdr_all[seq_len(ini)]   # fills the gamma window during warm-up
    theta_s <- theta[(ini + 1):N]

    # Adj-SAFFRON p-values from the alternative-distribution ECDF.
    p_value <- findInterval(z, z_ref_sorted) / n_ref
    p_value <- pmin(pmax(p_value, 1e-6), 1 - 1e-6)

    # ---- Labeled reference samples for the data-driven method ----
    z1_ref <- rgamma(n, shape = shape_alt, scale = scale_alt)
    z0_ref <- rexp(n, rate = 1)

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
    # FT is the offline fixed-threshold oracle at the SAME nominal budget alpha
    # as every other method. It solves the global MDR = alpha constraint over
    # the whole stream, so its TERMINAL cumulative MDR lands on alpha; but a
    # single frozen threshold cannot track the drifting prevalence, so its
    # RUNNING (cumulative) MDR trajectory transiently departs from alpha.
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
# Run the experiment
# ------------------------------------------------------------------------------
ini <- 500
m <- seq(from = 100, to = 1000, by = 50)
n <- 500           # labeled reference size (paper default n0 = n1 = 500)
alpha <- 0.1
reps <- 200
D <- 150           # legacy fixed-window fallback (UNUSED: DD now uses the
                   # paper growing window D_t = min{K0+t, max{10, (K0+t)^0.6}})
topk_window <- ini

# Cache of the intermediate numerical results, so plot/label tweaks do not
# require rerunning the experiment. Run with NFR_REUSE_CACHE=1 to skip the Monte
# Carlo and re-plot from the rds (folder convention, see the other scripts here).
fig_dir <- nfr_dir
results_file <- file.path(fig_dir, "NFR_Setting1_results.rds")
reuse_cache <- identical(Sys.getenv("NFR_REUSE_CACHE"), "1")

cat(sprintf("Starting Setting 1 + NFR credited barrier (d=%d)...\n", D_NFR))
if (reuse_cache && file.exists(results_file)) {
  cat("Reusing cached numerical results from", results_file, "\n")
  out_exp1 <- readRDS(results_file)
} else {
  out_exp1 <- exp1(m = m, n = n, ini = ini, alpha = alpha, reps = reps, D = D,
                   shape_alt = 3, scale_alt = 2,
                   pi_low = 0.03, pi_high = 0.35,
                   drop_center = 0.50, drop_width = 0.10,
                   pi_bounds = c(0.01, 0.99), topk_window = topk_window)
  saveRDS(out_exp1, results_file)
  cat("Saved numerical results to", results_file, "\n")
}

cat(sprintf("Static Lmdr threshold: mean %.3f (SD %.3f)\n",
            mean(out_exp1$static.threshold), sd(out_exp1$static.threshold)))
cat(sprintf("Rolling Top-k: mean k %.1f of %d (mean review fraction %.3f)\n",
            mean(out_exp1$topk.k), topk_window, mean(out_exp1$topk.rho)))

# Diagnostic: terminal (t = max m) MDR / FDR per method, flagging failures.
.diag_terminal <- function(out, m, alpha) {
  suf <- c("or", "dd", "off", "static", "topk", "lord")
  lab <- c("OMDRC.OR+NFR", "OMDRC.DD+NFR", "FT", "SCT",
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
.diag_terminal(out_exp1, m, alpha)

# Hard ledger check (the point of the NFR section): MDR_q in the TRUE q must be
# <= alpha on EVERY path for the oracle arm. The DD arm runs its ledger on q_hat
# so it can fail -- that failure is the documented cost of estimation error.
.diag_ledger <- function(out, alpha) {
  cat(sprintf("\n--- Barrier ledger check (d=%d, credited) ---\n", out$d_nfr))
  cat(sprintf("%-14s %10s %10s %12s %10s\n",
              "Arm", "mean MDR_q", "max MDR_q", "ledger", "blocked"))
  for (a in c("or", "dd")) {
    v <- out[[paste0("mdrq.", a)]]
    lab <- if (a == "or") "OMDRC.OR+NFR" else "OMDRC.DD+NFR"
    cat(sprintf("%-14s %10.4f %10.4f %12s %10.1f\n",
                lab, mean(v), max(v),
                if (max(v) <= alpha + 1e-9) "PASS" else "FAIL",
                mean(out[[paste0("blocked.", a)]])))
  }
}
.diag_ledger(out_exp1, alpha)

# ------------------------------------------------------------------------------
# Plotting (pointwise 95% Monte Carlo confidence bands + zoom inset)
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

df_mdr <- prepare_tidy_data(out_exp1, "mdr", m, method_suffixes, method_labels)
df_fdr <- prepare_tidy_data(out_exp1, "fdr", m, method_suffixes, method_labels)

g1.1_main <- ggplot(df_mdr,
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
  labs(subtitle = "(a.1)", x = "Time (t)", y = "MDR") +
  custom_theme +
  coord_cartesian(ylim = c(0, NA)) +
  scale_y_continuous(breaks = seq(0, 0.5, by = 0.1))

g1.1_final <- g1.1_main

g1.2 <- ggplot(df_fdr,
               aes(x = t, y = value, color = type, shape = type,
                   group = type)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = type),
              alpha = 0.16, colour = NA, show.legend = FALSE) +
  geom_line(size = 0.75) +
  geom_point(size = 1.4) +
  labs(subtitle = "(a.2)", x = "Time (t)", y = "FDR") +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  scale_shape_manual(values = my_shapes) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.10))) +
  custom_theme

# Always save the Setting 1 + NFR figure (PDF + PNG) into the nofreeride folder.
final_plot1 <- (g1.1_final + g1.2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
if (interactive()) print(final_plot1)
ggsave(file.path(fig_dir, "NFR_Setting1_fig.pdf"),
       plot = final_plot1, width = 10, height = 5)
ggsave(file.path(fig_dir, "NFR_Setting1_fig.png"),
       plot = final_plot1, width = 10, height = 5, dpi = 300)
cat("Figures written to", file.path(fig_dir, "NFR_Setting1_fig.pdf/png"), "\n")
