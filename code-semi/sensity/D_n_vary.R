# ==============================================================================
# SCRIPT: Sensitivity Analysis of Online Misdiscovery Rate Control (OMDRC.DD)
# Description: Evaluates the robustness of the data-driven OMDRC.DD (PC-DRE)
#              estimator against the Oracle baseline OMDRC.OR with respect to
#              its THREE governing estimation parameters:
#                n    - size of the labeled reference samples z0, z1
#                       (n0 = n1 = n), which controls the accuracy of the
#                       PC-DRE density ratio f1/f0 (Stage 1);
#                beta - growth exponent of the adaptive window
#                       D_t = max(D_min, floor((K0 + t)^beta)), which trades
#                       off the bias (tracking lag) and variance of the local
#                       prevalence estimate pi_t (Stage 2);
#                K0   - size of the initial unlabeled warm-up batch z_ini,
#                       which provides the FIRST local pi estimate and anchors
#                       the growing window through the offset (K0 + t)^beta.
#                       Under the PC-DRE default the density ratio is learned
#                       only from z0/z1, so K0 acts purely on the pi_t warm-up
#                       and the window construction.
#
#   D_min is NOT swept here: it is only a cold-start floor that binds when
#   floor((K0 + t)^beta) < D_min; with the paper defaults (K0 = 500,
#   beta = 0.6) the natural window is (K0+t)^0.6 in [41, 80] over the whole
#   horizon, so D_min = 10 NEVER binds -- it is a pure technical safeguard,
#   which is exactly why sweeping it would be uninformative. The GAM spline
#   dimension is left at its default (spline_k = 10).
#
#   DGP: **Setting 1 of the main text**. One-dimensional SKEWED components,
#   F0 = Exp(1), F1 = Gamma(3, 2); the exact density ratio gives the oracle
#   Lmdr. The signal
#   proportion pi_t follows Setting 1's smooth prior-DROP
#   (0.35 -> 0.03, drop_center 0.50, drop_width 0.10) over warm-up + stream;
#   the stream segment is anchored at k0_ref = 500 so it is IDENTICAL across
#   every K0 value. The monotone drop makes the trailing window lag pi_t
#   downward, so the beta- and K0-sweeps expose the bias-variance trade-off
#   of the window construction, while the n-sweep shows the graded
#   convergence of DD to OR as the labeled reference grows.
#
#   OUTPUT: a single-column 3 x 2 figure (rows: n-sweep / beta-sweep /
#   K0-sweep; columns: MDR / FDR) with pointwise 95% Monte-Carlo bands.
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. Environment Setup
# ------------------------------------------------------------------------------
library(Matrix)
library(foreach)
library(parallel)
library(doParallel)
library(ggplot2)
library(dplyr)
library(patchwork)
library(RColorBrewer)
library(ggpubr)

# Directory that CONTAINS OMDRC.R (with the current OMDRC_DD taking z0 & z1 and
# ratio_method = "gam" for the PC-DRE estimator).
code_dir <- file.path(Sys.getenv("OMDRC_ROOT", unset = getwd()), "code-semi")
setwd(code_dir)
source(file.path(code_dir, "OMDRC.R"))

# Setting 1's smooth prior-DROP prevalence path (high -> low): the monotone
# S-shaped decline makes the trailing window lag pi_t downward.
make_pi_drop <- function(N, pi_low = 0.03, pi_high = 0.35,
                         drop_center = 0.50, drop_width = 0.10) {
  make_smooth_prior_drop(N = N, pi_low = pi_low, pi_high = pi_high,
                         drop_center = drop_center, drop_width = drop_width)
}

# ------------------------------------------------------------------------------
# 1. Core Simulation Engine (PC-DRE version)
# ------------------------------------------------------------------------------
# MDR uses the Monte Carlo mean missed-signal count divided by the known
# expected signal count. FDR is the Monte Carlo mean of the pathwise FDP.
run_simulation <- function(m, ini, n, alpha, reps, D,
                           methods_to_run = c("OR", "DD"),
                           shape_alt = 3, scale_alt = 2,
                           pi_mode = c("vary", "fixed"),
                           pi_fixed = 0.1,
                           pi_low = 0.03, pi_high = 0.35,
                           drop_center = 0.50, drop_width = 0.10,
                           pi_bounds = c(0.01, 0.99),
                           D_mode = "growing", D_beta = 0.6, D_min = 10L,
                           k0_ref = NULL) {

  mdr.or <- if ("OR" %in% methods_to_run) matrix(NA_real_, reps, length(m)) else NULL
  fdr.or <- if ("OR" %in% methods_to_run) matrix(NA_real_, reps, length(m)) else NULL
  mdr.dd <- if ("DD" %in% methods_to_run) matrix(NA_real_, reps, length(m)) else NULL
  fdr.dd <- if ("DD" %in% methods_to_run) matrix(NA_real_, reps, length(m)) else NULL

  pi_mode <- match.arg(pi_mode)
  mx <- max(m)
  if (pi_mode == "vary") {
    # Setting 1's prior-drop over warm-up + stream. When k0_ref is given the
    # path is built over k0_ref + mx and the STREAM segment is anchored at its
    # last mx points, so pi_t over the evaluation stream is IDENTICAL for
    # every ini; the warm-up segment is the ini points immediately preceding
    # the stream (padded with pi_high when ini > k0_ref -- the path is flat at
    # pi_high there anyway).
    N_ref <- (if (is.null(k0_ref)) ini else k0_ref) + mx
    pi_path <- make_pi_drop(N_ref, pi_low = pi_low, pi_high = pi_high,
                            drop_center = drop_center, drop_width = drop_width)
    if (is.null(k0_ref)) {
      pi_ini    <- pi_path[seq_len(ini)]
      pi_stream <- pi_path[(ini + 1L):N_ref]
    } else {
      pi_stream <- pi_path[(k0_ref + 1L):N_ref]
      pi_ini <- if (ini <= k0_ref) {
        pi_path[(k0_ref - ini + 1L):k0_ref]
      } else {
        c(rep(pi_high, ini - k0_ref), pi_path[seq_len(k0_ref)])
      }
    }
  } else {
    # FIXED prevalence: isolates n / D from pi-estimation lag.
    pi_ini    <- rep(pi_fixed, ini)
    pi_stream <- rep(pi_fixed, mx)
  }

  n_cores <- max(1L, min(parallel::detectCores() - 8L, 40L))
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  # Guarantee worker cleanup even if foreach fails (avoids leaked PSOCK workers
  # that can poison subsequent makeCluster calls in the same script run).
  on.exit({ try(stopCluster(cl), silent = TRUE); registerDoSEQ() }, add = TRUE)
  clusterExport(cl, c("m", "ini", "n", "alpha", "reps", "D", "methods_to_run",
                      "shape_alt", "scale_alt", "pi_stream", "pi_ini",
                      "pi_bounds", "code_dir",
                      "D_mode", "D_beta", "D_min"),
                envir = environment())

  result <- foreach(r = 1:reps, .packages = c("mgcv")) %dopar% {
    source(file.path(code_dir, "OMDRC.R"))
    set.seed(r)
    mx <- max(m)

    # ---- Initial unlabeled warm-up window of size `ini` ----
    theta_ini <- rbinom(ini, size = 1, prob = pi_ini)
    z_ini <- rexp(ini, rate = 1)
    z_ini[theta_ini == 1] <- rgamma(sum(theta_ini == 1),
                                    shape = shape_alt, scale = scale_alt)

    # ---- Evaluation stream (1-D skewed components, time-varying pi) ----
    theta_stream <- rbinom(mx, size = 1, prob = pi_stream)
    z_stream <- rexp(mx, rate = 1)
    z_stream[theta_stream == 1] <- rgamma(sum(theta_stream == 1),
                                          shape = shape_alt, scale = scale_alt)
    # Exact Lmdr (oracle gold standard).
    r_true <- dgamma(z_stream, shape = shape_alt, scale = scale_alt) /
      dexp(z_stream, rate = 1)
    lmdr_stream <- (pi_stream * r_true) / ((1 - pi_stream) + pi_stream * r_true)
    r_ini <- dgamma(z_ini, shape = shape_alt, scale = scale_alt) /
      dexp(z_ini, rate = 1)
    lmdr_ini <- (pi_ini * r_ini) / ((1 - pi_ini) + pi_ini * r_ini)

    res_list <- list()

    # ---- Oracle OMDRC ----
    if ("OR" %in% methods_to_run) {
      de <- OMDRC_OR(lmdr_stream, alpha, w = 100L,
                     x.Lmdr_ini = lmdr_ini)$de
      res_list$mdr.or <- sapply(m, function(k)
        sum(theta_stream[1:k] * (1 - de[1:k])) /
          max(sum(pi_stream[1:k]), .Machine$double.eps))
      res_list$fdr.or <- sapply(m, function(k)
        sum((1 - theta_stream[1:k]) * de[1:k]) / max(sum(de[1:k]), 1))
    }

    # ---- Data-driven OMDRC.DD (PC-DRE, GAM classifier density ratio) ----
    if ("DD" %in% methods_to_run) {
      if (n < 3) {
        res_list$mdr.dd <- rep(NA_real_, length(m))
        res_list$fdr.dd <- rep(NA_real_, length(m))
      } else {
        # Labeled null / alternative reference samples of size n each (n0 = n1 = n).
        z0_ref <- rexp(n, rate = 1)
        z1_ref <- rgamma(n, shape = shape_alt, scale = scale_alt)
        # Guard the spline dimension against tiny reference samples; for
        # n >= 6 this equals the PC-DRE default spline_k = 10.
        spline_k <- max(3L, min(10L, as.integer(2 * n - 2)))

        de <- OMDRC_DD(z = z_stream, z_ini = z_ini,
                       z0 = z0_ref, z1 = z1_ref,
                       alpha = alpha, D = D,
                       ratio_method = "gam",
                       ratio_control = list(spline_k = spline_k),
                       pi_bounds = pi_bounds,
                       em_tol = 1e-6, em_max_iter = 50L,
                       D_mode = D_mode, D_beta = D_beta, D_min = D_min)$de
        res_list$mdr.dd <- sapply(m, function(k)
          sum(theta_stream[1:k] * (1 - de[1:k])) /
            max(sum(pi_stream[1:k]), .Machine$double.eps))
        res_list$fdr.dd <- sapply(m, function(k)
          sum((1 - theta_stream[1:k]) * de[1:k]) / max(sum(de[1:k]), 1))
      }
    }

    res_list
  }

  if ("OR" %in% methods_to_run) {
    for (r in 1:reps) {
      mdr.or[r, ] <- result[[r]]$mdr.or
      fdr.or[r, ] <- result[[r]]$fdr.or
    }
  }
  if ("DD" %in% methods_to_run) {
    for (r in 1:reps) {
      mdr.dd[r, ] <- result[[r]]$mdr.dd
      fdr.dd[r, ] <- result[[r]]$fdr.dd
    }
  }

  list(mdr.dd = mdr.dd, fdr.dd = fdr.dd, mdr.or = mdr.or, fdr.or = fdr.or)
}

# Transient PSOCK cluster failures have been observed to abort the long sweep
# midway with no R-level error message; retry each block up to 3 times.
run_sim_retry <- function(..., max_try = 3L) {
  for (a in seq_len(max_try)) {
    out <- tryCatch(run_simulation(...), error = function(e) {
      message("run_simulation failed (attempt ", a, "/", max_try, "): ",
              conditionMessage(e))
      NULL
    })
    if (!is.null(out)) return(out)
    Sys.sleep(5)
  }
  stop("run_simulation failed after ", max_try, " attempts")
}

# ------------------------------------------------------------------------------
# 2. Experimental Configuration
# ------------------------------------------------------------------------------
m_seq   <- seq(from = 100, to = 1000, by = 50)
alpha_v <- 0.1
reps_v  <- 200        # Monte-Carlo replicates reported in the manuscript.
                      # raising reps later just APPENDS new replicates.
# --- Setting 1 DGP parameters (1-D skewed components) ---
shape_alt_v <- 3      # F1 = Gamma(shape_alt, scale_alt) vs F0 = Exp(1)
scale_alt_v <- 2

# --- Prevalence path: Setting 1's smooth prior-drop ---
pi_mode_v     <- "vary"
pi_low_v      <- 0.03
pi_high_v     <- 0.35
drop_center_v <- 0.50
drop_width_v  <- 0.10
pi_bnds <- c(0.01, 0.99)   # paper default: broad bounds [pi_lower, pi_upper]

# --- Parameter grids for the three sensitivity sweeps ---
n_list    <- c(10, 20, 50, 100, 200, 500, 1000)      # labeled reference size (n0 = n1 = n)
beta_list <- c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)    # window growth exponent
k0_list   <- c(10, 20, 50, 100, 200, 500, 1000)      # initial unlabeled batch size K0
                                                   # (same grid as n_list)

# --- Default values: each sweep varies ONE parameter, holding the other two ---
default_n      <- 500    # paper default n0 = n1 = 500
default_ini    <- 500    # paper default K0; also the anchor (k0_ref) of the
                         # stream pi path, so the K0-sweep changes ONLY the
                         # warm-up length, never the stream prevalence
default_D      <- 150    # fallback for D_mode = "fixed" (unused in growing mode)
default_D_mode <- "growing"
default_D_beta <- 0.6    # paper default: D_t = min{K0+t, max(D_min, floor((K0+t)^0.6))}
default_D_min  <- 10     # cold-start floor (paper default): with K0 = 500 and
# beta = 0.6 the natural window (K0+t)^0.6 is in [41, 80] >> 10 over the whole
# horizon, so the floor NEVER binds -- a pure technical safeguard, which is
# why it is not swept.

# ------------------------------------------------------------------------------
# 3. Run the experiments (or reuse cached results for a fast re-plot)
# ------------------------------------------------------------------------------
results_file <- file.path(code_dir, "sensity", "D_n_vary_results.rds")
reuse_cache  <- TRUE   # set FALSE to rerun everything from scratch

sim_common <- function(...) {
  run_sim_retry(...,
                shape_alt = shape_alt_v, scale_alt = scale_alt_v,
                pi_mode = pi_mode_v,
                pi_low = pi_low_v, pi_high = pi_high_v,
                drop_center = drop_center_v, drop_width = drop_width_v,
                pi_bounds = pi_bnds)
}

run_or_baseline <- function() {
  cat("\n====== Running OMDRC.OR baseline (once) ======\n")
  sim_common(m = m_seq, ini = default_ini, n = default_n,
             alpha = alpha_v, reps = reps_v,
             D = default_D, methods_to_run = "OR")
}
run_n_sweep <- function() {
  lapply(n_list, function(nv) {
    cat(paste0("\n====== OMDRC.DD (PC-DRE) with n = ", nv, " ======\n"))
    sim_common(m = m_seq, ini = default_ini, n = nv,
               alpha = alpha_v, reps = reps_v, D = default_D,
               methods_to_run = "DD",
               D_mode = default_D_mode, D_beta = default_D_beta,
               D_min = default_D_min)
  })
}
run_beta_sweep <- function() {
  lapply(beta_list, function(bv) {
    cat(paste0("\n====== OMDRC.DD (PC-DRE) with beta = ", bv,
               " (D_t = floor((K0+t)^", bv, ")) ======\n"))
    sim_common(m = m_seq, ini = default_ini, n = default_n,
               alpha = alpha_v, reps = reps_v, D = default_D,
               methods_to_run = "DD",
               D_mode = "growing", D_beta = bv,
               D_min = default_D_min)
  })
}
run_k0_sweep <- function() {
  lapply(k0_list, function(kv) {
    cat(paste0("\n====== OMDRC.DD (PC-DRE) with K0 = ", kv,
               " (stream pi anchored at k0_ref = ", default_ini, ") ======\n"))
    sim_common(m = m_seq, ini = kv, n = default_n,
               alpha = alpha_v, reps = reps_v, D = default_D,
               methods_to_run = "DD",
               D_mode = "growing", D_beta = default_D_beta,
               D_min = default_D_min, k0_ref = default_ini)
  })
}

# Load the cache and reuse whichever components are present AND match the
# current grids; recompute only the missing/mismatched components, then merge
# everything back into the cache file. (Caches saved by the previous versions
# already contain res_n / res_beta with identical grids and seeds, so only the
# new K0 sweep is computed on the first run after the D_min -> K0 switch.)
cached <- if (reuse_cache && file.exists(results_file)) readRDS(results_file) else NULL

if (!is.null(cached$mdr.or.mat) && !is.null(cached$fdr.or.mat)) {
  cat("Reusing cached OMDRC.OR baseline.\n")
  mdr_or_mat <- cached$mdr.or.mat
  fdr_or_mat <- cached$fdr.or.mat
} else {
  base_res   <- run_or_baseline()
  mdr_or_mat <- base_res$mdr.or
  fdr_or_mat <- base_res$fdr.or
}
mdr_or_vec <- colMeans(mdr_or_mat, na.rm = TRUE)
fdr_or_vec <- colMeans(fdr_or_mat, na.rm = TRUE)

if (!is.null(cached$res_n) && isTRUE(all.equal(cached$n_list, n_list))) {
  cat("Reusing cached n-sweep.\n")
  res_n <- cached$res_n
} else {
  res_n <- run_n_sweep()
}

if (!is.null(cached$res_beta) && isTRUE(all.equal(cached$beta_list, beta_list))) {
  cat("Reusing cached beta-sweep.\n")
  res_beta <- cached$res_beta
} else {
  res_beta <- run_beta_sweep()
}

if (!is.null(cached$res_k0) && isTRUE(all.equal(cached$k0_list, k0_list))) {
  cat("Reusing cached K0-sweep.\n")
  res_k0 <- cached$res_k0
} else {
  res_k0 <- run_k0_sweep()
}

# Cache raw results so the figure can be re-plotted without rerunning.
saveRDS(list(m = m_seq, alpha = alpha_v,
             shape_alt = shape_alt_v, scale_alt = scale_alt_v, reps = reps_v,
             pi_mode = pi_mode_v,
             pi_low = pi_low_v, pi_high = pi_high_v,
             drop_center = drop_center_v, drop_width = drop_width_v,
             pi_bounds = pi_bnds,
             n_list = n_list, beta_list = beta_list, k0_list = k0_list,
             default_n = default_n, default_D_beta = default_D_beta,
             default_D_min = default_D_min, default_ini = default_ini,
             mdr_or = mdr_or_vec, fdr_or = fdr_or_vec,
             mdr.or.mat = mdr_or_mat, fdr.or.mat = fdr_or_mat,
             res_n = res_n, res_beta = res_beta, res_k0 = res_k0),
        file = results_file)

# Terminal MDR/FDR summary (OR vs growing-window DD at the default config).
idxT <- length(m_seq)
i_n_default    <- which(n_list == default_n)
i_beta_default <- which(beta_list == default_D_beta)
i_k0_default   <- which(k0_list == default_ini)
cat(sprintf("\n--- Terminal MDR/FDR at eval-t=%d (alpha=%.2f) ---\n", m_seq[idxT], alpha_v))
cat(sprintf("%-34s %8s %8s\n", "Method", "MDR", "FDR"))
cat(sprintf("%-34s %8.4f %8.4f\n", "OMDRC.OR", mdr_or_vec[idxT], fdr_or_vec[idxT]))
cat(sprintf("%-34s %8.4f %8.4f\n",
            sprintf("DD (n=%d, beta=%.1f, K0=%d)",
                    default_n, default_D_beta, default_ini),
            mean(res_k0[[i_k0_default]]$mdr.dd[, idxT], na.rm = TRUE),
            mean(res_k0[[i_k0_default]]$fdr.dd[, idxT], na.rm = TRUE)))
cat("\n--- Terminal MDR by sweep value ---\n")
cat("n   :", paste(sprintf("%g=%.3f", n_list,
    sapply(res_n, function(z) mean(z$mdr.dd[, idxT], na.rm = TRUE))), collapse = "  "), "\n")
cat("beta:", paste(sprintf("%g=%.3f", beta_list,
    sapply(res_beta, function(z) mean(z$mdr.dd[, idxT], na.rm = TRUE))), collapse = "  "), "\n")
cat("K0  :", paste(sprintf("%g=%.3f", k0_list,
    sapply(res_k0, function(z) mean(z$mdr.dd[, idxT], na.rm = TRUE))), collapse = "  "), "\n")

# ------------------------------------------------------------------------------
# 4. Data Processing
# ------------------------------------------------------------------------------
# Pointwise 95% Monte-Carlo band: estimate +/- 1.96 * SE over replicates.
band_cols <- function(mat) {
  est <- colMeans(mat, na.rm = TRUE)
  se  <- apply(mat, 2, stats::sd, na.rm = TRUE) / sqrt(nrow(mat))
  list(est = est, lower = pmax(0, est - 1.96 * se),
       upper = pmin(1, est + 1.96 * se))
}

get_plot_df <- function(res_list, param_name, param_values,
                        or_mdr_mat, or_fdr_mat) {
  bm_or <- band_cols(or_mdr_mat)
  bf_or <- band_cols(or_fdr_mat)
  df_or <- data.frame(t = m_seq,
                      MDR = bm_or$est, MDR_lo = bm_or$lower, MDR_hi = bm_or$upper,
                      FDR = bf_or$est, FDR_lo = bf_or$lower, FDR_hi = bf_or$upper,
                      Method = "OMDRC.OR")

  df_dd <- do.call(rbind, lapply(seq_along(param_values), function(i) {
    res <- res_list[[i]]
    bm <- band_cols(res$mdr.dd)
    bf <- band_cols(res$fdr.dd)
    data.frame(t = m_seq,
               MDR = bm$est, MDR_lo = bm$lower, MDR_hi = bm$upper,
               FDR = bf$est, FDR_lo = bf$lower, FDR_hi = bf$upper,
               Method = paste0("OMDRC.DD (", param_name, "=", param_values[i], ")"))
  }))

  df <- rbind(df_or, df_dd)
  lv <- c("OMDRC.OR", unique(df_dd$Method))
  df$Method <- factor(df$Method, levels = lv)
  df
}

# ------------------------------------------------------------------------------
# 5. Visualization (3 x 2 grid: rows n / beta / K0, columns MDR / FDR)
# ------------------------------------------------------------------------------
custom_theme <- theme_bw() +
  theme(
    plot.subtitle = element_text(hjust = 0.5, size = 22, face = "plain"),
    legend.position = "bottom",
    legend.title = element_blank(),
    # size 15: matches Figure2.R's legend.text so method labels render at the
    # same size across paper figures; 16 overflows the 10-inch device at ncol = 3.
    legend.text = element_text(size = 15),
    legend.key.width = unit(1.0, "cm"),
    axis.text = element_text(size = 19, colour = "black"),
    axis.title = element_text(size = 22),
    panel.grid.major = element_line(colour = "grey92", size = 0.4),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.8)
  )
my_shapes <- c(16, 17, 15, 3, 4, 2, 7, 1)
# Legend keys show line + marker. NOTE (ggplot2 3.3.5): there is no `linewidth`
# aesthetic pre-3.4, so putting `size` in override.aes also fattens the legend
# LINE into a thick bar that hides the marker. Therefore do NOT override size;
# keys inherit geom_line(size = 1.2) + geom_point(size = 2.8) directly.
legend_guides <- guides(
  color = guide_legend(ncol = 3, byrow = TRUE)
)

# Build one row (MDR + FDR panels sharing a collected legend) for a sweep.
# `labels` is a plotmath expression vector (same length/order as the factor
# levels of df$Method) so the legend renders n[0]==n[1], beta, K[0] as math.
make_row <- function(df, colors, tag, labels) {
  p_mdr <- ggplot(df, aes(t, MDR, color = Method, shape = Method)) +
    geom_ribbon(aes(ymin = MDR_lo, ymax = MDR_hi, fill = Method),
                alpha = 0.15, colour = NA, show.legend = FALSE) +
    geom_line(size = 1.2) + geom_point(size = 2.8) +
    geom_hline(yintercept = alpha_v, linetype = "dashed", size = 1.0) +
    scale_color_manual(values = colors, labels = labels) +
    scale_fill_manual(values = colors, labels = labels) +
    scale_shape_manual(values = my_shapes, labels = labels) +
    # extra right-side expansion so the "1000" tick label is not clipped
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.08))) +
    labs(x = "Time (t)", y = "MDR", subtitle = paste0("(", tag, ".1)")) +
    custom_theme + legend_guides
  p_fdr <- ggplot(df, aes(t, FDR, color = Method, shape = Method)) +
    geom_ribbon(aes(ymin = FDR_lo, ymax = FDR_hi, fill = Method),
                alpha = 0.15, colour = NA, show.legend = FALSE) +
    geom_line(size = 1.2) + geom_point(size = 2.8) +
    scale_color_manual(values = colors, labels = labels) +
    scale_fill_manual(values = colors, labels = labels) +
    scale_shape_manual(values = my_shapes, labels = labels) +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.08))) +
    labs(x = "Time (t)", y = "FDR", subtitle = paste0("(", tag, ".2)")) +
    custom_theme + legend_guides
  (p_mdr | p_fdr) + plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
}

# --- Row A: Sensitivity to n (Greens); legend shows n_ref = value ---
df_n <- get_plot_df(res_n, "n", n_list, mdr_or_mat, fdr_or_mat)
colors_n <- c("#F8766D", brewer.pal(9, "Greens")[3:9])
# NOTE plotmath: the quoted subscript n['ref'] renders "ref" UPRIGHT, matching
# the paper's n_{\mathrm{ref}}; a chained '==' is NOT valid R syntax, so the
# equals sign is inserted as a literal string via the '*' juxtaposition.
labels_n <- parse(text = c("OMDRC.OR",
                           paste0("OMDRC.DD~(n['ref']*'='*", n_list, ")")))
row_n <- make_row(df_n, colors_n, "a", labels_n)

# --- Row B: Sensitivity to beta (Purples); legend shows the Greek beta ---
df_beta <- get_plot_df(res_beta, "beta", beta_list, mdr_or_mat, fdr_or_mat)
colors_beta <- c("#F8766D", brewer.pal(9, "Purples")[3:9])
labels_beta <- parse(text = c("OMDRC.OR",
                              paste0("OMDRC.DD~(beta==", beta_list, ")")))
row_beta <- make_row(df_beta, colors_beta, "b", labels_beta)

# --- Row C: Sensitivity to K0 (Blues); legend shows the subscript K0 ---
df_k0 <- get_plot_df(res_k0, "K0", k0_list, mdr_or_mat, fdr_or_mat)
colors_k0 <- c("#F8766D", brewer.pal(9, "Blues")[3:9])
labels_k0 <- parse(text = c("OMDRC.OR",
                            paste0("OMDRC.DD~(K[0]==", k0_list, ")")))
row_k0 <- make_row(df_k0, colors_k0, "c", labels_k0)

# ------------------------------------------------------------------------------
# 6. Final Integration and Layout (3 x 2 grid)
# ------------------------------------------------------------------------------
final_plot <- ggarrange(row_n, row_beta, row_k0, ncol = 1, nrow = 3)

# Portrait-friendly aspect (width 10) so the figure reads well at single-column
# width, matching setting1,2,3/Setting123_3x2 (10 x 13, 3 rows); the extra
# height accommodates the three per-row legends.
if (interactive()) print(final_plot)
ggsave(file.path(code_dir, "sensity", "Rplot_D_n_sensitivity.pdf"),
       final_plot, width = 10, height = 16)
ggsave(file.path(code_dir, "sensity", "Rplot_D_n_sensitivity.png"),
       final_plot, width = 10, height = 16, dpi = 300)
