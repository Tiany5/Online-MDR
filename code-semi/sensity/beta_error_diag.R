# ==============================================================================
# SCRIPT: beta estimation-error diagnostics for OMDRC.DD (PC-DRE)
# Description: Companion diagnostic of the beta-sweeps in
#              sensity/D_n_vary.R (Setting 1) and
#              sensity/D_n_vary_setting2.R (Setting 2).
#              Instead of the downstream MDR/FDR, it looks DIRECTLY at the two
#              estimation targets governed by the window exponent beta of
#              D_t = min{K0+t, max(D_min, floor((K0+t)^beta))}:
#                (1) bias(hat(pi)_t)  = E[hat(pi)_t - pi_t]   (tracking lag),
#                (2) MAE(hat(pi)_t)   = E|hat(pi)_t - pi_t|,
#                (3) MAE(hat(L)_mdr,t)= E|hat(L)_mdr,t - L_mdr,t|,
#              where L_mdr,t uses the TRUE density ratio and TRUE pi_t.
#
#   Data generation, seeds (set.seed(r), r = 1..reps) and defaults
#   (n0 = n1 = 500, K0 = 500, D_min = 10, pi_bounds = c(0.01, 0.99),
#   em_tol = 1e-6, em_max_iter = 50) replicate the beta-sweep blocks of the
#   two sweep scripts EXACTLY, so the streams here are the very same streams
#   that produced the MDR/FDR rows there. The GAM density ratio is fitted
#   ONCE per replicate (beta does not enter Stage 1) and OMDRC_FROM_RATIO is
#   then run per beta, which is much cheaper than a full re-sweep.
#
#   OUTPUT: one 2 x 3 figure (rows: (a) Setting 1 prior-drop / (b) Setting 2
#   surge; columns: pi bias / pi MAE / Lmdr MAE) with pointwise 95%
#   Monte-Carlo bands, plus a compact rds cache of the per-t summaries.
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

code_dir <- file.path(Sys.getenv("OMDRC_ROOT", unset = getwd()), "code-semi")
setwd(code_dir)
source(file.path(code_dir, "OMDRC.R"))

# Setting 1's smooth prior-DROP prevalence path (high -> low).
make_pi_drop <- function(N, pi_low = 0.03, pi_high = 0.35,
                         drop_center = 0.50, drop_width = 0.10) {
  make_smooth_prior_drop(N = N, pi_low = pi_low, pi_high = pi_high,
                         drop_center = drop_center, drop_width = drop_width)
}
# Setting 2's smooth prevalence SURGE path (low -> high -> low).
make_pi_surge <- function(N, pi_low = 0.05, pi_high = 0.30,
                          center = 0.60, width = 0.15) {
  u <- seq(0, 1, length.out = as.integer(N))
  pmin(pmax(pi_low + (pi_high - pi_low) *
              exp(-((u - center)^2) / (2 * width^2)), pi_low), pi_high)
}

# ------------------------------------------------------------------------------
# 1. Configuration (matches the beta-sweep blocks of the two sweep scripts)
# ------------------------------------------------------------------------------
mx_v      <- 1000                                  # stream length max(m_seq)
alpha_v   <- 0.1
reps_v    <- 100                                   # deterministic seeds 1..reps
beta_list <- c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
default_n   <- 500                                 # n0 = n1 = 500
default_ini <- 500                                 # K0 (k0_ref = NULL in the
                                                   # beta sweep, so ini = 500)
default_D_min <- 10L
pi_bnds <- c(0.01, 0.99)

# Setting 1 DGP: F0 = Exp(1), F1 = Gamma(3, 2); prior-drop 0.35 -> 0.03.
shape_alt_v <- 3
scale_alt_v <- 2
# Setting 2 DGP: F0 = N(0,1), F1 = N(3,1); surge 0.05 -> 0.30 -> 0.05.
gs_mean_v <- 3

# ------------------------------------------------------------------------------
# 2. Diagnostic engine
# ------------------------------------------------------------------------------
# Per replicate: regenerate the beta-sweep stream (identical RNG order:
# theta_ini, z_ini, theta_stream, z_stream, z0_ref, z1_ref after set.seed(r)),
# fit the GAM ratio once, then run OMDRC_FROM_RATIO for every beta and record
# the per-t errors of hat(pi)_t and hat(L)_mdr,t against the truth.
run_beta_diag <- function(setting = c("setting1", "setting2"),
                          beta_list, mx, ini, n, alpha, reps,
                          pi_bounds, D_min) {
  setting <- match.arg(setting)

  # Prevalence path over warm-up + stream (k0_ref = NULL convention of the
  # beta sweep: path built over ini + mx, warm-up = first ini points).
  N_ref <- ini + mx
  pi_path <- if (setting == "setting1") {
    make_pi_drop(N_ref)
  } else {
    make_pi_surge(N_ref)
  }
  pi_ini    <- pi_path[seq_len(ini)]
  pi_stream <- pi_path[(ini + 1L):N_ref]

  n_cores <- max(1L, min(parallel::detectCores() - 8L, 40L))
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  on.exit({ try(stopCluster(cl), silent = TRUE); registerDoSEQ() }, add = TRUE)
  clusterExport(cl, c("setting", "beta_list", "mx", "ini", "n", "alpha",
                      "pi_ini", "pi_stream", "pi_bounds", "D_min",
                      "code_dir", "shape_alt_v", "scale_alt_v", "gs_mean_v"),
                envir = environment())

  result <- foreach(r = 1:reps, .packages = c("mgcv")) %dopar% {
    source(file.path(code_dir, "OMDRC.R"))
    set.seed(r)

    # ---- Warm-up + stream: identical draw order to the sweep scripts ----
    theta_ini <- rbinom(ini, size = 1, prob = pi_ini)
    if (setting == "setting1") {
      z_ini <- rexp(ini, rate = 1)
      z_ini[theta_ini == 1] <- rgamma(sum(theta_ini == 1),
                                      shape = shape_alt_v, scale = scale_alt_v)
      theta_stream <- rbinom(mx, size = 1, prob = pi_stream)
      z_stream <- rexp(mx, rate = 1)
      z_stream[theta_stream == 1] <- rgamma(sum(theta_stream == 1),
                                            shape = shape_alt_v,
                                            scale = scale_alt_v)
      r_true <- dgamma(z_stream, shape = shape_alt_v, scale = scale_alt_v) /
        dexp(z_stream, rate = 1)
      z0_ref <- rexp(n, rate = 1)
      z1_ref <- rgamma(n, shape = shape_alt_v, scale = scale_alt_v)
    } else {
      z_ini <- rnorm(ini)
      z_ini[theta_ini == 1] <- rnorm(sum(theta_ini == 1),
                                     mean = gs_mean_v, sd = 1)
      theta_stream <- rbinom(mx, size = 1, prob = pi_stream)
      z_stream <- rnorm(mx)
      z_stream[theta_stream == 1] <- rnorm(sum(theta_stream == 1),
                                           mean = gs_mean_v, sd = 1)
      r_true <- dnorm(z_stream, mean = gs_mean_v, sd = 1) /
        dnorm(z_stream, mean = 0, sd = 1)
      z0_ref <- rnorm(n)
      z1_ref <- rnorm(n, mean = gs_mean_v, sd = 1)
    }
    # True Lmdr_t built from the TRUE ratio and the TRUE pi_t.
    lmdr_true <- (pi_stream * r_true) / ((1 - pi_stream) + pi_stream * r_true)

    # ---- Stage 1 once per replicate: PC-DRE (GAM classifier) ratio ----
    ratio_model <- fit_ratio_gam(z0_ref, z1_ref, spline_k = 10L)
    LR_ini <- predict_ratio(ratio_model, matrix(z_ini, ncol = 1))
    LR     <- predict_ratio(ratio_model, matrix(z_stream, ncol = 1))

    # ---- Stage 2 per beta: growing-window EM + capacity rule ----
    pi_err <- matrix(NA_real_, length(beta_list), mx)   # hat(pi)_t - pi_t
    l_err  <- matrix(NA_real_, length(beta_list), mx)   # |hat(L)_t - L_t|
    for (b in seq_along(beta_list)) {
      fit <- OMDRC_FROM_RATIO(LR = LR, LR_ini = LR_ini, alpha = alpha,
                              pi_bounds = pi_bounds,
                              em_tol = 1e-6, em_max_iter = 50L,
                              D_mode = "growing", D_beta = beta_list[b],
                              D_min = D_min)
      pi_err[b, ] <- fit$pi_hat - pi_stream
      l_err[b, ]  <- abs(fit$Lmdr - lmdr_true)
    }
    list(pi_err = pi_err, l_err = l_err)
  }

  # Reassemble into per-beta (reps x mx) matrices and per-t summaries.
  summ <- lapply(seq_along(beta_list), function(b) {
    pe <- t(sapply(result, function(z) z$pi_err[b, ]))   # reps x mx
    le <- t(sapply(result, function(z) z$l_err[b, ]))
    n_r <- nrow(pe)
    bias_est <- colMeans(pe)
    bias_se  <- apply(pe, 2, stats::sd) / sqrt(n_r)
    mae_pi   <- colMeans(abs(pe))
    mae_pi_se <- apply(abs(pe), 2, stats::sd) / sqrt(n_r)
    mae_l    <- colMeans(le)
    mae_l_se <- apply(le, 2, stats::sd) / sqrt(n_r)
    data.frame(t = seq_len(ncol(pe)), beta = beta_list[b],
               bias = bias_est, bias_se = bias_se,
               mae_pi = mae_pi, mae_pi_se = mae_pi_se,
               mae_l = mae_l, mae_l_se = mae_l_se)
  })
  list(summary = do.call(rbind, summ), pi_stream = pi_stream)
}

run_diag_retry <- function(..., max_try = 3L) {
  for (a in seq_len(max_try)) {
    out <- tryCatch(run_beta_diag(...), error = function(e) {
      message("run_beta_diag failed (attempt ", a, "/", max_try, "): ",
              conditionMessage(e))
      NULL
    })
    if (!is.null(out)) return(out)
    Sys.sleep(5)
  }
  stop("run_beta_diag failed after ", max_try, " attempts")
}

# ------------------------------------------------------------------------------
# 3. Run (or reuse the compact per-t summary cache)
# ------------------------------------------------------------------------------
results_file <- file.path(code_dir, "sensity", "beta_error_diag_results.rds")
reuse_cache  <- TRUE   # set FALSE to rerun everything from scratch

cached <- if (reuse_cache && file.exists(results_file)) readRDS(results_file) else NULL
cache_ok <- !is.null(cached) &&
  isTRUE(all.equal(cached$beta_list, beta_list)) &&
  identical(cached$reps, reps_v) && identical(cached$mx, mx_v)

if (cache_ok) {
  cat("Reusing cached beta error diagnostics.\n")
  diag_s1 <- cached$diag_s1
  diag_s2 <- cached$diag_s2
} else {
  cat("\n====== beta error diagnostics: Setting 1 (Exp/Gamma, prior-drop) ======\n")
  diag_s1 <- run_diag_retry("setting1", beta_list = beta_list, mx = mx_v,
                            ini = default_ini, n = default_n, alpha = alpha_v,
                            reps = reps_v, pi_bounds = pi_bnds,
                            D_min = default_D_min)
  cat("\n====== beta error diagnostics: Setting 2 (Gaussian, surge) ======\n")
  diag_s2 <- run_diag_retry("setting2", beta_list = beta_list, mx = mx_v,
                            ini = default_ini, n = default_n, alpha = alpha_v,
                            reps = reps_v, pi_bounds = pi_bnds,
                            D_min = default_D_min)
  saveRDS(list(beta_list = beta_list, reps = reps_v, mx = mx_v,
               alpha = alpha_v, n = default_n, ini = default_ini,
               D_min = default_D_min, pi_bounds = pi_bnds,
               diag_s1 = diag_s1, diag_s2 = diag_s2),
          file = results_file)
}

# Terminal summary: time-averaged |bias| and MAEs per beta.
cat("\n--- Time-averaged estimation errors by beta ---\n")
for (nm in c("Setting 1", "Setting 2")) {
  dd <- if (nm == "Setting 1") diag_s1$summary else diag_s2$summary
  agg <- dd %>% group_by(beta) %>%
    summarise(mean_abs_bias = mean(abs(bias)), mean_mae_pi = mean(mae_pi),
              mean_mae_l = mean(mae_l), .groups = "drop")
  cat(nm, ":\n")
  for (i in seq_len(nrow(agg))) {
    cat(sprintf("  beta=%.1f  |bias(pi)|=%.4f  MAE(pi)=%.4f  MAE(Lmdr)=%.4f\n",
                agg$beta[i], agg$mean_abs_bias[i], agg$mean_mae_pi[i],
                agg$mean_mae_l[i]))
  }
}

# ------------------------------------------------------------------------------
# 4. Visualization (2 x 3: rows Setting 1 / Setting 2; cols bias / MAE / MAE)
# ------------------------------------------------------------------------------
custom_theme <- theme_bw() +
  theme(
    plot.subtitle = element_text(hjust = 0.5, size = 22, face = "bold"),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    legend.key.width = unit(1.0, "cm"),
    axis.text = element_text(size = 19, colour = "black"),
    axis.title = element_text(size = 22),
    panel.grid.major = element_line(colour = "grey92", size = 0.4),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.8)
  )
my_shapes <- c(16, 17, 15, 3, 4, 2, 7)
# Same beta colour convention as the sweep figures (Purples), but the oracle
# has zero estimation error by construction, so no OR curve is drawn.
colors_beta <- brewer.pal(9, "Purples")[3:9]
labels_beta <- parse(text = paste0("OMDRC.DD~(beta==", beta_list, ")"))
legend_guides <- guides(color = guide_legend(ncol = 4, byrow = TRUE))

# Markers are thinned to every 50th time point (the m_seq grid of the sweep
# figures); lines and ribbons use the full per-t resolution.
point_ts <- seq(100, mx_v, by = 50)

# Display-only smoothing: centred rolling mean over t (window w, partial at the
# edges) applied per beta to the per-t Monte-Carlo summaries. reps = 100 leaves
# visible per-t noise (especially in the Lmdr-MAE panels); the rolling mean
# reveals the systematic beta trend without touching the cached raw summaries.
roll_mean <- function(x, w = 25L) {
  kern <- rep(1, w)
  num <- stats::filter(x, kern, sides = 2)
  den <- stats::filter(rep(1, length(x)), kern, sides = 2)
  # stats::filter yields NA at the edges; recompute those with partial windows.
  out <- as.numeric(num / den)
  na_idx <- which(is.na(out))
  half <- w %/% 2
  for (i in na_idx) {
    j <- max(1L, i - half):min(length(x), i + half)
    out[i] <- mean(x[j])
  }
  out
}

make_err_panel <- function(df, yvar, sevar, ylab, tag, hline0 = FALSE,
                           clamp0 = TRUE) {
  df$est <- df[[yvar]]
  df$se  <- df[[sevar]]
  df <- df %>% group_by(beta) %>% arrange(t, .by_group = TRUE) %>%
    mutate(est = roll_mean(est), se = roll_mean(se)) %>% ungroup() %>%
    as.data.frame()
  df$lo <- df$est - 1.96 * df$se
  df$hi <- df$est + 1.96 * df$se
  if (clamp0) df$lo <- pmax(0, df$lo)
  df$Beta <- factor(df$beta, levels = beta_list)
  p <- ggplot(df, aes(t, est, color = Beta)) +
    geom_ribbon(aes(ymin = lo, ymax = hi, fill = Beta),
                alpha = 0.15, colour = NA, show.legend = FALSE)
  if (hline0) {
    p <- p + geom_hline(yintercept = 0, linetype = "dashed", size = 1.0)
  }
  p +
    geom_line(size = 1.0) +
    geom_point(data = df[df$t %in% point_ts, ],
               aes(shape = Beta), size = 2.4) +
    scale_color_manual(values = colors_beta, labels = labels_beta) +
    scale_fill_manual(values = colors_beta, labels = labels_beta) +
    scale_shape_manual(values = my_shapes, labels = labels_beta) +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.08))) +
    labs(x = "Time (t)", y = ylab, subtitle = paste0("(", tag, ")")) +
    custom_theme + legend_guides
}

make_diag_row <- function(diag, tag_prefix) {
  df <- diag$summary
  p1 <- make_err_panel(df, "bias", "bias_se",
                       expression(Bias~of~hat(pi)[t]),
                       paste0(tag_prefix, ".1"),
                       hline0 = TRUE, clamp0 = FALSE)
  p2 <- make_err_panel(df, "mae_pi", "mae_pi_se",
                       expression(MAE~of~hat(pi)[t]),
                       paste0(tag_prefix, ".2"))
  p3 <- make_err_panel(df, "mae_l", "mae_l_se",
                       expression(MAE~of~hat(L)[paste("mdr,", t)]),
                       paste0(tag_prefix, ".3"))
  (p1 | p2 | p3) + plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
}

row_s1 <- make_diag_row(diag_s1, "a")   # (a.*) Setting 1: prior-drop
row_s2 <- make_diag_row(diag_s2, "b")   # (b.*) Setting 2: surge

final_plot <- ggarrange(row_s1, row_s2, ncol = 1, nrow = 2)

if (interactive()) print(final_plot)
ggsave(file.path(code_dir, "sensity", "Rplot_beta_error_diag.pdf"),
       final_plot, width = 16, height = 11)
ggsave(file.path(code_dir, "sensity", "Rplot_beta_error_diag.png"),
       final_plot, width = 16, height = 11, dpi = 300)
