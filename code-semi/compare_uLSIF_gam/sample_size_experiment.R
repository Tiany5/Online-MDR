# Reference sample-size ablation of KDE vs uLSIF density-ratio estimation,
# reported over the whole online time line and across TWO data-generating
# settings, yielding a 2x2 comparison (rows = settings, columns = MDR / FDR).
#
# The labeled reference sample size is varied,
#     n0 = n1 in {500, 1000, 2000, 4000},
# to expose the causal chain  more reference data ==> better ratio ==> MDR down.
# uLSIF is nonparametric, so its capacity grows with n (sqrt(n) number of
# Gaussian centers, refined bandwidth grid, weaker ridge), giving a valid
# asymptotic setup rather than a fixed-capacity ceiling.
#
# Two settings (both under the same decreasing prior-drop pi_t path):
#   * "multimodal"     : F0 = N(0,1), F1 = 0.5 N(-3,0.7^2) + 0.5 N(2,0.7^2)
#                        (separated bimodal alternative; hard for uLSIF).
#   * "gaussian_shift" : F0 = N(0,1), F1 = N(gs_mean, 1) with gs_mean = 3
#                        (smooth monotone likelihood ratio; strengthened signal
#                        so uLSIF can control MDR at the largest sample size).
#
# KDE uses an orange colour family, uLSIF a blue family; within each family the
# shade goes from light (small n) to dark (large n) and the point shape also
# changes with n, so families stay separated and the progression is visible.
#
# Required shared file: code-semi/OMDRC.R
# Required packages: foreach, doParallel, ggplot2, dplyr, tidyr, patchwork

suppressPackageStartupMessages({
  library(foreach)
  library(doParallel)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0L) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1L]))))
  }
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    active_file <- rstudioapi::getActiveDocumentContext()$path
    if (nzchar(active_file)) return(dirname(normalizePath(active_file)))
  }
  getwd()
}

code_dir <- get_script_dir()
omdrc_file <- file.path(dirname(code_dir), "OMDRC.R")
if (!file.exists(omdrc_file)) {
  stop("Cannot find the shared code-semi/OMDRC.R library.")
}
source(omdrc_file)

# Fixed decreasing prior-drop path (shared by both settings).
make_pi_path <- function(N, pi_low = 0.08, pi_high = 0.30,
                         drop_center = 0.30, drop_width = 0.08) {
  make_smooth_prior_drop(N = N, pi_low = pi_low, pi_high = pi_high,
                         drop_center = drop_center, drop_width = drop_width)
}

# Nonparametric capacity schedule: number of Gaussian centers grows as sqrt(n).
n_centers_schedule <- function(s) {
  min(as.integer(s), as.integer(round(100 * sqrt(s / 200))))
}

# -----------------------------------------------------------------------------
# Monte Carlo experiment over reference sample sizes for one setting.
# -----------------------------------------------------------------------------
exp_sample_size <- function(setting, sizes, N, ini, alpha, reps, D,
                            time_grid = NULL,
                            pi_low = 0.08, pi_high = 0.30,
                            drop_center = 0.30, drop_width = 0.08,
                            pi_bounds = c(0.01, 0.99),
                            gs_mean = 3,
                            kde_control = list(),
                            ulsif_control = list(),
                            clf_control = list(),
                            sigma_multipliers = c(0.0625, 0.125, 0.25, 0.5, 1, 2),
                            lambda_grid = 10^seq(-8, -1, length.out = 8L),
                            n_cores = 1L,
                            seed = 20260723) {
  sizes <- as.integer(sizes)
  n_sizes <- length(sizes)
  if (is.null(time_grid)) time_grid <- seq(100L, N, by = 100L)
  time_grid <- as.integer(time_grid)
  n_time <- length(time_grid)
  # Each setting keeps its characteristic prevalence path:
  #   multimodal     -> monotone decreasing prior-drop,
  #   gaussian_shift -> Setting-3 non-monotone Gaussian surge (low->peak->low).
  pi_full <- if (setting == "multimodal") {
    make_pi_path(ini + N, pi_low, pi_high, drop_center, drop_width)
  } else {
    u <- seq(0, 1, length.out = as.integer(ini + N))
    pmin(pmax(0.05 + (0.30 - 0.05) *
                exp(-((u - 0.60)^2) / (2 * 0.15^2)), 0.05), 0.30)
  }
  true_pi_ini <- pi_full[seq_len(ini)]
  true_pi <- pi_full[ini + seq_len(N)]

  n_cores <- max(1L, as.integer(n_cores))
  if (n_cores > 1L) {
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
    on.exit({
      try(parallel::stopCluster(cl), silent = TRUE)
      foreach::registerDoSEQ()
    }, add = TRUE)
  } else {
    foreach::registerDoSEQ()
  }

  result <- foreach(
    r = seq_len(reps),
    .export = c("omdrc_file", "n_centers_schedule")
  ) %dopar% {
    source(omdrc_file)
    set.seed(seed + r)

    # Setting-specific data generators and oracle densities.
    if (setting == "multimodal") {
      .mu1 <- 2; .mu2 <- -3; .sd1 <- 0.7
      gen_null <- function(n) rnorm(n, 0, 1)
      gen_alt <- function(n) {
        cc <- rbinom(n, 1, 0.5)
        rnorm(n, ifelse(cc == 1L, .mu1, .mu2), .sd1)
      }
      dens0 <- function(x) dnorm(x, 0, 1)
      dens1 <- function(x) 0.5 * dnorm(x, .mu1, .sd1) + 0.5 * dnorm(x, .mu2, .sd1)
    } else {
      .am <- gs_mean; .asd <- 1
      gen_null <- function(n) rnorm(n, 0, 1)
      gen_alt <- function(n) rnorm(n, .am, .asd)
      dens0 <- function(x) dnorm(x, 0, 1)
      dens1 <- function(x) dnorm(x, .am, .asd)
    }

    # Evaluation stream (fixed within a replicate across all sample sizes).
    theta <- rbinom(N, size = 1, prob = true_pi)
    z <- ifelse(theta == 0L, gen_null(N), gen_alt(N))
    theta_ini <- rbinom(ini, size = 1, prob = true_pi_ini)
    z_ini <- ifelse(theta_ini == 0L, gen_null(ini), gen_alt(ini))

    f0 <- dens0(z); f1 <- dens1(z)
    lmdr_oracle <- true_pi * f1 / ((1 - true_pi) * f0 + true_pi * f1)
    sig_t <- cumsum(theta)[time_grid]

    fn_kde <- matrix(0, n_sizes, n_time); fn_ulsif <- matrix(0, n_sizes, n_time)
    fp_kde <- matrix(0, n_sizes, n_time); fp_ulsif <- matrix(0, n_sizes, n_time)
    disc_kde <- matrix(0, n_sizes, n_time); disc_ulsif <- matrix(0, n_sizes, n_time)
    mae_q_kde <- numeric(n_sizes); mae_q_ulsif <- numeric(n_sizes)
    rmse_pi_kde <- numeric(n_sizes); rmse_pi_ulsif <- numeric(n_sizes)
    fn_gam <- matrix(0, n_sizes, n_time); fp_gam <- matrix(0, n_sizes, n_time)
    disc_gam <- matrix(0, n_sizes, n_time)
    mae_q_gam <- numeric(n_sizes); rmse_pi_gam <- numeric(n_sizes)

    for (j in seq_along(sizes)) {
      s <- sizes[j]
      z0_ref <- gen_null(s)
      z1_ref <- gen_alt(s)

      fit_eval <- function(method, extra) {
        out <- OMDRC_DD(
          z = z, z_ini = z_ini, z0 = z0_ref, z1 = z1_ref,
          alpha = alpha,
          ratio_method = method, ratio_control = extra,
          pi_bounds = pi_bounds, em_tol = 1e-6, em_max_iter = 50L,
          D_mode = "growing", D_beta = 0.6, D_min = 10L
        )
        de <- out$de
        list(
          fn_t = cumsum(theta * (1 - de))[time_grid],
          fp_t = cumsum((1 - theta) * de)[time_grid],
          disc_t = cumsum(de)[time_grid],
          mae_q = mean(abs(out$Lmdr - lmdr_oracle)),
          rmse_pi = sqrt(mean((out$pi_hat - true_pi)^2))
        )
      }

      m_kde <- fit_eval("kde", kde_control)
      base_sigma <- .median_distance(rbind(z0_ref, z1_ref))
      u_ctrl <- modifyList(
        ulsif_control,
        list(
          n_centers = n_centers_schedule(s),
          sigma_grid = base_sigma * sigma_multipliers,
          lambda_grid = lambda_grid,
          seed = seed + 100000L + r + j
        )
      )
      m_uls <- fit_eval("ulsif", u_ctrl)

      # PC-DRE (probabilistic-classifier density-ratio estimator, mgcv GAM) uses
      # the SAME reference samples z0_ref, z1_ref and the SAME evaluation stream
      # z; the only difference from KDE/uLSIF is the density-ratio estimator.
      # Comparisons are therefore paired per replicate.
      m_clf <- fit_eval("gam", clf_control)

      fn_kde[j, ] <- m_kde$fn_t;   fn_ulsif[j, ] <- m_uls$fn_t
      fp_kde[j, ] <- m_kde$fp_t;   fp_ulsif[j, ] <- m_uls$fp_t
      disc_kde[j, ] <- m_kde$disc_t; disc_ulsif[j, ] <- m_uls$disc_t
      mae_q_kde[j] <- m_kde$mae_q; mae_q_ulsif[j] <- m_uls$mae_q
      rmse_pi_kde[j] <- m_kde$rmse_pi; rmse_pi_ulsif[j] <- m_uls$rmse_pi
      fn_gam[j, ] <- m_clf$fn_t; fp_gam[j, ] <- m_clf$fp_t
      disc_gam[j, ] <- m_clf$disc_t
      mae_q_gam[j] <- m_clf$mae_q; rmse_pi_gam[j] <- m_clf$rmse_pi
    }

    list(
      sig_t = sig_t,
      fn_kde = fn_kde, fn_ulsif = fn_ulsif, fn_gam = fn_gam,
      fp_kde = fp_kde, fp_ulsif = fp_ulsif, fp_gam = fp_gam,
      disc_kde = disc_kde, disc_ulsif = disc_ulsif, disc_gam = disc_gam,
      mae_q_kde = mae_q_kde, mae_q_ulsif = mae_q_ulsif, mae_q_gam = mae_q_gam,
      rmse_pi_kde = rmse_pi_kde, rmse_pi_ulsif = rmse_pi_ulsif,
      rmse_pi_gam = rmse_pi_gam
    )
  }

  if (n_cores > 1L) {
    parallel::stopCluster(cl)
    foreach::registerDoSEQ()
  }

  sig_all <- t(vapply(result, function(x) x$sig_t, numeric(n_time)))

  ratio_timeline <- function(num_key, den_all = NULL, den_key = NULL) {
    rows <- list()
    for (method in c("kde", "ulsif", "gam")) {
      nk <- paste0(num_key, "_", method)
      for (j in seq_along(sizes)) {
        num <- t(vapply(result, function(x) x[[nk]][j, ], numeric(n_time)))
        den <- if (is.null(den_key)) {
          den_all
        } else {
          t(vapply(result, function(x) x[[paste0(den_key, "_", method)]][j, ],
                   numeric(n_time)))
        }
        mean_num <- colMeans(num)
        mean_den <- pmax(colMeans(den), .Machine$double.eps)
        est <- mean_num / mean_den
        infl <- num - sweep(den, 2L, est, "*")
        se <- apply(infl, 2L, stats::sd) / (sqrt(reps) * mean_den)
        rows[[length(rows) + 1L]] <- data.frame(
          method = method, size = sizes[j], t = time_grid,
          value = est, lower = pmax(0, est - 1.96 * se),
          upper = pmin(1, est + 1.96 * se)
        )
      }
    }
    dplyr::bind_rows(rows)
  }

  timeline_df <- dplyr::bind_rows(
    ratio_timeline("fn", den_all = sig_all) |> dplyr::mutate(metric = "MDR"),
    ratio_timeline("fp", den_key = "disc") |> dplyr::mutate(metric = "FDR")
  ) |> dplyr::mutate(setting = setting)

  summarize_scalar <- function(method) {
    grab <- function(nm) {
      matrix(t(vapply(result, function(x) x[[nm]], numeric(n_sizes))),
             nrow = reps, ncol = n_sizes)
    }
    mae_q <- grab(paste0("mae_q_", method))
    rmse_pi <- grab(paste0("rmse_pi_", method))
    fn_mat <- matrix(
      t(vapply(result, function(x) x[[paste0("fn_", method)]][, n_time],
               numeric(n_sizes))), nrow = reps, ncol = n_sizes)
    sig_T <- sig_all[, n_time]
    mdr_T <- colMeans(fn_mat) / pmax(mean(sig_T), .Machine$double.eps)
    mean_se <- function(mat) {
      list(est = colMeans(mat), se = apply(mat, 2L, stats::sd) / sqrt(reps))
    }
    q_stat <- mean_se(mae_q); p_stat <- mean_se(rmse_pi)
    dplyr::bind_rows(
      data.frame(size = sizes, metric = "MDR_T", est = mdr_T, se = NA_real_),
      data.frame(size = sizes, metric = "MAE_q",
                 est = q_stat$est, se = q_stat$se),
      data.frame(size = sizes, metric = "RMSE_pi",
                 est = p_stat$est, se = p_stat$se)
    ) |> dplyr::mutate(method = method, setting = setting)
  }
  summary_df <- dplyr::bind_rows(
    summarize_scalar("kde"), summarize_scalar("ulsif"), summarize_scalar("gam")
  )

  list(timeline = timeline_df, summary = summary_df,
       sizes = sizes, time_grid = time_grid, setting = setting)
}

# -----------------------------------------------------------------------------
# Configuration and run (two settings)
# -----------------------------------------------------------------------------
settings <- c("multimodal", "gaussian_shift")
setting_labels <- c(multimodal = "Multimodal",
                    gaussian_shift = "Gaussian shift")
sizes <- c(500L, 1000L, 2000L, 4000L)
N <- 1000L
ini <- 500L
alpha <- 0.1
reps <- 150L
D <- 100L        # legacy fixed-window fallback (UNUSED: DD now uses the paper
                 # growing window D_t = min{K0+t, max{10, (K0+t)^0.6}})
pi_bounds <- c(0.01, 0.99)
gs_mean <- 2.5   # Gaussian-shift alt mean (signal strength): N(gs_mean, 1)
available_cores <- parallel::detectCores(logical = TRUE)
if (!is.finite(available_cores)) available_cores <- 1L
n_cores <- max(1L, min(10L, as.integer(available_cores) - 2L, reps))

kde_control <- list(
  grid_n = 2048L, pad = 5, bandwidth_multiplier = 1.5,
  density_floor = 1e-8, ratio_floor = 1e-8, ratio_cap = 1e6
)
ulsif_control <- list(
  n_folds = 5L, tune = TRUE, normalize = TRUE,
  calibration_fraction = 0.20, min_calibration = 50L,
  ratio_floor = 1e-8, ratio_cap = 1e6
)
clf_control <- list(
  spline_k = 10L, prob_floor = 1e-6,
  ratio_floor = 1e-8, ratio_cap = 1e6
)

# Only settings listed in `recompute` are (re)run; the others are reused from
# the cached results file (useful when one setting is unchanged). Empty vector
# = plot-only mode: both settings are read from the cache and the figure is
# rebuilt without any Monte-Carlo rerun.
recompute <- c()
cache_file <- file.path(code_dir, "samplesize_kde_ulsif_gam_results.rds")
cached <- if (file.exists(cache_file)) readRDS(cache_file) else NULL

run_one <- function(st) {
  reuse <- !is.null(cached) && !(st %in% recompute) &&
    st %in% unique(cached$timeline$setting)
  if (reuse) {
    cat(sprintf("Reusing cached setting: %s\n", st))
    list(
      timeline = dplyr::filter(cached$timeline, setting == st),
      summary  = dplyr::filter(cached$summary, setting == st)
    )
  } else {
    cat(sprintf("Running setting: %s\n", st))
    exp_sample_size(
      setting = st, sizes = sizes, N = N, ini = ini, alpha = alpha,
      reps = reps, D = D, pi_bounds = pi_bounds, gs_mean = gs_mean,
      kde_control = kde_control,
      ulsif_control = ulsif_control, clf_control = clf_control,
      n_cores = n_cores
    )
  }
}
out_list <- lapply(settings, run_one)
names(out_list) <- settings

timeline_all <- dplyr::bind_rows(lapply(out_list, `[[`, "timeline"))
summary_all <- dplyr::bind_rows(lapply(out_list, `[[`, "summary"))

# Console summary (terminal values).
method_label <- c(kde = "KDE", ulsif = "uLSIF", gam = "PC-DRE")
wide <- summary_all |>
  dplyr::mutate(method = method_label[method]) |>
  tidyr::pivot_wider(names_from = metric, values_from = c(est, se)) |>
  dplyr::arrange(setting, method, size)
cat(sprintf("\nreps = %d, N = %d, ini = %d, D_t = min{K0+t, max{10, (K0+t)^0.6}}, alpha = %.2f\n",
            reps, N, ini, alpha))
for (i in seq_len(nrow(wide))) {
  cat(sprintf(
    "%-14s %-6s n=%4d | MDR_T=%.3f  MAE(q)=%.4f  RMSE(pi)=%.4f\n",
    wide$setting[i], wide$method[i], wide$size[i],
    wide$est_MDR_T[i], wide$est_MAE_q[i], wide$est_RMSE_pi[i]
  ))
}

# -----------------------------------------------------------------------------
# 2x2 plot: rows = settings, columns = MDR / FDR
# -----------------------------------------------------------------------------
custom_theme <- theme_bw() +
  theme(
    plot.subtitle = element_text(size = 18, hjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    axis.text = element_text(size = 14, colour = "black"),
    axis.title = element_text(size = 16),
    panel.grid.minor = element_blank()
  )

orange_shades <- c("#fdbe85", "#fd8d3c", "#e6550d", "#a63603")
blue_shades   <- c("#9ecae1", "#4292c6", "#2171b5", "#08306b")
green_shades  <- c("#a1d99b", "#74c476", "#31a354", "#006d2c")
size_levels <- sort(unique(sizes))
group_levels <- c(paste0("KDE (n=", size_levels, ")"),
                  paste0("uLSIF (n=", size_levels, ")"),
                  paste0("PC-DRE (n=", size_levels, ")"))
group_colors <- setNames(
  c(orange_shades[seq_along(size_levels)],
    blue_shades[seq_along(size_levels)],
    green_shades[seq_along(size_levels)]),
  group_levels
)
group_shapes <- setNames(rep(c(16, 17, 15, 18)[seq_along(size_levels)], 3L),
                         group_levels)

# Legend labels are rendered with plotmath so that the reference sample size is
# displayed as "n_ref = <size>" matching the paper's n_{\mathrm{ref}}: the
# quoted subscript n["ref"] renders "ref" UPRIGHT.  Method names are quoted so
# that the hyphen in "PC-DRE" is not parsed as a minus sign; the "=" signs are
# quoted strings joined by `~` because chained `==` is not valid R syntax.
group_label_expr <- setNames(
  sprintf('"%s"~(n["ref"]~"="~%d)',
          rep(c("KDE", "uLSIF", "PC-DRE"), each = length(size_levels)),
          rep(as.integer(size_levels), 3L)),
  group_levels
)
legend_labeller <- function(v) parse(text = group_label_expr[as.character(v)])

tl_df <- timeline_all |>
  dplyr::mutate(
    method = method_label[method],
    grp = factor(paste0(method, " (n=", size, ")"), levels = group_levels),
    metric = factor(metric, levels = c("MDR", "FDR"))
  )

# Only display the time line up to plot_t_max (matches the article's T = 1000).
plot_t_max <- 1000

make_panel <- function(setting_key, metric_key, subtitle, hline = NULL) {
  d <- dplyr::filter(tl_df, setting == setting_key, metric == metric_key,
                     t <= plot_t_max)
  g <- ggplot(d, aes(x = t, y = value, color = grp, shape = grp, group = grp)) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = grp),
                alpha = 0.16, colour = NA, show.legend = FALSE) +
    geom_line(size = 0.8) +
    geom_point(size = 1.4) +
    scale_color_manual(values = group_colors, name = NULL,
                       labels = legend_labeller) +
    scale_fill_manual(values = group_colors, name = NULL,
                      labels = legend_labeller) +
    scale_shape_manual(values = group_shapes, name = NULL,
                       labels = legend_labeller) +
    labs(subtitle = subtitle, x = "Time (t)", y = metric_key) +
    custom_theme
  if (!is.null(hline)) {
    g <- g + geom_hline(yintercept = hline, linetype = "dashed",
                        size = 0.8, colour = "black")
  }
  g
}

g_2x2 <- (
  make_panel("multimodal", "MDR", "(a.1)", hline = alpha) +
  make_panel("multimodal", "FDR", "(a.2)")
) / (
  make_panel("gaussian_shift", "MDR", "(b.1)", hline = alpha) +
  make_panel("gaussian_shift", "FDR", "(b.2)")
) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom") &
  guides(color = guide_legend(nrow = 3, byrow = TRUE),
         shape = guide_legend(nrow = 3, byrow = TRUE))

if (interactive()) print(g_2x2)

# -----------------------------------------------------------------------------
# Save outputs (rds, csv, png, pdf).
# -----------------------------------------------------------------------------
saveRDS(list(timeline = timeline_all, summary = summary_all,
             sizes = sizes, settings = settings),
        file.path(code_dir, "samplesize_kde_ulsif_gam_results.rds"))
utils::write.csv(
  wide, file.path(code_dir, "samplesize_kde_ulsif_gam_summary.csv"),
  row.names = FALSE
)
for (ext in c("png", "pdf")) {
  ggsave(
    file.path(code_dir, paste0("samplesize_kde_ulsif_gam_2x2_timeline.", ext)),
    g_2x2, width = 13, height = 9, dpi = 300
  )
}
