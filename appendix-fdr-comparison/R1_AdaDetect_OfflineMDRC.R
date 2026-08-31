############################################################
## Offline MDRC vs AdaDetect
## Two offline comparisons:
##   (1) fixed pi, varying signal strength s
##   (2) fixed s, varying signal proportion pi
##
## Data model:
##   theta_i ~ Bernoulli(pi)
##   X_i | theta_i = 0 ~ N(0, 1)
##   X_i | theta_i = 1 ~ 0.5 N(2s, 0.7^2) + 0.5 N(-3s, 0.7^2)
##
## Important:
##   No oracle likelihood-ratio score is given to either method.
##
## Offline MDRC (revised, via code-semi/OMDRC.R):
##   estimates r = f1/f0 with the PC-DRE logistic-GAM classifier from an
##   independent labeled null sample and an independent labeled signal
##   sample, estimates a constant signal proportion from an independent
##   unlabeled calibration batch, forms Lmdr, and applies the offline
##   Lmdr rule at alpha = 0.1.
##
## AdaDetect:
##   uses split null samples and a data-adaptive KDE score,
##   then empirical conformal p-values + BH at q = 0.1.
############################################################

rm(list = ls())
gc()

suppressPackageStartupMessages({
  library(foreach)
  library(doParallel)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

############################################################
## 0. Locate this script, source the shared OMDRC.R engine,
##    and set the output directory to the script folder
############################################################

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

script_dir <- get_script_dir()

## This release stores the script one level below the repository root.
omdrc_file <- normalizePath(
  file.path(script_dir, "..", "code-semi", "OMDRC.R"),
  mustWork = FALSE
)
if (!file.exists(omdrc_file)) {
  stop("Cannot find code-semi/OMDRC.R relative to the script at: ", omdrc_file)
}
source(omdrc_file)

output_dir <- script_dir
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
cat("Script dir:", script_dir, "\n")
cat("OMDRC engine:", omdrc_file, "\n")
cat("Output directory:", output_dir, "\n")

############################################################
## 1. Robust KDE helper
############################################################

safe_bw <- function(x) {
  x <- x[is.finite(x)]
  if (length(unique(x)) < 3 || sd(x) <= 1e-12) {
    return(0.2)
  }
  h <- tryCatch(
    bw.SJ(x),
    error = function(e) NA_real_,
    warning = function(w) NA_real_
  )
  if (!is.finite(h) || h <= 0) {
    h <- bw.nrd0(x)
  }
  if (!is.finite(h) || h <= 0) {
    h <- 0.2
  }
  h
}

kde_density_on_grid <- function(x, from, to, grid_n = 1000, bw_mult = 1.2) {
  h <- safe_bw(x) * bw_mult
  d <- density(
    x,
    bw = h,
    from = from,
    to = to,
    n = grid_n,
    kernel = "gaussian"
  )
  d
}

eval_density <- function(den, xout, floor_val = 1e-8) {
  y <- approx(den$x, den$y, xout = xout, rule = 2)$y
  pmax(y, floor_val)
}

############################################################
## 2. Offline MDRC via the revised data-driven engine (OMDRC.R)
##
##   1. r_hat = f1 / f0 is fit with the paper's primary PC-DRE
##      estimator (logistic-GAM classifier, fit_ratio_gam) from an
##      independent labeled null sample x0 and labeled signal sample x1.
##   2. A single constant signal proportion pi_hat is estimated from an
##      independent unlabeled calibration batch x_cal by the local
##      ratio-based EM (.estimate_local_pi_ratio).
##   3. Lmdr(x) = pi_hat r_hat(x) / (1 - pi_hat + pi_hat r_hat(x)).
##   4. The offline Lmdr rule OMDRC_OFF applies the level-alpha
##      MDR-controlling threshold.
##
## fit_ratio_gam, predict_ratio, .estimate_local_pi_ratio,
## .ratio_to_lmdr and OMDRC_OFF are all provided by code-semi/OMDRC.R.
############################################################

MDRC_OFF_OMDRC <- function(x_test,
                           x0,
                           x1,
                           x_cal,
                           alpha,
                           pi_bounds = c(0.005, 0.5),
                           spline_k = 10L) {

  ratio_model <- fit_ratio_gam(
    x0 = matrix(as.numeric(x0), ncol = 1L),
    x1 = matrix(as.numeric(x1), ncol = 1L),
    spline_k = spline_k
  )

  ratio_cal <- predict_ratio(ratio_model, matrix(as.numeric(x_cal), ncol = 1L))
  pi_fit <- .estimate_local_pi_ratio(
    ratio = ratio_cal,
    init = mean(pi_bounds),
    pi_bounds = pi_bounds
  )
  pi_hat <- pi_fit$pi

  ratio_test <- predict_ratio(ratio_model, matrix(as.numeric(x_test), ncol = 1L))
  Lmdr_test <- .ratio_to_lmdr(ratio_test, pi_hat)

  decision <- OMDRC_OFF(Lmdr_test, alpha)

  list(
    de = decision$de,
    lambda = decision$lambda,
    pi_hat = pi_hat,
    Lmdr = Lmdr_test,
    ratio = ratio_test
  )
}

############################################################
## 3. BH procedure
############################################################

BH_reject <- function(p_value, alpha) {
  m <- length(p_value)
  o <- order(p_value)
  p_sorted <- p_value[o]
  
  k_set <- which(p_sorted <= alpha * seq_len(m) / m)
  k <- if (length(k_set) == 0) 0 else max(k_set)
  
  reject <- rep(0L, m)
  if (k > 0) {
    reject[o[seq_len(k)]] <- 1L
  }
  
  reject
}

############################################################
## 4. Empirical conformal p-values
############################################################

empirical_conformal_p <- function(score_test, score_cal) {
  score_cal_sorted <- sort(score_cal)
  ell <- length(score_cal_sorted)
  
  n_leq <- findInterval(
    score_test,
    score_cal_sorted,
    rightmost.closed = TRUE
  )
  
  n_greater <- ell - n_leq
  p_value <- (1 + n_greater) / (ell + 1)
  
  pmin(pmax(p_value, 0), 1)
}

############################################################
## 5. AdaDetect-style data-driven KDE score
############################################################

AdaDetect_KDE <- function(x_test,
                          x_null,
                          alpha,
                          split_prop = 0.5,
                          grid_n = 1000,
                          pad = 3) {
  
  n_null <- length(x_null)
  k <- floor(split_prop * n_null)
  k <- max(20, min(k, n_null - 20))
  
  x_null_train <- x_null[seq_len(k)]
  x_null_cal <- x_null[(k + 1):n_null]
  x_mixed <- c(x_null_cal, x_test)
  
  all_x <- c(x_null_train, x_mixed)
  from <- min(all_x) - pad
  to <- max(all_x) + pad
  
  den_null <- kde_density_on_grid(
    x = x_null_train,
    from = from,
    to = to,
    grid_n = grid_n,
    bw_mult = 1.2
  )
  
  den_mix <- kde_density_on_grid(
    x = x_mixed,
    from = from,
    to = to,
    grid_n = grid_n,
    bw_mult = 1.2
  )
  
  score_fun <- function(z) {
    f_mix <- eval_density(den_mix, z)
    f_null <- eval_density(den_null, z)
    log(f_mix) - log(f_null)
  }
  
  score_cal <- score_fun(x_null_cal)
  score_test <- score_fun(x_test)
  
  p_value <- empirical_conformal_p(score_test, score_cal)
  de <- BH_reject(p_value, alpha)
  
  list(
    de = de,
    p_value = p_value,
    score_test = score_test,
    score_cal = score_cal
  )
}

############################################################
## 6. Metrics
############################################################

compute_metrics <- function(theta, de) {
  theta <- as.integer(theta)
  de <- as.integer(de)
  
  n_signal <- sum(theta)
  n_rej <- sum(de)
  
  c(
    MDR = sum(theta * (1 - de)) / max(n_signal, 1),
    FDR = sum((1 - theta) * de) / max(n_rej, 1),
    TDR = sum(theta * de) / max(n_signal, 1),
    Discoveries = n_rej,
    False_Pos = sum((1 - theta) * de),
    Missed = sum(theta * (1 - de)),
    Signals = n_signal
  )
}

############################################################
## 7. Data generators
##
##   theta ~ Bernoulli(pi)
##   X | theta = 0 ~ N(0, 1)
##   X | theta = 1 ~ 0.5 N(2s, 0.7^2) + 0.5 N(-3s, 0.7^2)
############################################################

draw_signal <- function(n,
                        s,
                        mu1_base = 2,
                        mu2_base = -3,
                        sd1 = 0.7) {
  comp <- rbinom(n, size = 1, prob = 0.5)
  mu <- ifelse(comp == 1, mu1_base * s, mu2_base * s)
  rnorm(n, mean = mu, sd = sd1)
}

draw_mixture <- function(n,
                         pi_signal,
                         s,
                         mu1_base = 2,
                         mu2_base = -3,
                         sd1 = 0.7) {
  theta <- rbinom(n, size = 1, prob = pi_signal)
  z_null <- rnorm(n, mean = 0, sd = 1)
  z_signal <- draw_signal(n, s, mu1_base, mu2_base, sd1)
  x <- ifelse(theta == 0, z_null, z_signal)
  list(x = x, theta = theta)
}

############################################################
## 8. One experiment over a parameter grid
############################################################

run_grid_experiment <- function(grid_values,
                                grid_name = c("s", "pi"),
                                m = 2000,
                                alpha_mdr = 0.1,
                                q_adadetect = 0.1,
                                reps = 200,
                                n0 = 500,
                                n1 = 500,
                                K0 = 500,
                                ncores = 10,
                                fixed_pi = 0.1,
                                fixed_s = 1.4,
                                grid_n = 1000) {
  
  grid_name <- match.arg(grid_name)
  
  metric_names <- c("MDR", "FDR", "TDR", "Discoveries", "False_Pos",
                    "Missed", "Signals")
  G <- length(grid_values)
  
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  on.exit({
    try(stopCluster(cl), silent = TRUE)
    registerDoSEQ()
  }, add = TRUE)
  
  result <- foreach(
    r = seq_len(reps),
    .packages = c("stats", "mgcv"),
    .export = c(
      "omdrc_file",
      "safe_bw",
      "kde_density_on_grid",
      "eval_density",
      "MDRC_OFF_OMDRC",
      "BH_reject",
      "empirical_conformal_p",
      "AdaDetect_KDE",
      "compute_metrics",
      "draw_signal",
      "draw_mixture"
    )
  ) %dopar% {
    
    source(omdrc_file)
    set.seed(20260708 + r)
    
    off_mat <- matrix(NA_real_, nrow = G, ncol = length(metric_names))
    ada_mat <- matrix(NA_real_, nrow = G, ncol = length(metric_names))
    
    colnames(off_mat) <- metric_names
    colnames(ada_mat) <- metric_names
    
    for (g in seq_along(grid_values)) {
      
      if (grid_name == "s") {
        s_now <- grid_values[g]
        pi_now <- fixed_pi
      } else {
        s_now <- fixed_s
        pi_now <- grid_values[g]
      }
      
      ## Test stream: m observations from the mixture model.
      test <- draw_mixture(n = m, pi_signal = pi_now, s = s_now)
      
      ## Shared independent null reference sample (both methods), size n0.
      x0 <- rnorm(n0, mean = 0, sd = 1)
      
      ## Offline MDRC only: labeled signal reference (n1) and an
      ## independent unlabeled calibration batch (K0) drawn from the
      ## same mixture as the test stream.
      x1 <- draw_signal(n = n1, s = s_now)
      x_cal <- draw_mixture(n = K0, pi_signal = pi_now, s = s_now)$x
      
      decision_off <- MDRC_OFF_OMDRC(
        x_test = test$x,
        x0 = x0,
        x1 = x1,
        x_cal = x_cal,
        alpha = alpha_mdr
      )
      
      decision_ada <- AdaDetect_KDE(
        x_test = test$x,
        x_null = x0,
        alpha = q_adadetect,
        split_prop = 0.5,
        grid_n = grid_n
      )
      
      off_mat[g, ] <- compute_metrics(test$theta, decision_off$de)
      ada_mat[g, ] <- compute_metrics(test$theta, decision_ada$de)
    }
    
    list(off = off_mat, ada = ada_mat)
  }
  
  arr_off <- array(NA_real_, dim = c(reps, G, length(metric_names)))
  arr_ada <- array(NA_real_, dim = c(reps, G, length(metric_names)))
  
  dimnames(arr_off) <- list(NULL, NULL, metric_names)
  dimnames(arr_ada) <- list(NULL, NULL, metric_names)
  
  for (r in seq_len(reps)) {
    arr_off[r, , ] <- result[[r]]$off
    arr_ada[r, , ] <- result[[r]]$ada
  }
  
  list(
    off = arr_off,
    ada = arr_ada,
    grid_values = grid_values,
    grid_name = grid_name,
    m = m,
    alpha_mdr = alpha_mdr,
    q_adadetect = q_adadetect,
    reps = reps,
    n0 = n0,
    n1 = n1,
    K0 = K0,
    fixed_pi = fixed_pi,
    fixed_s = fixed_s
  )
}

############################################################
## 9. Run two experiments
############################################################

m <- 2000

alpha <- 0.1
q_adadetect <- 0.1

reps <- 200
ncores <- 10

## Shared null reference n0; offline-MDRC-only labeled signal n1 and
## unlabeled calibration batch K0.
n0 <- 500
n1 <- 500
K0 <- 500

s_grid <- seq(0.4, 2.2, by = 0.1)
pi_grid <- seq(0.02, 0.45, by = 0.03)

rds_s <- file.path(output_dir, "OfflineMDRC_vs_AdaDetect_vary_s.rds")
rds_pi <- file.path(output_dir, "OfflineMDRC_vs_AdaDetect_vary_pi_s1p4.rds")

## Set FORCE_RECOMPUTE <- TRUE to rerun the (~10 min) Monte Carlo simulation;
## otherwise reuse the cached .rds and only rebuild the summary + figure.
FORCE_RECOMPUTE <- FALSE

if (!FORCE_RECOMPUTE && file.exists(rds_s) && file.exists(rds_pi)) {
  cat("Reusing cached simulation results",
      "(set FORCE_RECOMPUTE <- TRUE to rerun).\n")
  out_s <- readRDS(rds_s)
  out_pi <- readRDS(rds_pi)
} else {
  system.time({
    out_s <- run_grid_experiment(
      grid_values = s_grid,
      grid_name = "s",
      m = m,
      alpha_mdr = alpha,
      q_adadetect = q_adadetect,
      reps = reps,
      n0 = n0,
      n1 = n1,
      K0 = K0,
      ncores = ncores,
      fixed_pi = 0.1,
      fixed_s = 1.4,
      grid_n = 1000
    )
  })

  system.time({
    out_pi <- run_grid_experiment(
      grid_values = pi_grid,
      grid_name = "pi",
      m = m,
      alpha_mdr = alpha,
      q_adadetect = q_adadetect,
      reps = reps,
      n0 = n0,
      n1 = n1,
      K0 = K0,
      ncores = ncores,
      fixed_pi = 0.1,
      fixed_s = 1.4,
      grid_n = 1000
    )
  })

  saveRDS(out_s, file = rds_s)
  saveRDS(out_pi, file = rds_pi)
}

cat("s experiment replications:", out_s$reps, "\n")
cat("pi experiment replications:", out_pi$reps, "\n")
cat("Offline array dim for s:", dim(out_s$off), "\n")
cat("AdaDetect array dim for s:", dim(out_s$ada), "\n")
cat("Offline array dim for pi:", dim(out_pi$off), "\n")
cat("AdaDetect array dim for pi:", dim(out_pi$ada), "\n")

############################################################
## 10. Summarise results
############################################################

summarise_array <- function(arr, method, metric, x_grid, experiment) {
  if (identical(metric, "MDR")) {
    num <- arr[, , "Missed"]
    den <- arr[, , "Signals"]
    mean_num <- colMeans(num, na.rm = TRUE)
    mean_den <- pmax(colMeans(den, na.rm = TRUE), .Machine$double.eps)
    mean_val <- mean_num / mean_den
    influence <- num - sweep(den, 2L, mean_val, "*")
    se_val <- apply(influence, 2, sd, na.rm = TRUE) /
      (sqrt(nrow(arr)) * mean_den)
  } else {
    mat <- arr[, , metric]
    mean_val <- colMeans(mat, na.rm = TRUE)
    se_val <- apply(mat, 2, sd, na.rm = TRUE) / sqrt(nrow(mat))
  }
  
  data.frame(
    x = x_grid,
    method = method,
    metric = metric,
    experiment = experiment,
    mean = mean_val,
    se = se_val,
    lower = mean_val - 1.96 * se_val,
    upper = mean_val + 1.96 * se_val
  )
}

df_s <- bind_rows(
  summarise_array(out_s$off, "Offline MDRC", "MDR", s_grid, "vary_s"),
  summarise_array(out_s$ada, "AdaDetect", "MDR", s_grid, "vary_s"),
  summarise_array(out_s$off, "Offline MDRC", "FDR", s_grid, "vary_s"),
  summarise_array(out_s$ada, "AdaDetect", "FDR", s_grid, "vary_s"),
  summarise_array(out_s$off, "Offline MDRC", "TDR", s_grid, "vary_s"),
  summarise_array(out_s$ada, "AdaDetect", "TDR", s_grid, "vary_s"),
  summarise_array(out_s$off, "Offline MDRC", "False_Pos", s_grid, "vary_s"),
  summarise_array(out_s$ada, "AdaDetect", "False_Pos", s_grid, "vary_s")
)

df_pi <- bind_rows(
  summarise_array(out_pi$off, "Offline MDRC", "MDR", pi_grid, "vary_pi"),
  summarise_array(out_pi$ada, "AdaDetect", "MDR", pi_grid, "vary_pi"),
  summarise_array(out_pi$off, "Offline MDRC", "FDR", pi_grid, "vary_pi"),
  summarise_array(out_pi$ada, "AdaDetect", "FDR", pi_grid, "vary_pi"),
  summarise_array(out_pi$off, "Offline MDRC", "TDR", pi_grid, "vary_pi"),
  summarise_array(out_pi$ada, "AdaDetect", "TDR", pi_grid, "vary_pi"),
  summarise_array(out_pi$off, "Offline MDRC", "False_Pos", pi_grid, "vary_pi"),
  summarise_array(out_pi$ada, "AdaDetect", "False_Pos", pi_grid, "vary_pi")
)

df_plot <- bind_rows(df_s, df_pi)
df_plot$method <- factor(df_plot$method, levels = c("Offline MDRC", "AdaDetect"))

write.csv(
  df_plot,
  file = file.path(output_dir, "OfflineMDRC_vs_AdaDetect_2x2_summary.csv"),
  row.names = FALSE
)

############################################################
## 11. Plot 2 x 2 figure
## Layout:
##   left column  = FDR
##   right column = MDR
##
## Confidence bands:
##   shaded region = Monte Carlo mean +/- 1.96 standard errors
############################################################

my_colors <- c(
  "Offline MDRC" = "#00BFC4",
  "AdaDetect" = "#F8766D"
)

my_shapes <- c(
  "Offline MDRC" = 15,
  "AdaDetect" = 16
)

## ggplot2 compatibility: versions < 3.4.0 (e.g. 3.3.5 on this server) do not
## support the `linewidth` argument and require `size` for line thickness.
.gg_new <- utils::packageVersion("ggplot2") >= "3.4.0"
line_lw <- function(x) if (.gg_new) list(linewidth = x) else list(size = x)

custom_theme <- theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 17),
    plot.subtitle = element_text(size = 17, hjust = 0.5, margin = margin(b = 5)),
    legend.title = element_blank(),
    legend.text = element_text(size = 15),
    axis.text = element_text(size = 14, colour = "black"),
    axis.title = element_text(size = 15),
    panel.grid.major = do.call(element_line,
                               c(list(colour = "grey92"), line_lw(0.4))),
    panel.grid.minor = element_blank(),
    panel.border = do.call(element_rect,
                           c(list(colour = "black", fill = NA), line_lw(0.8))),
    plot.margin = margin(8, 8, 8, 8)
  )

plot_panel <- function(df,
                       experiment_name,
                       metric_name,
                       xlab,
                       ylab,
                       subtitle,
                       hline = NULL,
                       ylim = c(0, 1),
                       y_breaks = NULL) {
  
  p <- ggplot(
    df %>% filter(experiment == experiment_name, metric == metric_name),
    aes(x = x, y = mean, color = method, shape = method, fill = method)
  ) +
    geom_ribbon(
      aes(ymin = lower, ymax = upper),
      alpha = 0.30,
      color = NA,
      show.legend = FALSE
    ) +
    do.call(geom_line, line_lw(0.85)) +
    geom_point(size = 2.1) +
    scale_color_manual(values = my_colors) +
    scale_fill_manual(values = my_colors) +
    scale_shape_manual(values = my_shapes) +
    labs(
      subtitle = subtitle,
      x = xlab,
      y = ylab
    ) +
    custom_theme +
    coord_cartesian(ylim = ylim)
  
  if (!is.null(y_breaks)) {
    p <- p + scale_y_continuous(breaks = y_breaks)
  }
  
  if (!is.null(hline)) {
    p <- p +
      do.call(geom_hline, c(
        list(
          yintercept = hline,
          linetype = "dashed",
          color = "black"
        ),
        line_lw(0.8)
      ))
  }
  
  p
}

p1 <- plot_panel(
  df = df_plot,
  experiment_name = "vary_s",
  metric_name = "FDR",
  xlab = "Signal strength scale (s)",
  ylab = "FDR",
  subtitle = "(a.1)",
  hline = q_adadetect,
  ylim = c(0, 1),
  y_breaks = seq(0, 1, by = 0.25)
)

p2 <- plot_panel(
  df = df_plot,
  experiment_name = "vary_s",
  metric_name = "MDR",
  xlab = "Signal strength scale (s)",
  ylab = "MDR",
  subtitle = "(a.2)",
  hline = alpha,
  ylim = c(0, 1),
  y_breaks = seq(0, 1, by = 0.25)
)

p3 <- plot_panel(
  df = df_plot,
  experiment_name = "vary_pi",
  metric_name = "FDR",
  xlab = expression("Signal proportion " * pi),
  ylab = "FDR",
  subtitle = "(b.1)",
  hline = q_adadetect,
  ylim = c(0, 0.75),
  y_breaks = seq(0, 0.75, by = 0.25)
)

p4 <- plot_panel(
  df = df_plot,
  experiment_name = "vary_pi",
  metric_name = "MDR",
  xlab = expression("Signal proportion " * pi),
  ylab = "MDR",
  subtitle = "(b.2)",
  hline = alpha,
  ylim = c(0, 1),
  y_breaks = seq(0, 1, by = 0.25)
)

final_plot <- (p1 + p2) / (p3 + p4) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

if (interactive()) print(final_plot)

ggsave(
  filename = file.path(output_dir, "OfflineMDRC_vs_AdaDetect_2x2_FDR_left_MDR_right.pdf"),
  plot = final_plot,
  width = 10.8,
  height = 8
)

ggsave(
  filename = file.path(output_dir, "OfflineMDRC_vs_AdaDetect_2x2_FDR_left_MDR_right.png"),
  plot = final_plot,
  width = 10.8,
  height = 8,
  dpi = 300
)
