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
#   DGP: Setting (a) of the main text, with one-dimensional skewed components:
#     F0 = Exp(1), F1 = Gamma(shape_alt = 3, scale_alt = 2)
#     pi_t : smooth prior-drop 0.35 -> 0.03 (drop_center 0.50, drop_width 0.10)
#   The exact density ratio
#   r(x) = dgamma(x, shape_alt, scale_alt) / dexp(x, 1) gives the oracle Lmdr.
#   The signal proportion pi_t follows a smooth prior-DROP
#   (0.35 -> 0.03, drop_center 0.50, drop_width 0.10) over warm-up + stream;
#   the stream segment is anchored at k0_ref = 500 so it is IDENTICAL across
#   every K0 value. The monotone drop makes the trailing window lag pi_t
#   downward, so the beta- and K0-sweeps expose the bias-variance trade-off
#   of the window construction, while the n-sweep shows the graded
#   convergence of DD to OR as the labeled reference grows.
#
#   OUTPUT: a single-column 4 x 2 figure (rows: n-sweep / beta-sweep /
#   K0-sweep / w-sweep; columns: MDR / FDR) with pointwise 95% Monte-Carlo bands.
#
# ==============================================================================
# MANUSCRIPT BARRIER VERSION
# ==============================================================================
# Both OMDRC arms use the credited local barrier from Algorithms 1 and 2:
#   OMDRC.OR -> OMDRC.OR+NFR = OMDRC_OR_NFR(true q,      alpha, d, true q_ini)
#   OMDRC.DD -> OMDRC.DD+NFR = OMDRC_OR_NFR(estimated q, alpha, d, q_hat_ini)
# The barrier gate is `q_t < gamma_t` with gamma_t the causal windowed FT
# threshold over a trailing window of length d (_nfr_core.R::nfr_gamma_path);
# "credited" means a barrier-forced rejection still credits +alpha*q.
#
# A FOURTH sweep row is added for the barrier's own window length d, which is
# the ONLY new free parameter the barrier introduces. Note the semantic
# difference that must be stated in the caption: n / beta / K0 are ESTIMATION
# parameters (Stage 1 accuracy, pi_t window bias-variance, warm-up), whereas d
# is a POLICY parameter of the constraint.
#
# WHY THE d ROW IS ESSENTIALLY FREE: d enters neither Stage 1 nor Stage 2 -- the
# estimated score path (GAM ratio + windowed EM pi_t) depends only on the
# trailing window of density RATIOS, never on d and never on the decisions. So
# OMDRC_DD is fitted ONCE per replicate and the barrier is replayed for every d
# on the same score path (same trick as beta_error_diag.R). The whole d sweep
# therefore costs about as much as a single configuration.
#
# The oracle baseline drawn in ALL FOUR rows is OMDRC.OR+NFR at d = default_d
# (100), so the DD-vs-OR gap is purely estimation error with the barrier present
# in both arms. The OR arm's own d-dependence is printed in the log instead of
# being drawn (it is pathwise controlled for every d, see the ledger table).
#
# Outputs: Rplot_NFR_D_n_sensitivity.{pdf,png} + nfr_D_n_vary_results.rds.
# The non-barrier originals (D_n_vary.R, Rplot_D_n_sensitivity.*) are untouched.
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

# Shared NFR core: nfr_gamma_path() + the credited barrier OMDRC_OR_NFR().
# DO NOT EDIT _nfr_core.R -- its sha256 is recorded in several application
# validation-lock files.
nfr_core_file <- file.path(code_dir, "nofreeride", "_nfr_core.R")
if (!file.exists(nfr_core_file)) stop("Cannot find _nfr_core.R at: ", nfr_core_file)
source(nfr_core_file)

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
                           k0_ref = NULL,
                           d_vec = 100L) {

  # d_vec may hold SEVERAL barrier windows; the expensive parts (GAM ratio, EM
  # pi_t) are computed once per replicate and the barrier is replayed per d.
  d_vec <- as.integer(d_vec)
  nd <- length(d_vec)
  run_or <- "OR" %in% methods_to_run
  run_dd <- "DD" %in% methods_to_run
  mk <- function() lapply(seq_len(nd), function(.) matrix(NA_real_, reps, length(m)))
  mdr.or <- if (run_or) mk() else NULL
  fdr.or <- if (run_or) mk() else NULL
  mdr.dd <- if (run_dd) mk() else NULL
  fdr.dd <- if (run_dd) mk() else NULL
  mdrq.or <- if (run_or) matrix(NA_real_, reps, nd) else NULL
  mdrq.dd <- if (run_dd) matrix(NA_real_, reps, nd) else NULL

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
                      "pi_bounds", "code_dir", "nfr_core_file",
                      "D_mode", "D_beta", "D_min", "d_vec", "nd"),
                envir = environment())

  result <- foreach(r = 1:reps, .packages = c("mgcv")) %dopar% {
    source(file.path(code_dir, "OMDRC.R"))
    source(nfr_core_file)
    set.seed(r)
    mx <- max(m)

    # ---- Initial unlabeled warm-up window of size `ini` ----
    # RNG draw ORDER is identical to D_n_vary.R, so the streams are the same.
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
    # True Lmdr of the warm-up batch: fills the barrier's gamma window while
    # t < d, mirroring x.Lmdr_ini in the nofreeride scripts.
    r_ini <- dgamma(z_ini, shape = shape_alt, scale = scale_alt) /
      dexp(z_ini, rate = 1)
    lmdr_ini <- (pi_ini * r_ini) / ((1 - pi_ini) + pi_ini * r_ini)

    # Cumulative MDR / FDR read out at the time grid m, plus the per-path
    # ledger quantity MDR_q = sum_skip q / sum q in the TRUE q.
    metrics <- function(de) {
      list(mdr = sapply(m, function(k)
             sum(theta_stream[1:k] * (1 - de[1:k])) /
               max(sum(pi_stream[1:k]), .Machine$double.eps)),
           fdr = sapply(m, function(k)
             sum((1 - theta_stream[1:k]) * de[1:k]) / max(sum(de[1:k]), 1)),
           mdrq = sum(lmdr_stream * (1 - de)) / sum(lmdr_stream))
    }

    res_list <- list()

    # ---- Oracle OMDRC + credited barrier, replayed for every d ----
    if ("OR" %in% methods_to_run) {
      res_list$mdr.or <- matrix(NA_real_, nd, length(m))
      res_list$fdr.or <- matrix(NA_real_, nd, length(m))
      res_list$mdrq.or <- rep(NA_real_, nd)
      for (j in seq_len(nd)) {
        de <- OMDRC_OR_NFR(lmdr_stream, alpha, d = d_vec[j],
                           x.Lmdr_ini = lmdr_ini)$de
        mm <- metrics(de)
        res_list$mdr.or[j, ] <- mm$mdr
        res_list$fdr.or[j, ] <- mm$fdr
        res_list$mdrq.or[j] <- mm$mdrq
      }
    }

    # ---- OMDRC.DD (PC-DRE) + credited barrier, replayed for every d ----
    if ("DD" %in% methods_to_run) {
      res_list$mdr.dd <- matrix(NA_real_, nd, length(m))
      res_list$fdr.dd <- matrix(NA_real_, nd, length(m))
      res_list$mdrq.dd <- rep(NA_real_, nd)
      if (n >= 3) {
        # Labeled null / alternative reference samples of size n each (n0 = n1 = n).
        z0_ref <- rexp(n, rate = 1)
        z1_ref <- rgamma(n, shape = shape_alt, scale = scale_alt)
        # Guard the spline dimension against tiny reference samples; for
        # n >= 6 this equals the PC-DRE default spline_k = 10.
        spline_k <- max(3L, min(10L, as.integer(2 * n - 2)))

        # ONE fit per replicate: the estimated score path is alpha- and
        # d-independent, so all barrier windows reuse DR / DR_ini.
        fit <- OMDRC_DD(z = z_stream, z_ini = z_ini,
                        z0 = z0_ref, z1 = z1_ref,
                        alpha = alpha, D = D,
                        ratio_method = "gam",
                        ratio_control = list(spline_k = spline_k),
                        pi_bounds = pi_bounds,
                        em_tol = 1e-6, em_max_iter = 50L,
                        D_mode = D_mode, D_beta = D_beta, D_min = D_min)
        for (j in seq_len(nd)) {
          de <- OMDRC_OR_NFR(fit$DR, alpha, d = d_vec[j],
                             x.Lmdr_ini = fit$DR_ini)$de
          mm <- metrics(de)
          res_list$mdr.dd[j, ] <- mm$mdr
          res_list$fdr.dd[j, ] <- mm$fdr
          res_list$mdrq.dd[j] <- mm$mdrq
        }
      }
    }

    res_list
  }

  for (r in 1:reps) {
    for (j in seq_len(nd)) {
      if (run_or) {
        mdr.or[[j]][r, ] <- result[[r]]$mdr.or[j, ]
        fdr.or[[j]][r, ] <- result[[r]]$fdr.or[j, ]
        mdrq.or[r, j] <- result[[r]]$mdrq.or[j]
      }
      if (run_dd) {
        mdr.dd[[j]][r, ] <- result[[r]]$mdr.dd[j, ]
        fdr.dd[[j]][r, ] <- result[[r]]$fdr.dd[j, ]
        mdrq.dd[r, j] <- result[[r]]$mdrq.dd[j]
      }
    }
  }

  # One element per d, each shaped exactly like the non-barrier script's output
  # so the downstream plotting helpers are unchanged.
  lapply(seq_len(nd), function(j) list(
    d = d_vec[j],
    mdr.dd = if (run_dd) mdr.dd[[j]] else NULL,
    fdr.dd = if (run_dd) fdr.dd[[j]] else NULL,
    mdr.or = if (run_or) mdr.or[[j]] else NULL,
    fdr.or = if (run_or) fdr.or[[j]] else NULL,
    mdrq.or = if (run_or) mdrq.or[, j] else NULL,
    mdrq.dd = if (run_dd) mdrq.dd[, j] else NULL))
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
shape_alt_v <- 3      # F1 = Gamma(shape_alt, scale_alt) vs F0 = Exp(1).
# scale_alt = 2 is the manuscript value used by nfr_setting1_fig.R. It may be
# overridden from the shell (NFR_SCALE_ALT) for a diagnostic run; any value != 2
# writes to its own rds / figure via out_tag below, so
# the aligned s = 2 outputs are never overwritten.
scale_alt_v <- as.numeric(Sys.getenv("NFR_SCALE_ALT", "2"))
out_tag <- if (isTRUE(all.equal(scale_alt_v, 2))) "" else
  paste0("_scale", sub("\\.", "p", format(scale_alt_v, trim = TRUE)))
cat(sprintf("DGP: F1 = Gamma(shape=%g, scale=%g); output tag = '%s'\n",
            shape_alt_v, scale_alt_v, out_tag))

# --- Prevalence path: main-text Setting-1 prior-drop (0.35 -> 0.03) ---
# This sensitivity study shares the main-text Setting (a) prevalence path.
pi_mode_v     <- "vary"
pi_low_v      <- 0.03
pi_high_v     <- 0.35
drop_center_v <- 0.50
drop_width_v  <- 0.10
pi_bnds <- c(0.01, 0.99)   # paper default: broad bounds [pi_lower, pi_upper]

# --- Parameter grids for the four sensitivity sweeps ---
n_list    <- c(10, 20, 50, 100, 200, 500, 1000)      # labeled reference size (n0 = n1 = n)
beta_list <- c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)    # window growth exponent
k0_list   <- c(10, 20, 50, 100, 200, 500, 1000)      # initial unlabeled batch size K0
                                                   # (same grid as n_list)
# Local-barrier window w: the grid spans two orders of magnitude and includes
# the manuscript default w = 100. With T = 1000 and K0 = 500 the value
# w = 1000 makes the gamma window cover essentially the whole available history,
# i.e. it degenerates to a GLOBAL fixed-threshold gate -- a deliberate extreme
# reference point rather than a practical setting.
d_list <- c(25, 50, 100, 200, 500, 1000)

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
default_d      <- 100L   # NFR barrier window held fixed in the n/beta/K0 rows,
                         # and the d value of the OR+NFR baseline in every row

# ------------------------------------------------------------------------------
# 3. Run the experiments (or reuse cached results for a fast re-plot)
# ------------------------------------------------------------------------------
results_file <- file.path(code_dir, "sensity",
                          paste0("nfr_D_n_vary_results", out_tag, ".rds"))
reuse_cache  <- TRUE   # set FALSE to rerun everything from scratch

sim_common <- function(...) {
  run_sim_retry(...,
                shape_alt = shape_alt_v, scale_alt = scale_alt_v,
                pi_mode = pi_mode_v,
                pi_low = pi_low_v, pi_high = pi_high_v,
                drop_center = drop_center_v, drop_width = drop_width_v,
                pi_bounds = pi_bnds)
}
# The engine always returns a LIST over d; the n/beta/K0 sweeps use a single
# barrier window, so unwrap the only element to keep their downstream shape.
sim_one_d <- function(...) sim_common(..., d_vec = default_d)[[1]]

run_or_baseline <- function() {
  cat("\n====== Running OMDRC.OR+NFR baseline (once, all d) ======\n")
  # The oracle arm is cheap, so it is evaluated on the FULL d grid in one pass:
  # the row baselines use d = default_d and the rest goes into the log table.
  sim_common(m = m_seq, ini = default_ini, n = default_n,
             alpha = alpha_v, reps = reps_v,
             D = default_D, methods_to_run = "OR",
             d_vec = sort(unique(c(default_d, d_list))))
}
run_n_sweep <- function() {
  lapply(n_list, function(nv) {
    cat(paste0("\n====== OMDRC.DD+NFR (PC-DRE) with n = ", nv, " ======\n"))
    sim_one_d(m = m_seq, ini = default_ini, n = nv,
              alpha = alpha_v, reps = reps_v, D = default_D,
              methods_to_run = "DD",
              D_mode = default_D_mode, D_beta = default_D_beta,
              D_min = default_D_min)
  })
}
run_beta_sweep <- function() {
  lapply(beta_list, function(bv) {
    cat(paste0("\n====== OMDRC.DD+NFR (PC-DRE) with beta = ", bv,
               " (D_t = floor((K0+t)^", bv, ")) ======\n"))
    sim_one_d(m = m_seq, ini = default_ini, n = default_n,
              alpha = alpha_v, reps = reps_v, D = default_D,
              methods_to_run = "DD",
              D_mode = "growing", D_beta = bv,
              D_min = default_D_min)
  })
}
run_k0_sweep <- function() {
  lapply(k0_list, function(kv) {
    cat(paste0("\n====== OMDRC.DD+NFR (PC-DRE) with K0 = ", kv,
               " (stream pi anchored at k0_ref = ", default_ini, ") ======\n"))
    sim_one_d(m = m_seq, ini = kv, n = default_n,
              alpha = alpha_v, reps = reps_v, D = default_D,
              methods_to_run = "DD",
              D_mode = "growing", D_beta = default_D_beta,
              D_min = default_D_min, k0_ref = default_ini)
  })
}
run_d_sweep <- function() {
  cat(paste0("\n====== OMDRC.DD+NFR barrier-window sweep, d in {",
             paste(d_list, collapse = ", "),
             "} (ONE GAM fit per replicate) ======\n"))
  # All defaults for n / beta / K0; only the barrier window varies. Because the
  # score path is d-independent this single call covers the whole row.
  sim_common(m = m_seq, ini = default_ini, n = default_n,
             alpha = alpha_v, reps = reps_v, D = default_D,
             methods_to_run = "DD",
             D_mode = "growing", D_beta = default_D_beta,
             D_min = default_D_min, d_vec = d_list)
}

# Load the cache and reuse whichever components are present AND match the
# current grids; recompute only the missing/mismatched components, then merge
# everything back into the cache file. (Caches saved by the previous versions
# already contain res_n / res_beta with identical grids and seeds, so only the
# new K0 sweep is computed on the first run after the D_min -> K0 switch.)
cached <- if (reuse_cache && file.exists(results_file)) readRDS(results_file) else NULL

# Cache guard: the per-component reuse below only checks that the *grids* match,
# NOT the DGP. A change to the alternative (shape_alt) or the prevalence path
# invalidates every stored matrix, so drop the whole cache when the DGP differs
# from what produced it.
if (!is.null(cached)) {
  dgp_same <- isTRUE(all.equal(cached$shape_alt, shape_alt_v)) &&
    isTRUE(all.equal(cached$scale_alt, scale_alt_v)) &&
    isTRUE(all.equal(cached$pi_low, pi_low_v)) &&
    isTRUE(all.equal(cached$pi_high, pi_high_v)) &&
    isTRUE(all.equal(cached$drop_center, drop_center_v)) &&
    isTRUE(all.equal(cached$drop_width, drop_width_v)) &&
    isTRUE(all.equal(cached$reps, reps_v))
  if (!dgp_same) {
    cat("Cache DGP/reps mismatch -> recomputing everything from scratch.\n")
    cached <- NULL
  }
}

or_d_grid <- sort(unique(c(default_d, d_list)))
if (!is.null(cached$or_all) && isTRUE(all.equal(cached$or_d_grid, or_d_grid))) {
  cat("Reusing cached OMDRC.OR+NFR baseline.\n")
  or_all <- cached$or_all
} else {
  or_all <- run_or_baseline()
}
i_or_default <- which(or_d_grid == default_d)
mdr_or_mat <- or_all[[i_or_default]]$mdr.or
fdr_or_mat <- or_all[[i_or_default]]$fdr.or
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

if (!is.null(cached$res_d) && isTRUE(all.equal(cached$d_list, d_list))) {
  cat("Reusing cached d-sweep.\n")
  res_d <- cached$res_d
} else {
  res_d <- run_d_sweep()
}

# Cache raw results so the figure can be re-plotted without rerunning.
saveRDS(list(m = m_seq, alpha = alpha_v,
             shape_alt = shape_alt_v, scale_alt = scale_alt_v, reps = reps_v,
             pi_mode = pi_mode_v,
             pi_low = pi_low_v, pi_high = pi_high_v,
             drop_center = drop_center_v, drop_width = drop_width_v,
             pi_bounds = pi_bnds,
             n_list = n_list, beta_list = beta_list, k0_list = k0_list,
             d_list = d_list, or_d_grid = or_d_grid, default_d = default_d,
             barrier_ledger = "credited",
             default_n = default_n, default_D_beta = default_D_beta,
             default_D_min = default_D_min, default_ini = default_ini,
             mdr_or = mdr_or_vec, fdr_or = fdr_or_vec,
             mdr.or.mat = mdr_or_mat, fdr.or.mat = fdr_or_mat,
             or_all = or_all,
             res_n = res_n, res_beta = res_beta, res_k0 = res_k0,
             res_d = res_d),
        file = results_file)

# Terminal MDR/FDR summary (OR+NFR vs growing-window DD+NFR at the defaults).
idxT <- length(m_seq)
i_n_default    <- which(n_list == default_n)
i_beta_default <- which(beta_list == default_D_beta)
i_k0_default   <- which(k0_list == default_ini)
i_d_default    <- which(d_list == default_d)
cat(sprintf("\n--- Terminal MDR/FDR at eval-t=%d (alpha=%.2f, credited barrier) ---\n",
            m_seq[idxT], alpha_v))
cat(sprintf("%-34s %8s %8s\n", "Method", "MDR", "FDR"))
cat(sprintf("%-34s %8.4f %8.4f\n",
            sprintf("OMDRC.OR+NFR (d=%d)", default_d),
            mdr_or_vec[idxT], fdr_or_vec[idxT]))
cat(sprintf("%-34s %8.4f %8.4f\n",
            sprintf("DD+NFR (n=%d, beta=%.1f, K0=%d, d=%d)",
                    default_n, default_D_beta, default_ini, default_d),
            mean(res_k0[[i_k0_default]]$mdr.dd[, idxT], na.rm = TRUE),
            mean(res_k0[[i_k0_default]]$fdr.dd[, idxT], na.rm = TRUE)))
cat("\n--- Terminal MDR by sweep value ---\n")
cat("n   :", paste(sprintf("%g=%.3f", n_list,
    sapply(res_n, function(z) mean(z$mdr.dd[, idxT], na.rm = TRUE))), collapse = "  "), "\n")
cat("beta:", paste(sprintf("%g=%.3f", beta_list,
    sapply(res_beta, function(z) mean(z$mdr.dd[, idxT], na.rm = TRUE))), collapse = "  "), "\n")
cat("K0  :", paste(sprintf("%g=%.3f", k0_list,
    sapply(res_k0, function(z) mean(z$mdr.dd[, idxT], na.rm = TRUE))), collapse = "  "), "\n")
cat("d   :", paste(sprintf("%g=%.3f", d_list,
    sapply(res_d, function(z) mean(z$mdr.dd[, idxT], na.rm = TRUE))), collapse = "  "), "\n")

# Per-path ledger check in the TRUE q: MDR_q = sum_skip q / sum q must be
# <= alpha on EVERY path for the ORACLE arm (capacity identity). The DD arm runs
# its ledger on q_hat, so it can and does fail -- that is the price of
# estimation error, and this table is the reps-robust evidence (the mean MDR
# curves at reps=100 cannot settle it).
cat(sprintf("\n--- Barrier ledger check, MDR_q in true q (alpha=%.2f) ---\n", alpha_v))
cat(sprintf("%-22s %10s %10s %8s\n", "Arm / config", "mean", "max", "ledger"))
ledger_row <- function(lab, v) {
  cat(sprintf("%-22s %10.4f %10.4f %8s\n", lab, mean(v, na.rm = TRUE),
              max(v, na.rm = TRUE),
              if (max(v, na.rm = TRUE) <= alpha_v + 1e-9) "PASS" else "FAIL"))
}
for (j in seq_along(or_d_grid))
  ledger_row(sprintf("OR+NFR  d=%g", or_d_grid[j]), or_all[[j]]$mdrq.or)
for (j in seq_along(d_list))
  ledger_row(sprintf("DD+NFR  d=%g", d_list[j]), res_d[[j]]$mdrq.dd)
cat(sprintf("%-22s %10s %10s %8s\n", "-- n sweep (d=100)", "", "", ""))
for (j in seq_along(n_list))
  ledger_row(sprintf("DD+NFR  n=%g", n_list[j]), res_n[[j]]$mdrq.dd)

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
# 5. Visualization (4 x 2 grid: rows n / beta / K0 / d, columns MDR / FDR)
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

# --- Row D: Sensitivity to the NFR barrier window d (Oranges) ---
# NOTE: unlike n / beta / K0 this is NOT an estimation parameter but the POLICY
# parameter of the anti-free-riding constraint; state that in the caption.
df_d <- get_plot_df(res_d, "d", d_list, mdr_or_mat, fdr_or_mat)
colors_d <- c("#F8766D", brewer.pal(9, "Oranges")[3:9][seq_along(d_list)])
labels_d <- parse(text = c("OMDRC.OR",
                           paste0("OMDRC.DD~(w==", d_list, ")")))
row_d <- make_row(df_d, colors_d, "d", labels_d)

# ------------------------------------------------------------------------------
# 6. Final Integration and Layout (4 x 2 grid)
# ------------------------------------------------------------------------------
final_plot <- ggarrange(row_n, row_beta, row_k0, row_d, ncol = 1, nrow = 4)

# Portrait-friendly aspect (width 10); the height is scaled from the 3-row
# original (10 x 16) to keep the per-row aspect identical with four rows.
if (interactive()) print(final_plot)
fig_base <- file.path(code_dir, "sensity",
                      paste0("Rplot_NFR_D_n_sensitivity", out_tag))
ggsave(paste0(fig_base, ".pdf"), final_plot, width = 10, height = 21)
ggsave(paste0(fig_base, ".png"), final_plot, width = 10, height = 21, dpi = 300)
cat("\nFigure written to", paste0(fig_base, ".{pdf,png}"), "\n")
