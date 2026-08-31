# =========================================================================
# SCRIPT FOR REPRODUCING THE CREDIT CARD FRAUD DETECTION EXPERIMENT
# (REVISED MANUSCRIPT VERSION, CAUSAL PI_BOUNDS)
# OMDRC.DD uses a 28D additive-GAM classifier for the density ratio and carries
# the credited local barrier in Algorithm 2, i.e. the reported OMDRC.DD rule.
# The manuscript output is the six-panel single-column Figure 3 written at the
# end of this script. OMDRC.DD prior bounds are set by PI_BOUNDS_MODE below:
# "causal"
# anchors them on a warm-up-only EM fit, "tuned" uses the post hoc interval.
# =========================================================================
project_root <- normalizePath(Sys.getenv("OMDRC_ROOT", unset = getwd()),
                              mustWork = TRUE)
setwd(file.path(project_root, "application/CCFD"))

# --- 1. SETUP: Install and Load Packages ---
required_packages <- c("data.table", "mgcv", "ggplot2", "dplyr", "patchwork", "onlineFDR", "scales")

install_if_missing <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}
cat("Loading required packages...\n")
sapply(required_packages, install_if_missing)

# OMDRC.R must sit in the working directory or its parent.
OMDRC_PATH <- if (file.exists("OMDRC.R")) "OMDRC.R" else "../OMDRC.R"
if (file.exists(OMDRC_PATH)) {
  cat("Loading OMDRC functions from", OMDRC_PATH, "...\n")
  source(OMDRC_PATH)
} else {
  stop("OMDRC.R not found in the working directory or its parent. Please place it there.")
}


# --- 2. DATA / RESULT-CACHE CONFIGURATION ---
features_to_use <- paste0("V", 1:28)
RESULT_CACHE <- "ccfd_results_cache.rds"
USE_RESULT_CACHE <- TRUE
# The release does not ship caches: regenerate automatically on the first run.
FORCE_RECOMPUTE <- !file.exists(RESULT_CACHE)
USE_ROBUST_PREPROCESSING <- TRUE
RUN_ADJ_SAFFRON <- TRUE
SAFFRON_CACHE <- "ccfd_adj_saffron_decisions.rds"
ALLOW_FAST_SAFFRON_FALLBACK <- FALSE
RESULT_CACHE_CONFIG <- list(
  features_to_use = features_to_use,
  USE_ROBUST_PREPROCESSING = USE_ROBUST_PREPROCESSING,
  RUN_ADJ_SAFFRON = RUN_ADJ_SAFFRON,
  ALLOW_FAST_SAFFRON_FALLBACK = ALLOW_FAST_SAFFRON_FALLBACK,
  SAFFRON_PVAL_REF = "z1_signal_cdf_swapped",
  ALPHA = 0.1,
  K0 = 2000,
  N_LABELED = 200,
  N0_LABELED = 3000,
  D_BETA = 0.6,
  D_MIN = 10L,
  W_BARRIER = 100L,
  GAM_SPLINE_K = 10L,
  PI_BOUNDS = c(0.0005, 0.0018)
)

cat(sprintf("Using %d PCA features with an additive GAM density-ratio estimator%s.\n",
            length(features_to_use),
            if (USE_ROBUST_PREPROCESSING) " and robust reference/calibration preprocessing" else ""))

fit_robust_preprocessor <- function(x0, x1, x_ini, winsor_probs = c(0.005, 0.995)) {
  anchor <- rbind(x0, x1, x_ini)
  lower <- apply(anchor, 2, quantile, probs = winsor_probs[1], na.rm = TRUE, names = FALSE)
  upper <- apply(anchor, 2, quantile, probs = winsor_probs[2], na.rm = TRUE, names = FALSE)
  anchor <- sweep(anchor, 2, lower, pmax)
  anchor <- sweep(anchor, 2, upper, pmin)

  center <- apply(anchor, 2, median, na.rm = TRUE)
  scale_value <- apply(anchor, 2, IQR, na.rm = TRUE)
  scale_value[!is.finite(scale_value) | scale_value <= 1e-12] <- 1

  function(x) {
    ans <- sweep(x, 2, lower, pmax)
    ans <- sweep(ans, 2, upper, pmin)
    ans <- sweep(ans, 2, center, "-")
    sweep(ans, 2, scale_value, "/")
  }
}


# --- 3. EXPERIMENT SETUP / COMPUTATION CACHE ---
if (USE_RESULT_CACHE && file.exists(RESULT_CACHE) && !FORCE_RECOMPUTE) {
  cat("Loading cached CCFD computation results from ", RESULT_CACHE, "...\n", sep = "")
  cache <- readRDS(RESULT_CACHE)
  if (is.null(cache$config) || !identical(cache$config, RESULT_CACHE_CONFIG)) {
    stop(sprintf(
      "Result cache '%s' was created with different data/model/parameter/decision settings. Set FORCE_RECOMPUTE <- TRUE only when you intentionally changed those settings.",
      RESULT_CACHE
    ))
  }
  list2env(cache[names(cache) != "config"], envir = environment())
  if (RUN_ADJ_SAFFRON && length(saffron_decisions) != n_stream) {
    stop("Cached Adj-SAFFRON decisions have the wrong length.")
  }
  cat(sprintf("Cached discoveries -- OMDRC.DD: %d | SCT: %d | RTK: %d%s\n",
              sum(omdrc_decisions), sum(sct_decisions), sum(rtk_decisions),
              if (RUN_ADJ_SAFFRON) sprintf(" | Adj-SAFFRON: %d", sum(saffron_decisions)) else ""))
} else {
  if (USE_RESULT_CACHE && !FORCE_RECOMPUTE && !file.exists(RESULT_CACHE)) {
    stop(sprintf(
      "Result cache '%s' is missing. To create it once, set FORCE_RECOMPUTE <- TRUE; otherwise keep the existing cache and only redraw figures.",
      RESULT_CACHE
    ))
  }

  set.seed(123)
  ALPHA      <- RESULT_CACHE_CONFIG$ALPHA
  K0         <- RESULT_CACHE_CONFIG$K0
  N_LABELED  <- RESULT_CACHE_CONFIG$N_LABELED
  N0_LABELED <- RESULT_CACHE_CONFIG$N0_LABELED
  D_BETA     <- RESULT_CACHE_CONFIG$D_BETA
  D_MIN      <- RESULT_CACHE_CONFIG$D_MIN
  W_BARRIER  <- RESULT_CACHE_CONFIG$W_BARRIER
  GAM_SPLINE_K <- RESULT_CACHE_CONFIG$GAM_SPLINE_K
  PI_BOUNDS  <- RESULT_CACHE_CONFIG$PI_BOUNDS

  cat("Loading creditcard.csv dataset...\n")
  DATA_PATH <- if (file.exists("creditcard.csv")) "creditcard.csv" else "../ccfd_data/creditcard.csv"
  credit_data <- fread(DATA_PATH)

  normal_tx <- credit_data[Class == 0]
  fraud_tx  <- credit_data[Class == 1]

  # Labeled reference samples for the density-ratio classifier.
  z1_samples <- sample_n(fraud_tx,  N_LABELED)
  z0_samples <- sample_n(normal_tx, N0_LABELED)
  z1 <- as.matrix(as.data.frame(z1_samples)[, features_to_use])
  z0 <- as.matrix(as.data.frame(z0_samples)[, features_to_use])

  reference_samples <- bind_rows(z1_samples, z0_samples)
  online_stream_df  <- anti_join(credit_data, reference_samples,
                                 by = names(credit_data)) %>%
    arrange(Time)

  X_full_stream     <- as.matrix(as.data.frame(online_stream_df)[, features_to_use])
  theta_full_stream <- online_stream_df$Class

  z_ini    <- X_full_stream[1:K0, , drop = FALSE]
  z        <- X_full_stream[(K0 + 1):nrow(X_full_stream), , drop = FALSE]
  theta    <- theta_full_stream[(K0 + 1):length(theta_full_stream)]
  n_stream <- nrow(z)

  if (USE_ROBUST_PREPROCESSING) {
    preprocess <- fit_robust_preprocessor(z0, z1, z_ini)
    z0 <- preprocess(z0)
    z1 <- preprocess(z1)
    z_ini <- preprocess(z_ini)
    z <- preprocess(z)
  }

  cat(sprintf(
    "Experiment Setup:\n - Target MDR (alpha): %.2f\n - Feature dimension: %d\n - Initial unlabeled (K0): %d\n - Labeled frauds (n1): %d\n - Labeled normals (n0): %d\n - Online stream size: %d\n",
    ALPHA, ncol(z), K0, N_LABELED, N0_LABELED, n_stream
  ))


  # --- 4. RUN ALGORITHMS ---
  cat("\nRunning OMDRC.DD (28D additive-GAM density ratio)...\n")
  omdrc_results <- OMDRC_DD(z = z, z_ini = z_ini, z0 = z0, z1 = z1,
                            alpha = ALPHA, ratio_method = "gam",
                            ratio_control = list(spline_k = GAM_SPLINE_K),
                            pi_bounds = PI_BOUNDS,
                            em_tol = 1e-6, em_max_iter = 50L,
                            D_mode = "growing", D_beta = D_BETA,
                            D_min = D_MIN, w = W_BARRIER)
  omdrc_decisions <- omdrc_results$de

  # Score-based baselines that reuse the SAME Lmdr scores as OMDRC.DD.
  cat("Running SCT (static calibration threshold) and RTK (rolling Top-k)...\n")
  sct_res <- STATIC_LMDR_DD(omdrc_results$DR, omdrc_results$DR_ini, alpha = ALPHA)
  rtk_res <- ROLLING_TOPK_DD(omdrc_results$DR, omdrc_results$DR_ini, alpha = ALPHA,
                             window_size = length(omdrc_results$DR_ini))
  sct_decisions <- sct_res$de
  rtk_decisions <- rtk_res$de

  if (RUN_ADJ_SAFFRON) {
    saffron_cache <- if (file.exists(SAFFRON_CACHE)) readRDS(SAFFRON_CACHE) else NULL
    saffron_cache_valid <- !is.null(saffron_cache) &&
      identical(saffron_cache$reference, "F1_ratio_ecdf_1_minus_R") &&
      identical(saffron_cache$alpha, ALPHA) &&
      identical(saffron_cache$n_stream, n_stream) &&
      identical(saffron_cache$fallback, ALLOW_FAST_SAFFRON_FALLBACK)
    if (saffron_cache_valid) {
      cat("Loading cached Adj-SAFFRON decisions from ", SAFFRON_CACHE, "...\n", sep = "")
      saffron_decisions <- saffron_cache$saffron_decisions
      p_values_saffron <- saffron_cache$p_values_saffron
    } else {
      cat("Running Adj-SAFFRON under the swapped FNR formulation...\n")
      # The swapped null is F1. Its empirical CDF supplies p-values, and the
      # anomaly decision is 1-R, exactly as in the simulation section.
      ratio_signal_ref <- predict_ratio(omdrc_results$ratio_model, z1)
      r1_sorted <- sort(ratio_signal_ref)
      p_values_saffron <- findInterval(omdrc_results$LR, r1_sorted) / length(r1_sorted)
      p_values_saffron <- pmin(pmax(p_values_saffron, 1e-6), 1 - 1e-6)
      if (ALLOW_FAST_SAFFRON_FALLBACK) {
        cat("Using fast Adj-SAFFRON-style fallback because ALLOW_FAST_SAFFRON_FALLBACK = TRUE.\n")
        spending_threshold <- ALPHA / (seq_len(n_stream) * log(seq_len(n_stream) + 1)^2)
        saffron_decisions <- 1L - as.integer(p_values_saffron <= spending_threshold)
      } else {
        saffron_fdr_results <- SAFFRON(p_values_saffron, alpha = ALPHA)
        saffron_decisions <- as.integer(1L - saffron_fdr_results$R)
      }
      saveRDS(list(
        saffron_decisions = saffron_decisions,
        p_values_saffron = p_values_saffron,
        alpha = ALPHA,
        n_stream = n_stream,
        fallback = ALLOW_FAST_SAFFRON_FALLBACK,
        reference = "F1_ratio_ecdf_1_minus_R"
      ), SAFFRON_CACHE)
    }
    if (length(saffron_decisions) != n_stream) {
      stop("Cached Adj-SAFFRON decisions have the wrong length.")
    }
  } else {
    cat("Skipping Adj-SAFFRON baseline (RUN_ADJ_SAFFRON = FALSE).\n")
    p_values_saffron <- NULL
    saffron_decisions <- rep(NA_integer_, n_stream)
  }

  cat(sprintf("Discoveries -- OMDRC.DD: %d | SCT: %d | RTK: %d%s\n",
              sum(omdrc_decisions), sum(sct_decisions), sum(rtk_decisions),
              if (RUN_ADJ_SAFFRON) sprintf(" | Adj-SAFFRON: %d", sum(saffron_decisions)) else ""))

  saveRDS(list(
    config = RESULT_CACHE_CONFIG,
    ALPHA = ALPHA,
    K0 = K0,
    N_LABELED = N_LABELED,
    N0_LABELED = N0_LABELED,
    D_BETA = D_BETA,
    D_MIN = D_MIN,
    W_BARRIER = W_BARRIER,
    GAM_SPLINE_K = GAM_SPLINE_K,
    PI_BOUNDS = PI_BOUNDS,
    features_to_use = features_to_use,
    USE_ROBUST_PREPROCESSING = USE_ROBUST_PREPROCESSING,
    RUN_ADJ_SAFFRON = RUN_ADJ_SAFFRON,
    SAFFRON_CACHE = SAFFRON_CACHE,
    ALLOW_FAST_SAFFRON_FALLBACK = ALLOW_FAST_SAFFRON_FALLBACK,
    theta = theta,
    n_stream = n_stream,
    omdrc_results = omdrc_results,
    omdrc_decisions = omdrc_decisions,
    sct_decisions = sct_decisions,
    rtk_decisions = rtk_decisions,
    saffron_decisions = saffron_decisions,
    p_values_saffron = p_values_saffron
  ), RESULT_CACHE)
  cat("Saved CCFD computation results to ", RESULT_CACHE, ".\n", sep = "")
}


# =========================================================================
# --- 4a. Adj-SAFFRON p-VALUES: alternative-reference, swapped (1-R) ---
# =========================================================================
# Adj-SAFFRON is defined (Settings 1-3) by SWAPPING the null/alternative so the
# procedure targets the FNR: the null becomes "the observation is a signal",
# p-values are uniform under that null (alternative-reference ECDF), and the
# ANOMALY decision is 1 - R (a point is flagged UNLESS SAFFRON rejects the
# "is-signal" null). We mirror that here so the CCFD baseline matches the
# simulation exactly (cf. setting3.R):
#   p_t = ECDF_{r(z1)}(r(X_t)) = P_hat_{F1}( r <= r(X_t) )   (increasing in score)
#   discovery_t = 1 - R_t
# The alternative reference r(z1) = ratio scored on the labeled fraud sample z1
# (held out from the stream, so this is causal). z1 is NOT in the main cache; it
# is scored with the fitted GAM density-ratio model (no GAM refit).
#
# Consequence on this 0.1%-prevalence stream: for a null-looking transaction the
# F1-side CDF is essentially 0, so SAFFRON rejects the "is-signal" null almost
# everywhere and the 1-R convention flags only the highest-scoring tail. It is
# therefore an under-rejecting baseline here, like SCT/RTK.
if (RUN_ADJ_SAFFRON) {
  cat(sprintf("Adj-SAFFRON (F1-reference, 1-R): discoveries = %d | rejection rate = %.2f%%\n",
              sum(saffron_decisions), 100 * mean(saffron_decisions)))
}


# =========================================================================
# --- 4b. OMDRC.DD PRIOR BOUNDS: "tuned" (oracle) vs "causal" (warm-up) ---
# =========================================================================
# pi_bounds clips the EM prevalence estimate, and because q = pi*r/(1-pi+pi*r)
# saturates for the fraud-like tail but is linear in pi for the nulls, raising
# the bounds inflates the null q's, makes each skip more expensive, and forces
# more rejections (lower MDR, more alerts).  Two ways to set it:
#
#   "tuned"  : c(0.0015, 0.0030), picked post hoc because it drove
#              max_{t>20000} MDR to 0.0970 <= alpha = 0.10.  This uses the TEST
#              labels, so it cannot be defended as a data-driven choice.
#   "causal" : anchor pi0 by a single EM fit on the K0 warm-up UNLABELED ratios
#              only (no stream, no labels), then take [pi0/c, pi0*c] as a loose
#              sanity guard.  On CCFD the warm-up says pi0 = 0.001735, i.e. the
#              hand-tuned interval is recovered without ever touching the test
#              labels -- this is the defensible version.
#
# NOTE (bug fix): OMDRC_FROM_RATIO expects RAW density ratios in `LR_ini`.  The
# previous code passed omdrc_results$DR_ini, which is the Lmdr SCORE (median
# 5.2e-06 vs 1.0e-02 for the ratio), so the warm-up EM collapsed to its lower
# bound and the first ~220 stream windows were fed wrong-scale history.  We now
# pass omdrc_results$LR_ini.  SCT / RTK correctly consume DR_ini (Lmdr) and are
# untouched, as are the GAM density ratio and Adj-SAFFRON.
PI_BOUNDS_MODE  <- "causal"            # "causal" | "tuned"
PI_BOUNDS_ORACLE <- c(0.0015, 0.0030)  # used when PI_BOUNDS_MODE == "tuned"
PI_CAUSAL_WIDTH <- 3                   # bounds = [pi0/c, pi0*c]
PI_CAUSAL_RANGE <- c(1e-5, 0.5)        # uninformative range for the warm-up EM

LR_INI_STREAM <- omdrc_results$LR_ini  # raw warm-up density ratios (see NOTE)

if (identical(PI_BOUNDS_MODE, "causal")) {
  pi0_causal <- .estimate_local_pi_ratio(
    ratio = LR_INI_STREAM, init = mean(PI_CAUSAL_RANGE),
    pi_bounds = PI_CAUSAL_RANGE, tol = 1e-10, max_iter = 500L
  )$pi
  PI_BOUNDS_TUNED <- c(max(pi0_causal / PI_CAUSAL_WIDTH, 1e-6),
                       min(pi0_causal * PI_CAUSAL_WIDTH, 0.5))
  cat(sprintf(paste0("Prior bounds (causal): warm-up EM on K0=%d unlabeled",
                     " ratios gives pi0 = %.6f -> pi_bounds = [%.5f, %.5f]\n"),
              length(LR_INI_STREAM), pi0_causal,
              PI_BOUNDS_TUNED[1], PI_BOUNDS_TUNED[2]))
} else {
  PI_BOUNDS_TUNED <- PI_BOUNDS_ORACLE
  cat(sprintf("Prior bounds (oracle-tuned): [%.5f, %.5f]\n",
              PI_BOUNDS_TUNED[1], PI_BOUNDS_TUNED[2]))
}
TUNED_DECISIONS_CACHE <- "ccfd_omdrc_tuned_decisions.rds"
PI_SIG <- list(mode = PI_BOUNDS_MODE, bounds = PI_BOUNDS_TUNED,
               lr_ini = "LR_ini", decision_rule = "local_barrier_w100")

omdrc_decisions_orig <- omdrc_decisions
tuned_cache <- if (file.exists(TUNED_DECISIONS_CACHE)) readRDS(TUNED_DECISIONS_CACHE) else NULL
if (!is.null(tuned_cache) &&
    identical(tuned_cache$pi_sig, PI_SIG) &&
    identical(tuned_cache$alpha, ALPHA) &&
    identical(tuned_cache$n_stream, n_stream)) {
  cat("Loading tuned OMDRC.DD decisions from ", TUNED_DECISIONS_CACHE, "...\n", sep = "")
  omdrc_decisions <- integer(n_stream)
  omdrc_decisions[tuned_cache$idx] <- 1L
} else {
  cat(sprintf("Recomputing OMDRC.DD decisions with pi_bounds = [%.5f, %.5f]...\n",
              PI_BOUNDS_TUNED[1], PI_BOUNDS_TUNED[2]))
  omdrc_decisions <- as.integer(OMDRC_FROM_RATIO(
    LR = omdrc_results$LR, LR_ini = LR_INI_STREAM,
    alpha = ALPHA, pi_bounds = PI_BOUNDS_TUNED,
    em_tol = 1e-6, em_max_iter = 50L,
    D_mode = "growing", D_beta = 0.6, D_min = 10L, w = 100L
  )$de)
  saveRDS(list(pi_sig = PI_SIG, pi_bounds = PI_BOUNDS_TUNED, alpha = ALPHA,
               n_stream = n_stream, idx = which(omdrc_decisions == 1L)),
          TUNED_DECISIONS_CACHE)
}
cat(sprintf("OMDRC.DD discoveries: %d tuned vs %d original.\n",
            sum(omdrc_decisions), sum(omdrc_decisions_orig)))


# =========================================================================
# --- 4c. LOCAL-BARRIER REPLAY AND ALPHA-SWEEP PRECOMPUTATION ---
# =========================================================================
# Apply the manuscript default barrier window w=100 to the OMDRC.DD score path
# across the main run and the alpha sweep. gamma_t is a causal windowed FT
# threshold; the Lmdr score path is alpha-independent, so one gamma matrix (all
# sweep alphas) is precomputed once.
source(file.path(project_root, "code-semi/nofreeride/_nfr_core.R"))
D_NFR <- 100L
BARRIER_LEDGER <- "credited"          # OMDRC_OR_NFR (credited ledger)
BARRIER_SIG <- list(ledger = BARRIER_LEDGER, d = D_NFR, pi_sig = PI_SIG)
GAMMA_CACHE <- "ccfd_nfr_gamma_cache.rds"

# Tuned OMDRC.DD score path (alpha-independent); its $de equals the cached
# tuned decisions (verified). Needed to build gamma and to feed the barrier.
base_tuned <- OMDRC_FROM_RATIO(LR = omdrc_results$LR, LR_ini = LR_INI_STREAM,
                               alpha = ALPHA, pi_bounds = PI_BOUNDS_TUNED,
                               em_tol = 1e-6, em_max_iter = 50L,
                               D_mode = "growing", D_beta = 0.6,
                               D_min = 10L, w = D_NFR)
qh_nfr <- base_tuned$DR; qhi_nfr <- base_tuned$DR_ini

# SCT and RTK use exactly the same tuned Lmdr path as the reported OMDRC.DD
# arm, as specified in the manuscript baseline definitions.
sct_res <- STATIC_LMDR_DD(qh_nfr, qhi_nfr, alpha = ALPHA)
rtk_res <- ROLLING_TOPK_DD(qh_nfr, qhi_nfr, alpha = ALPHA,
                           window_size = length(qhi_nfr))
sct_decisions <- as.integer(sct_res$de)
rtk_decisions <- as.integer(rtk_res$de)

# Causal windowed FT-gamma for ALL sweep alphas in one pass (sort window once).
gamma_paths_multi <- function(qh, qhi, alphas, d) {
  q <- pmin(pmax(as.numeric(qh), 0), 1)
  hist0 <- pmin(pmax(as.numeric(qhi), 0), 1)
  m <- length(q); off <- length(hist0); full <- c(hist0, q)
  out <- matrix(Inf, m, length(alphas))
  for (t in seq_len(m)) {
    hi <- off + t; lo <- max(1L, hi - d + 1L)
    s <- sort(full[lo:hi]); tot <- sum(s)
    if (!is.finite(tot) || tot <= .Machine$double.eps) next
    cs <- cumsum(s) / tot
    k <- findInterval(alphas, cs)             # #{cs <= alpha}, cs monotone up
    ok <- k < length(s)
    out[t, ok] <- s[k[ok] + 1L]
  }
  out
}

ALPHA_SWEEP_ALL <- c(0.05, 0.08, 0.10, 0.12, 0.15, 0.20, 0.25, 0.30, 0.35)
gamma_cache <- if (file.exists(GAMMA_CACHE)) readRDS(GAMMA_CACHE) else NULL
if (!is.null(gamma_cache) && identical(gamma_cache$sig, BARRIER_SIG) &&
    identical(gamma_cache$alphas, ALPHA_SWEEP_ALL) &&
    nrow(gamma_cache$gamma) == n_stream) {
  cat("Loading cached NFR gamma matrix from ", GAMMA_CACHE, "...\n", sep = "")
  gamma_mat <- gamma_cache$gamma
} else {
  cat(sprintf("Computing NFR gamma matrix (d=%d, %d alphas, n=%d)...\n",
              D_NFR, length(ALPHA_SWEEP_ALL), n_stream))
  t0 <- Sys.time()
  gamma_mat <- gamma_paths_multi(qh_nfr, qhi_nfr, ALPHA_SWEEP_ALL, D_NFR)
  cat(sprintf("  gamma matrix done in %.1f s\n",
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  # sanity: one-pass gamma matches the reference nfr_gamma_path on a prefix
  chk <- nfr_gamma_path(qh_nfr[1:3000], ALPHA, D_NFR, qhi_nfr)
  stopifnot(max(abs(pmin(chk, 1e6) -
                     pmin(gamma_mat[1:3000, which(ALPHA_SWEEP_ALL == ALPHA)], 1e6))) < 1e-9)
  saveRDS(list(sig = BARRIER_SIG, alphas = ALPHA_SWEEP_ALL, gamma = gamma_mat),
          GAMMA_CACHE)
}
gamma_col <- function(a) gamma_mat[, which(abs(ALPHA_SWEEP_ALL - a) < 1e-9)]

# Credited-barrier decisions from the alpha-independent score path qh_nfr.
barrier_decisions <- function(a) {
  as.integer(OMDRC_OR_NFR(qh_nfr, a, d = D_NFR,
                          x.Lmdr_ini = qhi_nfr, gamma = gamma_col(a))$de)
}

omdrc_decisions_public <- omdrc_decisions
omdrc_decisions <- barrier_decisions(ALPHA)
stopifnot(identical(omdrc_decisions, omdrc_decisions_public))
cat(sprintf("OMDRC.DD local-barrier replay verified (w=%d): %d discoveries.\n",
            D_NFR, sum(omdrc_decisions)))


# =========================================================================
# --- 5. EVALUATION & VISUALIZATION (FOUR-PANEL PLOT VERSION) ---
# =========================================================================

# --- Step 1: calculate_metrics() additionally returns the cumulative rejection rate ---
calculate_metrics <- function(decisions, true_labels) {
  len <- length(decisions)
  true_labels <- true_labels[1:len]
  
  cumulative_tp <- cumsum(true_labels * decisions)
  cumulative_fp <- cumsum((1 - true_labels) * decisions)
  cumulative_discoveries <- cumsum(decisions)
  cumulative_signals <- cumsum(true_labels)
  
  # Cumulative rejection rate
  cumulative_rejection_rate <- cumulative_discoveries / seq_along(decisions)
  
  safe_discoveries <- ifelse(cumulative_discoveries == 0, 1, cumulative_discoveries)
  safe_signals <- ifelse(cumulative_signals == 0, 1, cumulative_signals)
  
  fdr_empirical <- cumulative_fp / safe_discoveries
  missed_signals <- cumulative_signals - cumulative_tp
  mdr_empirical <- missed_signals / safe_signals
  
  return(data.frame(
    mdr = mdr_empirical, 
    precision = 1 - fdr_empirical,
    fdr = fdr_empirical,
    tp = cumulative_tp,
    rej_rate = cumulative_rejection_rate
  ))
}

omdrc_metrics <- calculate_metrics(omdrc_decisions, theta)
sct_metrics <- calculate_metrics(sct_decisions, theta)
rtk_metrics <- calculate_metrics(rtk_decisions, theta)
if (RUN_ADJ_SAFFRON) {
  saffron_metrics <- calculate_metrics(saffron_decisions, theta)
}

# --- Step 2: assemble plot_df with all four methods (rejection rate included) ---
plot_df <- data.frame(
  time = 0:n_stream,
  mdr_omdrc = c(0, omdrc_metrics$mdr),
  mdr_sct = c(0, sct_metrics$mdr),
  mdr_rtk = c(0, rtk_metrics$mdr),
  tp_omdrc = c(0, omdrc_metrics$tp),
  tp_sct = c(0, sct_metrics$tp),
  tp_rtk = c(0, rtk_metrics$tp),
  mdr_gap_omdrc = c(0, omdrc_metrics$mdr - ALPHA),
  mdr_gap_sct = c(0, sct_metrics$mdr - ALPHA),
  mdr_gap_rtk = c(0, rtk_metrics$mdr - ALPHA),
  precision_omdrc = c(NA_real_, omdrc_metrics$precision),
  precision_sct = c(NA_real_, sct_metrics$precision),
  precision_rtk = c(NA_real_, rtk_metrics$precision),
  fdr_omdrc = c(NA_real_, omdrc_metrics$fdr),
  fdr_sct = c(NA_real_, sct_metrics$fdr),
  fdr_rtk = c(NA_real_, rtk_metrics$fdr),
  rej_rate_omdrc = c(0, omdrc_metrics$rej_rate),
  rej_rate_sct = c(0, sct_metrics$rej_rate),
  rej_rate_rtk = c(0, rtk_metrics$rej_rate)
)
if (RUN_ADJ_SAFFRON) {
  plot_df$mdr_saffron <- c(0, saffron_metrics$mdr)
  plot_df$tp_saffron <- c(0, saffron_metrics$tp)
  plot_df$mdr_gap_saffron <- c(0, saffron_metrics$mdr - ALPHA)
  plot_df$precision_saffron <- c(NA_real_, saffron_metrics$precision)
  plot_df$fdr_saffron <- c(NA_real_, saffron_metrics$fdr)
  plot_df$rej_rate_saffron <- c(0, saffron_metrics$rej_rate)
}

# Finite-sample (Monte-Carlo) tolerance for the CUMULATIVE empirical MDR.  Only
# S_t = sum_{i<=t} theta_i signals have been seen by time t, so a single missed
# fraud moves MDR_t by 1/S_t (0.35 pp at the end of this stream) and the path
# fluctuates around the target with SE = sqrt(alpha(1-alpha)/S_t).  Figure 3
# panels (a) and (c) shade [alpha, alpha + 1.96 SE] so that excursions of a few
# frauds are not read as a violation of the constraint.  The band is drawn only
# once >= 25 signals have accrued (before that the SE exceeds alpha itself).
MC_TOL_Z <- 1.96
MC_TOL_MIN_SIG <- 25L
sig_cum <- cumsum(theta)
mdr_tol_hi <- ALPHA + MC_TOL_Z * sqrt(ALPHA * (1 - ALPHA) / pmax(sig_cum, 1))
mdr_tol_hi[sig_cum < MC_TOL_MIN_SIG] <- NA_real_
plot_df$mdr_tol_hi <- c(NA_real_, mdr_tol_hi)

# Method display names / colours consistent with the Setting 1-3 figures.
method_levels <- c("OMDRC.DD", "SCT", "RTK")
if (RUN_ADJ_SAFFRON) method_levels <- c(method_levels, "Adj-SAFFRON")
my_colors_all <- c("OMDRC.DD" = "#7CAE00", "SCT" = "#619CFF",
                   "RTK" = "#FF61C3", "Adj-SAFFRON" = "#C77CFF")
# Solid lines for all methods (colour-only encoding, matching the Setting 1-3
# figures); the legend colours carry the distinction.
my_linetypes_all <- c("OMDRC.DD" = "solid", "SCT" = "solid",
                      "RTK" = "solid", "Adj-SAFFRON" = "solid")
my_colors <- my_colors_all[method_levels]
my_linetypes <- my_linetypes_all[method_levels]
plot_linewidth <- 1.0

# Typography: the figures are placed in a two-column LaTeX layout, so they are
# down-scaled to the text width and every label must survive that reduction.
# All panel text derives from FONT_BASE and all annotate() text from the two
# ANN_SIZE_* constants, so the whole set can be rescaled in one place.
FONT_BASE <- 16
ANN_SIZE_MAIN <- 4.6   # in-panel reference labels (target MDR, total frauds)
ANN_SIZE_SMALL <- 3.8  # explanatory notes (clipping, band definition, k*)

# --- Build the four panels (four-method comparison) ---
plot_a <- ggplot(plot_df, aes(x = time)) +
  geom_ribbon(aes(ymin = ALPHA, ymax = mdr_tol_hi), fill = "red", alpha = 0.08) +
  geom_line(aes(y = mdr_omdrc, color = "OMDRC.DD", linetype = "OMDRC.DD"), size = plot_linewidth) +
  geom_line(aes(y = mdr_sct, color = "SCT", linetype = "SCT"), size = plot_linewidth) +
  geom_line(aes(y = mdr_rtk, color = "RTK", linetype = "RTK"), size = plot_linewidth) +
  geom_hline(yintercept = ALPHA, linetype = "dotted", color = "red", size = 1) +
  annotate("text", x = n_stream * 0.62, y = ALPHA + 0.055, label = paste("Target MDR =", ALPHA), color = "red", size = ANN_SIZE_MAIN) +
  labs(title = "(a) MDR Control", x = "Time (t)", y = "Empirical MDR") +
  theme_minimal(base_size = FONT_BASE) +
  coord_cartesian(ylim = c(0, NA))
if (RUN_ADJ_SAFFRON) {
  plot_a <- plot_a + geom_line(aes(y = mdr_saffron, color = "Adj-SAFFRON", linetype = "Adj-SAFFRON"), size = plot_linewidth)
}

total_fraud_count_in_stream <- sum(theta)
plot_b <- ggplot(plot_df, aes(x = time)) +
  geom_line(aes(y = tp_omdrc, color = "OMDRC.DD", linetype = "OMDRC.DD"), size = plot_linewidth) +
  geom_line(aes(y = tp_sct, color = "SCT", linetype = "SCT"), size = plot_linewidth) +
  geom_line(aes(y = tp_rtk, color = "RTK", linetype = "RTK"), size = plot_linewidth) +
  geom_hline(yintercept = total_fraud_count_in_stream, linetype = "dotted", color = "grey30", size = 1) +
  annotate("text", x = n_stream * 0.6, y = total_fraud_count_in_stream * 0.9, label = paste("Total Frauds:", total_fraud_count_in_stream), color = "grey30", size = ANN_SIZE_MAIN) +
  labs(title = "(b) True Discoveries", x = "Time (t)", y = "True Frauds Detected") +
  theme_minimal(base_size = FONT_BASE)
if (RUN_ADJ_SAFFRON) {
  plot_b <- plot_b + geom_line(aes(y = tp_saffron, color = "Adj-SAFFRON", linetype = "Adj-SAFFRON"), size = plot_linewidth)
}

plot_c <- ggplot(plot_df, aes(x = time)) +
  geom_line(aes(y = fdr_omdrc, color = "OMDRC.DD", linetype = "OMDRC.DD"), size = plot_linewidth) +
  geom_line(aes(y = fdr_sct, color = "SCT", linetype = "SCT"), size = plot_linewidth) +
  geom_line(aes(y = fdr_rtk, color = "RTK", linetype = "RTK"), size = plot_linewidth) +
  labs(title = "(c) Empirical FDR", x = "Time (t)", y = "Empirical FDR") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_minimal(base_size = FONT_BASE)
if (RUN_ADJ_SAFFRON) {
  plot_c <- plot_c + geom_line(aes(y = fdr_saffron, color = "Adj-SAFFRON", linetype = "Adj-SAFFRON"), size = plot_linewidth)
}

# --- Step 3: the rejection-rate panel (plot_d) ---
# This is the REJECTION-RATE panel of Figure 4.  The y-axis is clipped to
# REJ_RATE_YMAX because every method starts at a 100% cumulative rejection rate
# for the first few transactions (a burn-in artefact of the cumulative ratio);
# on the steady-state scale the curves live in [0, 5%].
REJ_RATE_YMAX <- 0.10
plot_d <- ggplot(plot_df, aes(x = time)) +
  geom_line(aes(y = rej_rate_omdrc, color = "OMDRC.DD", linetype = "OMDRC.DD"), size = plot_linewidth) +
  geom_line(aes(y = rej_rate_sct, color = "SCT", linetype = "SCT"), size = plot_linewidth) +
  geom_line(aes(y = rej_rate_rtk, color = "RTK", linetype = "RTK"), size = plot_linewidth) +
  labs(title = "(b) Rejection Rate", x = "Time (t)",
       y = "Cumulative Rejection Rate") +
  scale_y_continuous(labels = scales::percent) + # format the y axis as percentages
  annotate("text", x = n_stream * 0.02, y = REJ_RATE_YMAX * 0.97, hjust = 0,
           vjust = 1, size = ANN_SIZE_SMALL, color = "grey45",
           label = "initial transient clipped") +
  coord_cartesian(ylim = c(0, REJ_RATE_YMAX)) +
  theme_minimal(base_size = FONT_BASE)
if (RUN_ADJ_SAFFRON) {
  plot_d <- plot_d + geom_line(aes(y = rej_rate_saffron, color = "Adj-SAFFRON", linetype = "Adj-SAFFRON"), size = plot_linewidth)
}

# =========================================================================
# --- 6. ALPHA SWEEP FOR FIGURE 3 (c)/(d): max_T MDR_t & TOTAL DISCOVERIES ---
# =========================================================================
# For each target alpha in [0.05, 0.35], recompute the online decisions of
# all methods from the CACHED density-ratio scores / p-values (no GAM refit):
#   (e) maximum of the cumulative MDR(t) over the stream (after burn-in);
#   (f) total number of discoveries (rejections) over the whole stream.
# (these become Figure 3 panels (e) and (f) in the manuscript layout.)
# OMDRC.DD uses the causal prior bounds above. Decisions per (alpha, method)
# are cached in DECISIONS_CACHE (invalidated if pi_bounds change).

cat("\nRunning alpha sweep for panels (e)/(f)...\n")

# ALPHA_SWEEP_MIN controls the smallest target reported in panels (e)/(f).
# The manuscript plot begins at 0.08; set this to 0.05 for an additional
# diagnostic point (the cached score fit is reused).
ALPHA_SWEEP_MIN <- 0.08
alpha_sweep <- c(0.05, 0.08, 0.10, 0.12, 0.15, 0.20, 0.25, 0.30, 0.35)
alpha_sweep <- alpha_sweep[alpha_sweep >= ALPHA_SWEEP_MIN]
BURN_IN <- 20000
DECISIONS_CACHE <- "figure3_sweep_decisions_cache.rds"
SAFFRON_SWEEP_CACHE <- "saffron_alpha_sweep_f1ref_1minusR_cache.rds"
post_burn <- (BURN_IN + 1):n_stream

load_or_compute_saffron <- function(pvals, alpha_val) {
  if (file.exists(SAFFRON_SWEEP_CACHE)) {
    sc <- readRDS(SAFFRON_SWEEP_CACHE)
    key <- as.character(alpha_val)
    if (!is.null(sc[[key]])) return(sc[[key]])
  }
  cat(sprintf("  [compute] SAFFRON alpha=%.2f...\n", alpha_val))
  res <- SAFFRON(pvals, alpha = alpha_val)
  de <- as.integer(1L - res$R)   # swapped (FNR) convention: discovery = 1 - R
  sc <- if (file.exists(SAFFRON_SWEEP_CACHE)) readRDS(SAFFRON_SWEEP_CACHE) else list()
  sc[[as.character(alpha_val)]] <- de
  saveRDS(sc, SAFFRON_SWEEP_CACHE)
  de
}

dec_cache <- if (file.exists(DECISIONS_CACHE)) readRDS(DECISIONS_CACHE) else list()
if (!identical(dec_cache$pi_sig, PI_SIG) ||
    !identical(dec_cache$barrier, BARRIER_SIG)) {
  cat("Sweep decisions cache missing or stale (pi_bounds/barrier changed); recomputing.\n")
  dec_cache <- list(pi_sig = PI_SIG, barrier = BARRIER_SIG, decisions = list())
}

get_sweep_decisions <- function(a) {
  key <- as.character(a)
  if (!is.null(dec_cache$decisions[[key]])) {
    cat(sprintf("  [cache] decisions for alpha=%.2f loaded\n", a))
    return(lapply(dec_cache$decisions[[key]], function(idx) {
      de <- integer(n_stream); de[idx] <- 1L; de
    }))
  }
  cat(sprintf("  [compute] decisions for alpha=%.2f...\n", a))
  # The Lmdr score path is alpha-independent, so barrier_decisions() reuses the
  # precomputed gamma for the manuscript OMDRC.DD rule.
  de_omdrc_nfr <- barrier_decisions(a)
  res_sct <- STATIC_LMDR_DD(qh_nfr, qhi_nfr, alpha = a)
  res_rtk <- ROLLING_TOPK_DD(qh_nfr, qhi_nfr,
                             alpha = a, window_size = length(qhi_nfr))
  de_list <- list("OMDRC.DD" = de_omdrc_nfr,
                  "SCT" = as.integer(res_sct$de),
                  "RTK" = as.integer(res_rtk$de))
  if (RUN_ADJ_SAFFRON) {
    de_list[["Adj-SAFFRON"]] <- load_or_compute_saffron(p_values_saffron, a)
  }
  dec_cache$decisions[[key]] <<- lapply(de_list, function(de) which(de == 1L))
  saveRDS(dec_cache, DECISIONS_CACHE)
  de_list
}

alpha_summary <- list()
for (a in alpha_sweep) {
  de_list <- get_sweep_decisions(a)
  for (mn in method_levels) {
    if (is.null(de_list[[mn]])) next
    m_all <- calculate_metrics(de_list[[mn]], theta)
    alpha_summary[[length(alpha_summary) + 1]] <- data.frame(
      alpha   = a,
      method  = mn,
      max_mdr = max(m_all$mdr[post_burn]),
      max_fdr = max(m_all$fdr[post_burn]),
      n_disc  = sum(de_list[[mn]])
    )
  }
}
alpha_summary_df <- do.call(rbind, alpha_summary)
alpha_summary_df$method <- factor(alpha_summary_df$method, levels = method_levels)

cat(sprintf("\n--- Alpha sweep summary (max over t > %d, burn-in excluded) ---\n", BURN_IN))
cat("\n[Max MDR]\n")
print(reshape(alpha_summary_df[, c("alpha", "method", "max_mdr")],
              idvar = "alpha", timevar = "method", direction = "wide"),
      row.names = FALSE)
cat("\n[Total discoveries]\n")
print(reshape(alpha_summary_df[, c("alpha", "method", "n_disc")],
              idvar = "alpha", timevar = "method", direction = "wide"),
      row.names = FALSE)

# --- Alpha sweep panels: legends suppressed, the collected Method legend of
# Figure 3 panels (a)/(b) identifies the colours. ---
method_shapes_all <- c("OMDRC.DD" = 16, "SCT" = 17, "RTK" = 15, "Adj-SAFFRON" = 18)
method_shapes <- method_shapes_all[method_levels]
# Axis breaks follow the sweep actually plotted: the 0.05-grid restricted to the
# reported range, plus the smallest target itself when it is far enough from the
# next grid point to be labelled without the two labels colliding.
sweep_grid <- seq(0.05, 0.35, by = 0.05)
sweep_grid <- sweep_grid[sweep_grid >= min(alpha_sweep)]
sweep_breaks <- if (min(sweep_grid) - min(alpha_sweep) >= 0.03) {
  sort(c(min(alpha_sweep), sweep_grid))
} else {
  sweep_grid
}
# At the two-column font sizes the dense 0.05-grid labels run into each other,
# so panels (e)/(f) are labelled only at 0.1, 0.2, 0.3 (the plotted sweep is
# unchanged); the ticks are also drawn in a smaller font via axis_tick_theme.
sweep_breaks <- c(0.1, 0.2, 0.3)
sweep_breaks <- sweep_breaks[sweep_breaks >= min(alpha_sweep) &
                             sweep_breaks <= max(alpha_sweep)]
sweep_point_size <- 3.0
# Smaller axis-tick text just for the alpha-sweep panels, so the x labels stop
# overlapping while the axis titles keep the larger FONT_BASE size.
axis_tick_theme <- theme(axis.text = element_text(size = rel(0.75)))


# Same MC tolerance as Figure 3 panel (a), now as a function of the target alpha
# with the full signal count of the stream (the sweep maxima are attained near
# t = n, where S_t is within a few of total_fraud_count_in_stream).
n_sig_eval <- total_fraud_count_in_stream
tol_band_df <- data.frame(alpha = seq(min(alpha_sweep), max(alpha_sweep),
                                      length.out = 200))
tol_band_df$hi <- tol_band_df$alpha +
  MC_TOL_Z * sqrt(tol_band_df$alpha * (1 - tol_band_df$alpha) / n_sig_eval)

plot_e <- ggplot(alpha_summary_df,
                 aes(x = alpha, y = max_mdr, color = method,
                     shape = method, group = method)) +
  geom_ribbon(data = tol_band_df, aes(x = alpha, ymin = alpha, ymax = hi),
              inherit.aes = FALSE, fill = "grey55", alpha = 0.25) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey50", size = 0.8) +
  annotate("text", x = min(alpha_sweep) + 0.005, y = 0.048,
           label = paste0("MDR = alpha (constraint boundary)\n",
                          "band: alpha + 1.96 x MC SE (", n_sig_eval, " signals)"),
           color = "grey50", size = ANN_SIZE_SMALL, hjust = 0, vjust = 1, lineheight = 1.05) +
  geom_line(size = plot_linewidth, show.legend = FALSE) +
  geom_point(size = sweep_point_size, show.legend = FALSE) +
  labs(title = "(c) Max MDR vs alpha",
       x = "Target alpha", y = "Max Empirical MDR") +
  scale_shape_manual(values = method_shapes, guide = "none") +
  scale_x_continuous(breaks = sweep_breaks) +
  scale_y_continuous(
    sec.axis = sec_axis(~ . * n_sig_eval,
                        name = sprintf("Missed frauds (of %d)", n_sig_eval))) +
  # Extend the axis to cover the WHOLE tolerance band and the MDR = alpha
  # boundary over the entire sweep (both were previously clipped at large
  # alpha), so the constraint reference is drawn in full rather than running
  # off the top-right corner.
  coord_cartesian(ylim = c(0, max(max(alpha_summary_df$max_mdr),
                                  max(tol_band_df$hi)) * 1.05)) +
  theme_minimal(base_size = FONT_BASE) +
  theme(panel.grid.minor = element_blank()) +
  axis_tick_theme

plot_f <- ggplot(alpha_summary_df,
                 aes(x = alpha, y = n_disc, color = method,
                     shape = method, group = method)) +
  geom_line(size = plot_linewidth, show.legend = FALSE) +
  geom_point(size = sweep_point_size, show.legend = FALSE) +
  labs(title = "(d) Discoveries vs alpha",
       x = "Target alpha", y = "Total Discoveries (log scale)") +
  scale_shape_manual(values = method_shapes, guide = "none") +
  scale_x_continuous(breaks = sweep_breaks) +
  scale_y_log10(labels = scales::comma,
                breaks = c(100, 300, 1000, 3000, 10000, 30000, 100000)) +
  theme_minimal(base_size = FONT_BASE) +
  theme(panel.grid.minor = element_blank()) +
  axis_tick_theme


# =========================================================================
# --- 6b. MDR-vs-DISCOVERIES FRONTIER (FIGURE 4, PANEL (a)) ---
# =========================================================================
# Reviewer rebuttal to "SCT/RTK make far fewer discoveries than OMDRC": all four
# methods threshold the SAME estimated Lmdr score, so a fixed alert budget k
# pins down the achievable MDR.  Ranking the stream by that score and rejecting
# the top k traces the frontier MDR(k) = 1 - TP(k)/S; the target MDR = alpha is
# reachable only from k >= k*.  SCT/RTK sit at k ~ 300-500, where the frontier
# itself is at MDR ~ 0.22 -- their leaner alert set is an INFEASIBLE operating
# point, not a more efficient one.  NOTE: this is the top-k curve of the shared
# score under a single GLOBAL threshold; online procedures use time-varying
# thresholds and can therefore sit slightly ABOVE it at the same budget (RTK
# does), so the curve is a reference frontier, not a hard bound.
# The operating points are read from the LIVE decision vectors of this run, so
# the panel can never drift out of sync with panels (a)-(d).
front_ord <- order(omdrc_results$DR, decreasing = TRUE)
front_mdr <- 1 - cumsum(theta[front_ord]) / total_fraud_count_in_stream
kstar <- which(front_mdr <= ALPHA)[1]

op_point <- function(de) data.frame(n_disc = sum(de),
                                    mdr = 1 - sum(theta * de) / total_fraud_count_in_stream)
op_list <- list("OMDRC.DD" = omdrc_decisions, "SCT" = sct_decisions,
                "RTK" = rtk_decisions)
if (RUN_ADJ_SAFFRON) op_list[["Adj-SAFFRON"]] <- saffron_decisions
op_df <- do.call(rbind, lapply(method_levels, function(mn)
  cbind(method = mn, op_point(op_list[[mn]]))))
op_df$method <- factor(op_df$method, levels = method_levels)

cat(sprintf("\nFrontier: controlling MDR <= %.2f on this stream needs k* = %d rejections\n",
            ALPHA, kstar))
cat(sprintf("Operating points at alpha = %.2f:\n", ALPHA))
print(op_df, row.names = FALSE)

# Drop the first 50 points of the frontier: with k < 50 the curve is still above
# 0.8 and only stretches the axis.
front_df <- data.frame(k = 50:n_stream, mdr = front_mdr[50:n_stream])
front_ymax <- min(0.7, max(op_df$mdr) + 0.45)

plot_frontier <- ggplot() +
  # Same MC tolerance as Figure 3: with S = 282 signals the end-of-stream MDR has
  # SE = sqrt(alpha(1-alpha)/S), so the operating points are judged against
  # [alpha, alpha + 1.96 SE] rather than against alpha exactly.
  annotate("rect", xmin = 50, xmax = n_stream, ymin = ALPHA,
           ymax = ALPHA + MC_TOL_Z * sqrt(ALPHA * (1 - ALPHA) / n_sig_eval),
           fill = "red", alpha = 0.08) +
  geom_line(data = front_df, aes(x = k, y = mdr), color = "grey35",
            size = 0.9) +
  geom_hline(yintercept = ALPHA, linetype = "dotted", color = "red", size = 1) +
  geom_vline(xintercept = kstar, linetype = "dashed", color = "grey45",
             size = 0.7) +
  annotate("text", x = kstar * 1.15, y = front_ymax * 0.92, hjust = 0,
           size = ANN_SIZE_SMALL, color = "grey35",
           label = sprintf("MDR <= %.2f needs\n>= %s rejections", ALPHA,
                           format(kstar, big.mark = ","))) +
  # The label must clear BOTH the MC tolerance band (which sits immediately
  # above alpha) and the dotted alpha line, so it goes just below alpha: the
  # left part of the panel is empty there because the frontier only crosses
  # alpha at k*.
  annotate("text", x = 55, y = ALPHA - 0.045, hjust = 0, color = "red",
           size = ANN_SIZE_SMALL, label = sprintf("Target MDR = %.2f", ALPHA)) +
  geom_point(data = op_df, aes(x = n_disc, y = mdr, color = method,
                               shape = method),
             size = 4.0, show.legend = FALSE) +
  scale_shape_manual(values = method_shapes, guide = "none") +
  scale_x_log10(breaks = c(100, 1000, 10000, 100000),
                labels = scales::comma) +
  labs(title = "(a) MDR-Discoveries Frontier",
       x = "# discoveries (rejections, log scale)",
       y = "Empirical MDR at end of stream") +
  coord_cartesian(ylim = c(0, front_ymax)) +
  theme_minimal(base_size = FONT_BASE) +
  theme(panel.grid.minor = element_blank())


# =========================================================================
# --- 7. ASSEMBLE AND SAVE FIGURE 3 (CONTROL) AND FIGURE 4 (COST) ---
# =========================================================================
# The OMDRC.DD arm carries the local-barrier decisions from Algorithm 2,
# while the figure legend retains the concise method name used elsewhere.
#
# The old single 3x2 figure mixed two claims that a reader must not read against
# each other ("MDR is controlled" and "FDR is ~98%"), so it is split:
#   Figure 3 (2x2, control)  : (a) MDR(t)  (b) true discoveries
#                              (c) max MDR vs alpha  (d) discoveries vs alpha
#   Figure 4 (1x3, cost)     : (a) MDR-discoveries frontier
#                              (b) rejection rate  (c) empirical FDR
# Figure 4 must be read left to right: k* explains why FDR is structurally
# pinned near 1 at 0.1% prevalence.
method_labels <- c("OMDRC.DD" = "OMDRC.DD", "SCT" = "SCT",
                   "RTK" = "RTK", "Adj-SAFFRON" = "Adj-SAFFRON")[method_levels]

apply_method_scales <- function(p) {
  p + plot_layout(guides = 'collect') &
    scale_color_manual(name = "Method", values = my_colors,
                       breaks = method_levels, labels = method_labels) &
    scale_linetype_manual(name = "Method", values = my_linetypes,
                          breaks = method_levels, labels = method_labels) &
    theme(legend.position = 'bottom')
}

# Figure 3: MDR control.  Legend is collected from panels (a)/(b); the alpha
# sweep panels suppress their own guides.
fig3_plot <- apply_method_scales((plot_a + plot_b) / (plot_e + plot_f))
ggsave("Figure3_CCFD_mdr_control_revision_diagnostic.pdf", plot = fig3_plot,
       width = 10, height = 8, units = "in")
ggsave("Figure3_CCFD_mdr_control_revision_diagnostic.png", plot = fig3_plot,
       width = 10, height = 8, units = "in", dpi = 300)
cat("Saved Figure3_CCFD_mdr_control_revision_diagnostic.pdf/png.\n")

# Figure 4: operational cost and feasibility.  Legend is collected from panels
# (b)/(c); the frontier panel suppresses its own guides (same colour code).
fig4_plot <- apply_method_scales(plot_frontier + plot_d + plot_c)
ggsave("Figure4_CCFD_cost_feasibility_revision_diagnostic.pdf", plot = fig4_plot,
       width = 13, height = 5, units = "in")
ggsave("Figure4_CCFD_cost_feasibility_revision_diagnostic.png", plot = fig4_plot,
       width = 13, height = 5, units = "in", dpi = 300)
cat("Saved Figure4_CCFD_cost_feasibility_revision_diagnostic.pdf/png.\n")

# Legacy 3x2 layout retained for the manuscript's six-panel presentation.
# Reuse the clipped rejection-rate panel, but restore its position-specific
# label to (d); Figure 4 continues to use the same panel as (b).
SAVE_LEGACY_3x2 <- FALSE
if (SAVE_LEGACY_3x2) {
  plot_d_3x2 <- plot_d + labs(title = "(d) Rejection Rate")
  plot_e_3x2 <- plot_e + labs(title = "(e) Max MDR vs alpha")
  plot_f_3x2 <- plot_f + labs(title = "(f) Discoveries vs alpha")
  legacy_plot <- apply_method_scales(
    (plot_a + plot_b) / (plot_c + plot_d_3x2) / (plot_e_3x2 + plot_f_3x2))
  # Canvas kept slightly narrower than the panel content needs, so that after
  # LaTeX scales it to the two-column text width the labels stay legible.
  ggsave("Figure3_CCFD_online_comparison_legacy_copy_revision.pdf", plot = legacy_plot,
         width = 10, height = 11, units = "in")
  ggsave("Figure3_CCFD_online_comparison_legacy_copy_revision.png", plot = legacy_plot,
         width = 10, height = 11, units = "in", dpi = 300)
  cat("Saved Figure3_CCFD_online_comparison_legacy_copy_revision.pdf/png.\n")
}


# =========================================================================
# --- 8. REVISED SECTION 5 FIGURE: ORIGINAL PIPELINE, NEW OUTPUT NAMES ---
# =========================================================================
# This figure is the manuscript revision. It uses the same cached data split,
# fitted score path, and decisions as the original Figure 3. The original
# Figure3_CCFD_online_comparison.png/pdf are never written by this script.
# Panels (e)-(f) report the maximum empirical MDR and total discoveries as
# functions of alpha, matching the revised manuscript caption.

revision_summary <- list()
for (a in alpha_sweep) {
  de_list <- get_sweep_decisions(a)
  for (mn in method_levels) {
    if (is.null(de_list[[mn]])) next
    de <- as.integer(de_list[[mn]])
    tp <- sum(theta * de)
    flagged <- sum(de)
    mdr_path <- calculate_metrics(de, theta)$mdr
    revision_summary[[length(revision_summary) + 1L]] <- data.frame(
      alpha = a,
      method = mn,
      max_mdr = max(mdr_path),
      terminal_mdr = 1 - tp / sum(theta),
      detected_frauds = tp,
      flagged = flagged,
      false_alarms = flagged - tp,
      review_rate = flagged / length(theta),
      empirical_fdr = if (flagged > 0) (flagged - tp) / flagged else 0
    )
  }
}
revision_summary_df <- do.call(rbind, revision_summary)
revision_summary_df$method <- factor(revision_summary_df$method,
                                     levels = method_levels)
write.csv(revision_summary_df, "CCFD_section5_alpha_summary_revised.csv",
          row.names = FALSE)

plot_a_revision <- ggplot(plot_df, aes(x = time)) +
  geom_line(aes(y = mdr_omdrc, color = "OMDRC.DD"), size = plot_linewidth) +
  geom_line(aes(y = mdr_sct, color = "SCT"), size = plot_linewidth) +
  geom_line(aes(y = mdr_rtk, color = "RTK"), size = plot_linewidth) +
  geom_hline(yintercept = ALPHA, linetype = "dashed", color = "grey30",
             size = 0.8) +
  labs(title = "(a) Empirical MDR", x = "Transaction index",
       y = "Cumulative MDR") +
  coord_cartesian(ylim = c(0, max(plot_df$mdr_saffron, na.rm = TRUE) * 1.04)) +
  theme_minimal(base_size = FONT_BASE)
if (RUN_ADJ_SAFFRON) {
  plot_a_revision <- plot_a_revision +
    geom_line(aes(y = mdr_saffron, color = "Adj-SAFFRON"),
              size = plot_linewidth)
}

plot_b_revision <- ggplot(plot_df, aes(x = time)) +
  geom_line(aes(y = tp_omdrc, color = "OMDRC.DD"), size = plot_linewidth) +
  geom_line(aes(y = tp_sct, color = "SCT"), size = plot_linewidth) +
  geom_line(aes(y = tp_rtk, color = "RTK"), size = plot_linewidth) +
  geom_hline(yintercept = total_fraud_count_in_stream, linetype = "dashed",
             color = "grey30", size = 0.8) +
  labs(title = "(b) Detected frauds", x = "Transaction index",
       y = "Cumulative true positives") +
  theme_minimal(base_size = FONT_BASE)
if (RUN_ADJ_SAFFRON) {
  plot_b_revision <- plot_b_revision +
    geom_line(aes(y = tp_saffron, color = "Adj-SAFFRON"),
              size = plot_linewidth)
}

plot_c_revision <- ggplot(plot_df, aes(x = time)) +
  geom_line(aes(y = fdr_omdrc, color = "OMDRC.DD"), size = plot_linewidth,
            na.rm = TRUE) +
  geom_line(aes(y = fdr_sct, color = "SCT"), size = plot_linewidth,
            na.rm = TRUE) +
  geom_line(aes(y = fdr_rtk, color = "RTK"), size = plot_linewidth,
            na.rm = TRUE) +
  labs(title = "(c) Empirical FDR", x = "Transaction index",
       y = "Cumulative FDR") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1)) +
  theme_minimal(base_size = FONT_BASE)
if (RUN_ADJ_SAFFRON) {
  plot_c_revision <- plot_c_revision +
    geom_line(aes(y = fdr_saffron, color = "Adj-SAFFRON"),
              size = plot_linewidth, na.rm = TRUE)
}

plot_d_revision <- ggplot(plot_df, aes(x = time)) +
  geom_line(aes(y = rej_rate_omdrc, color = "OMDRC.DD"),
            size = plot_linewidth) +
  geom_line(aes(y = rej_rate_sct, color = "SCT"), size = plot_linewidth) +
  geom_line(aes(y = rej_rate_rtk, color = "RTK"), size = plot_linewidth) +
  labs(title = "(d) Review rate", x = "Transaction index",
       y = "Cumulative fraction flagged") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  coord_cartesian(ylim = c(0, REJ_RATE_YMAX)) +
  theme_minimal(base_size = FONT_BASE)
if (RUN_ADJ_SAFFRON) {
  plot_d_revision <- plot_d_revision +
    geom_line(aes(y = rej_rate_saffron, color = "Adj-SAFFRON"),
              size = plot_linewidth)
}

plot_e_revision <- ggplot(
  revision_summary_df,
  aes(x = alpha, y = max_mdr, color = method,
      shape = method, group = method)
) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey35", size = 0.8) +
  geom_line(size = plot_linewidth, show.legend = FALSE) +
  geom_point(size = sweep_point_size, show.legend = FALSE) +
  scale_shape_manual(values = method_shapes, guide = "none") +
  scale_x_continuous(breaks = c(0.1, 0.2, 0.3)) +
  labs(title = "(e) Max MDR vs target",
       x = "Target MDR level", y = "Maximum empirical MDR") +
  theme_minimal(base_size = FONT_BASE) +
  theme(panel.grid.minor = element_blank())

plot_f_revision <- ggplot(
  revision_summary_df,
  aes(x = alpha, y = flagged, color = method,
      shape = method, group = method)
) +
  geom_line(size = plot_linewidth, show.legend = FALSE) +
  geom_point(size = sweep_point_size, show.legend = FALSE) +
  scale_shape_manual(values = method_shapes, guide = "none") +
  scale_x_continuous(breaks = c(0.1, 0.2, 0.3)) +
  scale_y_log10(labels = scales::comma,
                breaks = c(50, 100, 300, 1000, 3000, 10000, 30000)) +
  labs(title = "(f) Discoveries vs target",
       x = "Target MDR level", y = "Total discoveries (log scale)") +
  theme_minimal(base_size = FONT_BASE) +
  theme(panel.grid.minor = element_blank())

figure3_section5_revision <-
  ((plot_a_revision + plot_b_revision) /
   (plot_c_revision + plot_d_revision) /
   (plot_e_revision + plot_f_revision)) +
  plot_layout(guides = "collect") &
  scale_color_manual(name = "Method", values = my_colors,
                     breaks = method_levels, labels = method_labels) &
  theme(legend.position = "bottom")

ggsave("Figure3_CCFD_online_comparison_revised.pdf",
       plot = figure3_section5_revision,
       width = 10, height = 11, units = "in")
ggsave("Figure3_CCFD_online_comparison_revised.png",
       plot = figure3_section5_revision,
       width = 10, height = 11, units = "in", dpi = 300)
cat("Saved Figure3_CCFD_online_comparison_revised.pdf/png.\n")

# --- Single-column manuscript version -------------------------------------
# The journal uses a two-column layout.  This version is designed to be placed
# at \columnwidth (about one half of the page width), rather than shrinking the
# 10-inch-wide version above.  A 7 x 7.8 inch source canvas gives a final height
# of about 3.8 inches at column width, while the 14-point source typography
# scales to approximately 7 points in the manuscript.
SC_FONT_BASE <- 14
SC_ANN_MAIN <- 4.0
SC_ANN_SMALL <- 3.3
sc_time_breaks <- c(0, 100000, 200000)
sc_time_labels <- c("0", "100k", "200k")
sc_theme <- theme_minimal(base_size = SC_FONT_BASE) +
  theme(
    plot.title = element_text(size = rel(1.04), margin = margin(b = 2)),
    axis.title = element_text(size = rel(0.92)),
    axis.text = element_text(size = rel(0.80)),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t = 3, r = 4, b = 2, l = 3)
  )

plot_a_singlecolumn <- ggplot(plot_df, aes(x = time)) +
  geom_ribbon(aes(ymin = ALPHA, ymax = mdr_tol_hi),
              fill = "red", alpha = 0.08) +
  geom_line(aes(y = mdr_omdrc, color = "OMDRC.DD"),
            size = plot_linewidth) +
  geom_line(aes(y = mdr_sct, color = "SCT"), size = plot_linewidth) +
  geom_line(aes(y = mdr_rtk, color = "RTK"), size = plot_linewidth) +
  geom_hline(yintercept = ALPHA, linetype = "dotted", color = "red",
             size = 0.9) +
  annotate("text", x = n_stream * 0.98, y = ALPHA + 0.025,
           label = sprintf("Target MDR = %.1f", ALPHA), color = "red",
           size = SC_ANN_SMALL, hjust = 1, vjust = 0) +
  scale_x_continuous(breaks = sc_time_breaks, labels = sc_time_labels) +
  labs(title = "(a) MDR Control", x = "Time (t)", y = "Empirical MDR") +
  coord_cartesian(ylim = c(0, max(plot_df$mdr_saffron, na.rm = TRUE) * 1.04)) +
  sc_theme
if (RUN_ADJ_SAFFRON) {
  plot_a_singlecolumn <- plot_a_singlecolumn +
    geom_line(aes(y = mdr_saffron, color = "Adj-SAFFRON"),
              size = plot_linewidth)
}

plot_b_singlecolumn <- ggplot(plot_df, aes(x = time)) +
  geom_line(aes(y = tp_omdrc, color = "OMDRC.DD"),
            size = plot_linewidth) +
  geom_line(aes(y = tp_sct, color = "SCT"), size = plot_linewidth) +
  geom_line(aes(y = tp_rtk, color = "RTK"), size = plot_linewidth) +
  geom_hline(yintercept = total_fraud_count_in_stream, linetype = "dotted",
             color = "grey30", size = 0.9) +
  annotate("text", x = n_stream * 0.53,
           y = total_fraud_count_in_stream * 0.90,
           label = paste("Total Frauds:", total_fraud_count_in_stream),
           color = "grey30", size = SC_ANN_MAIN) +
  scale_x_continuous(breaks = sc_time_breaks, labels = sc_time_labels) +
  labs(title = "(b) True Discoveries", x = "Time (t)",
       y = "True Frauds Detected") +
  sc_theme
if (RUN_ADJ_SAFFRON) {
  plot_b_singlecolumn <- plot_b_singlecolumn +
    geom_line(aes(y = tp_saffron, color = "Adj-SAFFRON"),
              size = plot_linewidth)
}

plot_c_singlecolumn <- ggplot(plot_df, aes(x = time)) +
  geom_line(aes(y = fdr_omdrc, color = "OMDRC.DD"),
            size = plot_linewidth, na.rm = TRUE) +
  geom_line(aes(y = fdr_sct, color = "SCT"),
            size = plot_linewidth, na.rm = TRUE) +
  geom_line(aes(y = fdr_rtk, color = "RTK"),
            size = plot_linewidth, na.rm = TRUE) +
  scale_x_continuous(breaks = sc_time_breaks, labels = sc_time_labels) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1)) +
  labs(title = "(c) Empirical FDR", x = "Time (t)",
       y = "Empirical FDR") +
  sc_theme
if (RUN_ADJ_SAFFRON) {
  plot_c_singlecolumn <- plot_c_singlecolumn +
    geom_line(aes(y = fdr_saffron, color = "Adj-SAFFRON"),
              size = plot_linewidth, na.rm = TRUE)
}

plot_d_singlecolumn <- ggplot(plot_df, aes(x = time)) +
  geom_line(aes(y = rej_rate_omdrc, color = "OMDRC.DD"),
            size = plot_linewidth) +
  geom_line(aes(y = rej_rate_sct, color = "SCT"), size = plot_linewidth) +
  geom_line(aes(y = rej_rate_rtk, color = "RTK"), size = plot_linewidth) +
  annotate("text", x = n_stream * 0.025, y = REJ_RATE_YMAX * 0.965,
           hjust = 0, vjust = 1, size = SC_ANN_SMALL, color = "grey45",
           label = "initial transient clipped") +
  scale_x_continuous(breaks = sc_time_breaks, labels = sc_time_labels) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(title = "(d) Rejection Rate", x = "Time (t)",
       y = "Cumulative Rejection Rate") +
  coord_cartesian(ylim = c(0, REJ_RATE_YMAX)) +
  sc_theme
if (RUN_ADJ_SAFFRON) {
  plot_d_singlecolumn <- plot_d_singlecolumn +
    geom_line(aes(y = rej_rate_saffron, color = "Adj-SAFFRON"),
              size = plot_linewidth)
}

plot_e_singlecolumn <- ggplot(
  revision_summary_df,
  aes(x = alpha, y = max_mdr, color = method,
      shape = method, group = method)
) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey35", size = 0.8) +
  annotate("text", x = 0.085, y = 0.052, label = "MDR = target",
           hjust = 0, color = "grey40", size = SC_ANN_SMALL) +
  geom_line(size = plot_linewidth, show.legend = FALSE) +
  geom_point(size = sweep_point_size, show.legend = FALSE) +
  scale_shape_manual(values = method_shapes, guide = "none") +
  scale_x_continuous(breaks = c(0.1, 0.2, 0.3)) +
  labs(title = "(e) Maximum MDR", x = "Target MDR level",
       y = "Maximum Empirical MDR") +
  sc_theme

plot_f_singlecolumn <- ggplot(
  revision_summary_df,
  aes(x = alpha, y = flagged, color = method,
      shape = method, group = method)
) +
  geom_line(size = plot_linewidth, show.legend = FALSE) +
  geom_point(size = sweep_point_size, show.legend = FALSE) +
  scale_shape_manual(values = method_shapes, guide = "none") +
  scale_x_continuous(breaks = c(0.1, 0.2, 0.3)) +
  scale_y_log10(labels = scales::comma,
                breaks = c(100, 300, 1000, 3000, 10000)) +
  labs(title = "(f) Total Discoveries", x = "Target MDR level",
       y = "Total Discoveries (log scale)") +
  sc_theme

figure3_section5_singlecolumn <-
  ((plot_a_singlecolumn + plot_b_singlecolumn) /
   (plot_c_singlecolumn + plot_d_singlecolumn) /
   (plot_e_singlecolumn + plot_f_singlecolumn)) +
  plot_layout(guides = "collect") &
  scale_color_manual(name = "Method", values = my_colors,
                     breaks = method_levels, labels = method_labels) &
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(size = rel(0.82)),
    legend.text = element_text(size = rel(0.76)),
    legend.key.width = grid::unit(8, "mm"),
    legend.margin = margin(t = -2, r = 0, b = 0, l = 0)
  )

ggsave("Figure3_CCFD_online_comparison_revised_singlecolumn.pdf",
       plot = figure3_section5_singlecolumn,
       width = 7, height = 7.8, units = "in")
ggsave("Figure3_CCFD_online_comparison_revised_singlecolumn.png",
       plot = figure3_section5_singlecolumn,
       width = 7, height = 7.8, units = "in", dpi = 300)
cat("Saved Figure3_CCFD_online_comparison_revised_singlecolumn.pdf/png.\n")
cat("Original Figure3_CCFD_online_comparison.pdf/png were not modified.\n")


# =========================================================================
# END OF SCRIPT
# =========================================================================
