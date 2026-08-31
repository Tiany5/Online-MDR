# ==============================================================================
# NFR VARIANT of code-semi/setting1,2,3/compare_gam_dim.R
#
# ABLATION: high-dimensional density-ratio estimation for OMDRC.DD in Setting 3,
# now run under the credited anti-free-riding (NFR) barrier.
#
# The data-generating process is copied VERBATIM from compare_gam_dim.R
# (F0 = N_d(0, Sigma), F1 = N_d(mu, Sigma), Sigma_jk = 0.5^|j-k|,
# mu = (1.5 x 5, 0 x 15), multi-regime staircase pi_t with levels
# (0.35, 0.18, 0.04)). Three ways of scoring the stream are compared:
#   OMDRC.OR      - oracle: exact 20-dim linear log density ratio x'Sigma^{-1}mu.
#   DD-iForest    - reduce 20-dim -> 1-dim Isolation-Forest anomaly score, then
#                   fit a GAM density ratio on that 1-D score (ablation arm).
#   DD-GAM20      - NO dimension reduction: additive GAM classifier density ratio
#                   directly on the 20-dim features (y ~ s(V1)+...+s(V20)).
#
# The ONLY change w.r.t. compare_gam_dim.R is that ALL THREE arms now carry the
# credited local barrier at w = 100 (Algorithms 1 and 2):
#   OMDRC.OR    -> OMDRC_OR_NFR(true 20-dim q,  alpha, d)
#   DD-iForest  -> OMDRC_OR_NFR(estimated q_iso, alpha, d)
#   DD-GAM20    -> OMDRC_OR_NFR(estimated q_g20, alpha, d)
# "credited" = a barrier-forced rejection credits +alpha*q to the live ledger
# (_nfr_core.R::OMDRC_OR_NFR). OMDRC_DD is still run for both estimated arms,
# but only its estimated score path (DR / DR_ini) is used; the decisions come
# from the barrier ledger.
#
# Goal (unchanged): show DD-GAM20 recovers near-oracle FDR (the true log-ratio is
# additive/linear here, so the additive GAM is well specified), whereas
# DD-iForest pays a large FDR penalty from the lossy 20->1 compression. Under the
# barrier we additionally print the pathwise ledger check MDR_q = sum_skip q/sum q
# in the TRUE q: the oracle arm must PASS on every path, while the two estimated
# arms run their ledger on q_hat and may FAIL -- the documented cost of
# estimation error under the credited ledger.
# ==============================================================================

suppressPackageStartupMessages({
  library(foreach)
  library(parallel)
  library(doParallel)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

code_dir <- file.path(Sys.getenv("OMDRC_ROOT", unset = getwd()), "code-semi")
omdrc_file <- file.path(code_dir, "OMDRC.R")
if (!file.exists(omdrc_file)) stop("Cannot find OMDRC.R in: ", code_dir)
source(omdrc_file)

# Shared NFR core (gamma path + the barrier ledgers).
nfr_dir <- file.path(code_dir, "nofreeride")
nfr_core_file <- file.path(nfr_dir, "_nfr_core.R")
if (!file.exists(nfr_core_file)) stop("Cannot find _nfr_core.R in: ", nfr_dir)
source(nfr_core_file)

# Barrier window: d = 100, identical to nfr_setting{1,2,3}_fig.R and the probes.
D_NFR <- 100L

# ---- Setting-3 configuration (kept identical to compare_gam_dim.R) ----
ini   <- 500
m     <- seq(from = 100, to = 1000, by = 50)
n     <- 1000            # labeled reference size per class
alpha <- 0.1
reps  <- 100
D     <- 150
d            <- 20
mu_val       <- 1.5
n_signal_dims<- 5
rho          <- 0.5
n_train      <- 2000
pi_levels <- c(0.35, 0.18, 0.04)
pi_breaks <- c(0.30, 0.60)
pi_width  <- 0.03
pi_bounds <- c(0.01, 0.99)
iso_ntrees <- 100L
iso_sample <- 256L

make_pi_staircase <- function(N, levels = c(0.35, 0.18, 0.04),
                              breaks = c(0.30, 0.60), width = 0.03) {
  u <- seq(0, 1, length.out = as.integer(N))
  pi <- rep(levels[1], length(u))
  for (i in seq_along(breaks)) {
    pi <- pi - (levels[i] - levels[i + 1]) * stats::plogis((u - breaks[i]) / width)
  }
  pmin(pmax(pi, min(levels)), max(levels))
}

# Correlated-Gaussian components and exact linear log-ratio (oracle).
idx <- seq_len(d)
Sigma <- rho^abs(outer(idx, idx, "-"))
mu <- c(rep(mu_val, n_signal_dims), rep(0, d - n_signal_dims))
Sigma_chol <- chol(Sigma)
w_lin <- solve(Sigma, mu)
cc_lin <- 0.5 * sum(mu * w_lin)

N <- max(m) + ini
pi_path <- make_pi_staircase(N, levels = pi_levels, breaks = pi_breaks,
                             width = pi_width)

fig_dir <- nfr_dir
results_file <- file.path(fig_dir, "NFR_compare_gam_dim_results.rds")
reuse_cache <- identical(Sys.getenv("NFR_REUSE_CACHE"), "1")

if (reuse_cache && file.exists(results_file)) {
  cat("Reusing cached numerical results from", results_file, "\n")
  result <- readRDS(results_file)
} else {

cat("Starting NFR high-dimensional density-ratio ablation (credited barrier, d=",
    D_NFR, ")...\n", sep = "")
n_cores <- max(1L, min(parallel::detectCores(logical = TRUE), reps))
cl <- makeCluster(n_cores)
registerDoParallel(cl)
clusterExport(cl, c("m", "n", "ini", "alpha", "D", "N", "d", "mu", "Sigma_chol",
                    "w_lin", "cc_lin", "n_train", "pi_path", "pi_bounds",
                    "iso_ntrees", "iso_sample", "omdrc_file", "nfr_core_file",
                    "D_NFR"),
              envir = environment())

result <- foreach(r = seq_len(reps), .packages = "isotree") %dopar% {
  source(omdrc_file)
  source(nfr_core_file)
  set.seed(r)

  rmv <- function(nrow, mean_vec) {
    if (nrow <= 0L) return(matrix(0, 0, d))
    Z <- matrix(rnorm(nrow * d), nrow, d)
    sweep(Z %*% Sigma_chol, 2L, mean_vec, "+")
  }

  # ---- Evaluation stream (d-dimensional, time-varying pi) ----
  theta <- rbinom(N, size = 1, prob = pi_path)
  X <- matrix(0, N, d)
  if (any(theta == 0)) X[theta == 0, ] <- rmv(sum(theta == 0), rep(0, d))
  if (any(theta == 1)) X[theta == 1, ] <- rmv(sum(theta == 1), mu)

  # Exact 20-dim oracle Lmdr via the linear log-ratio.
  R_all <- exp(as.numeric(X %*% w_lin) - cc_lin)
  lmdr_all <- (pi_path * R_all) / ((1 - pi_path) + pi_path * R_all)

  # Labeled references + pure-null iForest training set.
  X0_ref <- rmv(n, rep(0, d))
  X1_ref <- rmv(n, mu)
  X_train <- rmv(n_train, rep(0, d))

  # Split into initial calibration batch and stream.
  X_ini <- X[seq_len(ini), , drop = FALSE]
  X_str <- X[(ini + 1):N, , drop = FALSE]
  lmdr  <- lmdr_all[(ini + 1):N]
  lmdr_ini <- lmdr_all[seq_len(ini)]   # fills the gamma window during warm-up
  theta_s <- theta[(ini + 1):N]

  # ---- (1) Oracle + credited barrier ----
  d_or <- OMDRC_OR_NFR(lmdr, alpha, d = D_NFR, x.Lmdr_ini = lmdr_ini)

  # ---- (2) DD-iForest: 20-dim -> 1-dim score -> GAM ratio -> barrier ----
  iso <- isotree::isolation.forest(
    as.data.frame(X_train), ntrees = iso_ntrees,
    sample_size = min(iso_sample, n_train), nthreads = 1L, seed = r
  )
  score <- function(mat) as.numeric(predict(iso, as.data.frame(mat)))
  s_ini  <- score(X_ini)
  s_str  <- score(X_str)
  s0_ref <- score(X0_ref)
  s1_ref <- score(X1_ref)
  d_iso <- OMDRC_DD(z = s_str, z_ini = s_ini, z0 = s0_ref, z1 = s1_ref,
                    alpha = alpha, ratio_method = "gam",
                    pi_bounds = pi_bounds, em_tol = 1e-6, em_max_iter = 50L,
                    D_mode = "growing", D_beta = 0.6, D_min = 10L)
  d_iso_nfr <- OMDRC_OR_NFR(d_iso$DR, alpha, d = D_NFR,
                            x.Lmdr_ini = d_iso$DR_ini)

  # ---- (3) DD-GAM20: additive GAM ratio on the 20-dim features -> barrier ----
  d_g20 <- OMDRC_DD(z = X_str, z_ini = X_ini, z0 = X0_ref, z1 = X1_ref,
                    alpha = alpha, ratio_method = "gam",
                    pi_bounds = pi_bounds, em_tol = 1e-6, em_max_iter = 50L,
                    D_mode = "growing", D_beta = 0.6, D_min = 10L)
  d_g20_nfr <- OMDRC_OR_NFR(d_g20$DR, alpha, d = D_NFR,
                            x.Lmdr_ini = d_g20$DR_ini)

  expected_signals <- cumsum(pi_path[(ini + 1):N])[m]
  metrics <- function(de) {
    cmiss <- cumsum(theta_s * (1 - de)); cfp <- cumsum((1 - theta_s) * de)
    cde <- cumsum(de)
    list(mdr = cmiss[m] / pmax(expected_signals, .Machine$double.eps),
         fdr = cfp[m] / pmax(cde[m], 1))
  }
  # Ledger check: MDR_q = sum_skip q / sum q in the TRUE q. For the ORACLE arm
  # the capacity identity pins this at <= alpha on EVERY path; the two estimated
  # arms run their ledger on q_hat, so the identity is broken and this is exactly
  # where the estimation error shows up.
  mdrq <- function(de) sum(lmdr * (1 - de)) / sum(lmdr)

  list(or = metrics(d_or$de), iso = metrics(d_iso_nfr$de),
       g20 = metrics(d_g20_nfr$de),
       mdrq = c(or = mdrq(d_or$de), iso = mdrq(d_iso_nfr$de),
                g20 = mdrq(d_g20_nfr$de)),
       blocked = c(or = sum(d_or$blocked), iso = sum(d_iso_nfr$blocked),
                   g20 = sum(d_g20_nfr$blocked)))
}

stopCluster(cl); registerDoSEQ()
saveRDS(result, results_file)
cat("Saved numerical results to", results_file, "\n")
}

# ---- Aggregate terminal (t = 1000) metrics ----
idxT <- length(m)
agg <- function(key, metric) {
  vals <- sapply(result, function(x) x[[key]][[metric]][idxT])
  c(mean = mean(vals, na.rm = TRUE),
    se = sd(vals, na.rm = TRUE) / sqrt(sum(is.finite(vals))))
}
tab <- data.frame(
  Method = c("OMDRC.OR+NFR (20-dim exact)",
             "DD-iForest+NFR (1-dim score + GAM)",
             "DD-GAM20+NFR (direct 20-dim additive GAM)"),
  MDR = c(agg("or", "mdr")["mean"], agg("iso", "mdr")["mean"], agg("g20", "mdr")["mean"]),
  FDR = c(agg("or", "fdr")["mean"], agg("iso", "fdr")["mean"], agg("g20", "fdr")["mean"]),
  MDR_se = c(agg("or", "mdr")["se"], agg("iso", "mdr")["se"], agg("g20", "mdr")["se"]),
  FDR_se = c(agg("or", "fdr")["se"], agg("iso", "fdr")["se"], agg("g20", "fdr")["se"])
)

cat(sprintf("\n=== Setting 3 density-ratio ablation UNDER NFR (t=%d, alpha=%.2f, reps=%d, d=%d, credited) ===\n",
            m[idxT], alpha, length(result), D_NFR))
cat(sprintf("%-44s %8s %8s\n", "Method", "MDR", "FDR"))
for (i in seq_len(nrow(tab))) {
  tag <- if (tab$MDR[i] > alpha + 1e-9) " [MDR>alpha]" else ""
  cat(sprintf("%-44s %6.3f(%.3f) %6.3f(%.3f)%s\n",
              tab$Method[i], tab$MDR[i], tab$MDR_se[i],
              tab$FDR[i], tab$FDR_se[i], tag))
}
cat("\nFDR gap DD-iForest+NFR - OMDRC.OR+NFR : ",
    sprintf("%.3f", tab$FDR[2] - tab$FDR[1]), "\n")
cat("FDR gap DD-GAM20+NFR   - OMDRC.OR+NFR : ",
    sprintf("%.3f", tab$FDR[3] - tab$FDR[1]), "\n")

saveRDS(tab, file.path(fig_dir, "NFR_compare_gam_dim_summary.rds"))

# Hard ledger check (the point of the NFR section): MDR_q in the TRUE q must be
# <= alpha on EVERY path for the oracle arm. The estimated arms run their ledger
# on q_hat so they can fail -- that failure is the cost of estimation error.
.diag_ledger <- function(res, alpha, d_nfr) {
  cat(sprintf("\n--- Barrier ledger check (d=%d, credited) ---\n", d_nfr))
  cat(sprintf("%-18s %10s %10s %12s %10s\n",
              "Arm", "mean MDR_q", "max MDR_q", "ledger", "blocked"))
  labs <- c(or = "OMDRC.OR+NFR", iso = "DD-iForest+NFR", g20 = "DD-GAM20+NFR")
  for (a in names(labs)) {
    v <- vapply(res, function(x) as.numeric(x$mdrq[[a]]), numeric(1))
    b <- vapply(res, function(x) as.numeric(x$blocked[[a]]), numeric(1))
    cat(sprintf("%-18s %10.4f %10.4f %12s %10.1f\n",
                labs[[a]], mean(v), max(v),
                if (max(v) <= alpha + 1e-9) "PASS" else "FAIL", mean(b)))
  }
}
.diag_ledger(result, alpha, D_NFR)

# ---- Timeline figure (MDR / FDR over the full time grid) ----
# Plotted legend labels are deliberately plain (no "+NFR" suffix), matching
# nfr_setting{1,2,3}_fig.R; the console tables above keep the explicit suffix.
method_keys   <- c("or", "iso", "g20")
method_labels <- c("OMDRC.OR", "OMDRC.DD (iForest)", "OMDRC.DD (PC-DRE)")
my_colors <- c("OMDRC.OR" = "#984ea3", "OMDRC.DD (iForest)" = "#E69F00",
               "OMDRC.DD (PC-DRE)" = "#009E73")
my_shapes <- c("OMDRC.OR" = 16, "OMDRC.DD (iForest)" = 17, "OMDRC.DD (PC-DRE)" = 15)

build_df <- function(metric) {
  do.call(rbind, lapply(seq_along(method_keys), function(j) {
    mat <- t(sapply(result, function(x) x[[method_keys[j]]][[metric]]))  # reps x |m|
    est <- colMeans(mat, na.rm = TRUE)
    se  <- apply(mat, 2, sd, na.rm = TRUE) / sqrt(nrow(mat))
    data.frame(t = m, value = est,
               lower = pmax(0, est - 1.96 * se),
               upper = pmin(1, est + 1.96 * se),
               type = factor(method_labels[j], levels = method_labels))
  }))
}
df_mdr <- build_df("mdr"); df_fdr <- build_df("fdr")

thm <- theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 18),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 13),
        panel.grid.minor = element_blank())

g_mdr <- ggplot(df_mdr, aes(t, value, color = type, shape = type, group = type)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = type),
              alpha = 0.16, colour = NA, show.legend = FALSE) +
  geom_line(size = 0.8) + geom_point(size = 1.8) +
  geom_hline(yintercept = alpha, linetype = "dashed", colour = "black") +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  scale_shape_manual(values = my_shapes) +
  labs(x = "Time (t)", y = "MDR") + thm

g_fdr <- ggplot(df_fdr, aes(t, value, color = type, shape = type, group = type)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = type),
              alpha = 0.16, colour = NA, show.legend = FALSE) +
  geom_line(size = 0.8) + geom_point(size = 1.8) +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  scale_shape_manual(values = my_shapes) +
  labs(x = "Time (t)", y = "FDR") + thm

g <- (g_mdr + g_fdr) + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(file.path(fig_dir, "compare_gam_dim_timeline.pdf"), g,
       width = 10, height = 5)
ggsave(file.path(fig_dir, "compare_gam_dim_timeline.png"), g,
       width = 10, height = 5, dpi = 300)
cat("Saved figure: compare_gam_dim_timeline.{pdf,png}\n")
