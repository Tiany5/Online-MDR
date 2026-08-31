# =============================================================================
# Figure 2: Oracle mechanism analysis under fixed and Gaussian-surge pi_t.
# =============================================================================
# The manuscript's OMDRC.OR uses the local barrier; OMDRC.OR_nob is the
# capacity-only comparison obtained by removing it. FT is the offline fixed-
# threshold oracle. The shared helper calls the barrier window `d`; this is the
# manuscript parameter w and is set to 100 below.
#
# The barrier adds, on top of the capacity ledger, a causal local FT-type
# threshold gamma_t computed on a trailing window of length d: a skip is only
# allowed when BOTH (1-alpha) q_t <= C_t (capacity) AND q_t < gamma_t (barrier).
# A "credited" barrier-forced rejection credits +alpha*q_t to the live ledger,
# so C_t >= 0 is preserved and the capacity identity (hence MDR_q <= alpha on
# every path) still goes through for the ORACLE arm.
#
# Figure contents:
#   * capacity panels (a.1/b.1): the capacity ledger C_t is the barrier ledger;
#     barrier line (1-alpha) gamma_t is overlaid, and the free-riding attempts
#     the barrier converts into rejections are marked as points on the cost
#     curve;
#   * MDR/FDR panels (a.2--b.3): OMDRC.OR, OMDRC.OR_nob, and FT are compared.
#
# Main output:
#   Figure2_oracle_mechanism_setting2_surge_pi_NFR.pdf
#   Figure2_oracle_mechanism_setting2_surge_pi_NFR.png
# =============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# Shared NFR core (gamma path + the credited barrier OMDRC_OR_NFR).
code_dir <- file.path(Sys.getenv("OMDRC_ROOT", unset = getwd()), "code-semi")
nfr_core_file <- file.path(code_dir, "nofreeride", "_nfr_core.R")
if (!file.exists(nfr_core_file)) stop("Cannot find _nfr_core.R in: ",
                                      dirname(nfr_core_file))
source(nfr_core_file)

# -----------------------------------------------------------------------------
# 1. User-adjustable settings
# -----------------------------------------------------------------------------

SIM_PARAMS <- list(
  alpha        = 0.10,
  n_time       = 1000L,
  n_reps       = 200L,
  # Realization shown in the capacity panels (see Figure2.R for the rationale).
  seed_path    = 1826L,
  gs_mean      = 3,      # Setting 2: F1 = N(gs_mean, 1), F0 = N(0, 1)
  pi_fixed     = 0.05,   # fixed regime uses the surge baseline pi_low
  pi_low       = 0.05,   # Setting 2 surge path parameters
  pi_high      = 0.30,
  surge_center = 0.60,
  surge_width  = 0.15,
  # NFR barrier window (trailing length d for gamma_t). Matches the credited
  # barrier used in nofreeride/nfr_setting2_fig.R and application/CCFD.
  d_nfr        = 100L
)

SAVE_FIGURE <- TRUE
OUTPUT_STEM <- "Figure2_oracle_mechanism_setting2_surge_pi_NFR"

# Cache of the Monte Carlo MDR curves. Set REUSE_CACHE <- TRUE (or export
# NFR_REUSE_CACHE=1) to reload and re-plot only.
RESULTS_FILE <- paste0(OUTPUT_STEM, "_results.rds")
REUSE_CACHE <- identical(Sys.getenv("NFR_REUSE_CACHE"), "1")

# -----------------------------------------------------------------------------
# 2. Core oracle procedures
# -----------------------------------------------------------------------------

# Capacity-only oracle used for the OMDRC.OR_nob comparison.
OMDRC_OR_capacity <- function(lmdr, alpha) {
  lmdr <- pmin(pmax(as.numeric(lmdr), 0), 1)
  n <- length(lmdr)

  decision <- integer(n)
  capacity_before <- numeric(n)
  capacity_after <- numeric(n)
  capacity <- 0

  for (t in seq_len(n)) {
    capacity_before[t] <- capacity
    local_cost <- (1 - alpha) * lmdr[t]

    if (local_cost <= capacity) {
      decision[t] <- 0L
      capacity <- capacity - local_cost
    } else {
      decision[t] <- 1L
      capacity <- capacity + alpha * lmdr[t]
    }

    if (capacity < 0 && capacity > -1e-12) {
      capacity <- 0
    }
    capacity_after[t] <- capacity
  }

  list(
    de = decision,
    capacity_before = capacity_before,
    capacity_after = capacity_after,
    lmdr = lmdr
  )
}

OMDRC_OR <- function(lmdr, alpha) {
  out <- OMDRC_OR_capacity(lmdr, alpha)
  list(de = out$de)
}

# NFR barrier oracle: capacity ledger PLUS the trailing-window gamma_t barrier.
# Returns the same fields as OMDRC_OR_capacity (capacity_before) so the capacity
# panel can reuse the plotting code, plus the barrier line and blocked flags.
OMDRC_OR_NFR_capacity <- function(lmdr, alpha, d) {
  out <- OMDRC_OR_NFR(lmdr, alpha, d = d)   # from _nfr_core.R
  m <- length(out$de)
  list(
    de = out$de,
    capacity_before = out$capacity[seq_len(m)],   # capacity BEFORE decision t
    gamma = out$gamma,                            # barrier threshold on q_t
    blocked = out$blocked,                        # free-riding attempts blocked
    lmdr = pmin(pmax(as.numeric(lmdr), 0), 1)
  )
}

# FT: a noncausal offline oracle. It observes the full realized Lmdr sequence
# and selects one fixed threshold to satisfy the terminal MDR constraint.
OMDRC_OFF <- function(lmdr, alpha) {
  x <- pmin(pmax(as.numeric(lmdr), 0), 1)
  n <- length(x)

  if (n == 0L) {
    stop("lmdr must be non-empty.")
  }

  total_mass <- sum(x)
  if (!is.finite(total_mass) || total_mass <= .Machine$double.eps) {
    return(list(lambda = Inf, de = integer(n)))
  }

  xs <- sort(x, decreasing = FALSE)
  missed_fraction <- cumsum(xs) / total_mass
  feasible <- which(missed_fraction <= alpha)

  if (length(feasible) == 0L) {
    lambda <- -Inf
    decision <- rep(1L, n)
  } else {
    k <- max(feasible)
    lambda <- xs[k]
    decision <- as.integer(x > lambda)
  }

  list(lambda = lambda, de = decision)
}

# -----------------------------------------------------------------------------
# 3. Data-generating mechanism (identical to Figure2.R)
# -----------------------------------------------------------------------------

make_pi_surge <- function(n_time, pi_low, pi_high, center, width) {
  if (width <= 0) {
    stop("surge_width must be positive.")
  }
  if (pi_high < pi_low) {
    stop("pi_high must be at least pi_low.")
  }

  u <- seq(0, 1, length.out = as.integer(n_time))
  pi_t <- pi_low + (pi_high - pi_low) *
    exp(-((u - center)^2) / (2 * width^2))
  pmin(pmax(pi_t, pi_low), pi_high)
}

simulate_oracle_stream <- function(pi_t, seed, gs_mean) {
  n <- length(pi_t)
  set.seed(seed)

  theta <- rbinom(n, size = 1, prob = pi_t)
  z0 <- rnorm(n, mean = 0, sd = 1)
  z1 <- rnorm(n, mean = gs_mean, sd = 1)

  z <- ifelse(theta == 1L, z1, z0)

  f0 <- dnorm(z, mean = 0, sd = 1)
  f1 <- dnorm(z, mean = gs_mean, sd = 1)

  denominator <- (1 - pi_t) * f0 + pi_t * f1
  lmdr <- pi_t * f1 / pmax(denominator, .Machine$double.xmin)
  lmdr <- pmin(pmax(lmdr, 0), 1)

  list(theta = theta, z = z, lmdr = lmdr, pi_t = pi_t)
}

# -----------------------------------------------------------------------------
# 4. One-realization capacity panels (NFR barrier)
# -----------------------------------------------------------------------------

create_capacity_data <- function(pi_t, alpha, seed, gs_mean, d) {
  stream <- simulate_oracle_stream(
    pi_t = pi_t,
    seed = seed,
    gs_mean = gs_mean
  )

  online <- OMDRC_OR_NFR_capacity(stream$lmdr, alpha, d = d)
  plain <- OMDRC_OR_capacity(stream$lmdr, alpha)   # capacity-only (no barrier)
  offline <- OMDRC_OFF(stream$lmdr, alpha)

  ft_cost_threshold <- if (is.finite(offline$lambda)) {
    (1 - alpha) * offline$lambda
  } else {
    NA_real_
  }

  # Put the barrier on the SAME (cost) axis as Cost and Capacity: the skip gate
  # q_t < gamma_t is equivalent to (1-alpha) q_t < (1-alpha) gamma_t. Inactive
  # barrier (gamma_t = Inf) becomes NA so the line simply breaks there.
  barrier_cost_raw <- (1 - alpha) * online$gamma          # Inf where inactive
  barrier_cost <- barrier_cost_raw
  barrier_cost[!is.finite(barrier_cost)] <- NA_real_

  # Effective skip (acceptance) threshold for the NFR rule. A point is
  # NON-rejected iff BOTH (1-alpha) q_t <= C_t and (1-alpha) q_t < (1-alpha)
  # gamma_t, i.e. iff its cost falls below min(C_t, (1-alpha) gamma_t). When the
  # barrier is inactive (gamma_t = Inf) the minimum reduces to the capacity C_t.
  eff_skip <- pmin(online$capacity_before, barrier_cost_raw)

  plot_data <- data.frame(
    t = seq_along(pi_t),
    Capacity = online$capacity_before,
    CapacityPlain = plain$capacity_before,
    Cost = (1 - alpha) * online$lmdr,
    Barrier = barrier_cost,
    EffSkip = eff_skip
  )

  # Free-riding attempts the barrier converted into rejections, drawn on the
  # cost curve.
  blocked_idx <- which(online$blocked == 1L)
  blocked_data <- if (length(blocked_idx) > 0L) {
    data.frame(t = blocked_idx,
               Cost = (1 - alpha) * online$lmdr[blocked_idx])
  } else {
    data.frame(t = numeric(0), Cost = numeric(0))
  }

  list(
    plot_data = plot_data,
    blocked_data = blocked_data,
    ft_cost_threshold = ft_cost_threshold,
    stream = stream,
    online = online,
    offline = offline
  )
}

# -----------------------------------------------------------------------------
# 5. Monte Carlo MDR panels
# -----------------------------------------------------------------------------

# The manuscript defines MDR as a ratio of expectations:
#   sum_r missed_{r,t} / sum_r signals_{r,t}.
# Three arms: capacity-only OMDRC.OR_nob, barrier OMDRC.OR, and offline FT.
# A per-path ledger check on the TRUE q is accumulated
# for the barrier arm: MDR_q = sum(q (1 - de)) / sum(q) must be <= alpha on
# every path for the credited oracle barrier.
run_nfr_vs_ft_sim <- function(pi_t, alpha, n_reps, gs_mean, d,
                              progress_every = 100L) {
  n_time <- length(pi_t)

  missed_or_sum <- numeric(n_time)
  missed_nfr_sum <- numeric(n_time)
  missed_ft_sum <- numeric(n_time)
  signal_sum <- numeric(n_time)

  # FDR is accumulated as the pathwise empirical proportion V_t/max(R_t,1),
  # averaged over reps (the conventional E[V/R] estimator), in contrast to the
  # ratio-of-expectations MDR estimator above.
  fdr_or_sum <- numeric(n_time)
  fdr_nfr_sum <- numeric(n_time)
  fdr_ft_sum <- numeric(n_time)

  mdrq_nfr <- numeric(n_reps)   # terminal MDR_q (true q) for the barrier arm
  mdrq_or <- numeric(n_reps)
  blocked_nfr <- numeric(n_reps)

  for (r in seq_len(n_reps)) {
    if (progress_every > 0L && r %% progress_every == 0L) {
      message("  repetition ", r, " / ", n_reps)
    }

    stream <- simulate_oracle_stream(
      pi_t = pi_t,
      seed = r,
      gs_mean = gs_mean
    )

    de_or <- OMDRC_OR(stream$lmdr, alpha)$de
    nfr <- OMDRC_OR_NFR(stream$lmdr, alpha, d = d)
    de_nfr <- nfr$de
    de_ft <- OMDRC_OFF(stream$lmdr, alpha)$de

    cumulative_signals <- cumsum(stream$theta)
    signal_sum <- signal_sum + cumulative_signals
    missed_or_sum <- missed_or_sum + cumsum(stream$theta * (1 - de_or))
    missed_nfr_sum <- missed_nfr_sum + cumsum(stream$theta * (1 - de_nfr))
    missed_ft_sum <- missed_ft_sum + cumsum(stream$theta * (1 - de_ft))

    false_pos <- 1 - stream$theta   # a rejection of a null observation
    rej_or <- pmax(cumsum(de_or), 1)
    rej_nfr <- pmax(cumsum(de_nfr), 1)
    rej_ft <- pmax(cumsum(de_ft), 1)
    fdr_or_sum <- fdr_or_sum + cumsum(false_pos * de_or) / rej_or
    fdr_nfr_sum <- fdr_nfr_sum + cumsum(false_pos * de_nfr) / rej_nfr
    fdr_ft_sum <- fdr_ft_sum + cumsum(false_pos * de_ft) / rej_ft

    q_tot <- sum(stream$lmdr)
    if (q_tot > .Machine$double.eps) {
      mdrq_nfr[r] <- sum(stream$lmdr * (1 - de_nfr)) / q_tot
      mdrq_or[r] <- sum(stream$lmdr * (1 - de_or)) / q_tot
    }
    blocked_nfr[r] <- sum(nfr$blocked)
  }

  denominator <- pmax(signal_sum, 1)

  list(
    mdr = data.frame(
      t = seq_len(n_time),
      OMDRC.OR = missed_or_sum / denominator,
      OMDRC.OR.NFR = missed_nfr_sum / denominator,
      FT = missed_ft_sum / denominator
    ),
    fdr = data.frame(
      t = seq_len(n_time),
      OMDRC.OR = fdr_or_sum / n_reps,
      OMDRC.OR.NFR = fdr_nfr_sum / n_reps,
      FT = fdr_ft_sum / n_reps
    ),
    ledger = list(
      mdrq_nfr = mdrq_nfr, mdrq_or = mdrq_or, blocked_nfr = blocked_nfr
    )
  )
}

# -----------------------------------------------------------------------------
# 6. Shared plotting style
# -----------------------------------------------------------------------------

COLORS_CAPACITY <- c(
  "EffSkip" = "#D55E00",   # matches OMDRC.OR+NFR in the MDR panels
  "CapPlain" = "#466300",  # matches OMDRC.OR (plain, no barrier)
  "Cost" = "#F8766D",
  "FT" = "#00BFC4"         # matches FT
)

LINE_TYPES_CAPACITY <- c(
  "EffSkip" = "solid",
  "CapPlain" = "solid",
  "Cost" = "solid",
  "FT" = "dashed"
)

# The two ledgers come from DIFFERENT decision paths, so they need different
# symbols: C_t is the barrier arm's own capacity (the one entering the min), and
# C_t^nob is the capacity-only arm's ledger. The barrier (1-alpha) gamma_t is NOT
# drawn separately: it is the binding arm of the min at ~96% of the steps, so its
# line would sit exactly under the thick EffSkip curve.
LEGEND_LABELS_CAPACITY <- expression(
  paste(plain(min), group("(", list(C[t], (1 - alpha) * gamma[t]), ")"),
        ": OMDRC.OR Threshold"),
  C[t]^{"nob"} * ": OMDRC.OR_nob Threshold",
  (1 - alpha) * Lmdr[t] * ": Local Cost",
  (1 - alpha) * lambda * ": FT Threshold"
)
names(LEGEND_LABELS_CAPACITY) <- c("EffSkip", "CapPlain", "Cost", "FT")

COLORS_MDR <- c(
  "OMDRC.OR" = "#466300",
  "OMDRC.OR.NFR" = "#D55E00",
  "FT" = "#00BFC4"
)

SHAPES_MDR <- c(
  "OMDRC.OR" = 16,
  "OMDRC.OR.NFR" = 17,
  "FT" = 15
)

# Display names: the barrier oracle is OMDRC.OR in the paper, and the plain
# capacity-only arm is relabelled OMDRC.OR_nob ("no barrier").
LABELS_MDR <- c(
  "OMDRC.OR" = "OMDRC.OR_nob",
  "OMDRC.OR.NFR" = "OMDRC.OR",
  "FT" = "FT"
)

CUSTOM_THEME <- theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    plot.subtitle = element_text(
      size = 17, hjust = 0.5, margin = margin(b = 5)
    ),
    legend.title = element_blank(),
    legend.text = element_text(size = 17),
    axis.text = element_text(size = 16, colour = "black"),
    axis.title = element_text(size = 17),
    panel.grid.major = element_line(colour = "grey92", size = 0.4),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.8)
  )

create_capacity_plot <- function(capacity_result, subtitle_text) {
  # The panel now centres on the effective skip (acceptance) threshold
  # min(C_t, (1-alpha) gamma_t): its shaded region below is where an observation
  # is NON-rejected; a Local Cost spike poking above it is rejected. The plain
  # capacity-only ledger C_t is kept as the no-barrier reference -- during the
  # surge it balloons and would free-ride high-Lmdr signals into skips, whereas
  # the NFR skip threshold stays capped by the barrier. Ledgers are cumulative
  # and not bounded by 1, so use a panel-specific honest upper limit.
  finite_values <- unlist(capacity_result$plot_data[
    c("EffSkip", "CapacityPlain", "Cost")
  ], use.names = FALSE)
  finite_values <- finite_values[is.finite(finite_values)]
  y_upper <- max(1, 1.05 * max(finite_values),
                 1.05 * capacity_result$ft_cost_threshold, na.rm = TRUE)

  p <- ggplot(capacity_result$plot_data, aes(x = t)) +
    # Acceptance (skip) region: cost below min(C_t, (1-alpha) gamma_t).
    geom_ribbon(aes(ymin = 0, ymax = EffSkip),
                fill = COLORS_CAPACITY[["EffSkip"]], alpha = 0.10) +
    geom_line(
      aes(y = Cost, colour = "Cost", linetype = "Cost"),
      alpha = 0.40,
      size = 0.30
    ) +
    # Plain capacity-only skip threshold (no barrier): free-riding reference,
    # kept faint so it does not clutter the panel.
    geom_line(
      aes(y = CapacityPlain, colour = "CapPlain", linetype = "CapPlain"),
      size = 0.50,
      alpha = 0.65
    ) +
    # NFR effective skip threshold = min(capacity, (1-alpha) barrier).
    geom_line(
      aes(y = EffSkip, colour = "EffSkip", linetype = "EffSkip"),
      size = 0.95
    )

  if (is.finite(capacity_result$ft_cost_threshold)) {
    ft_line <- data.frame(
      threshold = capacity_result$ft_cost_threshold,
      key = "FT"
    )
    p <- p + geom_hline(
      data = ft_line,
      aes(yintercept = threshold, colour = key, linetype = key),
      inherit.aes = FALSE,
      size = 0.80
    )
  }

  p +
    scale_colour_manual(
      values = COLORS_CAPACITY,
      breaks = c("EffSkip", "CapPlain", "Cost", "FT"),
      labels = LEGEND_LABELS_CAPACITY
    ) +
    scale_linetype_manual(
      values = LINE_TYPES_CAPACITY,
      breaks = c("EffSkip", "CapPlain", "Cost", "FT"),
      labels = LEGEND_LABELS_CAPACITY
    ) +
    labs(subtitle = subtitle_text, x = "Time (t)", y = "Value") +
    coord_cartesian(ylim = c(0, y_upper)) +
    guides(
      linetype = "none",
      colour = guide_legend(
        ncol = 1,
        override.aes = list(
          linetype = unname(LINE_TYPES_CAPACITY[c("EffSkip", "CapPlain",
                                                  "Cost", "FT")]),
          size = c(0.95, 0.55, 0.70, 0.80),
          alpha = 1
        )
      )
    ) +
    CUSTOM_THEME
}

create_mdr_plot <- function(mdr_wide, subtitle_text, alpha, y_upper) {
  method_levels <- c("FT", "OMDRC.OR", "OMDRC.OR.NFR")
  mdr_long <- mdr_wide %>%
    pivot_longer(
      cols = c("OMDRC.OR", "OMDRC.OR.NFR", "FT"),
      names_to = "Method",
      values_to = "MDR"
    ) %>%
    mutate(Method = factor(Method, levels = method_levels))

  point_indices <- seq(from = 50L, to = max(mdr_long$t), by = 50L)
  point_data <- mdr_long %>% filter(t %in% point_indices)

  ggplot(mdr_long, aes(x = t, y = MDR, colour = Method, shape = Method)) +
    geom_line(size = 0.60) +
    geom_point(data = point_data, size = 1.80) +
    geom_hline(
      yintercept = alpha,
      linetype = "dashed",
      colour = "black",
      size = 0.70
    ) +
    scale_colour_manual(values = COLORS_MDR, labels = LABELS_MDR) +
    scale_shape_manual(values = SHAPES_MDR, labels = LABELS_MDR) +
    labs(subtitle = subtitle_text, x = "Time (t)", y = "MDR") +
    coord_cartesian(ylim = c(0, y_upper)) +
    CUSTOM_THEME
}

# The first FDR_WARMUP steps are a cold start (a handful of decisions, FDR near
# 1) and are excluded BOTH from the y-range and from the drawn curves, so no line
# runs off the top of the panel.
FDR_WARMUP <- 20L

create_fdr_plot <- function(fdr_wide, subtitle_text, y_lower, y_upper) {
  method_levels <- c("FT", "OMDRC.OR", "OMDRC.OR.NFR")
  fdr_long <- fdr_wide %>%
    filter(t >= FDR_WARMUP) %>%
    pivot_longer(
      cols = c("OMDRC.OR", "OMDRC.OR.NFR", "FT"),
      names_to = "Method",
      values_to = "FDR"
    ) %>%
    mutate(Method = factor(Method, levels = method_levels))

  point_indices <- seq(from = 50L, to = max(fdr_long$t), by = 50L)
  point_data <- fdr_long %>% filter(t %in% point_indices)

  ggplot(fdr_long, aes(x = t, y = FDR, colour = Method, shape = Method)) +
    geom_line(size = 0.60) +
    geom_point(data = point_data, size = 1.80) +
    scale_colour_manual(values = COLORS_MDR, labels = LABELS_MDR) +
    scale_shape_manual(values = SHAPES_MDR, labels = LABELS_MDR) +
    labs(subtitle = subtitle_text, x = "Time (t)", y = "FDR") +
    coord_cartesian(ylim = c(y_lower, y_upper)) +
    CUSTOM_THEME
}

# -----------------------------------------------------------------------------
# 7. Construct fixed and surge prevalence settings (Setting 2 DGP)
# -----------------------------------------------------------------------------

pi_fixed <- rep(SIM_PARAMS$pi_fixed, SIM_PARAMS$n_time)

pi_smooth <- make_pi_surge(
  n_time = SIM_PARAMS$n_time,
  pi_low = SIM_PARAMS$pi_low,
  pi_high = SIM_PARAMS$pi_high,
  center = SIM_PARAMS$surge_center,
  width = SIM_PARAMS$surge_width
)

message(sprintf(
  "Setting 2 surge pi path: min %.3f, max %.3f (NFR barrier d = %d)",
  min(pi_smooth), max(pi_smooth), SIM_PARAMS$d_nfr
))

# -----------------------------------------------------------------------------
# 8. Capacity plots from one representative realization
# -----------------------------------------------------------------------------

capacity_fixed <- create_capacity_data(
  pi_t = pi_fixed,
  alpha = SIM_PARAMS$alpha,
  seed = SIM_PARAMS$seed_path,
  gs_mean = SIM_PARAMS$gs_mean,
  d = SIM_PARAMS$d_nfr
)

capacity_smooth <- create_capacity_data(
  pi_t = pi_smooth,
  alpha = SIM_PARAMS$alpha,
  seed = SIM_PARAMS$seed_path,
  gs_mean = SIM_PARAMS$gs_mean,
  d = SIM_PARAMS$d_nfr
)

p_capacity_fixed <- create_capacity_plot(
  capacity_fixed,
  subtitle_text = bquote(
    "(a.1) Fixed signal proportion: " * pi[t] == .(SIM_PARAMS$pi_fixed)
  )
)

p_capacity_smooth <- create_capacity_plot(
  capacity_smooth,
  subtitle_text = bquote(
    "(b.1) Surge signal proportion: " *
      pi[t] ~ "rises from" ~ .(SIM_PARAMS$pi_low) ~ "to" ~
      .(SIM_PARAMS$pi_high) ~ "and returns"
  )
)

capacity_column <- (p_capacity_fixed / p_capacity_smooth) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

# -----------------------------------------------------------------------------
# 9. Monte Carlo MDR plots
# -----------------------------------------------------------------------------

MC_PARAMS <- SIM_PARAMS[setdiff(names(SIM_PARAMS), "seed_path")]

message("Running fixed-prevalence oracle comparison (NFR)...")
cached <- if (REUSE_CACHE && file.exists(RESULTS_FILE)) {
  readRDS(RESULTS_FILE)
} else {
  NULL
}

if (!is.null(cached) && !identical(cached$params, MC_PARAMS)) {
  message("Cached results were produced with different SIM_PARAMS - rerunning.")
  cached <- NULL
}

# Older caches predate the FDR panels; force a rerun if the FDR curves are
# missing so the 2x3 figure can be built.
if (!is.null(cached) && (is.null(cached$fdr_fixed) || is.null(cached$fdr_smooth))) {
  message("Cached results lack FDR curves - rerunning.")
  cached <- NULL
}

if (!is.null(cached)) {
  message("Reusing cached Monte Carlo results from ", RESULTS_FILE)
  mdr_fixed <- cached$mdr_fixed
  mdr_smooth <- cached$mdr_smooth
  fdr_fixed <- cached$fdr_fixed
  fdr_smooth <- cached$fdr_smooth
  ledger_fixed <- cached$ledger_fixed
  ledger_smooth <- cached$ledger_smooth
} else {
  sim_fixed <- run_nfr_vs_ft_sim(
    pi_t = pi_fixed,
    alpha = SIM_PARAMS$alpha,
    n_reps = SIM_PARAMS$n_reps,
    gs_mean = SIM_PARAMS$gs_mean,
    d = SIM_PARAMS$d_nfr
  )
  mdr_fixed <- sim_fixed$mdr
  fdr_fixed <- sim_fixed$fdr
  ledger_fixed <- sim_fixed$ledger

  message("Running surge-prevalence oracle comparison (NFR)...")
  sim_smooth <- run_nfr_vs_ft_sim(
    pi_t = pi_smooth,
    alpha = SIM_PARAMS$alpha,
    n_reps = SIM_PARAMS$n_reps,
    gs_mean = SIM_PARAMS$gs_mean,
    d = SIM_PARAMS$d_nfr
  )
  mdr_smooth <- sim_smooth$mdr
  fdr_smooth <- sim_smooth$fdr
  ledger_smooth <- sim_smooth$ledger

  saveRDS(
    list(mdr_fixed = mdr_fixed, mdr_smooth = mdr_smooth,
         fdr_fixed = fdr_fixed, fdr_smooth = fdr_smooth,
         ledger_fixed = ledger_fixed, ledger_smooth = ledger_smooth,
         params = MC_PARAMS),
    RESULTS_FILE
  )
  message("Saved Monte Carlo results to ", RESULTS_FILE)
}

# Ledger PASS/FAIL check on the TRUE q for the credited oracle barrier: the
# capacity identity pins MDR_q <= alpha on EVERY path.
.report_ledger <- function(tag, ledger, alpha) {
  v <- ledger$mdrq_nfr
  message(sprintf(
    "[%s] OMDRC.OR+NFR ledger MDR_q: mean %.4f, max %.4f -> %s | mean blocked %.1f",
    tag, mean(v), max(v),
    if (max(v) <= alpha + 1e-9) "PASS" else "FAIL",
    mean(ledger$blocked_nfr)
  ))
}
.report_ledger("fixed", ledger_fixed, SIM_PARAMS$alpha)
.report_ledger("surge", ledger_smooth, SIM_PARAMS$alpha)

p_mdr_fixed <- create_mdr_plot(
  mdr_fixed,
  subtitle_text = "(a.2) Fixed signal proportion",
  alpha = SIM_PARAMS$alpha,
  y_upper = max(
    0.125,
    1.08 * max(mdr_fixed$FT, mdr_fixed$OMDRC.OR, mdr_fixed$OMDRC.OR.NFR,
               na.rm = TRUE)
  )
)

y_upper_smooth <- max(
  0.15,
  1.08 * max(mdr_smooth$FT, mdr_smooth$OMDRC.OR, mdr_smooth$OMDRC.OR.NFR,
             na.rm = TRUE)
)

message(sprintf(
  "Fixed pi: max FT MDR = %.3f; terminal FT MDR = %.3f; max OMDRC.OR+NFR MDR = %.3f",
  max(mdr_fixed$FT, na.rm = TRUE), tail(mdr_fixed$FT, 1),
  max(mdr_fixed$OMDRC.OR.NFR, na.rm = TRUE)
))
message(sprintf(
  "Surge pi: max FT MDR = %.3f; terminal FT MDR = %.3f; max OMDRC.OR+NFR MDR = %.3f",
  max(mdr_smooth$FT, na.rm = TRUE),
  tail(mdr_smooth$FT, 1),
  max(mdr_smooth$OMDRC.OR.NFR, na.rm = TRUE)
))

p_mdr_smooth <- create_mdr_plot(
  mdr_smooth,
  subtitle_text = "(b.2) Surge signal proportion",
  alpha = SIM_PARAMS$alpha,
  y_upper = y_upper_smooth
)

mdr_column <- (p_mdr_fixed / p_mdr_smooth) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom") &
  # The right block is only part of the figure width, so stack the three-method
  # legend vertically to avoid clipping the FT entry.
  guides(colour = guide_legend(ncol = 1), shape = guide_legend(ncol = 1))

# -----------------------------------------------------------------------------
# 9b. Monte Carlo FDR plots (companion of the MDR panels)
# -----------------------------------------------------------------------------

# Shared honest y-range per regime so the three arms remain distinguishable
# (FDR sits near 1 in these sparse settings).
.fdr_ylim <- function(fdr_wide, floor_span = 0.05) {
  vals <- unlist(fdr_wide[c("OMDRC.OR", "OMDRC.OR.NFR", "FT")], use.names = FALSE)
  vals <- vals[is.finite(vals) & fdr_wide$t >= FDR_WARMUP]   # drop the cold start
  lo <- min(vals, na.rm = TRUE)
  hi <- max(vals, na.rm = TRUE)
  pad <- max(0.01, 0.15 * (hi - lo))
  c(max(0, lo - pad), min(1, hi + pad))
}

fdr_ylim_fixed <- .fdr_ylim(fdr_fixed)
fdr_ylim_smooth <- .fdr_ylim(fdr_smooth)

message(sprintf(
  "Fixed pi: terminal FDR  OMDRC.OR = %.3f, +NFR = %.3f, FT = %.3f",
  tail(fdr_fixed$OMDRC.OR, 1), tail(fdr_fixed$OMDRC.OR.NFR, 1),
  tail(fdr_fixed$FT, 1)
))
message(sprintf(
  "Surge pi: terminal FDR  OMDRC.OR = %.3f, +NFR = %.3f, FT = %.3f",
  tail(fdr_smooth$OMDRC.OR, 1), tail(fdr_smooth$OMDRC.OR.NFR, 1),
  tail(fdr_smooth$FT, 1)
))

p_fdr_fixed <- create_fdr_plot(
  fdr_fixed,
  subtitle_text = "(a.3) Fixed signal proportion",
  y_lower = fdr_ylim_fixed[1],
  y_upper = fdr_ylim_fixed[2]
)

p_fdr_smooth <- create_fdr_plot(
  fdr_smooth,
  subtitle_text = "(b.3) Surge signal proportion",
  y_lower = fdr_ylim_smooth[1],
  y_upper = fdr_ylim_smooth[2]
)

# MDR and FDR share one three-arm legend collected across both columns.
metric_block <- ((p_mdr_fixed / p_mdr_smooth) | (p_fdr_fixed / p_fdr_smooth)) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom") &
  guides(colour = guide_legend(nrow = 1), shape = guide_legend(nrow = 1))

# -----------------------------------------------------------------------------
# 10. Final Figure 2 (NFR): 2 rows x 3 columns (capacity | MDR | FDR)
# -----------------------------------------------------------------------------

final_plot <- wrap_elements(capacity_column) |
  wrap_elements(metric_block)
# capacity : (MDR + FDR) so the three visual columns are ~1.8 : 1 : 1.
final_plot <- final_plot + plot_layout(widths = c(1.8, 2))

if (interactive()) print(final_plot)

if (isTRUE(SAVE_FIGURE)) {
  ggsave(
    filename = paste0(OUTPUT_STEM, ".pdf"),
    plot = final_plot,
    width = 19,
    height = 8.5,
    units = "in"
  )

  ggsave(
    filename = paste0(OUTPUT_STEM, ".png"),
    plot = final_plot,
    width = 19,
    height = 8.5,
    units = "in",
    dpi = 300
  )
  message("Saved figure: ", OUTPUT_STEM, ".{pdf,png}")
}

# Suggested manuscript description:
# The data-generating mechanism follows Setting 2: F0 = N(0,1) and F1 = N(3,1),
# with the constant pi_t = 0.05 (panels a) and the Gaussian surge path
#   pi_t = 0.05 + 0.25 exp{-(u-0.6)^2 / (2 * 0.15^2)},  u = (t-1)/(T-1)
# (panels b). OMDRC.OR uses the credited local barrier
# OMDRC_OR_NFR with window d = 100: a skip requires BOTH the capacity condition
# (1-alpha)Lmdr_t <= C_t AND the barrier condition Lmdr_t < gamma_t, where
# gamma_t is a causal trailing-window FT-type threshold. Panels (a.1)/(b.1) show
# the NFR capacity ledger C_t, the local cost (1-alpha)Lmdr_t, the barrier
# (1-alpha)gamma_t, and (as crosses) the free-riding skips the barrier converts
# into rejections. Panels (a.2)/(b.2) compare OMDRC.OR_nob, OMDRC.OR, and the
# offline oracle FT; both online arms
# keep MDR <= alpha on every path.
