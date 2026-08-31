# ==============================================================================
# Anti-free-riding (NFR) capacity rule -- shared core, sourced by the scripts in
# this folder (and inside parallel workers).
#
# Capacity-only comparison (called OMDRC.OR_nob in Figure 2):
#     delta_t = 0  <=>  (1 - alpha) q_t <= C_t
# Manuscript OMDRC rule with the local barrier:
#     within a causal trailing window of length d_t (current point included)
#         q_(1) <= ... <= q_(d_t),
#         k_t     = max{ j : (sum_{l<=j} q_(l)) / (sum_{l<=d_t} q_(l)) <= alpha }
#         gamma_t = q_(k_t + 1)
#     delta_t = 0  <=>  q_t < gamma_t  AND  (1 - alpha) q_t <= C_t
#
# Ledger is unchanged: skip -> C - (1-alpha) q_t, reject -> C + alpha q_t.
# NFR only converts skips into rejections, hence C_t >= 0 is preserved and the
# capacity argument (sum_skip q_t <= alpha * sum_all q_t) still goes through.
# ==============================================================================

#' Local FT-type threshold gamma_t on a window of Lmdr scores
#'
#' k_t = 0    -> gamma_t = min(window); `q_t < gamma_t` is then impossible and
#'               the point must be rejected.
#' k_t = d_t  -> gamma_t = +Inf, the constraint is inactive.
.local_gamma <- function(w, alpha) {
  s <- sort(as.numeric(w))
  tot <- sum(s)
  if (!is.finite(tot) || tot <= .Machine$double.eps) return(Inf)
  k <- sum(cumsum(s) / tot <= alpha)   # share is monotone -> the set is a prefix
  if (k >= length(s)) return(Inf)
  s[k + 1L]
}

#' Causal path of the local threshold gamma_t
#'
#' gamma_t depends ONLY on the trailing window of scores, never on the
#' decisions, so it can be precomputed once per (stream, d) and reused both by
#' the NFR rule and as a *diagnostic* applied to any other rule's decisions.
#' @export
nfr_gamma_path <- function(x.Lmdr, alpha, d, x.Lmdr_ini = NULL) {
  q <- pmin(pmax(as.numeric(x.Lmdr), 0), 1)
  m <- length(q)
  d <- as.integer(d)
  hist0 <- if (is.null(x.Lmdr_ini)) numeric(0) else
    pmin(pmax(as.numeric(x.Lmdr_ini), 0), 1)
  vapply(seq_len(m), function(t) {
    .local_gamma(utils::tail(c(hist0, q[seq_len(t)]), d), alpha)
  }, numeric(1))
}

#' Oracle online MDR control WITH the anti-free-riding constraint
#'
#' @param x.Lmdr      oracle Lmdr scores over the stream
#' @param alpha       MDR budget
#' @param d           trailing window length used for gamma_t
#' @param x.Lmdr_ini  Lmdr of the calibration batch, used to fill the window
#'                    during warm-up (t < d). Optional.
#' @param gamma       optional precomputed gamma path (see nfr_gamma_path)
#' @return list(de, gamma, capacity, blocked, d)
#'         `blocked[t] = 1` marks the free-riding attempts: capacity alone would
#'         have allowed a skip, but gamma_t forced a rejection.
OMDRC_OR_NFR <- function(x.Lmdr, alpha, d, x.Lmdr_ini = NULL, gamma = NULL) {
  q <- pmin(pmax(as.numeric(x.Lmdr), 0), 1)
  m <- length(q)
  if (!is.numeric(d) || length(d) != 1L || !is.finite(d) || d < 2)
    stop("d must be a scalar integer >= 2.")
  d <- as.integer(d)
  hist0 <- if (is.null(x.Lmdr_ini)) numeric(0) else
    pmin(pmax(as.numeric(x.Lmdr_ini), 0), 1)
  gamma_path <- if (is.null(gamma))
    nfr_gamma_path(q, alpha, d, hist0) else as.numeric(gamma)
  if (length(gamma_path) != m) stop("gamma must have length length(x.Lmdr).")

  delta <- integer(m)
  cap_path <- numeric(m + 1L)
  blocked <- integer(m)
  capacity <- 0

  for (t in seq_len(m)) {
    g <- gamma_path[t]
    cap_ok <- ((1 - alpha) * q[t] <= capacity)
    small <- (q[t] < g)
    if (cap_ok && small) {
      delta[t] <- 0L
      capacity <- capacity - (1 - alpha) * q[t]
    } else {
      delta[t] <- 1L
      capacity <- capacity + alpha * q[t]
      if (cap_ok) blocked[t] <- 1L
    }
    cap_path[t + 1L] <- capacity
  }

  list(de = delta, gamma = gamma_path, capacity = cap_path,
       blocked = blocked, d = d)
}

#' Anti-free-riding overlay with a shadow capacity ledger
#'
#' The plain capacity rule is run on a virtual (shadow) ledger to produce the
#' base decisions.  The NFR gate may only flip a base skip (0) to a discovery
#' (1); that flip is deliberately NOT credited to the shadow ledger.  Hence the
#' final skip set is a subset of the plain rule's skip set on every path.
#'
#' `capacity` is the controlling shadow ledger.  `actual_capacity` replays the
#' final decisions and is returned only as a diagnostic; it is not spendable.
#'
#' @return list(de, base_de, gamma, capacity, actual_capacity, blocked, d)
#' @export
OMDRC_NFR_SHADOW <- function(x.Lmdr, alpha, d, x.Lmdr_ini = NULL,
                             gamma = NULL) {
  q <- pmin(pmax(as.numeric(x.Lmdr), 0), 1)
  m <- length(q)
  if (!is.numeric(d) || length(d) != 1L || !is.finite(d) || d < 2)
    stop("d must be a scalar integer >= 2.")
  d <- as.integer(d)
  hist0 <- if (is.null(x.Lmdr_ini)) numeric(0) else
    pmin(pmax(as.numeric(x.Lmdr_ini), 0), 1)
  gamma_path <- if (is.null(gamma))
    nfr_gamma_path(q, alpha, d, hist0) else as.numeric(gamma)
  if (length(gamma_path) != m) stop("gamma must have length length(x.Lmdr).")

  base_de <- integer(m)
  delta <- integer(m)
  blocked <- integer(m)
  shadow_path <- numeric(m + 1L)
  actual_path <- numeric(m + 1L)
  shadow_capacity <- 0
  actual_capacity <- 0

  for (t in seq_len(m)) {
    base_skip <- ((1 - alpha) * q[t] <= shadow_capacity)
    if (base_skip) {
      base_de[t] <- 0L
      shadow_capacity <- shadow_capacity - (1 - alpha) * q[t]
    } else {
      base_de[t] <- 1L
      shadow_capacity <- shadow_capacity + alpha * q[t]
    }

    # The overlay can only add discoveries.  Its forced discovery never mints
    # spendable capacity; future base decisions continue on the shadow ledger.
    delta[t] <- base_de[t]
    if (base_skip && q[t] >= gamma_path[t]) {
      delta[t] <- 1L
      blocked[t] <- 1L
    }

    actual_capacity <- if (delta[t] == 0L)
      actual_capacity - (1 - alpha) * q[t] else
      actual_capacity + alpha * q[t]
    shadow_path[t + 1L] <- shadow_capacity
    actual_path[t + 1L] <- actual_capacity
  }

  stopifnot(all(delta >= base_de),
            all(shadow_path >= -100 * .Machine$double.eps),
            all(actual_path + 100 * .Machine$double.eps >= shadow_path))

  list(de = delta, base_de = base_de, gamma = gamma_path,
       capacity = shadow_path, actual_capacity = actual_path,
       blocked = blocked, d = d)
}

#' Anti-free-riding overlay with a NO-CREDIT (budget-saving) ledger
#'
#' SAST-like treatment of a barrier-forced rejection: when the plain capacity
#' rule would skip (`(1-alpha) q <= C`) but the gate fails (`q >= gamma`), the
#' point is rejected AND the capacity is left UNCHANGED -- the (1-alpha)q that a
#' skip would have spent is saved for later low-q skips, and no phantom
#' +alpha*q reward is minted (unlike OMDRC_OR_NFR).
#'
#' Validity: the spendable ledger `capacity` only debits (1-alpha)q on an
#' ORDINARY skip (guarded by the base-skip test, so capacity stays >= 0) and
#' credits +alpha*q on a natural rejection; forced rejections touch neither.
#' The FORMAL estimated ledger that credits EVERY rejection is
#'   formal_t = capacity_t + alpha * sum_{i<=t forced} q_i  >=  capacity_t >= 0,
#' so the estimated-capacity inequality (1-alpha) sum_skip q <= alpha sum_rej q
#' still holds and the perturbation transfer to the true q is unchanged.
#'
#' @return list(de, gamma, capacity, formal, blocked, d). `capacity` is the
#'   spendable ledger; `formal` is the estimated ledger crediting all rejections.
#' @export
OMDRC_NFR_NOCREDIT <- function(x.Lmdr, alpha, d, x.Lmdr_ini = NULL,
                               gamma = NULL, tol = 1e-12) {
  q <- pmin(pmax(as.numeric(x.Lmdr), 0), 1)
  m <- length(q)
  if (!is.numeric(d) || length(d) != 1L || !is.finite(d) || d < 2)
    stop("d must be a scalar integer >= 2.")
  d <- as.integer(d)
  hist0 <- if (is.null(x.Lmdr_ini)) numeric(0) else
    pmin(pmax(as.numeric(x.Lmdr_ini), 0), 1)
  gamma_path <- if (is.null(gamma))
    nfr_gamma_path(q, alpha, d, hist0) else as.numeric(gamma)
  if (length(gamma_path) != m) stop("gamma must have length length(x.Lmdr).")

  delta <- integer(m)
  blocked <- integer(m)
  cap_path <- numeric(m + 1L)     # spendable ledger (governs the base-skip test)
  formal_path <- numeric(m + 1L)  # formal estimated ledger (credit all rejects)
  capacity <- 0
  formal <- 0

  for (t in seq_len(m)) {
    cost <- (1 - alpha) * q[t]
    base_skip <- (cost <= capacity + tol)
    gate_pass <- (q[t] < gamma_path[t])
    if (base_skip && gate_pass) {           # ordinary skip
      delta[t] <- 0L
      capacity <- capacity - cost
      formal <- formal - cost
    } else if (base_skip && !gate_pass) {   # NFR forced rejection: save budget
      delta[t] <- 1L
      blocked[t] <- 1L                      # capacity unchanged; formal credited
      formal <- formal + alpha * q[t]
    } else {                                # natural rejection
      delta[t] <- 1L
      capacity <- capacity + alpha * q[t]
      formal <- formal + alpha * q[t]
    }
    cap_path[t + 1L] <- capacity
    formal_path[t + 1L] <- formal
  }

  stopifnot(all(cap_path >= -100 * .Machine$double.eps),
            all(formal_path >= -100 * .Machine$double.eps))

  list(de = delta, gamma = gamma_path, capacity = cap_path,
       formal = formal_path, blocked = blocked, d = d)
}

#' Replay the capacity ledger of an arbitrary decision sequence
#' (used to record C_t for the plain OMDRC_OR rule).
replay_capacity <- function(x.Lmdr, de, alpha) {
  q <- pmin(pmax(as.numeric(x.Lmdr), 0), 1)
  cap <- numeric(length(q) + 1L)
  cc <- 0
  for (t in seq_along(q)) {
    cc <- if (de[t] == 0L) cc - (1 - alpha) * q[t] else cc + alpha * q[t]
    cap[t + 1L] <- cc
  }
  cap
}
