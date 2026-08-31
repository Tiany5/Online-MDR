# Online MDR control with smoothly time-varying signal prevalence.
#
# The component distributions F0 and F1 are stable and learned from labeled
# reference samples. The local signal prevalence pi_t is estimated from a
# trailing window of unlabeled observations (including the current point).
# The density-ratio estimator is modular: the PRIMARY implementation is
# PC-DRE via a logistic GAM classifier (paper Stage 1); KDE, uLSIF, and a
# glmnet classifier wrapper are provided as alternatives. The helper
# `make_smooth_prior_drop()` reproduces the fixed simulation path used in Setting (a).

# -----------------------------------------------------------------------------
# Generic helpers
# -----------------------------------------------------------------------------

.as_numeric_matrix <- function(x, name = "x") {
  if (is.data.frame(x)) x <- as.matrix(x)
  if (is.vector(x) && !is.list(x)) x <- matrix(as.numeric(x), ncol = 1L)
  if (!is.matrix(x) || !is.numeric(x) || nrow(x) < 1L || ncol(x) < 1L) {
    stop(name, " must be a non-empty numeric vector, matrix, or data frame.")
  }
  storage.mode(x) <- "double"
  if (any(!is.finite(x))) stop(name, " contains non-finite values.")
  x
}

.validate_pi_bounds <- function(pi_bounds) {
  if (!is.numeric(pi_bounds) || length(pi_bounds) != 2L ||
      any(!is.finite(pi_bounds)) || pi_bounds[1L] <= 0 ||
      pi_bounds[2L] >= 1 || pi_bounds[1L] >= pi_bounds[2L]) {
    stop("pi_bounds must be c(pi_lower, pi_upper) with 0 < lower < upper < 1.")
  }
  as.numeric(pi_bounds)
}

#' Smoothly decreasing signal-prevalence path used in Setting (a)
#'
#' The path starts near `pi_high` and decreases smoothly toward `pi_low`:
#' pi_t = pi_low + (pi_high-pi_low) * logistic(-(t/T-center)/width).
#' @export
make_smooth_prior_drop <- function(N,
                                   pi_low = 0.03,
                                   pi_high = 0.35,
                                   drop_center = 0.50,
                                   drop_width = 0.10) {
  if (!is.numeric(N) || length(N) != 1L || !is.finite(N) || N < 2) {
    stop("N must be an integer at least 2.")
  }
  if (!is.numeric(pi_low) || !is.numeric(pi_high) ||
      length(pi_low) != 1L || length(pi_high) != 1L ||
      !is.finite(pi_low) || !is.finite(pi_high) ||
      pi_low <= 0 || pi_high >= 1 || pi_low >= pi_high) {
    stop("Require 0 < pi_low < pi_high < 1.")
  }
  if (!is.numeric(drop_center) || length(drop_center) != 1L ||
      !is.finite(drop_center) || drop_center <= 0 || drop_center >= 1) {
    stop("drop_center must lie in (0, 1).")
  }
  if (!is.numeric(drop_width) || length(drop_width) != 1L ||
      !is.finite(drop_width) || drop_width <= 0) {
    stop("drop_width must be positive.")
  }

  u <- seq(0, 1, length.out = as.integer(N))
  pi_t <- pi_low + (pi_high - pi_low) *
    stats::plogis(-(u - drop_center) / drop_width)
  pmin(pmax(pi_t, pi_low), pi_high)
}

.clip_ratio <- function(r, floor = 1e-8, cap = 1e6) {
  r <- as.numeric(r)
  if (!is.numeric(floor) || length(floor) != 1L || !is.finite(floor) || floor <= 0) {
    stop("ratio_floor must be a finite positive scalar.")
  }
  if (!is.numeric(cap) || length(cap) != 1L || is.na(cap) || cap <= floor) {
    stop("ratio_cap must exceed ratio_floor; Inf is allowed.")
  }
  r[!is.finite(r)] <- cap
  pmin(pmax(r, floor), cap)
}

# -----------------------------------------------------------------------------
# KDE density-ratio estimator (one-dimensional)
# -----------------------------------------------------------------------------

.safe_bw_sj <- function(x, multiplier = 1.5) {
  x <- as.numeric(x)
  sx <- stats::sd(x)
  if (length(x) < 2L || !is.finite(sx) || sx == 0) {
    return(max(abs(mean(x)), 1) * 1e-2)
  }

  bw <- tryCatch(stats::bw.SJ(x), error = function(e) stats::bw.nrd0(x))
  bw <- bw * multiplier
  if (!is.finite(bw) || bw <= 0) bw <- stats::bw.nrd0(x)
  max(bw, 1e-6)
}

.get_kde <- function(x_data, h, from, to, n, reflect = FALSE) {
  x_data <- as.numeric(x_data)
  if (reflect) {
    d <- stats::density(c(x_data, -x_data), bw = h, from = from, to = to, n = n)
    d$y <- d$y * 2
    d$y[d$x < 0] <- 0
  } else {
    d <- stats::density(x_data, bw = h, from = from, to = to, n = n)
  }
  d
}

.eval_kde <- function(den, x, density_floor = 1e-8) {
  val <- stats::approx(den$x, den$y, xout = as.numeric(x), rule = 2)$y
  pmax(as.numeric(val), density_floor)
}

#' Fit a KDE estimator of r(x) = f1(x) / f0(x)
#'
#' This estimator is restricted to one-dimensional observations.
#' @export
fit_ratio_kde <- function(x0, x1, grid_data = NULL,
                          grid_n = 2048L, pad = 5,
                          bandwidth_multiplier = 1.5,
                          density_floor = 1e-8,
                          ratio_floor = 1e-8,
                          ratio_cap = 1e6) {
  x0 <- .as_numeric_matrix(x0, "x0")
  x1 <- .as_numeric_matrix(x1, "x1")
  if (ncol(x0) != 1L || ncol(x1) != 1L) {
    stop("fit_ratio_kde supports one-dimensional observations only.")
  }
  if (ncol(x0) != ncol(x1)) stop("x0 and x1 must have the same dimension.")

  x0v <- as.numeric(x0[, 1L])
  x1v <- as.numeric(x1[, 1L])
  grid_values <- c(x0v, x1v)
  if (!is.null(grid_data)) {
    gd <- .as_numeric_matrix(grid_data, "grid_data")
    if (ncol(gd) != 1L) stop("grid_data must be one-dimensional for KDE.")
    grid_values <- c(grid_values, gd[, 1L])
  }

  positive_support <- all(grid_values >= 0)
  from0 <- if (positive_support) 0 else min(grid_values) - pad
  to0 <- max(grid_values) + pad
  h0 <- .safe_bw_sj(x0v, bandwidth_multiplier)
  h1 <- .safe_bw_sj(x1v, bandwidth_multiplier)
  den0 <- .get_kde(x0v, h0, from0, to0, as.integer(grid_n), positive_support)
  den1 <- .get_kde(x1v, h1, from0, to0, as.integer(grid_n), positive_support)

  structure(
    list(
      method = "kde",
      den0 = den0,
      den1 = den1,
      h0 = h0,
      h1 = h1,
      density_floor = density_floor,
      ratio_floor = ratio_floor,
      ratio_cap = ratio_cap,
      dimension = 1L,
      from = from0,
      to = to0
    ),
    class = "omdrc_ratio_model"
  )
}

# -----------------------------------------------------------------------------
# uLSIF direct density-ratio estimator
# -----------------------------------------------------------------------------

.squared_distance_matrix <- function(x, centers) {
  x <- .as_numeric_matrix(x, "x")
  centers <- .as_numeric_matrix(centers, "centers")
  if (ncol(x) != ncol(centers)) stop("x and centers must have the same dimension.")
  x2 <- rowSums(x^2)
  c2 <- rowSums(centers^2)
  pmax(outer(x2, c2, "+") - 2 * tcrossprod(x, centers), 0)
}

.gaussian_basis <- function(x, centers, sigma) {
  if (!is.finite(sigma) || sigma <= 0) stop("sigma must be positive.")
  exp(-.squared_distance_matrix(x, centers) / (2 * sigma^2))
}

.median_distance <- function(x, max_points = 300L) {
  x <- .as_numeric_matrix(x, "x")
  if (nrow(x) > max_points) x <- x[sample.int(nrow(x), max_points), , drop = FALSE]
  if (nrow(x) < 2L) return(1)
  dd <- stats::dist(x)
  med <- stats::median(as.numeric(dd))
  if (!is.finite(med) || med <= 0) {
    med <- sqrt(sum(apply(x, 2L, stats::var)))
  }
  if (!is.finite(med) || med <= 0) med <- 1
  med
}

.solve_ulsif <- function(phi0, phi1, lambda) {
  n0 <- nrow(phi0)
  H <- crossprod(phi0) / n0
  h <- colMeans(phi1)
  ridge <- H + diag(lambda, ncol(H))
  alpha <- tryCatch(
    solve(ridge, h),
    error = function(e) qr.solve(ridge, h, tol = 1e-10)
  )

  # uLSIF is unconstrained in the basis coefficients. Some fitted
  # coefficients may therefore be negative. Non-negativity is imposed on the
  # final predicted density ratio, not coefficient-by-coefficient.
  as.numeric(alpha)
}

.make_balanced_folds <- function(n, k) {
  if (n < 2L) return(rep(1L, n))
  k <- min(as.integer(k), n)
  sample(rep(seq_len(k), length.out = n))
}

#' Fit an unconstrained Least-Squares Importance Fitting estimator
#'
#' Estimates r(x)=f1(x)/f0(x) directly with Gaussian basis functions. Hyper-
#' parameters are selected by the standard uLSIF validation objective.
#' @export
fit_ratio_ulsif <- function(x0, x1,
                            n_centers = 100L,
                            sigma_grid = NULL,
                            lambda_grid = 10^seq(-6, -1, length.out = 6L),
                            n_folds = 5L,
                            tune = TRUE,
                            normalize = TRUE,
                            calibration_fraction = 0.20,
                            min_calibration = 50L,
                            ratio_floor = 1e-8,
                            ratio_cap = 1e6,
                            seed = NULL) {
  x0 <- .as_numeric_matrix(x0, "x0")
  x1 <- .as_numeric_matrix(x1, "x1")
  if (ncol(x0) != ncol(x1)) stop("x0 and x1 must have the same dimension.")
  if (!is.null(seed)) set.seed(seed)

  if (!is.numeric(calibration_fraction) || length(calibration_fraction) != 1L ||
      !is.finite(calibration_fraction) || calibration_fraction < 0 ||
      calibration_fraction >= 0.5) {
    stop("calibration_fraction must lie in [0, 0.5).")
  }

  # Reserve an independent part of the null reference sample for absolute
  # ratio calibration. This matters because local estimation of pi_t uses the
  # numerical scale of r(x), not merely its ranking.
  n0 <- nrow(x0)
  min_fit <- max(20L, as.integer(n_folds) + 2L)
  requested_cal <- max(as.integer(min_calibration),
                       floor(calibration_fraction * n0))
  n_cal <- if (isTRUE(normalize) && n0 > min_fit) {
    min(requested_cal, n0 - min_fit)
  } else {
    0L
  }

  if (n_cal > 0L) {
    calibration_idx <- sample.int(n0, n_cal, replace = FALSE)
    x0_cal <- x0[calibration_idx, , drop = FALSE]
    x0_fit <- x0[-calibration_idx, , drop = FALSE]
  } else {
    x0_fit <- x0
    x0_cal <- x0
  }

  n_centers <- max(1L, min(as.integer(n_centers), nrow(x1)))
  center_idx <- sample.int(nrow(x1), n_centers, replace = FALSE)
  centers <- x1[center_idx, , drop = FALSE]

  if (is.null(sigma_grid)) {
    base_sigma <- .median_distance(rbind(x0_fit, x1))
    # The smaller bandwidths are important for the separated modes in Setting (a).
    sigma_grid <- base_sigma * c(0.125, 0.25, 0.5, 1, 2)
  }
  sigma_grid <- unique(as.numeric(sigma_grid))
  sigma_grid <- sigma_grid[is.finite(sigma_grid) & sigma_grid > 0]
  lambda_grid <- unique(as.numeric(lambda_grid))
  lambda_grid <- lambda_grid[is.finite(lambda_grid) & lambda_grid > 0]
  if (length(sigma_grid) == 0L || length(lambda_grid) == 0L) {
    stop("sigma_grid and lambda_grid must contain positive values.")
  }

  best_score <- Inf
  best_sigma <- sigma_grid[1L]
  best_lambda <- lambda_grid[1L]
  cv_table <- expand.grid(sigma = sigma_grid, lambda = lambda_grid)
  cv_table$score <- NA_real_

  if (isTRUE(tune) && min(nrow(x0_fit), nrow(x1)) >= 4L) {
    k <- min(as.integer(n_folds), nrow(x0_fit), nrow(x1))
    folds0 <- .make_balanced_folds(nrow(x0_fit), k)
    folds1 <- .make_balanced_folds(nrow(x1), k)

    row_id <- 0L
    for (sigma in sigma_grid) {
      phi0_all <- .gaussian_basis(x0_fit, centers, sigma)
      phi1_all <- .gaussian_basis(x1, centers, sigma)
      for (lambda in lambda_grid) {
        row_id <- row_id + 1L
        fold_scores <- numeric(k)
        for (fold in seq_len(k)) {
          train0 <- folds0 != fold
          train1 <- folds1 != fold
          valid0 <- !train0
          valid1 <- !train1
          alpha <- .solve_ulsif(
            phi0_all[train0, , drop = FALSE],
            phi1_all[train1, , drop = FALSE],
            lambda
          )
          r0 <- as.numeric(phi0_all[valid0, , drop = FALSE] %*% alpha)
          r1 <- as.numeric(phi1_all[valid1, , drop = FALSE] %*% alpha)
          fold_scores[fold] <- 0.5 * mean(r0^2) - mean(r1)
        }
        score <- mean(fold_scores)
        cv_table$score[row_id] <- score
        if (is.finite(score) && score < best_score) {
          best_score <- score
          best_sigma <- sigma
          best_lambda <- lambda
        }
      }
    }
  } else {
    cv_table <- cv_table[1L, , drop = FALSE]
    cv_table$score <- NA_real_
  }

  phi0_fit <- .gaussian_basis(x0_fit, centers, best_sigma)
  phi1 <- .gaussian_basis(x1, centers, best_sigma)
  alpha <- .solve_ulsif(phi0_fit, phi1, best_lambda)

  fit_ratio0_raw <- as.numeric(phi0_fit %*% alpha)
  fit_ratio0 <- pmax(fit_ratio0_raw, 0)
  cal_ratio0_raw <- as.numeric(
    .gaussian_basis(x0_cal, centers, best_sigma) %*% alpha
  )
  cal_ratio0 <- pmax(cal_ratio0_raw, 0)

  normalizer <- 1
  if (isTRUE(normalize)) {
    # Since E_{F0}{r(X)} = 1, normalize with null observations not used to fit
    # the uLSIF coefficients. This avoids optimistic in-sample calibration.
    normalizer <- mean(cal_ratio0)
    if (!is.finite(normalizer) || normalizer <= .Machine$double.eps) {
      normalizer <- 1
    }
  }

  structure(
    list(
      method = "ulsif",
      centers = centers,
      sigma = best_sigma,
      lambda = best_lambda,
      alpha = alpha,
      normalizer = normalizer,
      normalize = isTRUE(normalize),
      calibration_fraction = calibration_fraction,
      calibration_size = nrow(x0_cal),
      fit_null_size = nrow(x0_fit),
      negative_train_fraction = mean(fit_ratio0_raw < 0),
      negative_calibration_fraction = mean(cal_ratio0_raw < 0),
      mean_fit_ratio = mean(fit_ratio0 / normalizer),
      mean_calibration_ratio = mean(cal_ratio0 / normalizer),
      cv_table = cv_table,
      ratio_floor = ratio_floor,
      ratio_cap = ratio_cap,
      dimension = ncol(x0)
    ),
    class = "omdrc_ratio_model"
  )
}

# -----------------------------------------------------------------------------
# Optional classifier density-ratio estimator for later high-dimensional use
# -----------------------------------------------------------------------------

#' Fit a classifier-based density-ratio estimator
#'
#' Requires package glmnet. This wrapper is not used by the current one-
#' dimensional simulation, but it can be used directly in a future high-
#' dimensional setting.
#' @export
fit_ratio_classifier <- function(x0, x1, alpha = 1, nfolds = 5L,
                                 type_measure = "deviance",
                                 ratio_floor = 1e-8,
                                 ratio_cap = 1e6,
                                 seed = NULL) {
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("Package 'glmnet' is required for fit_ratio_classifier().")
  }
  x0 <- .as_numeric_matrix(x0, "x0")
  x1 <- .as_numeric_matrix(x1, "x1")
  if (ncol(x0) != ncol(x1)) stop("x0 and x1 must have the same dimension.")
  if (!is.null(seed)) set.seed(seed)

  x <- rbind(x0, x1)
  y <- c(rep(0, nrow(x0)), rep(1, nrow(x1)))
  rho <- mean(y)
  fit <- glmnet::cv.glmnet(
    x = x,
    y = y,
    family = "binomial",
    alpha = alpha,
    nfolds = min(as.integer(nfolds), length(y)),
    type.measure = type_measure
  )

  structure(
    list(
      method = "classifier",
      fit = fit,
      rho = rho,
      ratio_floor = ratio_floor,
      ratio_cap = ratio_cap,
      dimension = ncol(x0)
    ),
    class = "omdrc_ratio_model"
  )
}

#' Fit a density ratio via an additive mgcv GAM classifier (PC-DRE)
#'
#' Classifier trick: pool F0 (y=0) and F1 (y=1), fit a cubic-regression-spline
#' logistic GAM, then convert the posterior to a density ratio with the prior
#' correction log((1 - rho)/rho). This is the paper's primary PC-DRE estimator
#' for both one-dimensional and multivariate inputs (`ratio_method = "gam"`).
#' @export
fit_ratio_gam <- function(x0, x1,
                          spline_k = 10L,
                          prob_floor = 1e-6,
                          ratio_floor = 1e-8,
                          ratio_cap = 1e6,
                          gamma = 1,
                          select = FALSE,
                          seed = NULL) {
  if (!requireNamespace("mgcv", quietly = TRUE)) {
    stop("Package 'mgcv' is required for fit_ratio_gam().")
  }
  x0 <- .as_numeric_matrix(x0, "x0")
  x1 <- .as_numeric_matrix(x1, "x1")
  if (ncol(x0) != ncol(x1)) stop("x0 and x1 must have the same dimension.")
  if (!is.null(seed)) set.seed(seed)

  d <- ncol(x0)
  n0 <- nrow(x0)
  n1 <- nrow(x1)
  if (n0 < 1L || n1 < 1L) stop("x0 and x1 must both be non-empty.")

  # Binary-classification pool. Labels are NOT reversed: y = 1 marks F1.
  rho_train <- n1 / (n0 + n1)
  if (!is.finite(rho_train) || rho_train <= 0 || rho_train >= 1) {
    stop("rho_train = n1 / (n0 + n1) must lie strictly in (0, 1).")
  }
  prior_correction <- log((1 - rho_train) / rho_train)
  if (!is.finite(prior_correction)) {
    stop("Prior correction log((1 - rho_train)/rho_train) is not finite.")
  }

  # Feature columns V1..Vd. Multi-dimensional inputs use an ADDITIVE spline
  # model y ~ s(V1) + ... + s(Vd) (no interaction terms), so the fit scales to
  # moderate d without the tensor-product curse of dimensionality. In 1-D this
  # reduces to the original single-smooth classifier.
  feature_names <- paste0("V", seq_len(d))
  X_all <- rbind(x0, x1)
  if (any(!is.finite(X_all))) {
    stop("Reference samples contain non-finite values.")
  }
  train_df <- as.data.frame(X_all)
  names(train_df) <- feature_names
  train_df$y <- c(rep(0L, n0), rep(1L, n1))

  rhs <- paste(sprintf("s(%s, bs = \"cr\", k = %d)",
                       feature_names, as.integer(spline_k)),
               collapse = " + ")
  gam_formula <- stats::as.formula(paste("y ~", rhs))

  # `gamma` inflates the effective degrees-of-freedom cost in the REML criterion
  # (gamma = 1 is plain REML, the default used by the simulations) and `select`
  # adds a shrinkage penalty on each smooth's null space.  Both default to the
  # historical behaviour, so existing experiments are bit-identical; raising
  # gamma / enabling select is what prevents the coefficient blow-up when the
  # two reference blocks happen to be (almost) perfectly separable.
  fit <- tryCatch(
    mgcv::gam(
      gam_formula,
      family = stats::binomial(link = "logit"),
      method = "REML",
      gamma = gamma,
      select = isTRUE(select),
      data = train_df
    ),
    error = function(e) {
      stop("GAM fitting failed in fit_ratio_gam(): ", conditionMessage(e))
    }
  )

  structure(
    list(
      method = "gam",
      fit = fit,
      feature_names = feature_names,
      rho_train = rho_train,
      prior_correction = prior_correction,
      prob_floor = prob_floor,
      ratio_floor = ratio_floor,
      ratio_cap = ratio_cap,
      gamma = gamma,
      select = isTRUE(select),
      dimension = d
    ),
    class = "omdrc_ratio_model"
  )
}

#' Predict an estimated density ratio from a fitted ratio model
#' @export
predict_ratio <- function(model, newdata) {
  if (!inherits(model, "omdrc_ratio_model")) {
    stop("model must be produced by a fit_ratio_*() function.")
  }
  x <- .as_numeric_matrix(newdata, "newdata")
  if (ncol(x) != model$dimension) stop("newdata has the wrong dimension.")

  ratio <- switch(
    model$method,
    kde = {
      f0 <- .eval_kde(model$den0, x[, 1L], model$density_floor)
      f1 <- .eval_kde(model$den1, x[, 1L], model$density_floor)
      f1 / pmax(f0, model$density_floor)
    },
    ulsif = {
      phi <- .gaussian_basis(x, model$centers, model$sigma)
      raw_ratio <- as.numeric(phi %*% model$alpha)
      # Standard uLSIF non-negative correction: truncate the estimated ratio,
      # not the basis coefficients.
      pmax(raw_ratio, 0) / model$normalizer
    },
    classifier = {
      eta <- as.numeric(stats::predict(
        model$fit,
        newx = x,
        s = "lambda.min",
        type = "response"
      ))
      eta <- pmin(pmax(eta, 1e-8), 1 - 1e-8)
      ((1 - model$rho) / model$rho) * eta / (1 - eta)
    },
    gam = {
      newdf <- as.data.frame(x)
      names(newdf) <- model$feature_names
      eta <- as.numeric(stats::predict(
        model$fit, newdata = newdf, type = "response"
      ))
      # Prevent probabilities exactly at 0 or 1 before the logit transform.
      floor_p <- model$prob_floor
      eta <- pmin(pmax(eta, floor_p), 1 - floor_p)
      # log r_hat(x) = qlogis(eta) + log((1 - rho_train) / rho_train).
      log_r_hat <- stats::qlogis(eta) + model$prior_correction
      exp(log_r_hat)
    },
    stop("Unknown ratio model method: ", model$method)
  )

  .clip_ratio(ratio, model$ratio_floor, model$ratio_cap)
}

# -----------------------------------------------------------------------------
# Local prevalence estimation and the generic OMDRC engine
# -----------------------------------------------------------------------------

# Local MLE/EM update based only on ratios r_s=f1(X_s)/f0(X_s).
.estimate_local_pi_ratio <- function(ratio, init, pi_bounds,
                                     tol = 1e-8, max_iter = 100L) {
  pi_bounds <- .validate_pi_bounds(pi_bounds)
  ratio <- .clip_ratio(ratio, 1e-12, 1e12)
  p <- min(max(as.numeric(init), pi_bounds[1L]), pi_bounds[2L])
  converged <- FALSE
  iterations <- 0L

  for (iter in seq_len(as.integer(max_iter))) {
    posterior <- p * ratio / pmax(1 - p + p * ratio, .Machine$double.xmin)
    p_new <- min(max(mean(posterior), pi_bounds[1L]), pi_bounds[2L])
    iterations <- iter
    if (abs(p_new - p) <= tol * (1 + abs(p))) {
      p <- p_new
      converged <- TRUE
      break
    }
    p <- p_new
  }

  list(pi = p, converged = converged, iterations = iterations)
}

.ratio_to_lmdr <- function(ratio, pi_value) {
  ratio <- .clip_ratio(ratio, 1e-12, 1e12)
  ans <- pi_value * ratio / pmax(1 - pi_value + pi_value * ratio,
                                .Machine$double.xmin)
  pmin(pmax(ans, 0), 1)
}

# Local no-free-riding threshold used by Algorithms 1 and 2 in the paper.
# The current score is included in the causal trailing window.
.omdrc_barrier_gamma <- function(scores, alpha) {
  scores <- sort(pmin(pmax(as.numeric(scores), 0), 1))
  total_mass <- sum(scores)
  if (!is.finite(total_mass) || total_mass <= .Machine$double.eps) return(Inf)
  k <- sum(cumsum(scores) / total_mass <= alpha)
  if (k >= length(scores)) Inf else scores[k + 1L]
}

.validate_barrier_window <- function(w) {
  if (!is.numeric(w) || length(w) != 1L || !is.finite(w) || w < 2) {
    stop("w must be a scalar integer at least 2.")
  }
  as.integer(w)
}

#' Run data-driven OMDRC from precomputed density-ratio scores
#'
#' This is the estimator-agnostic Algorithm 2 core. It estimates the local
#' prevalence, forms Lmdr scores, and applies both the capacity and local
#' no-free-riding conditions.
#' @export
OMDRC_FROM_RATIO <- function(LR, LR_ini, alpha, D = NULL,
                             pi_bounds = c(0.01, 0.99),
                             em_tol = 1e-8, em_max_iter = 100L,
                             D_mode = c("growing", "fixed"),
                             D_beta = 0.6, D_min = 10L,
                             w = 100L) {
  LR <- .clip_ratio(LR, 1e-12, 1e12)
  LR_ini <- .clip_ratio(LR_ini, 1e-12, 1e12)
  if (length(LR) < 1L || length(LR_ini) < 1L) stop("LR and LR_ini must be non-empty.")
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0 || alpha >= 1) {
    stop("alpha must be a scalar in (0, 1).")
  }
  pi_bounds <- .validate_pi_bounds(pi_bounds)
  w <- .validate_barrier_window(w)
  D_mode <- match.arg(D_mode)
  if (D_mode == "fixed") {
    if (!is.numeric(D) || length(D) != 1L || !is.finite(D) || D < 2) {
      stop("D must be a scalar integer at least 2 when D_mode = 'fixed'.")
    }
    D <- as.integer(D)
  } else {
    if (!is.numeric(D_beta) || length(D_beta) != 1L || !is.finite(D_beta) ||
        D_beta <= 0 || D_beta >= 1)
      stop("D_beta must be in (0, 1) when D_mode = 'growing'.")
    D_min <- as.integer(max(2L, D_min))
  }

  n_ini <- length(LR_ini)
  # Warm-up window at t = 0 follows the same paper formula with t = 0:
  # D_0 = min{K0, max{D_min, floor(K0^beta)}} in growing mode.
  D_init <- if (D_mode == "growing") {
    min(n_ini, max(D_min, as.integer(floor(n_ini^D_beta))))
  } else {
    min(D, n_ini)
  }
  initial_idx <- tail(seq_along(LR_ini), D_init)
  initial_fit <- .estimate_local_pi_ratio(
    ratio = LR_ini[initial_idx],
    init = mean(pi_bounds),
    pi_bounds = pi_bounds,
    tol = em_tol,
    max_iter = em_max_iter
  )
  pi_ini_hat <- initial_fit$pi
  Lmdr_ini <- .ratio_to_lmdr(LR_ini, pi_ini_hat)

  m <- length(LR)
  decision <- integer(m)
  Lmdr <- rep(NA_real_, m)
  pi_hat <- rep(NA_real_, m)
  capacity_path <- rep(NA_real_, m + 1L)
  gamma_path <- rep(NA_real_, m)
  em_iterations <- integer(m)
  em_converged <- logical(m)
  capacity <- 0
  capacity_path[1L] <- 0

  all_ratio <- c(LR_ini, LR)
  pi_start <- pi_ini_hat

  for (t in seq_len(m)) {
    current_position <- n_ini + t  # = K0 + t points available in X^{<= t}
    # Paper growing window: D_t = min{K0 + t, max{D_min, floor((K0 + t)^beta)}}.
    D_t <- if (D_mode == "growing")
      max(D_min, as.integer(floor(current_position^D_beta))) else D
    D_t <- min(D_t, current_position)
    # pi_t is estimated from the D_t most recent unlabeled observations in
    # X^{<= t}, INCLUDING the current point X_t (paper Stage 2 estimator).
    window_idx <- seq.int(current_position - D_t + 1L, current_position)

    local_fit <- .estimate_local_pi_ratio(
      ratio = all_ratio[window_idx],
      init = pi_start,
      pi_bounds = pi_bounds,
      tol = em_tol,
      max_iter = em_max_iter
    )
    current_pi <- local_fit$pi
    current_lmdr <- .ratio_to_lmdr(LR[t], current_pi)

    pi_hat[t] <- current_pi
    Lmdr[t] <- current_lmdr
    em_iterations[t] <- local_fit$iterations
    em_converged[t] <- local_fit$converged

    recent_scores <- utils::tail(c(Lmdr_ini, Lmdr[seq_len(t)]), w)
    current_gamma <- .omdrc_barrier_gamma(recent_scores, alpha)
    gamma_path[t] <- current_gamma

    if (current_lmdr < current_gamma &&
        (1 - alpha) * current_lmdr <= capacity) {
      decision[t] <- 0L
      capacity <- capacity - (1 - alpha) * current_lmdr
    } else {
      decision[t] <- 1L
      capacity <- capacity + alpha * current_lmdr
    }
    capacity_path[t + 1L] <- capacity
    pi_start <- current_pi
  }

  list(
    de = decision,
    Lmdr = Lmdr,
    Lmdr_ini = Lmdr_ini,
    pi_hat = pi_hat,
    pi_ini_hat = pi_ini_hat,
    LR = LR,
    LR_ini = LR_ini,
    capacity = capacity_path,
    gamma = gamma_path,
    barrier_window = w,
    em_iterations = em_iterations,
    em_converged = em_converged,
    DR = Lmdr,
    DR_ini = Lmdr_ini
  )
}

# -----------------------------------------------------------------------------
# Batch smooth-prevalence estimator (EM + penalized spline) -- OFFLINE benchmark
# -----------------------------------------------------------------------------
#
# Estimates pi_t = logit^{-1}(f(t)) as a smooth function of time from the
# unlabeled stream, given the density ratio r_s = f1(z_s)/f0(z_s). Uses ALL
# time points (batch, two-sided), so it is CONSISTENT at interior points for a
# smooth pi_t as T -> Inf with effective dof -> Inf and edf/T -> 0 (standard
# nonparametric regression asymptotics: at a fixed interior t the design
# density grows as data accrue on both sides). This is an OFFLINE benchmark --
# it uses future data relative to each t and is NOT a valid online estimator.
# It establishes the consistency ceiling for pi_t and the downstream MDR/FDR
# gain over the fixed-window online EM (which carries an irreducible tracking-
# error floor under drift).
#
# EM with GAM M-step (varying mixing proportion / "mixture-of-experts" EM):
#   E-step: w_s = pi(s) r_s / (1 - pi(s) + pi(s) r_s)
#   M-step: fit pi(t) = logit^{-1}(f(t)) by mgcv::gam(w ~ s(t), binomial)

#' Batch smooth prevalence estimator (offline benchmark)
#'
#' @param ratio density ratio r_s = f1(z_s)/f0(z_s) over the stream
#' @param t_idx time coordinate (default 1..length(ratio))
#' @param k basis dimension for s(t); choose large, the penalty controls edf
#' @param init initial pi (scalar or vector); default = mean(pi_bounds)
#' @param pi_bounds clipping bounds for pi
#' @param max_iter EM iterations
#' @param tol convergence: max|pi_new - pi_old| <= tol*(1+|pi|)
#' @param bs mgcv spline basis ("tp", "cr", "bs", ...)
#' @param method mgcv smoothing-parameter method ("REML" recommended)
#' @return list(pi_hat, fit, iterations, converged, edf)
#' @export
estimate_pi_gam_batch <- function(ratio, t_idx = seq_along(ratio),
                                  k = 40L, init = NULL,
                                  pi_bounds = c(0.005, 0.5),
                                  max_iter = 50L, tol = 1e-6,
                                  bs = "tp", method = "REML", fx = FALSE) {
  if (!requireNamespace("mgcv", quietly = TRUE)) {
    stop("estimate_pi_gam_batch requires the 'mgcv' package.")
  }
  ratio <- .clip_ratio(as.numeric(ratio), 1e-12, 1e12)
  Tn <- length(ratio)
  if (Tn < 5L) stop("Need at least 5 stream points for a smooth pi fit.")
  t_idx <- as.numeric(t_idx)
  if (length(t_idx) != Tn) stop("t_idx must have the same length as ratio.")
  pi_bounds <- .validate_pi_bounds(pi_bounds)
  # Aggregate by unique time index. The binomial M-step log-likelihood
  #   Sum_r [w_r log p(t_r) + (1-w_r) log(1-p(t_r))]
  # is EXACTLY a weighted binomial fit on (unique_t, mean_w) with weights =
  # counts. This makes the fit O(#unique t) regardless of replicates, so pooled
  # / repeated-t designs (e.g. consistency sweeps over R replicates) stay fast.
  t_uniq <- sort(unique(t_idx))
  n_u <- length(t_uniq)
  t_map <- match(t_idx, t_uniq)
  cnt <- as.numeric(tabulate(t_map, nbins = n_u))
  k <- max(4L, min(as.integer(k), max(4L, n_u - 1L)))
  df_u <- data.frame(t = t_uniq)

  if (is.null(init)) {
    p_u <- rep(mean(pi_bounds), n_u)
  } else if (length(init) == 1L) {
    p_u <- rep(init, n_u)
  } else if (length(init) == Tn) {
    p_u <- as.numeric(rowsum(as.numeric(init), group = t_map,
                             reorder = FALSE)) / cnt
  } else {
    p_u <- rep(mean(pi_bounds), n_u)
  }
  p_u <- pmin(pmax(p_u, pi_bounds[1L]), pi_bounds[2L])

  converged <- FALSE
  it <- 0L
  fit <- NULL
  for (iter in seq_len(as.integer(max_iter))) {
    it <- iter
    # E-step: per-observation responsibilities under current pi(t)
    p_obs <- p_u[t_map]
    w <- p_obs * ratio / pmax(1 - p_obs + p_obs * ratio, .Machine$double.xmin)
    w <- pmin(pmax(w, 0), 1)
    # aggregate to a weighted binomial response per unique t
    w_sum <- as.numeric(rowsum(w, group = t_map, reorder = FALSE))
    df_u$w <- w_sum / cnt
    # M-step: weighted smooth binomial GAM on unique t
    fit <- suppressWarnings(mgcv::gam(
      w ~ s(t, k = k, bs = bs, fx = fx), family = stats::binomial(),
      data = df_u, weights = cnt, method = method
    ))
    p_new_u <- as.numeric(stats::predict(fit, type = "response"))
    p_new_u <- pmin(pmax(p_new_u, pi_bounds[1L]), pi_bounds[2L])
    if (max(abs(p_new_u - p_u)) <= tol * (1 + max(abs(p_u)))) {
      p_u <- p_new_u
      converged <- TRUE
      break
    }
    p_u <- p_new_u
  }

  list(pi_hat = p_u[t_map], fit = fit, iterations = it, converged = converged,
       edf = if (!is.null(fit)) sum(fit$edf) else NA_real_,
       t_uniq = t_uniq, pi_hat_unique = p_u)
}

# -----------------------------------------------------------------------------
# Online one-sided local-linear EM estimator of pi_t (causal, best online tracker)
# -----------------------------------------------------------------------------
#
# At each t, fits a kernel-weighted local-LINEAR binomial mixture (EM + 2-param
# IRLS) over the PAST window [t-h, t-1] and reads pi_hat_t = logistic(a) at the
# boundary s=t. Local-LINEAR removes the O(h) one-sided lag bias of the local-
# constant windowed EM (leaving O(h^2)), giving a LOWER irreducible tracking
# floor -- the best online (causal) tracker. Still NOT consistent under sustained
# drift (Kalman steady-state floor > 0); consistent only where pi_t stops
# drifting. This is the honest "online best" -- NOT a consistent estimator.
#
#   E-step: w_s = pi(s) r_s / (1 - pi(s) + pi(s) r_s),  pi(s)=logistic(a + b(s-t))
#   M-step: weighted binomial IRLS of w_s on (s-t) with kernel weights K_h(t-s)

#' Online local-linear pi_t estimator (causal, best online tracker)
#'
#' @param ratio density ratio r_s = f1(z_s)/f0(z_s) over the stream
#' @param h kernel half-bandwidth (effective lookback ~ h points)
#' @param pi_bounds clipping bounds for pi
#' @param n_em EM iterations per time point
#' @param n_irls IRLS iterations per M-step
#' @return list(pi_hat = causal pi_t estimate, h = bandwidth)
#' @export
estimate_pi_online_ll <- function(ratio, h, pi_bounds = c(0.01, 0.99),
                                  n_em = 2L, n_irls = 6L) {
  ratio <- .clip_ratio(as.numeric(ratio), 1e-12, 1e12)
  Tn <- length(ratio)
  pi_bounds <- .validate_pi_bounds(pi_bounds)
  h <- as.numeric(h)
  if (!is.finite(h) || h < 2) stop("h must be >= 2.")
  logit_ <- stats::qlogis
  sig_   <- stats::plogis
  pi_hat <- rep(mean(pi_bounds), Tn)
  h2 <- as.integer(round(h))
  for (t in seq_len(Tn)) {
    if (t < 4L) next
    s <- max(1L, t - h2):(t - 1L)
    if (length(s) < 4L) { pi_hat[t] <- pi_hat[t - 1L]; next }
    dt <- as.numeric(t - s)            # in [1, h]
    u  <- dt / h
    K  <- pmax(0, (1 - u^2))^2         # biweight kernel (one-sided past only)
    sw <- sum(K)
    if (!is.finite(sw) || sw <= 0) { pi_hat[t] <- pi_hat[t - 1L]; next }
    K  <- K / sw
    rs <- ratio[s]
    xs <- as.numeric(s - t)            # < 0; local-linear covariate (0 = at t)
    a  <- logit_(pi_hat[t - 1L]); b <- 0
    for (em in seq_len(n_em)) {
      ps <- pmin(pmax(sig_(a + b * xs), 1e-6), 1 - 1e-6)
      w  <- ps * rs / pmax(1 - ps + ps * rs, .Machine$double.xmin)
      w  <- pmin(pmax(w, 1e-8), 1 - 1e-8)
      for (it in seq_len(n_irls)) {
        eta <- a + b * xs
        mu  <- pmin(pmax(sig_(eta), 1e-6), 1 - 1e-6)
        mu1 <- 1 - mu
        Wg  <- K * mu * mu1
        z   <- eta + (w - mu) / pmax(mu * mu1, 1e-8)
        S0  <- sum(Wg); Sx <- sum(Wg * xs); Sxx <- sum(Wg * xs * xs)
        Sy  <- sum(Wg * z);  Sxy <- sum(Wg * xs * z)
        det <- S0 * Sxx - Sx * Sx
        if (!is.finite(det) || abs(det) < 1e-14) break
        a <- (Sxx * Sy - Sx * Sxy) / det
        b <- (S0 * Sxy - Sx * Sy) / det
      }
    }
    pi_hat[t] <- pmin(pmax(sig_(a), pi_bounds[1L]), pi_bounds[2L])
  }
  list(pi_hat = pi_hat, h = h)
}

#' Data-Driven Online MDR Control with Time-Varying Prevalence
#'
#' Fits a chosen density-ratio estimator to labeled F0/F1 reference samples
#' (default: PC-DRE via a logistic GAM classifier, the paper's primary
#' implementation), estimates pi_t from the D_t most recent unlabeled
#' observations (growing window by default), and applies the capacity plus local
#' no-free-riding barrier rule from Algorithm 2.
#' @export
OMDRC_DD <- function(z, z_ini, z0, z1, alpha, D = NULL,
                     ratio_method = c("gam", "kde", "ulsif", "classifier"),
                     ratio_control = list(),
                     pi_bounds = c(0.01, 0.99),
                     em_tol = 1e-8, em_max_iter = 100L,
                     D_mode = c("growing", "fixed"),
                     D_beta = 0.6, D_min = 10L,
                     w = 100L) {
  ratio_method <- match.arg(ratio_method)
  z_mat <- .as_numeric_matrix(z, "z")
  z_ini_mat <- .as_numeric_matrix(z_ini, "z_ini")
  z0_mat <- .as_numeric_matrix(z0, "z0")
  z1_mat <- .as_numeric_matrix(z1, "z1")
  dims <- c(ncol(z_mat), ncol(z_ini_mat), ncol(z0_mat), ncol(z1_mat))
  if (length(unique(dims)) != 1L) stop("z, z_ini, z0, and z1 must have the same dimension.")

  fit_args <- c(
    list(x0 = z0_mat, x1 = z1_mat),
    ratio_control
  )
  if (ratio_method == "kde" && is.null(fit_args$grid_data)) {
    fit_args$grid_data <- z_ini_mat
  }
  ratio_model <- switch(
    ratio_method,
    kde = do.call(fit_ratio_kde, fit_args),
    ulsif = do.call(fit_ratio_ulsif, fit_args),
    classifier = do.call(fit_ratio_classifier, fit_args),
    gam = do.call(fit_ratio_gam, fit_args)
  )

  LR_ini <- predict_ratio(ratio_model, z_ini_mat)
  LR <- predict_ratio(ratio_model, z_mat)
  out <- OMDRC_FROM_RATIO(
    LR = LR,
    LR_ini = LR_ini,
    alpha = alpha,
    D = D,
    pi_bounds = pi_bounds,
    em_tol = em_tol,
    em_max_iter = em_max_iter,
    D_mode = D_mode,
    D_beta = D_beta,
    D_min = D_min,
    w = w
  )
  out$ratio_method <- ratio_method
  out$ratio_model <- ratio_model
  out
}

# -----------------------------------------------------------------------------
# Practical score-based baselines
# -----------------------------------------------------------------------------

.calibrate_initial_ranking <- function(scores, alpha) {
  scores <- pmin(pmax(as.numeric(scores), 0), 1)
  if (length(scores) == 0L || any(!is.finite(scores))) {
    stop("Initial Lmdr scores must be finite and non-empty.")
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0 || alpha >= 1) {
    stop("alpha must be a scalar in (0, 1).")
  }

  total_mass <- sum(scores)
  n_cal <- length(scores)
  if (!is.finite(total_mass) || total_mass <= .Machine$double.eps) {
    return(list(threshold = -Inf, k = n_cal, rho = 1, estimated_mdr = 0))
  }

  sorted_scores <- sort(scores, decreasing = TRUE)
  estimated_mdr_by_k <- pmax(0, 1 - cumsum(sorted_scores) / total_mass)
  k <- which(estimated_mdr_by_k <= alpha)[1L]
  if (is.na(k)) k <- n_cal
  threshold <- sorted_scores[k]
  initial_decision <- as.integer(scores >= threshold)

  list(
    threshold = threshold,
    k = as.integer(k),
    rho = k / n_cal,
    estimated_mdr = sum(scores * (1 - initial_decision)) / total_mass
  )
}

#' Static calibration threshold (SCT)
#' @export
STATIC_LMDR_DD <- function(Lmdr, Lmdr_ini, alpha) {
  calibration <- .calibrate_initial_ranking(Lmdr_ini, alpha)
  list(
    de = as.integer(Lmdr >= calibration$threshold),
    threshold = calibration$threshold,
    k_ini = calibration$k,
    rho_ini = calibration$rho,
    estimated_mdr_ini = calibration$estimated_mdr
  )
}

#' Causal rolling Top-k policy (RTK)
#' @export
ROLLING_TOPK_DD <- function(Lmdr, Lmdr_ini, alpha,
                            window_size = length(Lmdr_ini), k = NULL) {
  if (!is.numeric(window_size) || length(window_size) != 1L ||
      !is.finite(window_size) || window_size < 1) {
    stop("window_size must be a positive scalar integer.")
  }
  window_size <- min(as.integer(window_size), length(Lmdr_ini))
  calibration <- .calibrate_initial_ranking(tail(Lmdr_ini, window_size), alpha)
  if (is.null(k)) k <- calibration$k
  if (!is.numeric(k) || length(k) != 1L || !is.finite(k) || k < 1) {
    stop("k must be a positive scalar integer.")
  }
  k <- min(as.integer(k), window_size)

  decision <- integer(length(Lmdr))
  history <- if (window_size > 1L) tail(Lmdr_ini, window_size - 1L) else numeric(0)
  for (t in seq_along(Lmdr)) {
    score_window <- tail(c(history, Lmdr[t]), window_size)
    current_rank <- rank(-score_window, ties.method = "last")[length(score_window)]
    decision[t] <- as.integer(current_rank <= min(k, length(score_window)))
    history <- score_window
  }

  list(
    de = decision,
    k = k,
    rho = k / window_size,
    window_size = window_size,
    estimated_mdr_ini = calibration$estimated_mdr
  )
}

# -----------------------------------------------------------------------------
# Oracle procedures
# -----------------------------------------------------------------------------

#' Oracle online MDR control (Algorithm 1, including the local barrier)
#' @export
OMDRC_OR <- function(x.Lmdr, alpha, w = 100L, x.Lmdr_ini = NULL) {
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop("alpha must be a scalar in (0, 1).")
  }
  x.Lmdr <- pmin(pmax(as.numeric(x.Lmdr), 0), 1)
  x.Lmdr_ini <- if (is.null(x.Lmdr_ini)) numeric(0) else
    pmin(pmax(as.numeric(x.Lmdr_ini), 0), 1)
  w <- .validate_barrier_window(w)
  delta <- integer(length(x.Lmdr))
  gamma <- numeric(length(x.Lmdr))
  capacity_path <- numeric(length(x.Lmdr) + 1L)
  capacity <- 0
  for (t in seq_along(x.Lmdr)) {
    gamma[t] <- .omdrc_barrier_gamma(
      utils::tail(c(x.Lmdr_ini, x.Lmdr[seq_len(t)]), w), alpha
    )
    if (x.Lmdr[t] < gamma[t] &&
        (1 - alpha) * x.Lmdr[t] <= capacity) {
      delta[t] <- 0L
      capacity <- capacity - (1 - alpha) * x.Lmdr[t]
    } else {
      delta[t] <- 1L
      capacity <- capacity + alpha * x.Lmdr[t]
    }
    capacity_path[t + 1L] <- capacity
  }
  list(de = delta, gamma = gamma, capacity = capacity_path,
       barrier_window = w)
}

#' Offline oracle fixed threshold (FT)
#' @export
OMDRC_OFF <- function(x.Lmdr, alpha) {
  x <- pmin(pmax(as.numeric(x.Lmdr), 0), 1)
  total_mass <- sum(x)
  if (!is.finite(total_mass) || total_mass <= .Machine$double.eps) {
    return(list(lambda = -Inf, de = rep(1L, length(x))))
  }
  xs <- sort(x)
  feasible <- which(cumsum(xs) / total_mass <= alpha)
  k <- if (length(feasible) == 0L) 0L else max(feasible)
  if (k == 0L) {
    lambda <- -Inf
    de <- rep(1L, length(x))
  } else {
    lambda <- xs[k]
    de <- as.integer(x > lambda)
  }
  list(lambda = lambda, de = de)
}
