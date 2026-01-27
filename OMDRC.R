# ==============================================================================
# MODULE: Online Missed Discovery Rate Control (OMDRC) Functions
# Description: This file contains implementation for Data-Driven, Oracle, 
#              and Offline algorithms for controlling the Online MDR.
# ==============================================================================

library(kedd)

#' Data-Driven Online MDR Control (OMDRC-DD)
#'
#' Implements the data-driven online MDR control algorithm (Algorithm 2 in the paper).
#' This procedure operates on a stream of observations, making real-time decisions
#' by estimating the density ratio using Kernel Density Estimation (KDE) on a 
#' sliding window of recent data.
#'
#' @param z A numeric vector representing the incoming data stream $\{X_t\}$.
#' @param z_ini A numeric vector of initial unlabeled observations used to bootstrap 
#'   the mixture density estimation.
#' @param z1 A numeric vector of labeled samples from the alternative distribution 
#'   ($F_1$), representing prior knowledge of anomalies.
#' @param alpha A numeric value in (0, 1) specifying the target MDR level.
#' @param D An integer specifying the size of the sliding window used for estimating 
#'   the mixture density $f$.
#' @param M A numeric value for truncating the estimated density ratio to ensure
#'   numerical stability. Defaults to 10.
#' @param grid_n An integer for the number of grid points used in the KDE computation.
#'   Defaults to 1000.
#' @param pad A numeric value to extend the range of the KDE grid, preventing boundary
#'   effects. Defaults to 5.
#'
#' @return A list containing:
#' \item{de}{A binary vector of decisions for each point in `z` (1 for discovery, 0 for non-discovery).}
#' \item{DR}{A numeric vector of the estimated density ratios ($\widehat{\mathrm{DR}}_t$) for each observation.}
#' \item{h1}{The bandwidth calculated for the alternative density ($f_1$) estimation.}
#' \item{from0}{The lower bound of the grid used for density estimation.}
#' \item{to0}{The upper bound of the grid used for density estimation.}
#' @export
OMDRC_DD <- function(z, z_ini, z1, alpha, D,M=10, 
                     grid_n = 1000, pad = 5) {
  
  m <- length(z)
  decision <- integer(m)   
  DR <- rep(NA_real_, m)
  capacity <- 0 # Initial testing budget
  
  # Check if data is non-negative to apply reflection at boundary
  is_positive_only <- all(c(z_ini, z1) >= 0)
  
  all_data <- c(z, z_ini, z1)
  from0 <- if(is_positive_only) 0 else min(all_data) - pad
  to0   <- max(all_data) + pad
  
  # Bandwidth selection for the alternative density (f1)
  h1 <- bw.SJ(z1) * 1.5
  
  # Helper: KDE with optional boundary reflection
  get_kde <- function(x_data, h, from, to, n, reflect = FALSE) {
    if (reflect) {
      d <- density(c(x_data, -x_data), bw = h, from = from, to = to, n = n)
      d$y <- d$y * 2
      d$y[d$x < 0] <- 0
    } else {
      d <- density(x_data, bw = h, from = from, to = to, n = n)
    }
    return(d)
  }
  
  # Pre-calculate alternative density estimation (f1_hat)
  den1 <- get_kde(z1, h1, from0, to0, grid_n, reflect = is_positive_only)
  
  for (t in seq_len(m)) {
    # Extract current window for mixture density (f_hat) estimation
    hist_t <- if (t > 1) z[1:(t-1)] else numeric(0)
    tmp_x <- tail(c(z_ini, hist_t), D)
    
    # Bandwidth selection for the mixture density
    h0_t <- bw.nrd0(tmp_x) 
    den0 <- get_kde(tmp_x, h0_t, from0, to0, grid_n, reflect = is_positive_only)
    
    # Linear interpolation of densities at the current observation z[t]
    f1hat <- approx(den1$x, den1$y, xout = z[t], rule = 2)$y
    fhat  <- approx(den0$x, den0$y, xout = z[t], rule = 2)$y
    
    # Numerical stability
    fhat <- max(fhat, 1e-6)
    
    # Estimate the local density ratio (Local MDR)
    current_DR <- min(f1hat / fhat, M)
    DR[t] <- current_DR
    
    # Testing Budget (Capacity) Update Rule
    if ((1 - alpha) * current_DR <= capacity) {
      decision[t] <- 0L # Non-discovery
      capacity <- capacity - (1 - alpha) * current_DR
    } else {
      decision[t] <- 1L # Discovery (Signal detected)
      capacity <- capacity + alpha * current_DR
    }
  }
  
  list(de = decision, DR = DR, h1 = h1, from0 = from0, to0 = to0)
}

#' Linear Interpolation Helper
#'
#' Performs linear interpolation for density estimation on a discrete grid.
#'
#' @param x Coordinates where the density needs to be interpolated.
#' @param X Grid coordinates of the estimated densities.
#' @param Y Values of the estimated densities on the grid.
#' @return Numeric vector of interpolated density values.
#' @export
lin.itp <- function(x, X, Y){
  x.N <- length(x)
  X.N <- length(X)
  y <- rep(0, x.N)
  for (k in 1:x.N){
    i <- max(which((x[k]-X)>=0))
    if (i < X.N)
      y[k] <- Y[i] + (Y[i+1]-Y[i])/(X[i+1]-X[i])*(x[k]-X[i])
    else
      y[k] <- Y[i]
  }
  return(y)
}

#' Oracle Online MDR Control (OMDRC-OR)
#'
#' Implementation of the Oracle Online MDR procedure. This algorithm assumes 
#' the true Local MDR values (posterior probabilities) are known a priori.
#'
#' @param x_lmdr Numeric vector. True Local MDR values (posterior probability of being null).
#' @param alpha Numeric. Target MDR control level.
#'
#' @return A list containing:
#' \item{de}{Binary decisions (1 for discovery, 0 otherwise).}
#' \item{Chat}{Historical trace of the testing budget (capacity).}
#' @export
OMDRC_OR <- function(x_lmdr, alpha) {
  m <- length(x_lmdr)
  delta <- integer(m)   
  capacity_history <- numeric(m)
  capacity <- 0
  
  for (t in seq_len(m)) {
    capacity_history[t] <- capacity
    # Decision rule based on current testing budget
    if ((1 - alpha) * x_lmdr[t] <= capacity) {
      delta[t] <- 0L
      capacity <- capacity - (1 - alpha) * x_lmdr[t]
    } else {
      delta[t] <- 1L
      capacity <- capacity + alpha * x_lmdr[t]
    }
  }
  list(de = delta, Chat = capacity_history)
}

#' Offline MDR Control (OMDRC-OFF)
#'
#' Implementation of the offline MDR control procedure. Computes an optimal 
#' fixed threshold based on the complete batch of observations.
#'
#' @param x.Lmdr Numeric vector. Local MDR values for the batch.
#' @param alpha Numeric. Target MDR control level.
#'
#' @return A list containing:
#' \item{lambda}{The optimal fixed threshold.}
#' \item{de}{Binary decisions based on the threshold.}
#' @export
OMDRC_OFF <- function(x.Lmdr, alpha){
  x <- pmin(pmax(x.Lmdr, 0), 1)
  m <- length(x.Lmdr)
  
  # Sort Local MDR values and compute empirical MDR
  o <- order(x)
  xs <- x[o]
  cs <- cumsum(xs)
  mmdr <- cs / sum(x)
  
  # Find the largest index satisfying the MDR constraint
  k <- max(which(mmdr <= alpha), 0)
  
  if (k == 0) {
    lambda <- -Inf
    de <- rep(1L, m)
  } else {
    lambda <- xs[k]
    de <- as.integer(x > lambda)
  }
  
  list(lambda = lambda, de = de)
}
