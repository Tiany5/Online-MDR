library(kedd)

#' Data-Driven Online MDR Control (OMDRC-DD)
#'
#' This function implements the data-driven online Missed Discovery Rate control 
#' using Kernel Density Estimation (KDE) with adaptive bandwidth selection.
#'
#' @param z A numeric vector of the observed data stream.
#' @param z_ini A numeric vector of the initial unlabeled batch samples.
#' @param z1 A numeric vector of the labeled alternative samples (from f1).
#' @param alpha The target MDR control level.
#' @param D The sliding window size for the mixture density estimation.
#' @param grid_n Number of points in the KDE grid. Default is 1000.
#' @param pad Numerical padding for the density estimation range. Default is 5.
#'
#' @return A list containing:
#' \item{de}{Binary decisions (1 for rejection, 0 for non-rejection).}
#' \item{DR}{Estimated density ratios (Lmdr) for each observation.}
#' \item{h1}{Bandwidth used for the alternative density estimation.}
#' \item{from0}{The starting point of the estimation grid.}
#' \item{to0}{The end point of the estimation grid.}
#' @export
OMDRC_DD <- function(z, z_ini, z1, alpha, D, 
                     grid_n = 1000, pad = 5) {
  
  m <- length(z)
  decision <- integer(m)   
  DR <- rep(NA_real_, m)
  capacity <- 0
  
  is_positive_only <- all(c(z_ini, z1) >= 0)
  
  all_data <- c(z, z_ini, z1)
  from0 <- if(is_positive_only) 0 else min(all_data) - pad
  to0   <- max(all_data) + pad
  
  h1 <- bw.SJ(z1) * 1.5
  
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
  
  den1 <- get_kde(z1, h1, from0, to0, grid_n, reflect = is_positive_only)
  
  for (t in seq_len(m)) {
    hist_t <- if (t > 1) z[1:(t-1)] else numeric(0)
    tmp_x <- tail(c(z_ini, hist_t), D)
    
    h0_t <- bw.nrd0(tmp_x) 
    den0 <- get_kde(tmp_x, h0_t, from0, to0, grid_n, reflect = is_positive_only)
    
    f1hat <- approx(den1$x, den1$y, xout = z[t], rule = 2)$y
    fhat  <- approx(den0$x, den0$y, xout = z[t], rule = 2)$y
    
    fhat <- max(fhat, 1e-6)
    
    current_DR <- min(f1hat / fhat, 1/alpha)
    DR[t] <- current_DR
    
    if ((1 - alpha) * current_DR <= capacity) {
      decision[t] <- 0L
      capacity <- capacity - (1 - alpha) * current_DR
    } else {
      decision[t] <- 1L
      capacity <- capacity + alpha * current_DR
    }
  }
  
  list(de = decision, DR = DR, h1 = h1, from0 = from0, to0 = to0)
}

#' Linear Interpolation
#'
#' This function returns the linearly interpolated densities.
#'
#' @param x the coordinates of points where the density needs to be interpolated
#' @param X the coordinates of the estimated densities
#' @param Y the values of the estimated densities
#' @return A numeric vector of interpolated values.
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
#' Implements the oracle online MDR control procedure using true Local MDR values.
#'
#' @param x_lmdr Numeric vector. The true Local MDR values (posterior probabilities).
#' @param alpha Numeric. The target MDR control level.
#'
#' @return A list containing binary decisions (\code{de}) and capacity history (\code{Chat}).
#' @export
OMDRC_OR <- function(x_lmdr, alpha) {
  m <- length(x_lmdr)
  delta <- integer(m)   
  capacity_history <- numeric(m)
  capacity <- 0
  
  for (t in seq_len(m)) {
    capacity_history[t] <- capacity
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
#' This function implements the offline MDR control by finding an optimal fixed threshold.
#'
#' @param x.Lmdr A numeric vector of Lmdr values.
#' @param alpha The target MDR control level.
#'
#' @return A list containing:
#' \item{lambda}{The calculated fixed threshold.}
#' \item{de}{Binary decisions based on the threshold.}
#' @export
OMDRC_OFF <- function(x.Lmdr, alpha){
  x <- pmin(pmax(x.Lmdr, 0), 1)
  m <- length(x.Lmdr)
  
  o <- order(x)
  xs <- x[o]
  cs <- cumsum(xs)
  mmdr <- cs / sum(x)
  
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
