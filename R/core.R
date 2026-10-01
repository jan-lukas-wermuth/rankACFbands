# Internal validation deliberately rejects missing values: dropping a row changes
# the time distance between observations and hence the estimand.
.check_x <- function(x) {
  if (is.ordered(x)) x <- as.integer(x)
  if (!is.numeric(x) || !is.null(dim(x)) || length(x) < 3L ||
      anyNA(x) || any(!is.finite(x)))
    stop("x must be a finite univariate numeric series or ordered factor with at least 3 observations; missing values cannot be dropped without changing lags.", call. = FALSE)
  as.numeric(x)
}
.integer_scalar <- function(x, name, lower, upper = Inf) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x != floor(x) || x < lower || x > upper)
    stop(name, " must be an integer in [", lower, ", ", upper, "].", call. = FALSE)
  as.integer(x)
}
.check_H <- function(lag.max, n) {
  if (is.null(lag.max)) lag.max <- min(n - 2L, max(1L, floor(10 * log10(n))))
  .integer_scalar(lag.max, "lag.max", 1L, n - 2L)
}
.methods <- c("rho_b", "tau_b", "gamma", "rho", "tau")
.check_method <- function(method) match.arg(method, .methods)

# Mid-CDF and masses at the observations, using global order categories.
.marginal <- function(x) {
  id <- match(x, sort(unique(x)))
  count <- tabulate(id)
  grade <- (cumsum(count) - count / 2) / length(x)
  list(id = id, count = count, p = count[id] / length(x),
       grade = grade[id], K = length(count))
}

# Bivariate mid-distribution at the observed pairs. Fenwick prefix sums
# give O(n log K) time and O(n + K) memory, including tied coordinates.
.pair_mid <- function(xid, yid, K) {
  n <- length(xid)
  tree <- numeric(K)
  out <- numeric(n)
  ord <- order(xid)
  groups <- split(ord, xid[ord])
  for (ii in groups) {
    for (q in ii) {
      j <- yid[q]; total <- 0
      while (j > 0L) { total <- total + tree[j]; j <- j - bitwAnd(j, -j) }
      j <- yid[q] - 1L; below <- 0
      while (j > 0L) { below <- below + tree[j]; j <- j - bitwAnd(j, -j) }
      out[q] <- (total + below) / 2
    }
    # Equal first-coordinate observations contribute half their y mid-CDF.
    uy <- sort(unique(yid[ii]))
    cty <- tabulate(match(yid[ii], uy))
    out[ii] <- out[ii] + (cumsum(cty) - cty / 2)[match(yid[ii], uy)] / 2
    for (q in ii) {
      j <- yid[q]
      while (j <= K) { tree[j] <- tree[j] + 1; j <- j + bitwAnd(j, -j) }
    }
  }
  out / n
}
.weight_mid <- function(id, weight, K) {
  totals <- numeric(K)
  agg <- rowsum(weight, id, reorder = FALSE)
  totals[as.integer(rownames(agg))] <- agg[, 1L]
  (cumsum(totals) - totals / 2) / length(id)
}
.pair_stats <- function(xid, yid, K) {
  n <- length(xid)
  counts <- DoCount(xid, yid)
  total <- choose(n, 2)
  tau <- (counts$C - counts$D) / total
  nu  <- (total - counts$C - counts$D) / total
  # Still needed for the estimated first-order projections.
  joint <- .pair_mid(xid, yid, K)
  keys <- paste(xid, yid, sep = ":")
  mass <- as.numeric(table(keys)[keys]) / n
  list(tau = tau, nu = nu, joint = joint, mass = mass)
}
.compute <- function(x, H, projections = FALSE) {
  n <- length(x); marginal <- .marginal(x)
  rho <- tau <- nu <- numeric(H + 1L)
  rho[1L] <- 12 * mean((marginal$grade - 0.5)^2)
  tau[1L] <- 1 - sum(marginal$count * (marginal$count - 1)) / (n * (n - 1))
  nu[1L] <- 1 - tau[1L]
  pr <- if (projections) list() else NULL
  if (projections) {
    pr[["rho.0"]] <- 1 - marginal$p^2 - rho[1L]
    pr[["tau.0"]] <- 1 - marginal$p - tau[1L]
  }
  for (h in seq_len(H)) {
    i <- seq_len(n - h); j <- i + h
    ps <- .pair_stats(marginal$id[i], marginal$id[j], marginal$K)
    rho[h + 1L] <- 12 / n * sum((marginal$grade[i] - 0.5) * (marginal$grade[j] - 0.5))
    tau[h + 1L] <- ps$tau; nu[h + 1L] <- ps$nu
    if (projections) {
      # E[step_mid(Y-v)] = 1-G(v), with Y drawn from the full empirical marginal.
      g1 <- .weight_mid(marginal$id[i], 1 - marginal$grade[j], marginal$K)
      g2 <- .weight_mid(marginal$id[j], 1 - marginal$grade[i], marginal$K)
      pr[[paste0("rho.", h)]] <- 4 * (g1[marginal$id[i]] + g2[marginal$id[j]] +
        marginal$grade[i] * marginal$grade[j] - marginal$grade[i] - marginal$grade[j]) + 1 - rho[h + 1L]
      pr[[paste0("tau.", h)]] <- 4 * ps$joint - 2 * (marginal$grade[i] + marginal$grade[j]) + 1 - tau[h + 1L]
      pr[[paste0("nu.", h)]] <- marginal$p[i] + marginal$p[j] - ps$mass - nu[h + 1L]
    }
  }
  list(rho = rho, tau = tau, nu = nu, projections = pr, marginal = marginal)
}
.estimates <- function(raw, method) {
  if (method == "rho") return(raw$rho)
  if (method == "tau") return(raw$tau)
  if (method == "rho_b") {
    if (raw$rho[1L] <= 0) stop("Grade autocorrelation is undefined for a constant series.", call. = FALSE)
    return(raw$rho / raw$rho[1L])
  }
  if (method == "tau_b") {
    if (raw$tau[1L] <= 0) stop("Normalized Kendall autocorrelation is undefined for a constant series.", call. = FALSE)
    return(raw$tau / raw$tau[1L])
  }
  if (any(1 - raw$nu <= 1e-14))
    stop("Gamma is undefined at a requested lag: there are no pairs untied in both coordinates.", call. = FALSE)
  raw$tau / (1 - raw$nu)
}

#' Estimate a Rank Autocorrelation Function
#'
#' Computes the five estimators in Definitions 4 and 5 of the manuscript.
#' @param x Finite numeric vector, univariate time series, or ordered factor.
#'   Missing observations are rejected to preserve time spacing.
#' @param lag.max Largest lag, between 1 and length(x)-2. The default is
#'   min(length(x)-2, floor(10*log10(length(x)))).
#' @param method One of "rho_b" (grade), "tau_b" (normalized Kendall),
#'   "gamma", "rho" (classical Spearman), or "tau" (classical Kendall).
#' @param plot Logical; draw the correlogram.
#' @param ... Additional arguments passed to the plot method.
#' @return An object of class \code{rank_acf} with \code{acf} (lags 0 through H),
#'   \code{estimate} (lags 1 through H), \code{lag}, \code{n}, \code{method},
#'   \code{raw} (rho, tau and nu estimates), \code{ties}, and \code{call}.
#' @details Rho uses full-sample midranks and denominator T, not T-h.
#' Tau and nu use U-statistics on T-h lagged pairs. The denominators of
#' rho_b and tau_b are the full-series lag-zero estimates; therefore tau_b
#' is not cor(x[1:(T-h)], x[(h+1):T], method="kendall"). Finite-sample tau_b
#' values can exceed [-1,1] because the denominator uses a different sample.
#' Gamma excludes pairs tied in either coordinate. Ties use exact equality.
#' Computation uses O(H*T*log(T)) time and O(T) working memory for estimates.
#' @references Pohle, Weiss and Wermuth (2026), Rank Autocorrelations,
#'   September 29 draft, Definitions 4 and 5.
#' @examples
#' set.seed(11)
#' x <- rpois(150, 3)
#' fit <- rank_acf(x, lag.max = 8, method = "gamma")
#' coef(fit)
#' plot(fit)
#' @export
rank_acf <- function(x, lag.max = NULL, method = c("rho_b", "tau_b", "gamma", "rho", "tau"),
                     plot = FALSE, ...) {
  cl <- match.call(); x <- .check_x(x); H <- .check_H(lag.max, length(x))
  method <- .check_method(method)
  raw <- .compute(x, H)
  estimates <- .estimates(raw, method)
  names(estimates) <- as.character(0:H)
  result <- structure(list(acf = estimates, estimate = estimates[-1L], lag = seq_len(H),
    n = length(x), method = method, raw = raw[c("rho", "tau", "nu")],
    ties = rank_ties(x), call = cl), class = "rank_acf")
  if (plot) graphics::plot(result, ...)
  result
}

#' Empirical Tie Probabilities
#' @inheritParams rank_acf
#' @return Named vector: zeta is the sum of squared relative frequencies;
#'   zeta2 is the sum of cubed relative frequencies. These are plug-in estimates,
#'   not the unbiased U-statistic estimates used in tau(0) and nu(0).
#' @examples
#' rank_ties(c(1, 1, 2, 3))
#' @export
rank_ties <- function(x) {
  x <- .check_x(x)
  p <- .marginal(x)$count / length(x)
  c(zeta = sum(p^2), zeta2 = sum(p^3))
}

#' Estimated Hoeffding Projections for Rank Autocorrelations
#' @inheritParams rank_acf
#' @return Named list containing rho.h, tau.h and nu.h for h=0,...,H.
#'   Each numeric vector has length T-h and is evaluated at (X_t,X_(t+h)).
#'   Values are the literal plug-in projections in Section 4.2; they are not
#'   additionally sample-centered and do not include the U-statistic orders.
#' @examples
#' rank_projections(c(1, 2, 1, 3, 4, 2), lag.max = 2)
#' @export
rank_projections <- function(x, lag.max = NULL) {
  x <- .check_x(x); H <- .check_H(lag.max, length(x))
  .compute(x, H, projections = TRUE)$projections
}

#' Independence Asymptotic Variance
#' @inheritParams rank_acf
#' @return Scalar variance of sqrt(T) times the estimator under iid sampling.
#'   Divide by T for the sampling variance. Generalized measures require a
#'   nonconstant series.
#' @details Implements Corollaries 1 and 2 using empirical tie probabilities.
#' @examples
#' rank_iid_variance(c(1, 2, 1, 3, 2, 4), "rho_b")
#' @export
rank_iid_variance <- function(x, method = c("rho_b", "tau_b", "gamma", "rho", "tau")) {
  method <- .check_method(method); ties <- rank_ties(x)
  if (ties["zeta"] >= 1) stop("Inference requires a nonconstant series.", call. = FALSE)
  r <- unname(1 - ties["zeta2"]); t <- unname(1 - ties["zeta"])
  switch(method, rho = r^2, tau = 4 / 9 * r^2, rho_b = 1,
         tau_b = 4 / 9 * r^2 / t^2, gamma = 4 / 9 * r^2 / t^4)
}

#' @keywords internal
DoCount <- function (y, x, wts)
{
  if (missing(wts))
    wts <- rep_len(1L, length(x))
  ord <- order(y)
  ux <- sort(unique(x))
  n2 <- length(ux)
  idx <- DescTools::BinTree(n2)[match(x[ord], ux)] - 1L
  y <- cbind(y, 1)
  res <- .Call("conc", y[ord, ], as.double(wts[ord]),
               as.integer(idx), as.integer(n2), PACKAGE = "DescTools")
  return(list(pi.c = NA, pi.d = NA, C = res[2], D = res[1],
              T = res[3], N = res[4]))
}
