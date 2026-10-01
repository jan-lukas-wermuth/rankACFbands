.psd <- function(x, policy = "error") {
  x <- (x + t(x)) / 2
  e <- eigen(x, symmetric = TRUE)
  scale <- max(abs(e$values), .Machine$double.eps)
  bad <- min(e$values) < -1e-10 * scale
  if (bad && policy == "error") stop("Estimated covariance is not positive semidefinite; use psd='clip' to request eigenvalue truncation.", call. = FALSE)
  if (policy == "none") return(list(matrix = x, min.eigenvalue = min(e$values),
    adjusted = FALSE, relative.adjustment = 0))
  if (any(e$values < 0)) {
    fixed <- tcrossprod(sweep(e$vectors, 2L, sqrt(pmax(e$values, 0)), "*"))
    dimnames(fixed) <- dimnames(x)
    change <- sqrt(sum((fixed - x)^2)) / max(sqrt(sum(x^2)), .Machine$double.eps)
    if (bad) warning("HAC covariance has negative eigenvalues; truncated to zero. Inspect raw_covariance and psd diagnostics.", call. = FALSE)
    return(list(matrix = fixed, min.eigenvalue = min(e$values),
                adjusted = TRUE, relative.adjustment = change))
  }
  list(matrix = x, min.eigenvalue = min(e$values), adjusted = FALSE, relative.adjustment = 0)
}
.hac_compute <- function(raw, n, H, method, bandwidth, kernel, center, psd) {
  if (is.null(bandwidth)) bandwidth <- floor(2 * n^(1 / 3))
  bandwidth <- .integer_scalar(bandwidth, "bandwidth", 0L, n - 1L)
  kernel <- match.arg(kernel, c("bartlett", "parzen", "tukey-hanning"))
  weights <- if (bandwidth == 0) numeric() else {
    u <- seq_len(bandwidth) / (bandwidth + 1)
    switch(kernel, bartlett = 1 - u,
      parzen = ifelse(u <= 0.5, 1 - 6 * u^2 + 6 * u^3, 2 * (1 - u)^3),
      `tukey-hanning` = (1 + cos(pi * u)) / 2)
  }
  if (method %in% c("rho", "tau")) {
    kinds <- rep(method, H); lags <- seq_len(H)
    jac <- diag(H)
  } else if (method %in% c("rho_b", "tau_b")) {
    base <- if (method == "rho_b") "rho" else "tau"
    kinds <- rep(base, H + 1L); lags <- 0:H
    denominator <- raw[[base]][1L]
    jac <- cbind(-raw[[base]][-1L] / denominator^2, diag(H) / denominator)
  } else {
    kinds <- rep(c("tau", "nu"), H); lags <- rep(seq_len(H), each = 2)
    jac <- matrix(0, H, 2L * H)
    for (h in seq_len(H)) {
      jac[h, 2L * h - 1L] <- 1 / (1 - raw$nu[h + 1L])
      jac[h, 2L * h] <- raw$tau[h + 1L] / (1 - raw$nu[h + 1L])^2
    }
  }
  labels <- paste(kinds, lags, sep = ".")
  scores <- raw$projections[labels]
  if (center) scores <- lapply(scores, function(v) v - mean(v))
  orders <- ifelse(kinds == "rho", 3, 2)
  covariance <- matrix(0, length(scores), length(scores), dimnames = list(labels, labels))
  for (a in seq_along(scores)) for (b in seq_len(a)) {
    hstar <- max(lags[a], lags[b])
    ii <- seq_len(n - hstar)
    val <- sum(scores[[a]][ii] * scores[[b]][ii])
    if (bandwidth > 0) for (s in seq_len(min(bandwidth, n - hstar - 1L))) {
      ii <- seq_len(n - s - hstar)
      val <- val + weights[s] * (sum(scores[[a]][ii] * scores[[b]][ii + s]) +
                                  sum(scores[[b]][ii] * scores[[a]][ii + s]))
    }
    covariance[a, b] <- covariance[b, a] <- orders[a] * orders[b] * val / n
  }
  dimnames(jac) <- list(as.character(seq_len(H)), labels)
  transformed <- jac %*% covariance %*% t(jac)
  adjusted <- .psd(transformed, psd)
  list(covariance = adjusted$matrix, raw_covariance = transformed,
       primitive_covariance = covariance, jacobian = jac, bandwidth = bandwidth,
       kernel = kernel, center = center, psd = adjusted[-1L], projections = scores)
}

#' HAC Covariance Matrix for Rank Autocorrelations
#' @inheritParams rank_acf
#' @param bandwidth Nonnegative integer HAC bandwidth; defaults to
#'   min(T-1, floor(2*T^(1/3))). Zero retains contemporaneous products only;
#'   it is not the iid null covariance.
#' @param kernel HAC window: "bartlett", "parzen", or "tukey-hanning".
#' @param center Whether to additionally center each estimated projection.
#'   Default FALSE follows the manuscript literally.
#' @param psd Treatment of negative eigenvalues after the delta method:
#'   "clip" truncates them to zero (warning for material changes), "error"
#'   rejects material indefiniteness, "none" returns the unmodified estimate.
#' @return A list of class \code{rank_hac} containing \code{covariance} (the
#'   asymptotic covariance), \code{vcov} (covariance/T), \code{raw_covariance},
#'   \code{primitive_covariance}, \code{jacobian}, \code{projections},
#'   \code{psd} diagnostics, and estimator settings.
#' @details The empirical cross-products use denominator T and exactly
#'   T-s-max(h,h') observations for each primitive pair, as in Section 4.2.
#'   Both positive and negative HAC lags are included. U-statistic orders
#'   (3 for rho, 2 for tau and nu) multiply the entire covariance sum.
#'   Pair-specific truncation may yield a finite-sample indefinite matrix,
#'   even with a Bartlett window; the raw matrix is always retained.
#'   The asymptotic guarantee assumes fixed lag.max and a bandwidth growing
#'   more slowly than sqrt(T). Covariance computation costs O(H^2*T*bandwidth).
#' @references Pohle, Weiss and Wermuth (2026), Section 4.2 and Proposition 2.
#' @examples
#' set.seed(2)
#' x <- as.numeric(arima.sim(list(ar = 0.4), n = 150))
#' v <- rank_hac(x, lag.max = 4, method = "rho_b")
#' v$covariance
#' @export
rank_hac <- function(x, lag.max = NULL, method = c("rho_b", "tau_b", "gamma", "rho", "tau"),
                     bandwidth = NULL, kernel = c("bartlett", "parzen", "tukey-hanning"),
                     center = FALSE, psd = c("clip", "error", "none")) {
  x <- .check_x(x); H <- .check_H(lag.max, length(x)); method <- .check_method(method)
  kernel <- match.arg(kernel); psd <- match.arg(psd)
  if (!is.logical(center) || length(center) != 1L || is.na(center)) stop("center must be TRUE or FALSE.", call. = FALSE)
  raw <- .compute(x, H, projections = TRUE)
  estimate <- .estimates(raw, method)[-1L]
  result <- .hac_compute(raw, length(x), H, method, bandwidth, kernel, center, psd)
  result$vcov <- result$covariance / length(x)
  result$estimate <- stats::setNames(estimate, seq_len(H))
  result$n <- length(x); result$lag <- seq_len(H); result$method <- method
  class(result) <- "rank_hac"
  result
}
