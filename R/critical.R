.check_alpha <- function(alpha) {
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1)
    stop("alpha must be between 0 and 1.", call. = FALSE)
}
.with_seed <- function(seed, fun) {
  if (is.null(seed)) return(fun())
  seed <- .integer_scalar(seed, "seed", 0, .Machine$integer.max)
  existed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (existed) old <- get(".Random.seed", envir = .GlobalEnv)
  on.exit(if (existed) assign(".Random.seed", old, envir = .GlobalEnv) else
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
  set.seed(seed)
  fun()
}
.correlation <- function(covariance) {
  if (!is.matrix(covariance) || !is.numeric(covariance) || nrow(covariance) < 1L ||
      nrow(covariance) != ncol(covariance) || any(!is.finite(covariance)))
    stop("covariance must be a finite nonempty numeric square matrix.", call. = FALSE)
  tol <- 1e-10 * max(abs(covariance), .Machine$double.eps)
  if (max(abs(covariance - t(covariance))) > tol)
    stop("covariance must be symmetric.", call. = FALSE)
  if (any(diag(covariance) <= 0))
    stop("Inference requires positive asymptotic variances at all requested lags; one or more estimated variances are zero or negative.", call. = FALSE)
  z <- covariance / outer(sqrt(diag(covariance)), sqrt(diag(covariance)))
  e <- eigen((z + t(z)) / 2, symmetric = TRUE)
  if (min(e$values) < -1e-8) stop("covariance is not positive semidefinite.", call. = FALSE)
  z
}
.integration_budget <- function(nsim, maxpts) {
  if (!is.null(nsim)) {
    warning("nsim is deprecated; use maxpts for the integration budget (no Gaussian draws are generated).", call. = FALSE)
    maxpts <- nsim
  }
  .integer_scalar(maxpts, "maxpts", 100L, 1e9)
}
.calibrate <- function(covariance, alpha, type, alternative, nsim, seed, iid = FALSE,
                       abseps = 1e-5) {
  .check_alpha(alpha)
  alternative <- match.arg(alternative, c("two.sided", "greater", "less"))
  type <- match.arg(type, c("sup-t", "bonferroni", "pointwise", "sidak"))
  correlation <- .correlation(covariance); H <- nrow(correlation)
  two <- alternative == "two.sided"
  independent <- max(abs(correlation - diag(H))) < 1e-12
  if (type == "sidak" && !iid && !independent)
    stop("Sidak is implemented only for independent Gaussian components; use sup-t for general HAC inference.", call. = FALSE)
  if (type == "pointwise" || H == 1L)
    return(list(critical = stats::qnorm(1 - alpha / if (two) 2 else 1), maxima = NULL,
                calibration = "pointwise-normal", H = H, two = two))
  if (type == "bonferroni")
    return(list(critical = stats::qnorm(1 - alpha / (H * if (two) 2 else 1)), maxima = NULL,
                calibration = "bonferroni", H = H, two = two))
  if (independent || type == "sidak") {
    prob <- exp(log1p(-alpha) / H)
    return(list(critical = stats::qnorm(if (two) (1 + prob) / 2 else prob), maxima = NULL,
                calibration = "exact-independent-Gaussian", H = H, two = two))
  }
  maxpts <- .integer_scalar(nsim, "maxpts", 100L, 1e9)
  if (!is.numeric(abseps) || length(abseps) != 1L || !is.finite(abseps) ||
      abseps <= 0 || abseps >= 1)
    stop("abseps must be between 0 and 1.", call. = FALSE)
  # Reuse one integration seed throughout root finding AND p-value evaluation.
  # NULL advances the caller's RNG once to choose that seed.
  integration_seed <- if (is.null(seed)) sample.int(.Machine$integer.max, 1L) else
    .integer_scalar(seed, "seed", 0, .Machine$integer.max)
  diagnostics <- new.env(parent = emptyenv())
  diagnostics$max_error <- 0
  diagnostics$warned <- FALSE
  probability <- function(c) {
    if (two && c <= 0) return(0)
    if (is.infinite(c)) return(if (c > 0) 1 else 0)
    answer <- .with_seed(integration_seed, function()
      mvtnorm::pmvnorm(lower = rep(if (two) -c else -Inf, H),
                      upper = rep(c, H), corr = correlation,
                      algorithm = mvtnorm::GenzBretz(maxpts = maxpts,
                                                    abseps = abseps, releps = 0)))
    error <- attr(answer, "error")
    message <- attr(answer, "msg")
    if (!is.finite(answer[1L]) || is.null(error) || !is.finite(error) ||
        !(message %in% c("Normal Completion", "Completion with error > abseps")))
      stop("Multivariate normal integration failed: ", message, call. = FALSE)
    diagnostics$max_error <- max(diagnostics$max_error, error)
    if (error > abseps && !diagnostics$warned) {
      warning("Multivariate normal integration did not reach abseps; increase maxpts or relax abseps. Results are approximate.", call. = FALSE)
      diagnostics$warned <- TRUE
    }
    min(1, max(0, as.numeric(answer)))
  }
  target <- 1 - alpha
  lower <- if (two) 0 else stats::qnorm(alpha, lower.tail = FALSE) - 1
  upper <- max(1, stats::qnorm(alpha / (H * if (two) 2 else 1), lower.tail = FALSE))
  # The Bonferroni endpoint brackets the exact root. Enlarge if numerical
  # integration error puts its computed coverage slightly below the target.
  while (probability(upper) < target) upper <- 2 * upper + 1
  critical <- stats::uniroot(function(c) probability(c) - target,
                             interval = c(lower, upper), tol = 1e-7)$root
  list(critical = critical, maxima = NULL, probability = probability,
       calibration = "Gaussian-MVN-root", H = H, two = two,
       maxpts = maxpts, abseps = abseps, integration_seed = integration_seed,
       diagnostics = diagnostics)
}

.cal_pvalue <- function(statistic, calibration) {
  if (!is.null(calibration$probability))
    return(1 - calibration$probability(statistic))
  tail <- if (calibration$two) min(1, 2 * stats::pnorm(-statistic)) else stats::pnorm(statistic, lower.tail = FALSE)
  if (calibration$calibration == "exact-independent-Gaussian")
    return(-expm1(calibration$H * log1p(-tail)))
  if (calibration$calibration == "bonferroni") return(min(1, calibration$H * tail))
  tail
}

#' Gaussian Sup-t Critical Value
#' @param covariance Asymptotic covariance matrix with positive diagonal.
#' @param alpha Significance level in (0,1).
#' @param alternative "two.sided", "greater", or "less". One-sided calibration
#'   uses max(Z), not max(abs(Z)); Gaussian symmetry covers either direction.
#' @param nsim Deprecated alias for maxpts, retained for old calls.
#' @param maxpts Integration budget per multivariate normal probability
#'   evaluation (at least 100), default 1000000.
#' @param abseps Requested absolute integration error, default 1e-5.
#' @param seed Local random seed, default 1; the caller's random state is restored.
#'   NULL draws one integration seed from the caller's random state. The same
#'   integration seed is reused for root finding and p-values.
#' @return Numeric critical value. Independent components and dimension one
#'   use exact normal quantiles; otherwise uniroot inverts mvtnorm::pmvnorm.
#' @details Standardizes by the diagonal, preserving cross-lag dependence.
#'   Singular positive semidefinite correlation matrices are supported. A
#'   zero diagonal is rejected. Genz-Bretz integration is randomized; no sample
#'   of Gaussian maxima or empirical quantile is used. Increase maxpts and reduce
#'   abseps to check numerical sensitivity. A warning reports unmet integration
#'   tolerance. Root tolerance is 1e-7; probability accuracy is governed by abseps.
#' @references Montiel Olea and Plagborg-Moller (2019), Simultaneous confidence
#'   bands: Theory, implementation, and an application to SVARs,
#'   Journal of Applied Econometrics 34, 1-17, doi:10.1002/jae.2656.
#' @examples
#' supt_critical(diag(5), alpha = 0.05)
#' supt_critical(matrix(c(1, 0.6, 0.6, 1), 2), maxpts = 1000000)
#' @export
supt_critical <- function(covariance, alpha = 0.05,
                          alternative = c("two.sided", "greater", "less"),
                          nsim = NULL, seed = 1, maxpts = 1000000, abseps = 1e-5) {
  alternative <- match.arg(alternative)
  maxpts <- .integration_budget(nsim, maxpts)
  .calibrate(covariance, alpha, "sup-t", alternative, maxpts, seed, abseps = abseps)$critical
}
