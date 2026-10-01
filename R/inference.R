.null_vector <- function(null, H) {
  if (!is.numeric(null) || any(!is.finite(null)) || !(length(null) %in% c(1L, H)))
    stop("null must be a finite scalar or a vector with lag.max entries.", call. = FALSE)
  rep(null, length.out = H)
}
.band_result <- function(estimate, covariance, n, method, alpha, type,
                         alternative, nsim, seed, null, significance, iid, hac = NULL, abseps = 1e-5) {
  H <- length(estimate); null <- .null_vector(null, H)
  cal <- .calibrate(covariance, alpha, type, alternative, nsim, seed, iid, abseps)
  se <- sqrt(diag(covariance) / n)
  critical <- cal$critical
  at <- if (significance) null else estimate
  lower <- at - critical * se; upper <- at + critical * se
  if (alternative == "greater") {
    if (significance) lower[] <- -Inf else upper[] <- Inf
  }
  if (alternative == "less") {
    if (significance) upper[] <- Inf else lower[] <- -Inf
  }
  z <- (estimate - null) / se
  directed <- switch(alternative, two.sided = abs(z), greater = z, less = -z)
  point_p <- if (alternative == "two.sided") 2 * stats::pnorm(-abs(z)) else stats::pnorm(directed, lower.tail = FALSE)
  adjusted <- vapply(directed, .cal_pvalue, numeric(1), calibration = cal)
  global <- if (type == "pointwise" && H > 1L) NA_real_ else .cal_pvalue(max(directed), cal)
  tab <- data.frame(lag = seq_len(H), estimate = as.numeric(estimate), se = as.numeric(se),
                    lower = as.numeric(lower), upper = as.numeric(upper), null = null,
                    statistic = as.numeric(z), p.value = as.numeric(point_p),
                    adjusted.p.value = adjusted, reject = adjusted < alpha)
  structure(list(table = tab, estimate = stats::setNames(as.numeric(estimate), seq_len(H)),
    covariance = covariance, vcov = covariance / n, critical = critical,
    alpha = alpha, type = type, alternative = alternative,
    band = if (significance) "significance" else "confidence", null = null,
    p.value = global, reject = if (is.na(global)) NA else global < alpha,
    calibration = cal$calibration, nsim = 0L,
    maxpts = if (is.null(cal$maxpts)) 0L else cal$maxpts,
    abseps = cal$abseps, integration.seed = cal$integration_seed,
    integration.error = if (is.null(cal$diagnostics)) 0 else cal$diagnostics$max_error,
    seed = seed, method = method, n = n, lag = seq_len(H), hac = hac,
    simultaneous = type != "pointwise" || H == 1L), class = "rank_bands")
}

#' Independence Significance Bands for Rank Autocorrelations
#' @inheritParams rank_acf
#' @inheritParams supt_critical
#' @param type "sup-t" (default), "sidak", "bonferroni", or "pointwise".
#'   Under iid sampling sup-t equals the exact Gaussian Sidak calibration.
#' @return A \code{rank_bands} object. Its \code{table} contains lag, estimate,
#'   standard error, lower/upper bounds, null value, pointwise z statistic,
#'   pointwise and adjusted p-values, and rejection indicators. The object also
#'   contains the covariance, critical value and global test result.
#' @details These are null-centered non-rejection bands for serial independence,
#'   not confidence intervals for a possibly nonzero autocorrelation. For
#'   alternative="greater", rejection occurs above the finite upper threshold;
#'   for "less", below the finite lower threshold. Pointwise bands do not
#'   control familywise error; their global p.value and reject are NA for H>1.
#' @examples
#' set.seed(3)
#' rank_acf_sigbands(rpois(200, 3), lag.max = 6, method = "gamma", plot = FALSE)
#' @export
rank_acf_sigbands <- function(x, lag.max = NULL,
    method = c("rho_b", "tau_b", "gamma", "rho", "tau"), alpha = 0.05,
    type = c("sup-t", "sidak", "bonferroni", "pointwise"),
    alternative = c("two.sided", "greater", "less"), plot = TRUE, ...) {
  cl <- match.call(); method <- .check_method(method)
  type <- match.arg(type); alternative <- match.arg(alternative)
  fit <- rank_acf(x, lag.max, method)
  variance <- rank_iid_variance(x, method)
  covariance <- diag(variance, length(fit$lag))
  dimnames(covariance) <- list(names(fit$estimate), names(fit$estimate))
  result <- .band_result(fit$estimate, covariance, fit$n, method, alpha, type,
                         alternative, 0, NULL, 0, TRUE, TRUE)
  result$call <- cl
  if (plot) graphics::plot(result, ...)
  result
}

#' HAC Confidence Bands and General Hypothesis Tests
#' @inheritParams rank_acf
#' @inheritParams rank_hac
#' @inheritParams supt_critical
#' @param type "sup-t", "bonferroni", or "pointwise".
#' @param null Hypothesized autocorrelations, a scalar recycled across lags or
#'   a vector of length lag.max; defaults to zero. Does not affect band location.
#' @return A \code{rank_bands} object; see \code{rank_acf_sigbands}. The
#'   \code{hac} component contains full covariance diagnostics. Greater
#'   alternatives produce lower confidence bounds; less alternatives produce
#'   upper bounds. The unbounded endpoint is represented by Inf or -Inf.
#' @details General confidence bands are centered on the estimates and use the
#'   HAC matrix without imposing independence. Off-diagonal dependence is
#'   retained when integrating the standardized Gaussian distribution. At a positive
#'   null estimate-minus-null is tested; one-sided alternatives are selected
#'   before examining the data. Bands are not clipped to [-1,1]. For fixed H,
#'   coverage is asymptotic under the manuscript assumptions, positive marginal
#'   variances and consistent covariance estimation. Critical values invert
#'   multivariate normal coverage using uniroot; maxpts and abseps control
#'   integration accuracy. Rejection flags use p.value < alpha; comparisons
#'   to band endpoints can differ within numerical tolerance. No residual-estimation
#'   adjustment for fitted time-series models is implemented.
#' @examples
#' set.seed(4)
#' x <- as.numeric(arima.sim(list(ar = 0.4), n = 200))
#' fit <- rank_acf_confbands(x, lag.max = 5, plot = FALSE)
#' as.data.frame(fit)
#' plot(fit)
#' @export
rank_acf_confbands <- function(x, lag.max = NULL,
    method = c("rho_b", "tau_b", "gamma", "rho", "tau"), alpha = 0.05,
    type = c("sup-t", "bonferroni", "pointwise"),
    alternative = c("two.sided", "greater", "less"), null = 0,
    bandwidth = NULL, kernel = c("bartlett", "parzen", "tukey-hanning"),
    center = FALSE, psd = c("clip", "error", "none"), nsim = NULL,
    seed = 1, plot = TRUE, maxpts = 1000000, abseps = 1e-5, ...) {
  cl <- match.call(); method <- .check_method(method); type <- match.arg(type)
  alternative <- match.arg(alternative); kernel <- match.arg(kernel); psd <- match.arg(psd)
  maxpts <- .integration_budget(nsim, maxpts)
  hac <- rank_hac(x, lag.max, method, bandwidth, kernel, center, psd)
  result <- .band_result(hac$estimate, hac$covariance, hac$n, method, alpha, type,
                         alternative, maxpts, seed, null, FALSE, FALSE, hac, abseps)
  result$call <- cl
  if (plot) graphics::plot(result, ...)
  result
}

#' Rank Box-Pierce and General HAC Wald Tests
#' @inheritParams rank_acf
#' @inheritParams rank_hac
#' @param covariance "iid" for the rank Box-Pierce independence test or "hac"
#'   for a general Wald test of the vector of rank autocorrelations.
#' @param null Scalar or lag.max-vector of hypothesized coefficients. Must be
#'   zero for covariance="iid".
#' @return An \code{htest} object with statistic, degrees of freedom, p.value,
#'   estimates, null values, and the asymptotic covariance matrix.
#' @details Under iid sampling Q=T*sum(estimate^2)/asymptotic_variance and
#'   Q has a limiting chi-square distribution with lag.max degrees of freedom.
#'   The divisor is a variance, not a standard deviation. The HAC extension
#'   uses T*(estimate-null)'*solve(covariance)*(estimate-null). It requires a
#'   positive definite matrix; singular matrices are rejected rather than
#'   silently assigned an unjustified chi-square reference distribution.
#'   Failure to reject zero rank autocorrelations does not establish independence.
#' @examples
#' set.seed(5)
#' rank_box_test(rpois(200, 4), lag.max = 5, method = "tau")
#' @export
rank_box_test <- function(x, lag.max = NULL,
    method = c("rho_b", "tau_b", "gamma", "rho", "tau"),
    covariance = c("iid", "hac"), null = 0, bandwidth = NULL,
    kernel = c("bartlett", "parzen", "tukey-hanning"), center = FALSE,
    psd = c("clip", "error", "none")) {
  data.name <- deparse(substitute(x)); method <- .check_method(method)
  covariance <- match.arg(covariance); kernel <- match.arg(kernel); psd <- match.arg(psd)
  if (covariance == "iid") {
    fit <- rank_acf(x, lag.max, method)
    null <- .null_vector(null, length(fit$lag))
    if (any(null != 0)) stop("The iid Box-Pierce test requires null=0; use covariance='hac' for general hypotheses.", call. = FALSE)
    v <- diag(rank_iid_variance(x, method), length(fit$lag))
    estimate <- fit$estimate; n <- fit$n
  } else {
    fit <- rank_hac(x, lag.max, method, bandwidth, kernel, center, psd)
    v <- fit$covariance; estimate <- fit$estimate; n <- fit$n
    null <- .null_vector(null, length(estimate))
  }
  e <- eigen(v, symmetric = TRUE)
  if (min(e$values) <= max(e$values) * 1e-10 || max(e$values) <= 0)
    stop("Wald inference requires a nonsingular positive definite covariance; reduce lag.max or use sup-t bands if marginal variances are positive.", call. = FALSE)
  difference <- as.numeric(estimate) - null
  statistic <- n * sum((crossprod(e$vectors, difference)[, 1L])^2 / e$values)
  structure(list(statistic = c(Q = statistic), parameter = c(df = length(estimate)),
    p.value = stats::pchisq(statistic, length(estimate), lower.tail = FALSE),
    method = paste(if (covariance == "iid") "Rank Box-Pierce independence test" else "HAC Wald rank-autocorrelation test", paste0("(", method, ")")),
    data.name = data.name, alternative = "two.sided", estimate = estimate,
    null.value = stats::setNames(null, names(estimate)), covariance = v,
    covariance.type = covariance), class = "htest")
}

#' Sup-t Test of Rank Autocorrelations
#' @inheritParams rank_acf
#' @inheritParams rank_acf_confbands
#' @param covariance "hac" for arbitrary null vectors or "iid" for independence.
#' @param ... For covariance="hac", arguments forwarded to
#'   \code{rank_acf_confbands}, e.g. bandwidth, maxpts, abseps and seed. For "iid", no
#'   additional arguments are accepted.
#' @return An \code{htest} object, with the associated bands in \code{bands}.
#' @examples
#' set.seed(6)
#' rank_acf_test(rnorm(150), lag.max = 4, covariance = "iid")
#' @export
rank_acf_test <- function(x, lag.max = NULL,
    method = c("rho_b", "tau_b", "gamma", "rho", "tau"), null = 0,
    covariance = c("hac", "iid"), alpha = 0.05,
    alternative = c("two.sided", "greater", "less"), ...) {
  data.name <- deparse(substitute(x)); method <- .check_method(method)
  covariance <- match.arg(covariance); alternative <- match.arg(alternative)
  if (covariance == "iid") {
    if (length(list(...))) stop("Additional HAC arguments are not used for covariance='iid'.", call. = FALSE)
    if (!is.numeric(null) || any(!is.finite(null)) || any(null != 0))
      stop("iid testing requires null=0.", call. = FALSE)
    bands <- rank_acf_sigbands(x, lag.max, method, alpha, "sup-t", alternative, plot = FALSE)
    .null_vector(null, length(bands$lag))
  } else {
    bands <- rank_acf_confbands(x, lag.max, method, alpha, "sup-t", alternative,
                                null = null, plot = FALSE, ...)
  }
  z <- bands$table$statistic
  statistic <- max(switch(alternative, two.sided = abs(z), greater = z, less = -z))
  structure(list(statistic = c(`sup-t` = statistic), p.value = bands$p.value,
    method = paste("Simultaneous rank-autocorrelation test", paste0("(", method, "; ", covariance, ")")),
    data.name = data.name, alternative = alternative,
    estimate = bands$estimate, null.value = stats::setNames(bands$null, names(bands$estimate)),
    critical.value = bands$critical, reject = bands$reject, bands = bands), class = "htest")
}
