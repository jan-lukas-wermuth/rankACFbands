#' Methods for Rank Autocorrelation Results
#' @param x,object A rank_acf, rank_bands, or rank_hac result as appropriate.
#' @param digits Number of printed digits.
#' @param ... Additional arguments; plotting arguments are passed to base plot.
#' @param row.names,optional Arguments for as.data.frame.
#' @param main,xlab,ylab,ylim Base graphics labels and limits.
#' @param col,border Colors of estimates and band boundaries.
#' @return Print methods return the input invisibly; plot methods return the
#'   input invisibly. coef returns estimates, vcov returns sampling covariance
#'   (asymptotic covariance divided by T), and as.data.frame returns a table.
#' @name rank_methods
NULL

#' @rdname rank_methods
#' @export
print.rank_acf <- function(x, digits = 4, ...) {
  cat("Rank autocorrelation:", x$method, "  T =", x$n, "\n")
  print(data.frame(lag = 0:length(x$lag), estimate = as.numeric(x$acf)), digits = digits, row.names = FALSE)
  invisible(x)
}
#' @rdname rank_methods
#' @export
coef.rank_acf <- function(object, ...) object$estimate
#' @rdname rank_methods
#' @export
as.data.frame.rank_acf <- function(x, row.names = NULL, optional = FALSE, ...) {
  data.frame(lag = 0:length(x$lag), estimate = as.numeric(x$acf), row.names = row.names)
}
#' @rdname rank_methods
#' @export
plot.rank_acf <- function(x, main = paste("Rank autocorrelation:", x$method),
                          xlab = "Lag", ylab = "Autocorrelation", ylim = NULL,
                          col = "#145A70", ...) {
  if (is.null(ylim)) ylim <- range(c(0, x$estimate))
  graphics::plot(x$lag, x$estimate, type = "h", main = main, xlab = xlab,
                 ylab = ylab, ylim = ylim, col = col, ...)
  graphics::points(x$lag, x$estimate, pch = 16, col = col)
  graphics::abline(h = 0, col = "grey60")
  invisible(x)
}
#' @rdname rank_methods
#' @export
print.rank_bands <- function(x, digits = 4, ...) {
  cat(sprintf("%s rank %s bands (%s, %s); T = %d\n",
              x$method, x$band, x$type, x$alternative, x$n))
  cat(sprintf("Nominal %s coverage: %.1f%%; critical value: %.4f\n",
              if (x$simultaneous) "simultaneous" else "pointwise", 100 * (1 - x$alpha), x$critical))
  print(x$table, digits = digits, row.names = FALSE)
  if (x$simultaneous) cat("Global p-value:", format.pval(x$p.value, digits = digits),
                          "; reject:", x$reject, "\n")
  else cat("Pointwise bands do not give a global rejection decision.\n")
  if (!is.null(x$maxpts) && x$maxpts > 0)
    cat("Integration budget:", x$maxpts, "; requested error:", x$abseps, "\n")
  invisible(x)
}
#' @rdname rank_methods
#' @export
coef.rank_bands <- function(object, ...) object$estimate
#' @rdname rank_methods
#' @export
vcov.rank_bands <- function(object, ...) object$vcov
#' @rdname rank_methods
#' @export
as.data.frame.rank_bands <- function(x, row.names = NULL, optional = FALSE, ...) {
  result <- x$table
  if (!is.null(row.names)) rownames(result) <- row.names
  result
}
#' @rdname rank_methods
#' @export
plot.rank_bands <- function(x,
    main = paste(x$method, x$type, x$band, "bands"),
    xlab = "Lag", ylab = "Autocorrelation", ylim = NULL,
    col = "#145A70", border = "#B65A36", ...) {
  tab <- x$table
  if (is.null(ylim)) {
    vals <- c(0, tab$estimate, tab$lower, tab$upper, tab$null)
    ylim <- range(vals[is.finite(vals)])
    if (diff(ylim) == 0) ylim <- ylim + c(-0.1, 0.1)
    ylim <- ylim + c(-0.03, 0.22) * diff(ylim)
  }
  graphics::plot(tab$lag, tab$estimate, type = "h", main = main,
                 xlab = xlab, ylab = ylab, ylim = ylim, col = col, ...)
  graphics::abline(h = 0, col = "grey70")
  graphics::points(tab$lag, tab$estimate, pch = 16, col = col)
  for (bound in c("lower", "upper")) {
    keep <- is.finite(tab[[bound]])
    if (any(keep)) {
      graphics::lines(tab$lag[keep], tab[[bound]][keep], col = border, lty = 2, lwd = 2)
      graphics::points(tab$lag[keep], tab[[bound]][keep], pch = 3, col = border)
    }
  }
  if (x$band == "confidence") {
    keep <- is.finite(tab$lower) & is.finite(tab$upper)
    if (any(keep)) graphics::segments(tab$lag[keep], tab$lower[keep], tab$lag[keep],
                                      tab$upper[keep], col = grDevices::adjustcolor(border, 0.5))
  }
  graphics::legend("topright", legend = c("Estimate", paste0(100 * (1 - x$alpha), "% ",
    if (x$simultaneous) "simultaneous" else "pointwise", " ",
    if (x$alternative == "two.sided") "band" else if (x$band == "confidence")
      if (x$alternative == "greater") "lower bound" else "upper bound"
    else if (x$alternative == "greater") "upper threshold" else "lower threshold")),
    col = c(col, border), lty = c(1, 2), bty = "n", cex = 0.8)
  invisible(x)
}
#' @rdname rank_methods
#' @export
vcov.rank_hac <- function(object, ...) object$vcov
#' @rdname rank_methods
#' @export
coef.rank_hac <- function(object, ...) object$estimate
#' @rdname rank_methods
#' @export
print.rank_hac <- function(x, digits = 4, ...) {
  cat("HAC asymptotic covariance:", x$method, "  T =", x$n,
      "  bandwidth =", x$bandwidth, "  kernel =", x$kernel, "\n")
  print(x$covariance, digits = digits)
  if (x$psd$adjusted) cat("Eigenvalue adjustment applied; relative change:", x$psd$relative.adjustment, "\n")
  invisible(x)
}
