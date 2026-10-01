#' Rank Autocorrelations and Simultaneous Inference
#'
#' Estimates five rank autocorrelation functions for continuous, discrete,
#' ordinal and mixed data. Provides iid significance bands and Box-Pierce
#' tests, HAC covariance estimation, and Gaussian simultaneous sup-t bands.
#' Start with \code{rank_acf}, \code{rank_acf_sigbands}, or
#' \code{rank_acf_confbands}. See the installed guide in
#' \code{system.file("doc", "guide.html", package="rankACFbands")}.
#' @importFrom stats coef vcov
#' @importFrom utils head
#' @keywords package
"_PACKAGE"
