# rankACFbands

Rank autocorrelations and simultaneous inference for continuous, discrete, ordinal and mixed time series. Implements the estimators and inferential procedures specified in **Pohle, Weiss and Wermuth, _Rank Autocorrelations_, September 29, 2026 draft**.

The `rank_acf_sigbands()` / `rank_acf_confbands()` workflow is inspired by [ACFbands](https://github.com/TanjaZahn/ACFbands). This is an independent implementation; it does not copy or depend on that package.

## Install

Download `rankACFbands_0.1.1.tar.gz` and run:

```r
install.packages("mvtnorm")
install.packages("~/Downloads/rankACFbands_0.1.1.tar.gz",
                 repos = NULL, type = "source")
library(rankACFbands)
```

Requires R >= 4.1.0. The runtime dependency is **mvtnorm**, available from CRAN. This package itself contains no compiled code; mvtnorm normally installs as a binary on macOS and Windows. The ZIP is an editable source project, not a Windows binary package.

## Start here

```r
set.seed(42)
x <- rpois(300, lambda = 3)

# Five choices: "rho", "tau", "rho_b", "tau_b", "gamma"
rank_acf(x, lag.max = 10, method = "gamma", plot = TRUE)

# Independence: simultaneous null-centered significance bands
sig <- rank_acf_sigbands(x, lag.max = 10, method = "gamma")
as.data.frame(sig)
sig$p.value

# Independence: rank Box-Pierce test
rank_box_test(x, lag.max = 10, method = "gamma")

# Dependent data: confidence bands, estimated HAC covariance
x <- as.numeric(arima.sim(list(ar = 0.5), n = 400))
cb <- rank_acf_confbands(x, lag.max = 10, method = "rho_b",
                         maxpts = 1000000, seed = 42)
coef(cb)
vcov(cb)                 # Sampling covariance, divided by T
cb$covariance            # Asymptotic covariance of sqrt(T) times the estimator
cb$hac$psd               # Any covariance adjustment is recorded

# Lower simultaneous confidence bounds (a greater alternative)
lower <- rank_acf_confbands(x, lag.max = 10, method = "rho_b",
                            alternative = "greater", maxpts = 1000000, seed = 42)

# Joint test of a specified vector
rank_acf_test(x, lag.max = 10, method = "rho_b", null = rep(0, 10),
              covariance = "hac", maxpts = 1000000, seed = 42)

# Pointwise tests and intervals, not a simultaneous claim
pointwise <- rank_acf_confbands(x, lag.max = 10, type = "pointwise", plot = FALSE)
as.data.frame(pointwise)[3, ]    # Lag 3, including its pointwise p-value
```

## Function map

| Function | Purpose |
|---|---|
| `rank_acf()` | All five empirical rank ACFs, including lag zero |
| `rank_ties()` | Plug-in tie and double-tie probabilities |
| `rank_iid_variance()` | Corollaries 1 and 2 |
| `rank_projections()` | Estimated Hoeffding projections from Section 4.2 |
| `rank_hac()` | Primitive HAC covariance and delta-method covariance |
| `rank_acf_sigbands()` | Iid pointwise, Sidak/sup-t, and Bonferroni significance bands |
| `rank_acf_confbands()` | General pointwise, sup-t, and Bonferroni confidence bands |
| `rank_box_test()` | Iid rank Box-Pierce test; optional general HAC Wald test |
| `rank_acf_test()` | Simultaneous sup-t test of a specified vector |
| `supt_critical()` | Standalone Gaussian maximum calibration |

`plot()`, `print()`, `coef()`, `as.data.frame()` and `vcov()` methods are provided where applicable. All inferential functions work for each of the five coefficients.

## Important conventions

* Spearman uses **full-sample midranks**, with denominator T in the lagged sum. Grade ACF equals the usual denominator-T Pearson ACF of the full ranks. It is not pairwise `cor(..., method="spearman")`.
* Kendall uses a tau-a U-statistic on T-h lagged pairs. Normalized Kendall divides by the **full-series** tau(0), not separate lagged marginal tie factors. Its finite-sample estimate can occasionally exceed [-1,1]; estimates and bands are not silently clipped.
* Gamma discards pairs tied in either coordinate. An undefined denominator produces an error.
* The default HAC window is Bartlett, with bandwidth floor(2*T^(1/3)). Cross-products use exactly T-s-max(h,h') observations and denominator T.
* Pointwise inference and simultaneous inference are explicitly distinguished. Pointwise bands do not return a global rejection decision for multiple lags.
* Iid sup-t critical values are exact Sidak normal quantiles. General sup-t values use `uniroot()` to invert `mvtnorm::pmvnorm()` with the entire estimated correlation matrix. `maxpts` and `abseps` control integration accuracy. A supplied `seed` preserves the caller's random-number state.
* `alternative="greater"` gives a **lower confidence bound**, but an **upper significance threshold**. `"less"` reverses these directions.
* Missing and infinite values are rejected: deleting observations changes time spacing. Ordered factors are supported; unordered factors are rejected. Ties are exact, not tolerance-based.
* Asymptotic guarantees assume fixed H, the paper's mixing conditions, nonzero denominators, and positive marginal asymptotic variances. This package does not certify finite-sample coverage or adjust for fitted regression residuals.

## Manuscript-to-code decisions

See [the implementation guide](inst/doc/guide.html) for formulas, derivations of the efficient computation, finite-sample PSD adjustment, and draft ambiguities. In particular, the Box-Pierce statistic divides by **asymptotic variance**, and single-lag tests divide by the **standard error**. Editorial notes and proposed simulation competitors in the draft are not treated as specified methods.

## Development

```sh
R CMD build rankACFbands
R CMD check rankACFbands_0.1.1.tar.gz
```

The base-R validation suite checks estimators and projections against direct pairwise definitions, HAC entries against direct sums, Jacobians by finite differences, and inference behavior with ties, one-sided alternatives and singular covariance matrices. See `inst/examples/` for a reproducible workflow and a separate simulation diagnostic.

Version 0.1.1 is an initial research implementation. It has been checked against the specified formulas; broad finite-sample validation and CRAN submission remain separate research/release steps. The maintainer metadata is drawn from the supplied manuscript and can be edited in DESCRIPTION before publication.
