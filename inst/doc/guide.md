# rankACFbands: implementation and usage guide

Version 0.1.1 — implementation based on the September 29, 2026 manuscript.

## 1. Scope and assumptions

The package implements the statistical methods specified in Sections 2–4 of *Rank Autocorrelations* by Marc-Oliver Pohle, Christian H. Weiss and Jan-Lukas Wermuth. The manuscript is a draft: its simulation/application sections contain proposals rather than completed methods. Those editorial suggestions are not executable specifications and are not implemented as additional estimators. The package does not bundle the private manuscript or claim endorsement by its authors.

All five methods accept a regularly spaced univariate numeric series or an ordered factor. Strictly monotone transformations preserve the rank coefficients, including simultaneous reversal of the ordering. Missing or infinite values and unordered factors are rejected. To handle gaps, choose and justify a time-series preprocessing strategy yourself rather than deleting rows inside the estimator. No moment assumption on the observations is imposed by the code.

The theoretical inference requires the stated stationarity, ergodicity and mixing assumptions, fixed H, nonzero normalization denominators, and positive marginal limiting variances. HAC consistency requires a growing bandwidth that is o(sqrt(T)). Small-sample performance is not guaranteed. These routines do not include adjustments for residual estimation in dynamic or static regressions.

## 2. Installation and function selection

Install the source tarball with `install.packages(path, repos=NULL, type="source")`. All runtime dependencies are standard R packages. No compiler is needed. The ZIP is an editable source directory; open its `.Rproj` in RStudio if you want to develop the package.

* `rank_acf`: estimation only, returning lag zero and lags 1,...,H.
* `rank_acf_sigbands`: iid null-centered non-rejection bands.
* `rank_acf_confbands`: general estimate-centered confidence bands.
* `rank_acf_test`: sup-t test, returning an `htest` object.
* `rank_box_test`: iid rank Box-Pierce or general HAC Wald test.
* `rank_hac`: covariance and diagnostics without band calibration.
* `rank_projections`, `rank_ties`, `rank_iid_variance`: lower-level quantities.
* `supt_critical`: critical value from an arbitrary positive-diagonal PSD covariance.

The API follows the workflow of ACFbands, but it is not a drop-in replacement and does not copy its code. Supply `method="rho"`, `"tau"`, `"rho_b"`, `"tau_b"`, or `"gamma"`; the default is grade (`"rho_b"`). To compare methods, loop over these five strings. `plot=FALSE` returns all results without drawing.

## 3. Exact finite-sample definitions

Let T be the series length and G-hat(X_t)=(midrank(X_t)-1/2)/T, using the full sample.

* **rho(h)** = (12/T) sum from t=1 to T-h of (G-hat(X_t)-1/2)(G-hat(X_(t+h))-1/2).
* **tau(h)** = the average of sign(X_s-X_t) sign(X_(s+h)-X_(t+h)) over unordered pairs 1 <= s < t <= T-h.
* **nu(h)** = the fraction of those unordered pairs tied in at least one coordinate.
* **rho_b(h)** = rho(h)/rho(0).
* **tau_b(h)** = tau(h)/tau(0).
* **gamma(h)** = tau(h)/(1-nu(h)).

These definitions match the displayed formulas in Definitions 4–5. Rho uses denominator T, whereas tau and nu are U-statistics on T-h pairs. Rho_b equals `acf(rank(x), plot=FALSE)` with the same maximum lag. Pairwise `cor(head(x,-h), tail(x,-h), method="spearman")` instead reranks within the two lagged samples and generally differs.

Similarly, base R's Kendall correlation is pairwise tau-b, not the paper's tau-a. The paper's tau_b normalizes by the full-series tau(0). Its finite-sample value can exceed [-1,1] because the numerator and denominator use different samples. We preserve the defined estimator, without clipping. Normalized estimators equal one at lag zero; raw rho(0) and tau(0) need not.

`rank_ties` uses plug-in probabilities zeta=sum(p_j^2), zeta2=sum(p_j^3). In contrast, tau(0) uses the off-diagonal U-statistic tie correction. These distinct finite-sample conventions are intentional.

## 4. Efficient projection estimation

The implementation computes every projection in Section 4.2, without an extra order factor. For h>0 it uses the global empirical marginal distribution and the empirical pair distribution of (X_t,X_(t+h)), t=1,...,T-h. It uses the manuscript's separate formulas for h=0.

The bivariate mid-distribution is the empirical average of the product of two mid-step functions. A Fenwick cumulative-count tree computes these quantities at the observed pairs in O(T log T) time per lag, including tied values, without allocating T by T matrices. Marginal and joint tie counts are exact.

For the additional rho terms, write the mid-step as 1{q>z} + (1/2)1{q=z}. For any observed value v, the average of mid-step(X_t-v) over the full sample is 1-G-hat(v). Thus the average of G-hat_h(x,X_t) can be evaluated by a weighted cumulative sum over the lagged pairs. This algebra is exact for the empirical distributions, including ties, and avoids a triple loop. The test suite separately evaluates the original nested empirical averages and compares every projection.

Estimated projections are not necessarily exactly sample-centered. The default `center=FALSE` retains the manuscript's literal plug-in formula. `center=TRUE` is an optional finite-sample modification and is recorded in the result.

## 5. HAC covariance and delta method

For primitive projections a(h), b(h'), let h-star=max(h,h'). The estimator uses:

1. Contemporaneous products summed over t=1,...,T-h-star and divided by T.
2. For each s=1,...,bandwidth, both directional products summed over t=1,...,T-s-h-star and divided by T. Empty ranges contribute zero.
3. Weights w(s/(bandwidth+1)).
4. The multiplier r(a)r(b) on the **whole sum**, with r(rho)=3 and r(tau)=r(nu)=2.

Default Bartlett weights are 1-s/(bandwidth+1), and bandwidth=floor(2*T^(1/3)), capped at T-1. Parzen and Tukey-Hanning windows are also available. A zero bandwidth is supported for inspection, but it is not the iid independence variance: projections for different autocorrelation lags share observations.

Primitive ordering is rho.0,...,rho.H or tau.0,...,tau.H for normalized estimators and tau.1,nu.1,tau.2,nu.2,... for gamma. `jacobian` exposes the exact analytic derivatives. The final matrix is J Sigma J'. The tau_b and rho_b cross term has a negative sign; the gamma cross term has a positive sign.

`covariance` is the asymptotic matrix for sqrt(T)(estimate-parameter). `vcov(result)` divides by T. Standard errors are sqrt(diag(covariance)/T). `primitive_covariance`, `raw_covariance`, `jacobian` and `projections` remain available for inspection.

### Finite-sample positive semidefiniteness

The manuscript's pair-specific end truncations do not form a single common-sample Gram matrix. As a result, a finite-sample covariance can be indefinite, even with Bartlett weights. Default `psd="clip"` replaces negative eigenvalues of the final transformed matrix by zero. Material adjustments produce a warning. The result retains the raw matrix, its minimum eigenvalue, an adjustment flag and relative Frobenius change. No positive ridge is added.

Use `psd="error"` to reject material indefiniteness or `psd="none"` to inspect the literal raw estimate. Gaussian calibration rejects indefinite matrices and zero/negative marginal variances. A singular PSD matrix with positive diagonal is supported for sup-t integration. The HAC Wald test requires a positive definite matrix and refuses to assign a chi-square law to an arbitrary pseudoinverse.

## 6. Iid significance bands and Box-Pierce

The scalar asymptotic variances under independence are:

| Method | Asymptotic variance |
|---|---|
| rho | (1-zeta2)^2 |
| tau | (4/9)(1-zeta2)^2 |
| rho_b | 1 |
| tau_b | (4/9)(1-zeta2)^2/(1-zeta)^2 |
| gamma | (4/9)(1-zeta2)^2/(1-zeta)^4 |

The joint covariance is this variance times I_H. Pointwise critical values are ordinary normal quantiles. Two-sided simultaneous sup-t/Sidak critical values are qnorm((1+(1-alpha)^(1/H))/2). One-sided values are qnorm((1-alpha)^(1/H)). Bonferroni is available as an alternative.

The rank Box-Pierce statistic is T*sum(estimated_acf^2)/estimated_asymptotic_variance, compared to chi-square(H). It is a test against departures visible in the selected rank autocorrelations, not a characterization of every possible dependence alternative. The optional HAC Wald test generalizes the quadratic form to arbitrary specified rank-correlation vectors and is an explicit extension of the draft's editorial proposal.

## 7. General and one-sided inference

Sup-t calibration standardizes the estimated covariance to a correlation matrix and solves P(-c <= Z_h <= c for every h) = 1-alpha using stats::uniroot and mvtnorm::pmvnorm. One-sided calibration solves P(Z_h <= c for every h) = 1-alpha. P-values are one minus the same joint coverage evaluated at the observed maximum. For independent components or H=1 it uses exact normal quantiles. The Genz-Bretz probability integrator is randomized and supports singular PSD matrices; one fixed integration seed is reused throughout the calculation. No empirical Gaussian-maximum quantile is used. Rejection flags use p.value < alpha, so they agree with reported p-values. Band crossings can differ within integration/root tolerance. Results record maxpts, abseps, integration.seed and the largest estimated integration error across evaluations; an unmet tolerance generates a warning. The legacy nsim argument is a deprecated alias for maxpts.

For one-sided confidence bands, the critical value is the quantile of max(Z_h). A greater alternative returns estimate-critical*SE as the lower bound, with upper=Inf. A less alternative returns estimate+critical*SE as the upper bound, with lower=-Inf. Centered Gaussian symmetry makes the critical value equal for these two choices, but each has a separate 1-alpha guarantee; their intersection is not a two-sided 1-alpha band. Select direction before inspecting the estimates.

For significance bands under independence, greater means an upper rejection threshold about zero, with lower=-Inf; less means a lower threshold, with upper=Inf. This difference from confidence-band direction is deliberate.

`table$p.value` is pointwise. `table$adjusted.p.value` uses the chosen simultaneous calibration (or is unadjusted for pointwise bands). The object-level p.value is the joint maximum test. For multi-lag pointwise bands it is NA, as is the object-level rejection decision. An arbitrary scalar or vector `null` can be supplied to general confidence bands/tests. The band itself remains centered on the estimates.

## 8. Draft ambiguities resolved in this version

1. Definition 4 mixes n and T in the rank representation; code uses T throughout and follows the mid-distribution expression. Rho uses denominator T, as displayed. The adjacent editorial question is not an instruction to change the formula.
2. Section 3.4 prints a standard-deviation symbol in the Box-Pierce denominator. The chi-square derivation requires its square, the variance. That is what the package uses.
3. Section 4.3 correctly uses a square root in the scalar test denominator. The earlier guide incorrectly stated otherwise. Code divides by sqrt(covariance/T), in agreement with the manuscript.
4. The long-run covariance display is interpreted as an expectation of a product of centered projections, and the entire weighted sum carries the U-statistic order factors.
5. The proposal to reject after more than alpha*H pointwise exceedances does not control size in general and is not implemented. Simultaneous bands use their calibrated maximum instead.
6. The manuscript's simulation competitors and hypothetical empirical applications are not implemented as if they were new methods. One-sided bands, Bonferroni, additional HAC windows and the general HAC Wald test are documented extensions.
7. The manuscript assumes valid nondegenerate ratios. Undefined sample denominators and zero marginal standard errors produce informative errors instead of silent NaNs or invented confidence intervals.

## 9. Reproducibility, speed and release status

The package code is R and uses mvtnorm for joint Gaussian probabilities. Estimation is O(H*T*log(T)); computing all projections stores O(H*T) values. HAC computation is O(H^2*T*bandwidth). Gaussian integration costs depend on maxpts, abseps and H. Very large H or T can take appreciable time; the asymptotic results concern fixed H rather than a growing-dimensional regime.

A non-NULL seed preserves the caller's random state. Numerical eigenvalue truncation is reported, never hidden. Defaults are maxpts=1000000 and abseps=1e-5, with root tolerance 1e-7. Increase maxpts and reduce abseps to check sensitivity near a rejection threshold. Genz-Bretz error estimates are numerical diagnostics, not rigorous deterministic error bounds.

The automated validation compares optimized code with direct definitions, exact HAC sums and numerical derivatives. It checks iid variance scaling, exact Sidak values, positive and negative monotone transformations, ties, one-sided direction, RNG restoration, invalid inputs, singular covariance support and plotting. `coverage-diagnostic.R` is a reproducible supplementary experiment, not a formal proof of coverage. This is version 0.1.1, a research implementation rather than a CRAN-reviewed release. Inspect the source and resolve any future manuscript changes before publication.

## References

* Pohle, M.-O., Weiss, C. H. and Wermuth, J.-L. (2026). *Rank Autocorrelations*. Supplied draft dated September 29.
* Montiel Olea, J. L. and Plagborg-Moller, M. (2019). Simultaneous confidence bands: Theory, implementation, and an application to SVARs. *Journal of Applied Econometrics* 34, 1–17. DOI: 10.1002/jae.2656.
* ACFbands interface reference: https://github.com/TanjaZahn/ACFbands. No source code was copied.
