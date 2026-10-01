# rankACFbands 0.1.1

* General sup-t calibration now inverts multivariate normal coverage with uniroot and mvtnorm, replacing empirical quantiles of simulated Gaussian maxima.
* Critical values and p-values use the same seeded probability integrator. Rejection flags derive directly from p-values.
* Added maxpts and abseps accuracy controls and integration diagnostics. nsim remains a deprecated alias for maxpts.
* Independent Gaussian calibration remains analytic. Singular covariance matrices and one-sided bands remain supported.
* Corrected the guide: the manuscript scalar test already contains the required square root.
