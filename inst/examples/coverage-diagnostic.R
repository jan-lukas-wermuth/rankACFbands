# Optional simulation, not run during installation or R CMD check.
# Monte Carlo diagnostics are not evidence of uniform finite-sample coverage.
library(rankACFbands)
set.seed(20260930)
repetitions <- 200
T <- 300
H <- 5
methods <- c("rho", "tau", "rho_b", "tau_b", "gamma")
rejected <- matrix(FALSE, repetitions, length(methods), dimnames = list(NULL, methods))
for (b in seq_len(repetitions)) {
  x <- rpois(T, 3)
  for (m in methods) rejected[b, m] <- rank_acf_sigbands(x, H, m, plot = FALSE)$reject
}
print(data.frame(method = methods, rejection_rate = colMeans(rejected),
                 Monte_Carlo_SE = sqrt(colMeans(rejected) * (1 - colMeans(rejected)) / repetitions)))
