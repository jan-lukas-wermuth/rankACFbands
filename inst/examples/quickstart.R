library(rankACFbands)
set.seed(42)
# Ties: iid count observations.
x <- rpois(300, 3)
for (method in c("rho", "tau", "rho_b", "tau_b", "gamma")) {
  print(rank_acf(x, lag.max = 8, method = method))
  print(rank_box_test(x, lag.max = 8, method = method))
}
sig <- rank_acf_sigbands(x, lag.max = 8, method = "gamma", plot = FALSE)
plot(sig)
# Weak dependence: continuous AR(1).
x <- as.numeric(arima.sim(list(ar = 0.5), n = 400))
cb <- rank_acf_confbands(x, lag.max = 8, method = "rho_b", maxpts = 1000000,
                         seed = 42, plot = FALSE)
plot(cb)
print(as.data.frame(cb))
print(vcov(cb))
print(cb$hac$psd)
lower <- rank_acf_confbands(x, lag.max = 8, method = "rho_b",
                            alternative = "greater", maxpts = 1000000, seed = 42,
                            plot = FALSE)
plot(lower)
# Ordered categorical observations: user-specified ordering is respected.
ordinal <- ordered(c("poor", "good", "fair", "fair", "good", "poor", "fair", "good"),
                   levels = c("poor", "fair", "good"))
print(rank_acf(ordinal, lag.max = 2, method = "tau_b"))
