library(rankACFbands)
calibrate <- getFromNamespace(".calibrate", "rankACFbands")
pvalue <- getFromNamespace(".cal_pvalue", "rankACFbands")
band_result <- getFromNamespace(".band_result", "rankACFbands")
checks <- 0L
check <- function(x, label) {
  if (!isTRUE(x)) stop(label)
  checks <<- checks + 1L
}
# Independent reference: one-dimensional conditional-normal quadrature.
coverage <- function(c, r, two) {
  if (two && c <= 0) return(0)
  integrate(function(x) dnorm(x) *
    (pnorm((c-r*x)/sqrt(1-r*r)) -
       if (two) pnorm((-c-r*x)/sqrt(1-r*r)) else 0),
    lower=if (two) -c else -Inf, upper=c, rel.tol=1e-10)$value
}
for (r in c(-.7, .5, .95)) for (two in c(TRUE,FALSE)) for (alpha in c(.05,.2,.8)) {
  v <- matrix(c(1,r,r,1),2)
  alternative <- if (two) "two.sided" else "greater"
  cal <- calibrate(v, alpha, "sup-t", alternative, 1000000, 11, abseps=1e-7)
  check(abs(coverage(cal$critical,r,two)-(1-alpha))<2e-6,"quadrature coverage")
  for (offset in c(-.01,.01)) {
    observed <- cal$critical+offset
    check(abs(pvalue(observed,cal)-(1-coverage(observed,r,two)))<2e-6,"quadrature p-value")
  }
  check(abs(pvalue(cal$critical,cal)-alpha)<2e-6,"root p-value")
}
# Singular limits are known analytically.
check(abs(supt_critical(matrix(1,3,3))-qnorm(.975))<1e-5,"positive singular")
check(abs(supt_critical(matrix(c(1,-1,-1,1),2),alternative="greater")-qnorm(.975))<1e-5,"negative singular one-sided")
# The former near-cutoff rejection/p-value contradiction cannot recur.
v <- matrix(c(1,.5,.5,1),2)
c <- supt_critical(v)
for (offset in c(-1e-8,1e-8,-.01,.01)) {
  b <- band_result(c(c+offset,0)/10,v,100,"rho",.05,"sup-t",
                   "two.sided",1000000,1,0,FALSE,FALSE)
  check(identical(b$reject,b$p.value<.05),"global rejection consistent")
  check(all(b$table$reject==(b$table$adjusted.p.value<.05)),"lag rejection consistent")
  check(b$nsim==0 && b$maxpts==1000000,"integration metadata")
}
# A fixed seed restores the RNG; NULL chooses one reusable seed.
set.seed(123); old <- .Random.seed
v <- matrix(.35,4,4); diag(v) <- 1
cal <- calibrate(v,.05,"sup-t","two.sided",1000000,7)
check(identical(old,.Random.seed),"RNG restored")
check(identical(pvalue(2.5,cal),pvalue(2.5,cal)),"repeatable probability")
cal <- calibrate(v,.05,"sup-t","two.sided",1000000,NULL)
check(!identical(old,.Random.seed),"NULL seed advances RNG")
check(identical(pvalue(2.5,cal),pvalue(2.5,cal)),"NULL seed reused internally")
check(abs(supt_critical(diag(4))-qnorm((1+.95^(1/4))/2))<1e-12,"analytic path unchanged")
check(inherits(tryCatch(supt_critical(v,abseps=-1),error=identity),"error"),"invalid accuracy rejected")
# Insufficient integration budgets are visible to the caller.
v <- matrix(.6,8,8); diag(v) <- 1
warned <- FALSE
withCallingHandlers(supt_critical(v,maxpts=100,abseps=1e-12),
 warning=function(w) { warned <<- TRUE; invokeRestart("muffleWarning") })
check(warned,"integration tolerance warning")
cat("Validated",checks,"root-finding and probability checks.\n")
